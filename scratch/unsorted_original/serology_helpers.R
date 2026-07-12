# =====================================================================
# R/serology_helpers.R
# ---------------------------------------------------------------------
# Shared machinery for the systems-serology analysis, sourced by every
# analysis part. Sourcing config/endpoints.R first is required (it
# defines VISIT_RECODE, the feature vocabulary, arms, colours, etc.).
#
# Contents
#   1.  Small utilities         (pvalue_to_label, ratio_col, recode_visits)
#   2.  Data loading + recoding (load_serology_data)
#   3.  Participant / N tracker (requirement #4)
#   4.  Core model runner       (run_lm) — tracker-aware
#   5.  Display helpers         (show_model, r2_row)
#   6.  Derived chain variables (add_chain_derived) — visit-name driven
#   7.  Chain block runner      (run_chain_blocks) — all 6 blocks, 1 stratum
#   8.  Strata builder          (build_maternal_strata) — requirement #5
#   9.  Clinical helpers        (make_tertile, run_clinical_screen)
#   10. Context-mapping callout (mapnote) — requirement #2
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(relaimpo)
  library(knitr)
})


# =====================================================================
# 1. SMALL UTILITIES
# =====================================================================

# Vectorised p-value -> significance label.
pvalue_to_label <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "\u2014",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# log10 pro-inflammatory / tolerogenic ratio from log10 subclass values.
# Pro-inflammatory = IgG1 + IgG3 ; tolerogenic = IgG2 + IgG4.
ratio_col <- function(g1, g2, g3, g4) {
  log10((10^g1 + 10^g3) / (10^g2 + 10^g4))
}

# Recode the raw P0x / M0x visit codes to the recipient-friendly names.
# Driven entirely by VISIT_RECODE in config/endpoints.R, so the mapping
# lives in exactly one place. Unknown codes are left unchanged (and
# warned about) rather than silently turned into NA.
recode_visits <- function(x) {
  x   <- as.character(x)
  new <- VISIT_RECODE[x]
  unmapped <- is.na(new) & !is.na(x)
  if (any(unmapped)) {
    warning("recode_visits(): unmapped visit codes left unchanged: ",
            paste(unique(x[unmapped]), collapse = ", "))
    new[unmapped] <- x[unmapped]
  }
  unname(new)
}


# =====================================================================
# 2. DATA LOADING + RECODING
# =====================================================================
# Loads the antibody assay object, applies the visit recode ONCE, builds
# the arm_name + feature columns, and returns a tidy long data frame.
# Every downstream part starts from this, so the recode happens in a
# single, auditable place.
load_serology_data <- function(path = here::here("./data/c_set.RData"),
                               object_name = "ebaa_extra") {
  load(path)
  raw <- get(object_name)

  raw$arm_name   <- factor(raw$maternal_arm)         # TdaP vs TT (maternal vaccine)
  raw$feature    <- factor(paste(raw$antigen, raw$analyte, sep = "_"))
  raw$visit_name <- recode_visits(raw$visit_name)    # <-- the global recode
  raw$visit_name <- factor(raw$visit_name,
                           levels = intersect(VISIT_ORDER, unique(raw$visit_name)))

  # MANDATORY casing guard (single call site for the whole project): the
  # chain/concurrent paths match feature names case-sensitively against the
  # post-load mixed-case features built above, so a config/data casing drift
  # would silently drop features. verify_predictor_casing() (config/endpoints.R)
  # errors loudly on any mismatch. Guarded so the loader still works if the
  # config has not been sourced (e.g. an isolated unit test).
  if (exists("verify_predictor_casing", mode = "function"))
    verify_predictor_casing(raw)

  raw
}

# Load the clinical-assessment object (used by the clinical part).
load_clinical_data <- function(path = here::here("./data/clin_assess.RData"),
                               object_name = "clin_assess") {
  load(path)
  get(object_name)
}

# ---------------------------------------------------------------------
# Descriptive study-population table (manuscript Table 1), overall and by
# maternal arm. Dependency-light (base R + dplyr only). Continuous variables
# are summarised as median [IQR]; categoricals as n (%). The clinical object
# carries no arm column, so arm membership is joined from arm_df (one row per
# subject: subject_accession + arm_name) derived from the serology data.
# Returns a character data frame ready for knitr::kable().
# ---------------------------------------------------------------------
demographics_by_arm <- function(
    clin, arm_df,
    cont = c("maternal_age", "maternal_bmi", "maternal_Hb",
             "gestational_age_vaccination", "gestational_age_birth",
             "vaccine_birth_interval_days", "birth_weight"),
    cat  = c("parity", "infant_sex", "delivery_mode"),
    arm_levels = c("TdaP", "TT")) {

  d <- dplyr::left_join(clin, arm_df, by = "subject_accession")
  d$arm_name <- as.character(d$arm_name)
  cont <- intersect(cont, names(d))
  cat  <- intersect(cat,  names(d))

  idx <- c(list(Overall = seq_len(nrow(d))),
           stats::setNames(lapply(arm_levels,
                                  function(a) which(d$arm_name == a)), arm_levels))

  fmt_cont <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (all(is.na(x))) return("\u2014")
    sprintf("%.1f [%.1f, %.1f]",
            stats::median(x, na.rm = TRUE),
            stats::quantile(x, .25, na.rm = TRUE, names = FALSE),
            stats::quantile(x, .75, na.rm = TRUE, names = FALSE))
  }
  fmt_cat <- function(x, lvl) {
    x <- as.character(x); tot <- sum(!is.na(x)); n <- sum(x == lvl, na.rm = TRUE)
    if (tot == 0) "\u2014" else sprintf("%d (%.0f%%)", n, 100 * n / tot)
  }

  rows <- list(c("n", vapply(idx, function(i) as.character(length(i)), "")))
  for (v in cont)
    rows[[length(rows) + 1]] <- c(paste0(v, ", median [IQR]"),
                                  vapply(idx, function(i) fmt_cont(d[[v]][i]), ""))
  for (v in cat) {
    lv <- sort(unique(as.character(d[[v]][!is.na(d[[v]])])))
    rows[[length(rows) + 1]] <- c(paste0(v, ", n (%)"),
                                  rep("", length(idx)))
    for (l in lv)
      rows[[length(rows) + 1]] <- c(paste0("\u2003", l),
        vapply(idx, function(i) fmt_cat(d[[v]][i], l), ""))
  }
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("Characteristic", names(idx))
  rownames(out) <- NULL
  out
}


# =====================================================================
# 3. PARTICIPANT / N TRACKER  (requirement #4)
# =====================================================================
# A lightweight ledger that records, for every model fitted, how many
# subjects entered, how many were actually used (complete cases), how
# many were dropped, and which predictor(s) caused the drop. Rendered as
# running tables so a reader can follow exactly why different sections of
# the analysis contain different numbers of subjects.

new_ntracker <- function() {
  e <- new.env(parent = emptyenv())
  e$rows <- list()
  e
}

# Internal: append one accounting row. `vars` is the full set of columns
# the model needs (outcome + predictors); the drop is computed by
# listwise complete-case deletion over exactly those columns.
nt_log <- function(tracker, part, section, stratum, label, vars, data) {
  if (is.null(tracker)) return(invisible())
  vars_present <- intersect(vars, names(data))
  n_subset     <- nrow(data)

  if (length(vars_present) > 0) {
    cc         <- stats::complete.cases(data[, vars_present, drop = FALSE])
    n_used     <- sum(cc)
    na_counts  <- vapply(vars_present, function(v) sum(is.na(data[[v]])), integer(1))
    drivers    <- sort(na_counts[na_counts > 0], decreasing = TRUE)
    drivers_str <- if (length(drivers))
      paste(sprintf("%s (%d NA)", names(drivers), drivers), collapse = "; ") else "\u2014"
    missing_cols <- setdiff(vars, vars_present)
    if (length(missing_cols))
      drivers_str <- paste0(drivers_str,
                            "; [cols absent: ", paste(missing_cols, collapse = ", "), "]")
  } else {
    n_used <- 0L
    drivers_str <- paste0("no required columns present: ",
                          paste(vars, collapse = ", "))
  }

  tracker$rows[[length(tracker$rows) + 1]] <- data.frame(
    part = part, section = section, stratum = stratum, model = label,
    n_entered = n_subset, n_used = n_used, n_dropped = n_subset - n_used,
    dropped_due_to = drivers_str, stringsAsFactors = FALSE
  )
  invisible()
}

# Collect the whole ledger as a data frame.
nt_collect <- function(tracker) {
  if (is.null(tracker) || length(tracker$rows) == 0) return(data.frame())
  dplyr::bind_rows(tracker$rows)
}

# Render the ledger (optionally filtered to one part) as a kable.
nt_render <- function(tracker, part = NULL, caption = NULL) {
  df <- nt_collect(tracker)
  if (nrow(df) == 0) { cat("*No models logged.*\n\n"); return(invisible()) }
  if (!is.null(part)) df <- df[df$part == part, , drop = FALSE]
  if (is.null(caption))
    caption <- "Participant accounting: subjects entered, used, and dropped per model."
  print(knitr::kable(
    df[, c("section", "stratum", "model",
           "n_entered", "n_used", "n_dropped", "dropped_due_to")],
    col.names = c("Section", "Stratum", "Model",
                  "N entered", "N used", "N dropped", "Dropped due to (NA driver)"),
    caption = caption, row.names = FALSE,
    align = c("l", "l", "l", "r", "r", "r", "l")
  ))
  cat("\n")
}

# A compact per-stratum headcount: how many subjects have a non-missing
# outcome in each stratum, before any predictor-driven attrition.
nt_outcome_counts <- function(data, strata, outcome_col, caption = NULL) {
  rows <- lapply(names(strata), function(nm) {
    d <- strata[[nm]]$data
    data.frame(
      Stratum        = strata[[nm]]$label,
      N_subjects     = nrow(d),
      N_with_outcome = sum(!is.na(d[[outcome_col]])),
      stringsAsFactors = FALSE
    )
  })
  df <- dplyr::bind_rows(rows)
  if (is.null(caption))
    caption <- paste0("Subjects per stratum and number with a non-missing outcome (",
                      outcome_col, ").")
  print(knitr::kable(df, caption = caption, row.names = FALSE,
                     align = c("l", "r", "r")))
  cat("\n")
}


# =====================================================================
# 4. CORE MODEL RUNNER (tracker-aware)
# =====================================================================
# Fits an lm, computes relative importance (lmg), and — when a tracker is
# supplied — logs the participant accounting for this model. Returns a
# structured result consumed by show_model() / r2_row().
run_lm <- function(outcome_col, predictors, data, label = "",
                   show_relimp_calc = TRUE,
                   tracker = NULL, part = NA, section = NA, stratum = NA) {

  # ---- participant accounting (before fitting) ----
  nt_log(tracker, part, section, stratum, label,
         vars = c(outcome_col, predictors), data = data)

  f <- reformulate(predictors, response = outcome_col)
  m <- tryCatch(lm(f, data = data, na.action = na.omit),
                error = function(e) NULL)
  if (is.null(m)) {
    return(invisible(list(model = NULL, relimp = NULL, label = label,
                          R2 = NA, adj_R2 = NA, n_obs = NA, p_model = NA,
                          tidy_coef = NULL, ri_table = NULL)))
  }
  s <- summary(m)

  beta_lab <- "\u03b2"   # unicode in a STRING is fine; inside backticks it is not
  tidy_coef <- broom::tidy(m) %>%
    mutate(sig = pvalue_to_label(p.value),
           across(c(estimate, std.error, statistic), ~round(., 4)),
           p.value = signif(p.value, 3)) %>%
    rename(Predictor = term, SE = std.error, t = statistic, p = p.value, ` ` = sig) %>%
    rename(!!beta_lab := estimate)

  ri <- NULL
  if (isTRUE(show_relimp_calc)) {
    ri <- tryCatch(calc.relimp(m, type = "lmg"), error = function(e) NULL)
  }
  ri_table <- if (!is.null(ri)) {
    data.frame(Predictor = names(ri@lmg),
               lmg_pct   = round(ri@lmg * 100, 2),
               stringsAsFactors = FALSE) %>% arrange(desc(lmg_pct))
  } else NULL

  p_model <- tryCatch(
    pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE),
    error = function(e) NA)

  invisible(list(
    model = m, relimp = ri, label = label,
    R2 = round(s$r.squared, 4), adj_R2 = round(s$adj.r.squared, 4),
    n_obs = nobs(m), p_model = p_model,
    tidy_coef = tidy_coef, ri_table = ri_table
  ))
}


# =====================================================================
# 5. DISPLAY HELPERS
# =====================================================================

show_model <- function(res, show_relimp = TRUE) {
  if (is.null(res) || is.null(res$tidy_coef)) {
    cat("*Model could not be estimated (insufficient data or singular fit).*\n\n")
    return(invisible(NULL))
  }
  n_pred <- nrow(res$tidy_coef) - 1L            # non-intercept coefficient terms
  if (!is.na(res$n_obs) && res$n_obs <= n_pred + 1L) {
    cat(sprintf("\n*%s \u2014 model saturated: n = %s \u2264 predictors + 1 (%d); too small to estimate (see participant accounting).*\n\n",
                res$label, res$n_obs, n_pred + 1L))
    return(invisible(NULL))
  }
  cat(sprintf(
    "\n**n = %s &nbsp;|&nbsp; R\u00b2 = %.3f &nbsp;|&nbsp; adj-R\u00b2 = %.3f &nbsp;|&nbsp; p~model~ = %.3g**\n\n",
    res$n_obs, res$R2, res$adj_R2, res$p_model))
  print(knitr::kable(res$tidy_coef,
                     caption = paste0("Coefficients \u2014 ", res$label),
                     align = c("l", rep("r", ncol(res$tidy_coef) - 1))))
  cat("\n")
  if (isTRUE(show_relimp) && !is.null(res$ri_table)) {
    print(knitr::kable(res$ri_table %>% rename(`lmg (%)` = lmg_pct),
                       caption = paste0("Relative importance (lmg) \u2014 ", res$label),
                       align = c("l", "r")))
    cat("\n")
  }
}

# One-row R2 summary; `stratum` carries whichever stratum level is active.
# `saturated` flags models with n_obs <= n_predictors + 1 (same rule as show_model),
# so the R2 summary tables can suppress meaningless perfect-fit R2 values.
r2_row <- function(res, stratum, block, model_name) {
  n_pred <- if (is.null(res$tidy_coef)) NA_integer_ else nrow(res$tidy_coef) - 1L
  saturated <- is.null(res$tidy_coef) ||
    (!is.na(res$n_obs) && !is.na(n_pred) && res$n_obs <= n_pred + 1L)
  data.frame(stratum = stratum, block = block, model = model_name,
             n_obs = res$n_obs, R2 = res$R2, adj_R2 = res$adj_R2,
             saturated = saturated, stringsAsFactors = FALSE)
}

# Render every stored model id for one stratum (used by display sections).
show_all_models <- function(br, ids, show_relimp = TRUE) {
  for (id in ids) if (!is.null(br[[id]])) show_model(br[[id]], show_relimp = show_relimp)
}


# =====================================================================
# 6. DERIVED CHAIN VARIABLES (visit-name driven)
# =====================================================================
# Builds, for the maternal->cord->infant chain, the per-visit ratios and
# the Block 2/3/4/5/6 delta variables. Column names are generated from
# the (recoded) visit names in cfg$chain_visits, so nothing here is
# hard-coded to P09/M00/M02 and the same code works for every endpoint.
#
# Naming scheme produced:
#   <ANT>_ratio_<visit>                         e.g. PT_ratio_MatBirth
#   delta_<ANT>_ratio_<from>_<to>               e.g. delta_PT_ratio_PregEarly_MatBirth
#   delta_<feature>_<from>_<to>  (qty)          e.g. delta_PT_IgG1_CordBlood_InfMon2
add_chain_derived <- function(dw, cfg) {
  visits    <- cfg$chain_visits           # c(PregEarly, MatBirth, CordBlood, InfMon2)
  antigens  <- cfg$antigens
  feats_qty <- c(cfg$antigen_totals, cfg$all_subclasses)

  # ---- ratios at each visit ----
  for (v in visits) for (ant in antigens) {
    need <- paste0(ant, c("_IgG1", "_IgG2", "_IgG3", "_IgG4"), "_", v)
    if (all(need %in% names(dw))) {
      dw[[paste0(ant, "_ratio_", v)]] <-
        ratio_col(dw[[need[1]]], dw[[need[2]]], dw[[need[3]]], dw[[need[4]]])
    }
  }

  # ---- consecutive stage pairs ----
  pairs <- list(c(visits[1], visits[2]),   # Block 2 : ratio switch (pregnancy)
                c(visits[2], visits[3]),    # Block 3 qty + Block 4 ratio (transfer)
                c(visits[3], visits[4]))    # Block 5 qty + Block 6 ratio (decay/change)

  # ratio deltas for all three pairs
  for (p in pairs) for (ant in antigens) {
    ra <- paste0(ant, "_ratio_", p[1]); rb <- paste0(ant, "_ratio_", p[2])
    if (all(c(ra, rb) %in% names(dw)))
      dw[[paste0("delta_", ant, "_ratio_", p[1], "_", p[2])]] <- dw[[rb]] - dw[[ra]]
  }

  # quantity deltas only for the transfer and decay pairs
  for (p in pairs[2:3]) for (f in feats_qty) {
    fa <- paste0(f, "_", p[1]); fb <- paste0(f, "_", p[2])
    if (all(c(fa, fb) %in% names(dw)))
      dw[[paste0("delta_", f, "_", p[1], "_", p[2])]] <- dw[[fb]] - dw[[fa]]
  }

  dw
}

# Convenience: the visit names in their chain roles, for prose/labels.
chain_visit_roles <- function(cfg) {
  v <- cfg$chain_visits
  list(early = v[1], matbirth = v[2], cord = v[3], infant = v[4])
}


# =====================================================================
# 7. CHAIN BLOCK RUNNER (all six blocks, one stratum)
# =====================================================================
# Fits every Block 1-6 model on ONE data subset (one stratum) and returns
# the structured storage the display sections read from. Visit names come
# from cfg, so this is endpoint- and recode-agnostic. Every model is
# logged to the supplied N tracker under the given part/stratum.
run_chain_blocks <- function(d, cfg, tracker = NULL,
                             part = "Maternal chain", stratum = "") {
  oc  <- cfg$outcome_col
  rl  <- chain_visit_roles(cfg)
  ve  <- rl$early; vm <- rl$matbirth; vc <- rl$cord; vi <- rl$infant
  ant <- cfg$antigen_totals
  pt  <- PT_SUBCLASSES; fha <- FHA_SUBCLASSES; prn <- PRN_SUBCLASSES
  alls <- cfg$all_subclasses

  br <- list(); r2 <- list(); ratio_models <- list()

  L <- function(label, preds, sec, relimp = TRUE)
    run_lm(oc, preds, d, label, show_relimp_calc = relimp,
           tracker = tracker, part = part, section = sec, stratum = stratum)

  # ---- BLOCK 1 : maternal levels at early-pregnancy and at birth ----
  br[["1a"]] <- L(paste0("Antigen totals @ ", ve),       paste0(ant, "_", ve),  "Block 1")
  br[["1b"]] <- L(paste0("Antigen totals @ ", vm),       paste0(ant, "_", vm),  "Block 1")
  br[["1c"]] <- L(paste0("Antigen totals ", ve, " + ", vm),
                  c(paste0(ant, "_", ve), paste0(ant, "_", vm)), "Block 1", relimp = FALSE)
  br[["1d"]] <- L(paste0("PT subclasses @ ", ve),        paste0(pt, "_", ve),   "Block 1")
  br[["1e"]] <- L(paste0("PT subclasses @ ", vm),        paste0(pt, "_", vm),   "Block 1")
  br[["1f"]] <- L(paste0("FHA subclasses @ ", ve),       paste0(fha, "_", ve),  "Block 1")
  br[["1g"]] <- L(paste0("FHA subclasses @ ", vm),       paste0(fha, "_", vm),  "Block 1")
  br[["1h"]] <- L(paste0("PRN subclasses @ ", ve),       paste0(prn, "_", ve),  "Block 1")
  br[["1i"]] <- L(paste0("PRN subclasses @ ", vm),       paste0(prn, "_", vm),  "Block 1")
  br[["1j"]] <- L(paste0("Full 12-subclass @ ", ve),     paste0(alls, "_", ve), "Block 1")
  br[["1k"]] <- L(paste0("Full 12-subclass @ ", vm),     paste0(alls, "_", vm), "Block 1")
  br[["1l"]] <- L(paste0("Full 24-predictor ", ve, " + ", vm),
                  c(paste0(alls, "_", ve), paste0(alls, "_", vm)), "Block 1", relimp = FALSE)
  ids1 <- paste0("1", letters[1:12])
  mods1 <- c("antigen_early", "antigen_birth", "antigen_early+birth",
             "PT_sub_early", "PT_sub_birth", "FHA_sub_early", "FHA_sub_birth",
             "PRN_sub_early", "PRN_sub_birth", "full12_early", "full12_birth",
             "full24_early+birth")
  for (i in seq_along(ids1))
    r2[[ids1[i]]] <- r2_row(br[[ids1[i]]], stratum, "B1_quantity", mods1[i])

  # ---- BLOCK 2 : ratio switch (early pregnancy -> birth) ----
  b2 <- paste0("delta_", cfg$antigens, "_ratio_", ve, "_", vm)
  br[["2"]] <- L(paste0("\u0394 ratio ", ve, "\u2192", vm), b2, "Block 2")
  r2[["2"]] <- r2_row(br[["2"]], stratum, "B2_ratio_switch", "3-ratio")
  ratio_models[["B2"]] <- br[["2"]]$model

  # ---- BLOCK 3 : transfer quantity (birth -> cord) ----
  b3_ant  <- paste0("delta_", ant, "_", vm, "_", vc)
  b3_pt   <- paste0("delta_", pt,  "_", vm, "_", vc)
  b3_fha  <- paste0("delta_", fha, "_", vm, "_", vc)
  b3_prn  <- paste0("delta_", prn, "_", vm, "_", vc)
  br[["3a"]] <- L(paste0("\u0394 antigen totals ", vm, "\u2192", vc), b3_ant, "Block 3")
  br[["3b"]] <- L(paste0("\u0394 PT subclasses ",  vm, "\u2192", vc), b3_pt,  "Block 3")
  br[["3c"]] <- L(paste0("\u0394 FHA subclasses ", vm, "\u2192", vc), b3_fha, "Block 3")
  br[["3d"]] <- L(paste0("\u0394 PRN subclasses ", vm, "\u2192", vc), b3_prn, "Block 3")
  br[["3e"]] <- L(paste0("\u0394 Full 12-subclass ", vm, "\u2192", vc), c(b3_pt, b3_fha, b3_prn), "Block 3")
  ids3 <- paste0("3", letters[1:5])
  mods3 <- c("antigen_delta", "PT_sub_delta", "FHA_sub_delta", "PRN_sub_delta", "full12_delta")
  for (i in seq_along(ids3))
    r2[[ids3[i]]] <- r2_row(br[[ids3[i]]], stratum, "B3_transfer_qty", mods3[i])

  # ---- BLOCK 4 : transfer selectivity (birth -> cord) ----
  b4 <- paste0("delta_", cfg$antigens, "_ratio_", vm, "_", vc)
  br[["4"]] <- L(paste0("\u0394 ratio ", vm, "\u2192", vc), b4, "Block 4")
  r2[["4"]] <- r2_row(br[["4"]], stratum, "B4_transfer_ratio", "3-ratio")
  ratio_models[["B4"]] <- br[["4"]]$model

  # ---- BLOCK 5 : neonatal change quantity (cord -> infant) ----
  b5_ant  <- paste0("delta_", ant, "_", vc, "_", vi)
  b5_pt   <- paste0("delta_", pt,  "_", vc, "_", vi)
  b5_fha  <- paste0("delta_", fha, "_", vc, "_", vi)
  b5_prn  <- paste0("delta_", prn, "_", vc, "_", vi)
  br[["5a"]] <- L(paste0("\u0394 antigen totals ", vc, "\u2192", vi), b5_ant, "Block 5")
  br[["5b"]] <- L(paste0("\u0394 PT subclasses ",  vc, "\u2192", vi), b5_pt,  "Block 5")
  br[["5c"]] <- L(paste0("\u0394 FHA subclasses ", vc, "\u2192", vi), b5_fha, "Block 5")
  br[["5d"]] <- L(paste0("\u0394 PRN subclasses ", vc, "\u2192", vi), b5_prn, "Block 5")
  br[["5e"]] <- L(paste0("\u0394 Full 12-subclass ", vc, "\u2192", vi), c(b5_pt, b5_fha, b5_prn), "Block 5")
  ids5 <- paste0("5", letters[1:5])
  for (i in seq_along(ids5))
    r2[[ids5[i]]] <- r2_row(br[[ids5[i]]], stratum, "B5_change_qty", mods3[i])

  # ---- BLOCK 6 : neonatal change selectivity (cord -> infant) ----
  b6 <- paste0("delta_", cfg$antigens, "_ratio_", vc, "_", vi)
  br[["6"]] <- L(paste0("\u0394 ratio ", vc, "\u2192", vi), b6, "Block 6")
  r2[["6"]] <- r2_row(br[["6"]], stratum, "B6_change_ratio", "3-ratio")
  ratio_models[["B6"]] <- br[["6"]]$model

  list(
    br           = br,
    r2           = dplyr::bind_rows(r2),
    ratio_models = ratio_models,
    relimp = list(B1_birth = br[["1k"]]$relimp,   # full-12 maternal at birth
                  B3_transfer = br[["3e"]]$relimp,
                  B5_change   = br[["5e"]]$relimp)
  )
}


# =====================================================================
# 8. STRATA BUILDER  (requirement #5: maternal-arm primary)
# =====================================================================
# Returns an ordered, named list of strata with the maternal arm as the
# PRIMARY axis. For each maternal arm we produce:
#   (a) all infant arms pooled within that maternal arm, then
#   (b) each infant level separately (the maternal x infant cross).
# Each element carries $data, a human $label, and $marm / $iarm tags.
#
# Setting cross = FALSE yields only the pooled maternal margins.
build_maternal_strata <- function(d, maternal_arms = MATERNAL_ARMS,
                                  infant_arms = INFANT_ARMS, cross = TRUE) {
  out <- list()
  for (m in maternal_arms) {
    out[[m]] <- list(
      data  = dplyr::filter(d, arm_name == m),
      label = paste0(m, " \u2014 all infant arms"),
      marm  = m, iarm = "all"
    )
    if (isTRUE(cross)) for (i in infant_arms) {
      key <- paste0(m, "_", i)
      out[[key]] <- list(
        data  = dplyr::filter(d, arm_name == m, infant_arm == i),
        label = paste0(m, " / ", i, " infants"),
        marm  = m, iarm = i
      )
    }
  }
  out
}

# Optional infant-only margins (aP, wP pooled over maternal arm), retained
# for backward comparison with the earlier infant-primary write-up.
build_infant_strata <- function(d, infant_arms = INFANT_ARMS) {
  out <- list()
  for (i in infant_arms) {
    out[[i]] <- list(
      data  = dplyr::filter(d, infant_arm == i),
      label = paste0(i, " \u2014 all maternal arms"),
      marm  = "all", iarm = i
    )
  }
  out
}


# =====================================================================
# 9. CLINICAL HELPERS (used by the clinical-covariate part)
# =====================================================================

# Ordered tertile factor T1(low)/T2(mid)/T3(high).
make_tertile <- function(x) {
  qs <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE)
  cut(x, breaks = c(-Inf, qs[1], qs[2], Inf),
      labels = c("T1 (low)", "T2 (mid)", "T3 (high)"),
      include.lowest = TRUE, ordered_result = FALSE)
}

# Merge the outcome with clinical covariates and prepare covariate forms:
# 2-level factors for categoricals, tertile columns (<var>_T) for continuous.
# Tertiles are computed on the FULL merged cohort so cutpoints are identical
# across strata. Returns the merged data frame (one row per subject).
prepare_clinical <- function(data_outcome, clin, cov_named = COVARIATES_NAMED,
                             cat_vars = COVARIATES_CATEGORICAL) {
  d <- merge(data_outcome, clin, by = "subject_accession", all.x = TRUE)

  if ("infant_sex" %in% names(d))    d$infant_sex <- factor(d$infant_sex)

  if ("delivery_mode" %in% names(d)) {
    lev <- levels(factor(d$delivery_mode))
    d$delivery_mode <- forcats::fct_collapse(
      factor(d$delivery_mode),
      vaginal  = grep("(?i)vag|normal|spontan|SVD", lev, value = TRUE, perl = TRUE),
      cesarean = grep("(?i)ces|csec|section|C/S",   lev, value = TRUE, perl = TRUE),
      other    = grep("(?i)vag|normal|spontan|SVD|ces|csec|section|C/S",
                      lev, value = TRUE, perl = TRUE, invert = TRUE)
    ) %>% forcats::fct_drop()
    if ("vaginal" %in% levels(d$delivery_mode))
      d$delivery_mode <- relevel(d$delivery_mode, ref = "vaginal")
  }

  if ("parity" %in% names(d))
    d$parity <- factor(ifelse(d$parity == 0, "nulliparous", "parous"),
                       levels = c("nulliparous", "parous"))

  cont_vars <- setdiff(unname(unlist(cov_named)), cat_vars)
  for (v in cont_vars) if (v %in% names(d)) d[[paste0(v, "_T")]] <- make_tertile(d[[v]])
  d
}

# Fit one clinical univariate lm; tracker-aware; returns a structured result.
run_lm_clin <- function(data, outcome_col, predictor_col, label = "",
                        tracker = NULL, part = NA, section = NA, stratum = NA) {
  nt_log(tracker, part, section, stratum, label,
         vars = c(outcome_col, predictor_col), data = data)
  f <- reformulate(predictor_col, response = outcome_col)
  m <- tryCatch(lm(f, data = data, na.action = na.omit), error = function(e) NULL)
  if (is.null(m))
    return(list(label = label, predictor = predictor_col, R2 = NA, adj_R2 = NA,
                n_obs = NA, p_model = NA, tidy_coef = NULL, model = NULL))
  s <- summary(m)
  p_mod <- tryCatch(pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                       lower.tail = FALSE), error = function(e) NA)
  beta_lab <- "\u03b2"
  tidy_coef <- broom::tidy(m) %>%
    mutate(sig = pvalue_to_label(p.value),
           across(c(estimate, std.error, statistic), ~round(., 4)),
           p.value = signif(p.value, 3)) %>%
    rename(Predictor = term, SE = std.error, t = statistic, p = p.value, ` ` = sig) %>%
    rename(!!beta_lab := estimate)
  list(label = label, predictor = predictor_col,
       R2 = round(s$r.squared, 4), adj_R2 = round(s$adj.r.squared, 4),
       n_obs = nobs(m), p_model = p_mod, tidy_coef = tidy_coef, model = m)
}

# Fit a clinical multivariate lm from >= 2 selected predictors.
run_multi_clin <- function(data, outcome_col, predictors, label = "Multivariate",
                           tracker = NULL, part = NA, section = NA, stratum = NA) {
  if (length(predictors) < 2) return(NULL)
  nt_log(tracker, part, section, stratum, label,
         vars = c(outcome_col, predictors), data = data)
  f <- reformulate(predictors, response = outcome_col)
  m <- tryCatch(lm(f, data = data, na.action = na.omit), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  s <- summary(m)
  p_mod <- tryCatch(pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                       lower.tail = FALSE), error = function(e) NA)
  beta_lab <- "\u03b2"
  tidy_coef <- broom::tidy(m) %>%
    mutate(sig = pvalue_to_label(p.value),
           across(c(estimate, std.error, statistic), ~round(., 4)),
           p.value = signif(p.value, 3)) %>%
    rename(Predictor = term, SE = std.error, t = statistic, p = p.value, ` ` = sig) %>%
    rename(!!beta_lab := estimate)
  list(label = label, R2 = round(s$r.squared, 4), adj_R2 = round(s$adj.r.squared, 4),
       n_obs = nobs(m), p_model = p_mod, tidy_coef = tidy_coef, model = m)
}

# Display a clinical model result (header + coefficient table).
show_lm_clin <- function(res) {
  if (is.null(res) || is.null(res$tidy_coef)) {
    cat("*Model could not be estimated.*\n\n"); return(invisible(NULL))
  }
  cat(sprintf("\n**n = %s &nbsp;|&nbsp; R\u00b2 = %.3f &nbsp;|&nbsp; adj-R\u00b2 = %.3f &nbsp;|&nbsp; p~model~ = %.3g**\n\n",
              res$n_obs, res$R2, res$adj_R2, res$p_model))
  print(knitr::kable(res$tidy_coef, caption = paste0("Coefficients \u2014 ", res$label),
                     align = c("l", rep("r", ncol(res$tidy_coef) - 1))))
  cat("\n")
}

# Univariate screen + form selection + multivariate, for ONE stratum.
#   - categorical covariates: single factor model
#   - continuous covariates: continuous AND tertile forms; the significant
#     form with higher adj-R2 is selected (or the only significant form)
#   - if >= 2 covariates are selected, a joint multivariate model is fit
run_clinical_screen <- function(d, outcome_col, cov_list = COVARIATES_NAMED,
                                cat_vars = COVARIATES_CATEGORICAL, alpha = 0.05,
                                tracker = NULL, part = "Clinical", stratum = "") {
  univ_rows <- list(); univ_res <- list()

  for (lbl in names(cov_list)) {
    v <- cov_list[[lbl]]
    if (v %in% cat_vars) {
      r <- run_lm_clin(d, outcome_col, v, label = lbl,
                       tracker = tracker, part = part, section = "Univariate", stratum = stratum)
      univ_rows[[paste0(v, "_cat")]] <- data.frame(
        Label = lbl, Variable = v, Form = "categorical",
        n = ifelse(is.na(r$n_obs), NA_integer_, r$n_obs),
        R2 = r$R2, adj_R2 = r$adj_R2, p_model = r$p_model,
        sig = ifelse(is.na(r$p_model), "\u2014", pvalue_to_label(r$p_model)),
        stringsAsFactors = FALSE)
      univ_res[[paste0(v, "_cat")]] <- r
    } else {
      vT <- paste0(v, "_T")
      r_cont <- run_lm_clin(d, outcome_col, v,  label = paste0(lbl, " (continuous)"),
                            tracker = tracker, part = part, section = "Univariate", stratum = stratum)
      r_tert <- run_lm_clin(d, outcome_col, vT, label = paste0(lbl, " (tertile)"),
                            tracker = tracker, part = part, section = "Univariate", stratum = stratum)
      for (tag in c("cont", "tert")) {
        r   <- if (tag == "cont") r_cont else r_tert
        frm <- if (tag == "cont") "continuous" else "tertile"
        var <- if (tag == "cont") v else vT
        univ_rows[[paste0(v, "_", tag)]] <- data.frame(
          Label = lbl, Variable = var, Form = frm,
          n = ifelse(is.na(r$n_obs), NA_integer_, r$n_obs),
          R2 = r$R2, adj_R2 = r$adj_R2, p_model = r$p_model,
          sig = ifelse(is.na(r$p_model), "\u2014", pvalue_to_label(r$p_model)),
          stringsAsFactors = FALSE)
      }
      univ_res[[paste0(v, "_cont")]] <- r_cont
      univ_res[[paste0(v, "_tert")]] <- r_tert
    }
  }
  univ_table <- dplyr::bind_rows(univ_rows)

  selected <- list(); form_notes <- list()
  for (lbl in names(cov_list)) {
    v <- cov_list[[lbl]]
    if (v %in% cat_vars) {
      r <- univ_res[[paste0(v, "_cat")]]
      if (!is.na(r$p_model) && r$p_model < alpha) {
        selected[[lbl]] <- v; form_notes[[lbl]] <- "categorical \u2014 significant"
      }
    } else {
      rc <- univ_res[[paste0(v, "_cont")]]; rt <- univ_res[[paste0(v, "_tert")]]
      sc <- !is.na(rc$p_model) && rc$p_model < alpha
      st <- !is.na(rt$p_model) && rt$p_model < alpha
      if (sc && st) {
        if (!is.na(rc$adj_R2) && !is.na(rt$adj_R2) && rc$adj_R2 >= rt$adj_R2) {
          selected[[lbl]] <- v
          form_notes[[lbl]] <- sprintf("continuous (adj-R\u00b2=%.3f) > tertile (adj-R\u00b2=%.3f) \u2014 both sig",
                                       rc$adj_R2, rt$adj_R2)
        } else {
          selected[[lbl]] <- paste0(v, "_T")
          form_notes[[lbl]] <- sprintf("tertile (adj-R\u00b2=%.3f) > continuous (adj-R\u00b2=%.3f) \u2014 both sig",
                                       rt$adj_R2, rc$adj_R2)
        }
      } else if (sc) { selected[[lbl]] <- v;            form_notes[[lbl]] <- "continuous \u2014 only form significant"
      } else if (st) { selected[[lbl]] <- paste0(v, "_T"); form_notes[[lbl]] <- "tertile \u2014 only form significant" }
    }
  }

  sel_predictors <- unname(unlist(selected))
  multi_result <- if (length(sel_predictors) > 1)
    run_multi_clin(d, outcome_col, sel_predictors,
                   label = "Multivariate (all significant predictors)",
                   tracker = tracker, part = part, section = "Multivariate", stratum = stratum) else NULL

  list(univ_table = univ_table, selected = selected, form_notes = form_notes,
       multi_result = multi_result, univ_res = univ_res)
}


# =====================================================================
# 10. CONTEXT-MAPPING CALLOUT  (requirement #2)
# =====================================================================
# Emits a visually distinct block tying an analysis section to the
# specific interpretation-summary section / synthesis question it feeds.
# The knitted output (HTML + Word) is what gets handed to the downstream
# interpretation engine, so these callouts make the analysis<->narrative
# mapping explicit and machine-followable.
mapnote <- function(..., feeds = NULL) {
  cat("\n> **\U0001F9ED Maps to interpretation:** ", paste0(...), "\n")
  if (!is.null(feeds)) {
    cat(">\n")
    for (f in feeds) cat(">", "-", f, "\n")
  }
  cat("\n")
}


# =====================================================================
# 11. CONCURRENT (within-visit) ANALYSIS HELPERS
# =====================================================================
# Minimum complete-case n required to estimate a GGM network / run NCT.
NETWORK_MIN <- 30L

# Fits the concurrent regression layer on ONE stratum's wide data (one row
# per subject x visit, predictors and outcome measured at the same visit;
# models pool across the concurrent visits). Uses the generic run_lm, so
# everything is logged to the N tracker. Returns stored models + relimp.
run_concurrent_models <- function(d, cfg, tracker = NULL,
                                  part = "Concurrent", stratum = "") {
  oc  <- cfg$outcome_feature                      # bare feature name, e.g. WHOLE_SBA
  ant <- cfg$antigen_totals
  pt  <- PT_SUBCLASSES; fha <- FHA_SUBCLASSES; prn <- PRN_SUBCLASSES
  alls <- cfg$all_subclasses
  L <- function(label, preds, sec, relimp = TRUE)
    run_lm(oc, preds, d, label, show_relimp_calc = relimp,
           tracker = tracker, part = part, section = sec, stratum = stratum)

  br <- list()
  br[["antigen"]] <- L("Antigen totals \u2192 outcome",        ant, "Antigen")
  br[["PT"]]      <- L("PT subclasses \u2192 outcome",         pt,  "Subclass")
  br[["FHA"]]     <- L("FHA subclasses \u2192 outcome",        fha, "Subclass")
  br[["PRN"]]     <- L("PRN subclasses \u2192 outcome",        prn, "Subclass")
  br[["full12"]]  <- L("Full 12-subclass \u2192 outcome",      alls, "Full12")
  br[["ratio"]]   <- L("Pro/tolerogenic ratios \u2192 outcome",
                       c("PT_ratio", "FHA_ratio", "PRN_ratio"), "Ratio", relimp = FALSE)

  list(br = br,
       relimp_full = br[["full12"]]$relimp,
       relimp_antigen = br[["antigen"]]$relimp,
       ratio_model = br[["ratio"]]$model)
}

# Add the three pro-/tolerogenic ratio columns to concurrent wide data
# (predictors carry no visit suffix here, since each row is one visit).
add_concurrent_ratios <- function(dw) {
  for (ant in ANTIGENS) {
    need <- paste0(ant, c("_IgG1", "_IgG2", "_IgG3", "_IgG4"))
    if (all(need %in% names(dw)))
      dw[[paste0(ant, "_ratio")]] <-
        ratio_col(dw[[need[1]]], dw[[need[2]]], dw[[need[3]]], dw[[need[4]]])
  }
  dw
}


# =====================================================================
# 12. PART 4 HELPERS — MATERNAL-VACCINATION-EFFECT "BUBBLES"
# =====================================================================

# Manual N-tracker row (for accounting units that aren't lm complete-cases,
# e.g. unique-subject counts per bubble grid). Same columns as nt_log.
nt_log_manual <- function(tracker, part, section, stratum, label,
                          n_entered, n_used, note = "\u2014") {
  if (is.null(tracker)) return(invisible())
  tracker$rows[[length(tracker$rows) + 1]] <- data.frame(
    part = part, section = section, stratum = stratum, model = label,
    n_entered = n_entered, n_used = n_used, n_dropped = n_entered - n_used,
    dropped_due_to = note, stringsAsFactors = FALSE)
  invisible()
}

# Compute the TdaP-vs-TT effect grid (antigen x analyte) for one subset, and
# build the bubble plot. effect_size = median(TT) - median(TdaP) on the active
# (untransformed or standardized) scale; Wilcoxon rank-sum per cell; Bonferroni
# across the grid. Returns list(plot, stats). Cells with < 2 obs per arm get
# no test (p = 1, "No effect").
bubble_arm_effect <- function(data,
                              analyte_order = BUBBLE_ANALYTES,
                              antigen_order = BUBBLE_ANTIGENS,
                              effect_colors = BUBBLE_EFFECT_COLORS,
                              standstat = "untransformed", title = NULL) {
  ua <- unique(as.character(data$arm_name))
  if (!all(c("TT", "TdaP") %in% ua))
    stop("bubble_arm_effect(): arm_name must contain both 'TT' and 'TdaP'")

  stats_df <- data %>%
    dplyr::group_by(antigen, analyte) %>%
    dplyr::summarize(mTT   = list(log_assay_value[arm_name == "TT"]),
                     mTdaP = list(log_assay_value[arm_name == "TdaP"]),
                     .groups = "drop") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      n_TT   = length(stats::na.omit(unlist(mTT))),
      n_TdaP = length(stats::na.omit(unlist(mTdaP))),
      p_value = if (n_TT >= 2 && n_TdaP >= 2)
        stats::wilcox.test(unlist(mTT), unlist(mTdaP), exact = FALSE)$p.value else NA_real_,
      median_TT   = stats::median(unlist(mTT),   na.rm = TRUE),
      median_TdaP = stats::median(unlist(mTdaP), na.rm = TRUE),
      effect_size = median_TT - median_TdaP,
      direction = dplyr::case_when(
        is.na(effect_size) ~ "No effect",
        effect_size > 0    ~ "Lower TdaP",
        effect_size < 0    ~ "Higher TdaP",
        TRUE               ~ "No effect")) %>%
    dplyr::ungroup() %>%
    dplyr::select(-mTT, -mTdaP)

  stats_df$p_value     <- ifelse(stats_df$direction == "No effect" | is.na(stats_df$p_value),
                                 1, stats_df$p_value)
  stats_df$p_adj       <- stats::p.adjust(stats_df$p_value, method = "bonferroni", n = nrow(stats_df))
  stats_df$significant <- stats_df$p_adj < 0.05
  stats_df$border_size <- ifelse(stats_df$significant, 1.5, 0.5)
  stats_df$analyte <- forcats::fct_relevel(factor(stats_df$analyte),
                        intersect(analyte_order, unique(as.character(stats_df$analyte))))
  stats_df$antigen <- forcats::fct_relevel(factor(stats_df$antigen),
                        intersect(antigen_order, unique(as.character(stats_df$antigen))))
  stats_df$standstat <- standstat

  p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = analyte, y = antigen)) +
    ggplot2::geom_point(ggplot2::aes(size = abs(effect_size), fill = direction,
                                     stroke = border_size),
                        color = "black", shape = 21, show.legend = TRUE) +
    ggplot2::scale_fill_manual(name = "Effect direction", values = effect_colors) +
    ggplot2::scale_size_continuous(name = "Effect size (|median diff|)", range = c(3, 10)) +
    ggplot2::scale_x_discrete(name = "Analyte") +
    ggplot2::scale_y_discrete(name = "Antigen") +
    ggplot2::labs(title = title,
      caption = "Wilcoxon rank-sum (TT vs TdaP) on the active scale; Bonferroni-corrected.\nBold border = adjusted p < 0.05. Effect = median(TT) \u2212 median(TdaP).") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "right") +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1), color = "none",
                    size = ggplot2::guide_legend(order = 2))
  list(plot = p, stats = stats_df)
}

# Enhanced per-grid table from a bubble_arm_effect() stats data frame.
bubble_stats_kable <- function(stats, caption) {
  tab <- stats %>%
    dplyr::arrange(antigen, analyte) %>%
    dplyr::transmute(antigen, analyte, n_TT, n_TdaP,
                     med_TT = round(median_TT, 3), med_TdaP = round(median_TdaP, 3),
                     delta = round(effect_size, 3), direction,
                     p = signif(p_value, 3), p_bonf = signif(p_adj, 3),
                     sig = pvalue_to_label(p_adj))
  knitr::kable(tab,
    col.names = c("Antigen","Analyte","n TT","n TdaP","Median TT","Median TdaP",
                  "\u0394 (TT\u2212TdaP)","Direction","p","p (Bonf)",""),
    caption = caption,
    align = c("l","l","r","r","r","r","r","l","r","r","c"))
}


# =====================================================================
# 13. PAIRED WITHIN-SUBJECT CHANGE BETWEEN TWO VISITS
# =====================================================================
# Builds the per-subject change (log10) for every feature measured at BOTH
# `from_visit` and `to_visit`, classifies each subject as increase/decrease/
# no-change, and supports the cord->infant directionality analysis.

# Per-subject paired change. Returns one row per subject x feature that has
# both visits, with delta = value(to) - value(from) and a direction label.
paired_change <- function(data_raw, from_visit, to_visit, eps = 0) {
  base <- data_raw %>%
    dplyr::filter(visit_name %in% c(from_visit, to_visit)) %>%
    dplyr::select(subject_accession, arm_name, infant_arm, antigen, analyte,
                  feature, visit_name, log_assay_value) %>%
    dplyr::mutate(visit_name = as.character(visit_name)) %>%
    distinct(subject_accession, feature, visit_name, .keep_all = TRUE)

  wide <- base %>%
    tidyr::pivot_wider(names_from = visit_name, values_from = log_assay_value)
  if (!all(c(from_visit, to_visit) %in% names(wide)))
    stop("paired_change(): one or both visits absent after pivot.")
  wide$v_from <- wide[[from_visit]]
  wide$v_to   <- wide[[to_visit]]

  wide %>%
    dplyr::filter(!is.na(v_from), !is.na(v_to)) %>%
    dplyr::mutate(delta = v_to - v_from,
                  dir = dplyr::case_when(delta >  eps ~ "increase",
                                         delta < -eps ~ "decrease",
                                         TRUE         ~ "no change"))
}

# Coverage / attrition per feature: how many have the from visit, the to
# visit, and BOTH (the paired n that the directional test actually uses).
paired_coverage <- function(data_raw, from_visit, to_visit) {
  d <- data_raw %>%
    dplyr::filter(visit_name %in% c(from_visit, to_visit)) %>%
    dplyr::select(subject_accession, feature, visit_name) %>%
    distinct()
  from_n <- d %>% dplyr::filter(visit_name == from_visit) %>% dplyr::count(feature, name = "n_from")
  to_n   <- d %>% dplyr::filter(visit_name == to_visit)   %>% dplyr::count(feature, name = "n_to")
  both <- d %>% tidyr::pivot_wider(names_from = visit_name, values_from = visit_name,
                                   values_fn = length, values_fill = 0)
  if (!all(c(from_visit, to_visit) %in% names(both)))
    return(dplyr::full_join(from_n, to_n, by = "feature"))
  both$n_paired <- as.integer(both[[from_visit]] > 0 & both[[to_visit]] > 0)
  paired_n <- both %>% dplyr::group_by(feature) %>%
    dplyr::summarize(n_paired = sum(n_paired), .groups = "drop")
  from_n %>% dplyr::full_join(to_n, by = "feature") %>%
    dplyr::full_join(paired_n, by = "feature") %>%
    dplyr::mutate(across(c(n_from, n_to, n_paired), ~tidyr::replace_na(., 0L)))
}

# Per-feature (optionally per group level) directional summary: counts and
# proportion increasing/decreasing, median/mean delta, signed-rank and sign
# tests against no change.
summarise_paired <- function(paired, by = NULL, min_n = 1) {
  grp <- c("feature", by)
  out <- paired %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
    dplyr::summarize(
      n_paired = dplyr::n(),
      n_inc = sum(dir == "increase"), n_dec = sum(dir == "decrease"),
      n_eq  = sum(dir == "no change"),
      pct_inc = round(100 * n_inc / n_paired, 1),
      pct_dec = round(100 * n_dec / n_paired, 1),
      median_delta = round(stats::median(delta), 3),
      mean_delta   = round(mean(delta), 3),
      p_signrank = tryCatch(signif(stats::wilcox.test(delta, mu = 0, exact = FALSE)$p.value, 3),
                            error = function(e) NA_real_),
      p_sign = tryCatch(signif(stats::binom.test(n_inc, n_inc + n_dec, 0.5)$p.value, 3),
                        error = function(e) NA_real_),
      .groups = "drop") %>%
    dplyr::filter(n_paired >= min_n) %>%
    dplyr::mutate(direction = dplyr::case_when(
      pct_inc >= 60 ~ "mostly increase",
      pct_inc <= 40 ~ "mostly decrease",
      TRUE          ~ "mixed"))
  out
}
