# =====================================================================
# R/subclass_balance_helpers.R
# ---------------------------------------------------------------------
# Functions supporting the subclass-polarisation re-analysis, refactored
# so that COMPLEMENT (antigen-specific ADCD) proceeds through exactly the
# same stage machinery as the subclass balances. Every stage function is
# generic over a `tag`:
#     tag = "bal"   -> the ILR pro/tolerogenic balance  (<ANT>_bal_<visit>)
#     tag = "ADCD"  -> antibody-dependent complement deposition
#                       (<ANT>_ADCD_<visit>, as pulled from the long data)
#
# The maternal->cord->infant chain is decomposed into three STAGES:
#     production = LEVEL at delivery (MatBirth)
#     transfer   = delta  MatBirth -> CordBlood
#     decay      = delta  CordBlood -> InfMon2
# and a complement WALK traces  balance -> ADCD -> SBA  at each node
# (maternal, cord, baseline).
#
# SCALE NOTE (arbitrary-unit assays -> ordinal):
#   A balance is a linear combination of log-subclasses, so a per-subclass
#   unit change only ADDS a constant (shifts a regression intercept, not the
#   slope); stage DELTAS cancel the constant entirely. ADCD is likewise an
#   arbitrary-unit measure, so a-/b-path coefficients and individual balances
#   are NOT concentration-interpretable, but SLOPES, CIs, p values, delta-R2
#   and PROPORTION MEDIATED (a dimensionless ratio of effects on the SBA
#   scale) are all valid on ordinal data.
#
# Dependencies: base R + stats. dplyr/tidyr only in build_wide_frame().
# Sits beside config/endpoints.R and R/serology_helpers.R; table fitting is
# self-contained (does not require run_lm()) for portability + testability.
# =====================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# ---- default stage definitions (transitions in the chain) ------------------
default_balance_stages <- function() list(
  production = c("PregEarly", "MatBirth"),   # (kept for reference; production uses the level)
  transfer   = c("MatBirth",  "CordBlood"),  # placental transfer selectivity
  decay      = c("CordBlood", "InfMon2")      # post-natal decay to 8-week baseline
)

# default nodes for the complement walk: label -> visit
default_walk_nodes <- function() c(Maternal = "MatBirth",
                                    Cord = "CordBlood",
                                    Baseline = "InfMon2")

# ---------------------------------------------------------------------
# 1. ILR balance columns  <ANT>_bal_<visit>
# ---------------------------------------------------------------------
add_subclass_balances <- function(dw, antigens, visits,
                                  pro = c("IgG1", "IgG3"),
                                  tol = c("IgG2", "IgG4")) {
  for (a in antigens) for (v in visits) {
    pcols <- paste0(a, "_", pro, "_", v)
    tcols <- paste0(a, "_", tol, "_", v)
    if (!all(c(pcols, tcols) %in% names(dw))) next
    pro_mean <- rowMeans(dw[, pcols, drop = FALSE])   # na.rm=FALSE: honest complete-case
    tol_mean <- rowMeans(dw[, tcols, drop = FALSE])
    dw[[paste0(a, "_bal_", v)]] <- pro_mean - tol_mean
  }
  dw
}

# ---------------------------------------------------------------------
# 2. Stage delta columns  <ANT>_<tag>delta_<stage>   (generic over tag)
# ---------------------------------------------------------------------
add_stage_deltas <- function(dw, antigens, tag = "bal",
                             stages = default_balance_stages()) {
  for (a in antigens) for (st in names(stages)) {
    fc <- paste0(a, "_", tag, "_", stages[[st]][1])
    tc <- paste0(a, "_", tag, "_", stages[[st]][2])
    if (!all(c(fc, tc) %in% names(dw))) next
    dw[[paste0(a, "_", tag, "delta_", st)]] <- dw[[tc]] - dw[[fc]]
  }
  dw
}
# back-compatible wrapper
add_balance_deltas <- function(dw, antigens, stages = default_balance_stages())
  add_stage_deltas(dw, antigens, tag = "bal", stages = stages)

# per-antigen stage quantities for a given measure `tag`:
#   production = level at production_visit; transfer/decay = stage deltas
stage_quantity_cols <- function(antigen, tag = "bal", production_visit = "MatBirth") {
  c(production = paste0(antigen, "_", tag, "_", production_visit),
    transfer   = paste0(antigen, "_", tag, "delta_transfer"),
    decay      = paste0(antigen, "_", tag, "delta_decay"))
}

# ---------------------------------------------------------------------
# 3. Self-contained tidy linear fit (beta, 95% CI, p) for chosen terms
# ---------------------------------------------------------------------
tidy_lm <- function(outcome, predictors, data, keep = NULL) {
  keep <- keep %||% predictors
  f  <- stats::reformulate(predictors, outcome)
  cc <- stats::complete.cases(data[, all.vars(f), drop = FALSE])
  d  <- data[cc, , drop = FALSE]
  if (nrow(d) < (length(predictors) + 2L))
    return(list(model = NULL, n = nrow(d), r2 = NA, adj = NA, tbl = NULL))
  m  <- stats::lm(f, d)
  s  <- summary(m)
  co <- s$coefficients
  ci <- suppressWarnings(stats::confint(m))
  rows <- lapply(keep, function(k) {
    if (!k %in% rownames(co)) return(NULL)
    data.frame(term = k, beta = unname(co[k, 1]),
               lo = unname(ci[k, 1]), hi = unname(ci[k, 2]),
               p = unname(co[k, 4]), stringsAsFactors = FALSE)
  })
  list(model = m, n = nrow(d), r2 = s$r.squared, adj = s$adj.r.squared,
       tbl = do.call(rbind, rows))
}

# ---------------------------------------------------------------------
# 4. WITHIN-STAGE univariate: SBA ~ <stage quantity> + arm  (generic tag)
# ---------------------------------------------------------------------
fit_stage_univariate <- function(dw, antigens, outcome, tag = "bal",
                                 arm = "arm_name", production_visit = "MatBirth",
                                 nt = NULL, part = "Stage analysis") {
  stage_order <- c("production", "transfer", "decay")
  out <- list()
  for (a in antigens) {
    qcols <- stage_quantity_cols(a, tag, production_visit)
    for (st in stage_order) {
      qc <- qcols[[st]]
      if (!qc %in% names(dw)) next
      fit <- tidy_lm(outcome, c(qc, arm), dw, keep = qc)
      if (!is.null(nt) && exists("nt_log_manual"))
        nt_log_manual(nt, part = part, section = paste0("Within-stage (", tag, ")"),
                      stratum = st, label = paste0(a, " ", st, " ", tag, " + arm"),
                      n_entered = nrow(dw), n_used = fit$n)
      if (is.null(fit$tbl)) next
      out[[length(out) + 1L]] <- data.frame(
        antigen = a, stage = factor(st, levels = stage_order),
        n = fit$n, beta = fit$tbl$beta, lo = fit$tbl$lo, hi = fit$tbl$hi,
        p = fit$tbl$p, stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, out)
  if (!is.null(res)) res <- res[order(res$stage, res$antigen), ]
  res
}

# ---------------------------------------------------------------------
# 5. SEQUENTIAL (hierarchical) models on one common sample (generic tag)
#    m1: SBA ~ production + arm
#    m2: + transfer
#    m3: + decay
# ---------------------------------------------------------------------
fit_sequential_stages <- function(dw, antigen, outcome, tag = "bal",
                                  arm = "arm_name", production_visit = "MatBirth",
                                  nt = NULL, part = "Stage analysis") {
  q <- stage_quantity_cols(antigen, tag, production_visit)
  needed <- c(outcome, q[["production"]], q[["transfer"]], q[["decay"]])
  if (!all(needed %in% names(dw)))
    return(list(added = NULL, steps = NULL, note = "missing measure columns"))
  cc <- stats::complete.cases(dw[, c(needed, arm), drop = FALSE])
  d  <- dw[cc, , drop = FALSE]
  if (!is.null(nt) && exists("nt_log_manual"))
    nt_log_manual(nt, part = part, section = paste0("Sequential (", tag, ")"),
                  stratum = antigen,
                  label = paste0(antigen, " ", tag, ": production/transfer/decay + arm"),
                  n_entered = nrow(dw), n_used = nrow(d))
  specs <- list(
    production = c(q[["production"]], arm),
    transfer   = c(q[["production"]], q[["transfer"]], arm),
    decay      = c(q[["production"]], q[["transfer"]], q[["decay"]], arm))
  added_term <- c(production = q[["production"]], transfer = q[["transfer"]],
                  decay = q[["decay"]])
  steps <- list(); added <- list(); prev_r2 <- 0
  for (st in names(specs)) {
    fit <- tidy_lm(outcome, specs[[st]], d, keep = added_term[[st]])
    steps[[st]] <- data.frame(step = st, n = fit$n, R2 = fit$r2,
                              adj_R2 = fit$adj, delta_R2 = fit$r2 - prev_r2,
                              stringsAsFactors = FALSE)
    prev_r2 <- fit$r2
    if (!is.null(fit$tbl))
      added[[st]] <- data.frame(step = st, term = added_term[[st]],
                                beta = fit$tbl$beta, lo = fit$tbl$lo,
                                hi = fit$tbl$hi, p = fit$tbl$p, stringsAsFactors = FALSE)
  }
  full <- tidy_lm(outcome, specs[["decay"]], d, keep = unname(added_term))
  list(added = do.call(rbind, added), steps = do.call(rbind, steps),
       full_model = full$tbl, n = nrow(d))
}

# ---------------------------------------------------------------------
# 6. Feature discovery + generic wide builder + ADCD tag derivation
# ---------------------------------------------------------------------
discover_features <- function(data_long, pattern, antigen = NULL,
                              feature_col = "feature") {
  feats <- as.character(unique(data_long[[feature_col]]))
  hit <- grepl(pattern, feats, ignore.case = TRUE)
  if (!is.null(antigen)) hit <- hit & grepl(paste0("^", antigen, "_"), feats)
  feats[hit]
}

feature_visits <- function(data_long, feature, feature_col = "feature",
                           visit_col = "visit_name", value_col = "log_assay_value") {
  sub <- data_long[data_long[[feature_col]] == feature &
                     !is.na(data_long[[value_col]]), ]
  sort(unique(as.character(sub[[visit_col]])))
}

# The analyte token used by ADCD features (e.g. "ADCD" from "PRN_ADCD").
# Returns list(tag, feats) or NULL if none found for the given antigens.
derive_adcd_tag <- function(data_long, antigens, pattern = "ADCD") {
  feats <- unique(unlist(lapply(antigens, function(a)
    discover_features(data_long, pattern, antigen = a))))
  if (!length(feats)) return(NULL)
  tags <- unique(sub("^[A-Za-z0-9]+_", "", feats))   # strip "<ANT>_"
  list(tag = tags[1], feats = feats, all_tags = tags)
}

build_wide_frame <- function(data_long, feats, visits, arms,
                             id_cols = c("subject_accession", "arm_name", "infant_arm")) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE),
            requireNamespace("tidyr", quietly = TRUE))
  data_long |>
    dplyr::filter(.data$feature %in% feats,
                  .data$visit_name %in% visits,
                  .data$arm_name %in% arms) |>
    dplyr::select(dplyr::all_of(c(id_cols, "feature", "visit_name", "log_assay_value"))) |>
    dplyr::distinct(.data$subject_accession, .data$feature, .data$visit_name, .keep_all = TRUE) |>
    tidyr::pivot_wider(id_cols = dplyr::all_of(id_cols),
                       names_from = c("feature", "visit_name"),
                       values_from = "log_assay_value", names_sep = "_") |>
    as.data.frame()
}

# ---------------------------------------------------------------------
# 7. Mediation X -> M -> Y  (product-of-coefficients + bootstrap; mediation pkg if present)
# ---------------------------------------------------------------------
run_balance_mediation <- function(df, X, M, Y, cov = "arm_name",
                                  R = 2000, seed = 20240101) {
  keep <- c(X, M, Y, if (cov %in% names(df)) cov)
  df <- df[stats::complete.cases(df[, keep, drop = FALSE]), , drop = FALSE]
  if (nrow(df) < 20L) return(list(ok = FALSE, n = nrow(df), note = "too few complete cases"))
  rhs_m <- if (cov %in% names(df)) c(X, cov) else X
  rhs_y <- if (cov %in% names(df)) c(X, M, cov) else c(X, M)
  f_m <- stats::reformulate(rhs_m, M); f_y <- stats::reformulate(rhs_y, Y)
  est <- function(d) {
    a  <- stats::coef(stats::lm(f_m, d))[[X]]
    cy <- stats::coef(stats::lm(f_y, d)); b <- cy[[M]]; cp <- cy[[X]]
    c(indirect = a * b, direct = cp, total = a * b + cp)
  }
  pt <- est(df)
  set.seed(seed)
  bs <- replicate(R, {
    d <- df[sample(nrow(df), replace = TRUE), , drop = FALSE]
    tryCatch(est(d), error = function(e) c(indirect = NA, direct = NA, total = NA))
  })
  ci <- apply(bs, 1, stats::quantile, c(0.025, 0.975), na.rm = TRUE)
  pkg <- NULL
  if (requireNamespace("mediation", quietly = TRUE)) {
    pkg <- tryCatch({
      mm <- stats::lm(f_m, df); my <- stats::lm(f_y, df)
      set.seed(seed)
      med <- mediation::mediate(mm, my, treat = X, mediator = M, sims = 1000)
      list(acme = med$d0, acme_ci = med$d0.ci, ade = med$z0, ade_ci = med$z0.ci,
           total = med$tau.coef, total_ci = med$tau.ci, prop = med$n0, prop_ci = med$n0.ci)
    }, error = function(e) NULL)
  }
  list(ok = TRUE, n = nrow(df), point = pt, ci = ci,
       prop_mediated = unname(pt["indirect"] / pt["total"]), pkg = pkg, X = X, M = M, Y = Y)
}

# ---------------------------------------------------------------------
# 8. COMPLEMENT WALK: trace balance -> ADCD -> SBA at each chain node
# ---------------------------------------------------------------------
# For each node (Maternal/Cord/Baseline) at its visit v:
#   a-path : ADCD_v  ~ balance_v + arm     (does polarisation drive complement here?)
#   b-path : SBA_out ~ ADCD_v   + arm      (does complement here predict 8-week SBA?)
#   mediation: balance_v -> ADCD_v -> SBA_out
# Returns one tidy row per antigen x node.
stage_complement_walk <- function(dw, antigens, outcome, nodes = default_walk_nodes(),
                                  bal_tag = "bal", adcd_tag = "ADCD",
                                  arm = "arm_name", R = 2000, nt = NULL,
                                  part = "Stage analysis") {
  rows <- list()
  for (a in antigens) for (nlab in names(nodes)) {
    v  <- nodes[[nlab]]
    xc <- paste0(a, "_", bal_tag, "_", v)
    mc <- paste0(a, "_", adcd_tag, "_", v)
    if (!all(c(xc, mc, outcome) %in% names(dw))) next
    apath <- tidy_lm(mc, c(xc, arm), dw, keep = xc)          # ADCD ~ balance + arm
    bpath <- tidy_lm(outcome, c(mc, arm), dw, keep = mc)     # SBA  ~ ADCD    + arm
    med   <- run_balance_mediation(dw, X = xc, M = mc, Y = outcome, cov = arm, R = R)
    if (!is.null(nt) && exists("nt_log_manual"))
      nt_log_manual(nt, part = part, section = "Complement walk",
                    stratum = paste0(a, " @ ", nlab),
                    label = paste0(a, " balance->ADCD->SBA at ", nlab),
                    n_entered = nrow(dw), n_used = if (isTRUE(med$ok)) med$n else 0L)
    rows[[length(rows) + 1L]] <- data.frame(
      antigen = a, node = factor(nlab, levels = names(nodes)), visit = v,
      n = if (isTRUE(med$ok)) med$n else NA_integer_,
      a_beta = if (!is.null(apath$tbl)) apath$tbl$beta else NA,
      a_p    = if (!is.null(apath$tbl)) apath$tbl$p    else NA,
      b_beta = if (!is.null(bpath$tbl)) bpath$tbl$beta else NA,
      b_p    = if (!is.null(bpath$tbl)) bpath$tbl$p    else NA,
      indirect    = if (isTRUE(med$ok)) unname(med$point["indirect"]) else NA,
      indirect_lo = if (isTRUE(med$ok)) unname(med$ci[1, "indirect"]) else NA,
      indirect_hi = if (isTRUE(med$ok)) unname(med$ci[2, "indirect"]) else NA,
      prop_med    = if (isTRUE(med$ok)) med$prop_mediated else NA,
      stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, rows)
  if (!is.null(res)) res <- res[order(res$antigen, res$node), ]
  res
}

# ---------------------------------------------------------------------
# 9. balance vs original linear-sum ratio (robustness)
# ---------------------------------------------------------------------
compare_balance_vs_ratio <- function(dw, antigen, visit, outcome, arm = "arm_name") {
  bal <- paste0(antigen, "_bal_", visit)
  rat <- paste0(antigen, "_ratio_", visit)   # from add_chain_derived()
  out <- list()
  for (nm in c(balance = bal, ratio = rat)) {
    if (!nm %in% names(dw)) next
    fit <- tidy_lm(outcome, c(nm, arm), dw, keep = nm)
    if (is.null(fit$tbl)) next
    out[[length(out) + 1L]] <- data.frame(
      representation = if (nm == bal) "ILR balance (geometric)" else "linear-sum ratio",
      column = nm, n = fit$n, beta = fit$tbl$beta, lo = fit$tbl$lo, hi = fit$tbl$hi,
      p = fit$tbl$p, R2 = fit$r2, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}
