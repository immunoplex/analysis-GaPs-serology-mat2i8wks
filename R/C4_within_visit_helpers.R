# =============================================================================
# R/within_visit_helpers.R
#
# Reusable helpers for:
#   07_within_visit_correlations_and_responder_subgroups.Rmd
#
# Design notes
# ------------
# These functions were factored out of the original monolithic analysis chunk
# so that they (a) take all inputs as explicit arguments instead of reaching
# into the knit environment for globals, and (b) can be unit-tested or reused
# by other modules. None of them attach packages; they reference dplyr/tidyr/
# ggplot2/cluster/ggcorrplot via `::` or rely on the caller having attached
# dplyr (for the `%>%` pipe and bare-name NSE, matching the rest of the
# codebase). Source this file AFTER dplyr is attached, but note that sourcing
# only *defines* these functions -- the pipe only needs to resolve when they
# are *called*.
#
# Conventions inherited from the existing pipeline:
#   * long-format columns: subject_accession, arm_name, feature,
#     visit_name, log_assay_value
#   * maternal arm reference level is "TT"
#   * titres are skewed -> Spearman / rank-based transforms by default
# =============================================================================


# ---- small numeric / lookup utilities --------------------------------------

#' Blom rank-based normal score, computed within the supplied vector.
#' NAs are preserved. Used to standardise a feature within a single timepoint
#' so the cross-plate offset cancels before cross-timepoint combination.
blom <- function(x) {
  r <- rank(x, na.last = "keep")
  n <- sum(!is.na(x))
  stats::qnorm((r - 3 / 8) / (n + 1 / 4))
}

#' Case-insensitive resolution of one requested feature against the panel.
#' Returns NA_character_ if absent. Mirrors resolve_feat() in 06_prediction_models.
resolve_feat <- function(w, avail_feats) {
  hit <- avail_feats[toupper(avail_feats) == toupper(w)]
  if (length(hit)) hit[[1]] else NA_character_
}

#' Vectorised resolve_feat that drops anything not in the panel.
resolve_feats <- function(ws, avail_feats) {
  out <- vapply(ws, resolve_feat, character(1), avail_feats = avail_feats)
  as.character(stats::na.omit(out))
}

#' Pick a single visit display label by matching any of several regex patterns
#' against the recoded visit levels. Warns (and returns NA) if nothing matches,
#' so an empty downstream cell is attributable to design, not a silent bug.
pick_visit <- function(visit_levels, patterns, label = "(unnamed)") {
  hit <- visit_levels[Reduce(`|`, lapply(patterns, function(p)
    grepl(p, visit_levels, ignore.case = TRUE)))]
  if (!length(hit)) {
    warning(sprintf("No visit matched '%s' (have: %s)",
                    label, paste(visit_levels, collapse = ", ")))
    return(NA_character_)
  }
  hit[[1]]
}


# ---- feature x visit availability -------------------------------------------

# Whole-bacterial functional / binding endpoints are NOT assayed on maternal
# serum. Per config/endpoints.R (each endpoint's `concurrent_visits` are exactly
# the visits at which its outcome_feature is measured):
#   WHOLE_PTNA   : InfMon2, InfMon5, InfMon9
#   WHOLE_SBA    : InfMon2, InfMon5, InfMon9
#   WHOLE_WT_IgG : CordBlood, InfMon2, InfMon5, InfMon9   (the one exception that
#                  is also measured on cord blood)
# i.e. NONE of the WHOLE_* assays exist at PregEarly or MatBirth.
# Antigen-specific totals are available at every visit (PT/FHA/PRN/DT/TT/FIM
# _IgG totals exist maternal+infant). NB: FIM *subclasses* are infant-only, but
# the FIM_IgG total used in this module exists everywhere, so it is not gated.
WHOLE_ASSAY_VISITS_DEFAULT <- list(
  WHOLE_PTNA   = c("InfMon2", "InfMon5", "InfMon9"),
  WHOLE_SBA    = c("InfMon2", "InfMon5", "InfMon9"),
  WHOLE_WT_IgG = c("CordBlood", "InfMon2", "InfMon5", "InfMon9")
)

#' Derive the whole-assay availability map from a sourced ENDPOINTS list so the
#' truth lives in config, not here: each endpoint's concurrent_visits are the
#' visits where its outcome_feature is measured. Falls back to the default above.
whole_assay_visit_map <- function(endpoints = NULL) {
  if (is.null(endpoints) || !length(endpoints)) return(WHOLE_ASSAY_VISITS_DEFAULT)
  out <- tryCatch(
    stats::setNames(
      lapply(endpoints, function(e) as.character(e$concurrent_visits)),
      vapply(endpoints, function(e) as.character(e$outcome_feature), character(1))),
    error = function(e) NULL)
  if (is.null(out) || !length(out) || any(!nzchar(names(out))))
    WHOLE_ASSAY_VISITS_DEFAULT else out
}

#' Is `feat` measured at `visit`? WHOLE_* features are gated by the map
#' (case-insensitive on the key); every other feature is assumed available.
feature_measured_at <- function(feat, visit, whole_map = WHOLE_ASSAY_VISITS_DEFAULT) {
  key <- names(whole_map)[toupper(names(whole_map)) == toupper(feat)]
  if (!length(key)) return(TRUE)            # not a gated whole-bacterial assay
  visit %in% whole_map[[key[1]]]
}

#' Validate that a recoded visit label is actually present in the data; warn
#' (and still return the label) if not, so a typo in VISIT_RECODE is visible.
require_visit <- function(label, visit_levels) {
  if (!label %in% visit_levels)
    warning(sprintf("Visit label '%s' not found in data (have: %s)",
                    label, paste(visit_levels, collapse = ", ")))
  label
}


# ---- per-visit wide frame ----------------------------------------------------

#' One subject x feature wide frame for a single visit.
#' Keeps arm_name (releveled to `ref_arm` when present) so Q4 can split on arm.
make_wide <- function(data_raw, visit, assays, ref_arm = "TT") {
  w <- dplyr::filter(data_raw, visit_name == visit, feature %in% assays)
  w <- dplyr::select(w, subject_accession, arm_name, feature, log_assay_value)
  w <- tidyr::pivot_wider(
    w,
    id_cols     = c(subject_accession, arm_name),
    names_from  = feature,
    values_from = log_assay_value,
    values_fn   = function(x) mean(x, na.rm = TRUE))
  arm <- factor(w$arm_name)
  if (ref_arm %in% levels(arm)) arm <- stats::relevel(arm, ref = ref_arm)
  w$arm_name <- arm
  w
}


# ---- within-subject correlation engine --------------------------------------

#' Pairwise correlation, across subjects, of the assays measured on the SAME
#' blood draw. Pairwise-complete; reports r, raw p, BH-adjusted p (across the
#' unique off-diagonal pairs) and the n behind every cell. Default Spearman.
within_visit_cor <- function(w, assays, method = "spearman", min_n = 4L) {
  cols <- intersect(assays, names(w))
  M <- as.matrix(w[, cols, drop = FALSE])
  storage.mode(M) <- "double"

  p <- length(cols)
  R <- matrix(NA_real_, p, p, dimnames = list(cols, cols))
  P <- matrix(NA_real_, p, p, dimnames = list(cols, cols))
  N <- matrix(0L,       p, p, dimnames = list(cols, cols))

  for (i in seq_len(p)) for (j in seq_len(p)) {
    ok <- is.finite(M[, i]) & is.finite(M[, j])
    N[i, j] <- sum(ok)
    if (sum(ok) >= min_n && i != j) {
      ct <- tryCatch(
        stats::cor.test(M[ok, i], M[ok, j], method = method, exact = FALSE),
        error = function(e) NULL)
      if (!is.null(ct)) {
        R[i, j] <- unname(ct$estimate)
        P[i, j] <- ct$p.value
      }
    } else if (i == j) {
      R[i, j] <- 1
    }
  }

  # Benjamini-Hochberg across the unique off-diagonal tests
  Padj <- P
  ut <- upper.tri(P)
  if (any(ut)) {
    Padj[ut] <- stats::p.adjust(P[ut], method = "BH")
    Padj[lower.tri(Padj)] <- t(Padj)[lower.tri(Padj)]
  }

  list(r = R, p = P, p_adj = Padj, n = N, cols = cols, method = method)
}

#' Flatten a within_visit_cor() result to one row per unique assay pair.
cor_long <- function(cc, timepoint) {
  cols <- cc$cols
  out <- list(); k <- 0
  for (i in seq_along(cols)) for (j in seq_along(cols)) if (j > i) {
    k <- k + 1
    out[[k]] <- data.frame(
      timepoint = timepoint, assay_x = cols[i], assay_y = cols[j],
      r = cc$r[i, j], n = cc$n[i, j], p = cc$p[i, j], p_BH = cc$p_adj[i, j])
  }
  if (!length(out)) {
    return(data.frame(
      timepoint = character(), assay_x = character(), assay_y = character(),
      r = numeric(), n = integer(), p = numeric(), p_BH = numeric()))
  }
  dplyr::bind_rows(out)
}

#' Heatmap of a correlation matrix. Uses ggcorrplot when available, otherwise a
#' base-ggplot tile fallback. Prints a note (and returns invisibly) when < 2
#' assays are available so the calling chunk degrades gracefully.
plot_cor <- function(cc, title,
                     has_ggcorr = requireNamespace("ggcorrplot", quietly = TRUE)) {
  if (length(cc$cols) < 2) {
    cat("*Fewer than two assays available - no matrix.*\n\n")
    return(invisible())
  }
  if (has_ggcorr) {
    print(ggcorrplot::ggcorrplot(
      cc$r, lab = TRUE, type = "lower", title = title,
      colors = c("#d01c8b", "white", "#2166ac")))
  } else {
    df <- expand.grid(x = cc$cols, y = cc$cols)
    df$r <- as.vector(cc$r)
    print(
      ggplot2::ggplot(df, ggplot2::aes(x, y, fill = r)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(
          ggplot2::aes(label = ifelse(is.na(r), "", sprintf("%.2f", r))),
          size = 3) +
        ggplot2::scale_fill_gradient2(
          low = "#d01c8b", mid = "white", high = "#2166ac",
          midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
        ggplot2::labs(title = title, x = NULL, y = NULL) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
  }
  invisible()
}


# ---- responder-panel assembly & subgrouping ---------------------------------

#' Long, within-timepoint-standardised responder table across an arbitrary set
#' of visits. `binding_feats` are included at every visit in `visits`;
#' `whole_feats` are included only at visits where the availability map
#' (feature_measured_at) says they are assayed -- so a maternal-visit set
#' silently carries no whole-bacterial columns, and an infant-visit set carries
#' them at every infant visit. One Blom z per subject x feature x visit,
#' standardised WITHIN each (feature, visit) so the cross-plate offset cancels.
build_responder_panel <- function(data_raw, visits, binding_feats, whole_feats,
                                   avail_feats,
                                   whole_map = WHOLE_ASSAY_VISITS_DEFAULT) {
  out <- list(); k <- 0
  for (v in visits) {
    if (is.na(v)) next
    fset <- binding_feats
    for (wf in whole_feats)
      if (feature_measured_at(wf, v, whole_map)) fset <- c(fset, wf)
    fset <- intersect(unique(fset), avail_feats)
    if (!length(fset)) next
    sub <- data_raw %>%
      dplyr::filter(visit_name == v, feature %in% fset) %>%
      dplyr::group_by(subject_accession, arm_name, feature) %>%
      dplyr::summarize(val = mean(log_assay_value, na.rm = TRUE),
                       .groups = "drop")
    if (!nrow(sub)) next
    sub <- sub %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(z = blom(val)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(visit = v, key = paste(feature, visit, sep = "@"))
    k <- k + 1; out[[k]] <- sub
  }
  if (!length(out)) return(dplyr::tibble())
  dplyr::bind_rows(out)
}

#' Subject x visit responder score: mean within-timepoint z across the panel
#' features available at that visit. Used to trace a subject (or subgroup)
#' trajectory across an ordered set of visits (e.g. InfMon2 -> InfMon5 -> InfMon9).
responder_score_by_visit <- function(resp_long) {
  if (!nrow(resp_long)) return(resp_long)
  resp_long %>%
    dplyr::group_by(subject_accession, arm_name, visit) %>%
    dplyr::summarize(score_visit = mean(z, na.rm = TRUE),
                     n_feat = dplyr::n(), .groups = "drop")
}

#' Add the transparent tertile bands and the breadth-gated "universal" flags to
#' a subject-level table from summarise_responder_subjects(). Shared by stages.
add_tertile_bands <- function(subject_df, breadth_thresh = 0.6) {
  q <- stats::quantile(subject_df$universal_z, c(1 / 3, 2 / 3), na.rm = TRUE)
  subject_df %>%
    dplyr::mutate(
      score_band = dplyr::case_when(
        universal_z >= q[2] ~ "High score",
        universal_z <= q[1] ~ "Low score",
        TRUE                ~ "Mid"),
      universal_high = score_band == "High score" & breadth_hi >= breadth_thresh,
      universal_low  = score_band == "Low score"  & breadth_lo >= breadth_thresh)
}

#' Subject-level universal-responder score + breadth (fraction of features in
#' the subject's top / bottom tertile).
summarise_responder_subjects <- function(resp_long) {
  resp_long %>%
    dplyr::group_by(subject_accession, arm_name) %>%
    dplyr::summarize(
      n_feat      = dplyr::n(),
      universal_z = mean(z, na.rm = TRUE),
      breadth_hi  = mean(z > stats::qnorm(2 / 3), na.rm = TRUE),
      breadth_lo  = mean(z < stats::qnorm(1 / 3), na.rm = TRUE),
      .groups = "drop")
}

#' Wide standardised matrix (subjects x feature@visit) for clustering.
build_resp_wide <- function(resp_long) {
  resp_long %>%
    dplyr::select(subject_accession, arm_name, key, z) %>%
    tidyr::pivot_wider(
      id_cols     = c(subject_accession, arm_name),
      names_from  = key, values_from = z,
      values_fn   = function(x) mean(x, na.rm = TRUE))
}

#' Data-driven subgroup detection by PAM (k-medoids) with silhouette-chosen k.
#' Returns a list; $ok is FALSE (with $reason) when prerequisites are unmet so
#' the calling chunk can fall back to the transparent tertile rule.
detect_subgroups <- function(resp_wide, resp_subject,
                             has_cluster = requireNamespace("cluster", quietly = TRUE),
                             k_range = 2:4, min_complete = 20L,
                             min_feat = NULL) {
  feat_cols <- setdiff(names(resp_wide), c("subject_accession", "arm_name"))
  if (!length(feat_cols)) {
    return(list(ok = FALSE, reason = "no responder feature columns",
                subgroup_tbl = NULL))
  }

  # Default completeness rule: half the panel (>= 4). Callers running a
  # single-assay panel (1-3 feature@visit columns) pass min_feat = 1 so PAM is
  # not force-skipped down to the tertile fallback.
  if (is.null(min_feat)) min_feat <- max(4, floor(0.5 * length(feat_cols)))
  panel_n  <- rowSums(is.finite(as.matrix(resp_wide[, feat_cols, drop = FALSE])))
  cc_idx   <- which(panel_n >= min_feat)

  if (!has_cluster || length(cc_idx) < min_complete) {
    return(list(
      ok = FALSE,
      reason = sprintf(
        "need package 'cluster' and >= %d sufficiently complete subjects (have %d)",
        min_complete, length(cc_idx)),
      subgroup_tbl = NULL))
  }

  X <- as.matrix(resp_wide[cc_idx, feat_cols, drop = FALSE])
  X[!is.finite(X)] <- 0  # column median z is 0 by construction (Blom)

  sil    <- sapply(k_range, function(k) cluster::pam(X, k = k)$silinfo$avg.width)
  best_k <- k_range[which.max(sil)]
  pamfit <- cluster::pam(X, k = best_k)

  sg <- data.frame(
    subject_accession = resp_wide$subject_accession[cc_idx],
    arm_name = resp_wide$arm_name[cc_idx],
    cluster  = factor(pamfit$clustering),
    universal_z = resp_subject$universal_z[
      match(resp_wide$subject_accession[cc_idx],
            resp_subject$subject_accession)])

  # label clusters by mean universal score (Low ... High)
  lab <- sg %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(mz = mean(universal_z, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(mz)
  lab$label <- if (best_k == 2) {
    c("Low responders", "High responders")
  } else {
    c("Low", paste0("Mid", seq_len(best_k - 2)), "High")[seq_len(best_k)]
  }
  sg$subgroup <- lab$label[match(sg$cluster, lab$cluster)]

  list(ok = TRUE, best_k = best_k, sil = sil, subgroup_tbl = sg)
}


# =============================================================================
# ADDED FOR THE REFACTOR (three features)
#   1. Infant-arm (wP vs aP) stratified within-visit correlations + scatter panel
#   2. Single-assay responder subgroup detection
#   3. Linkage-class collapsing + per-class serological profile tables
# All functions take inputs as explicit arguments (no knit globals) and use the
# same conventions as the rest of the file: long columns subject_accession /
# arm_name / feature / visit_name / log_assay_value; dplyr attached at call time.
# =============================================================================


# ---- (1) infant-arm resolution & stratified correlations --------------------

#' Normalise a free-text priming label to the canonical {wP, aP}. Anything that
#' is not recognisably whole-cell or acellular becomes NA.
.normalize_priming <- function(x) {
  v <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(v))
  out[grepl("^wp", v) | grepl("whole", v)]      <- "wP"
  out[grepl("^ap", v) | grepl("acell", v)]      <- "aP"
  out
}

#' Locate the infant-priming (wP vs aP) variable and return a per-subject
#' lookup data.frame(subject_accession, infant_arm). The maternal arm
#' (`arm_name`, TdaP/TT) is a different variable and is NOT used here.
#'
#' Search order: (a) named candidates in `data_raw`, (b) any other column in
#' `data_raw` whose values normalise cleanly to {wP, aP}, then the same in the
#' clinical table `clin`. Returns NULL (with a warning) if nothing qualifies, so
#' the calling chunk can degrade to the unstratified analysis gracefully.
resolve_infant_arm <- function(data_raw, clin = NULL,
                               candidates = c("infant_arm", "priming",
                                              "infant_priming", "inf_arm",
                                              "arm_infant", "priming_arm",
                                              "prime", "infant_vaccine",
                                              "vaccine_infant"),
                               id_col = "subject_accession") {
  try_frame <- function(df) {
    if (is.null(df) || !is.data.frame(df) || !id_col %in% names(df)) return(NULL)
    nm   <- names(df)
    cand <- nm[tolower(nm) %in% tolower(candidates)]
    others <- setdiff(nm, c(id_col, cand))
    pick <- function(col) {
      z <- .normalize_priming(df[[col]])
      lv <- unique(stats::na.omit(z))
      # accept only a column that is essentially a clean priming label
      if (mean(!is.na(z)) >= 0.5 && length(lv) >= 1 && all(lv %in% c("wP", "aP")))
        data.frame(subject_accession = df[[id_col]], infant_arm = z,
                   stringsAsFactors = FALSE)
      else NULL
    }
    for (col in c(cand, others)) {
      m <- pick(col)
      if (!is.null(m)) return(list(col = col, map = m))
    }
    NULL
  }

  hit <- try_frame(data_raw); src <- "data_raw"
  if (is.null(hit)) { hit <- try_frame(clin); src <- "clin" }
  if (is.null(hit)) {
    warning("resolve_infant_arm(): no infant-arm/priming column found; ",
            "wP/aP stratification and arm colouring will be skipped.")
    return(NULL)
  }
  out <- hit$map
  out <- out[!is.na(out$infant_arm), , drop = FALSE]
  out <- out[!duplicated(out$subject_accession), , drop = FALSE]  # constant/subject
  out$infant_arm <- factor(out$infant_arm, levels = c("wP", "aP"))
  attr(out, "source_col") <- hit$col
  attr(out, "source_tbl") <- src
  out
}

#' Attach the resolved infant arm to a per-visit wide frame (keyed on subject).
attach_infant_arm <- function(w, infant_arm_tbl, id_col = "subject_accession") {
  if (is.null(infant_arm_tbl)) {
    w$infant_arm <- factor(rep(NA_character_, nrow(w)), levels = c("wP", "aP"))
    return(w)
  }
  w$infant_arm <- infant_arm_tbl$infant_arm[
    match(w[[id_col]], infant_arm_tbl$subject_accession)]
  w
}

#' Within-visit correlations computed OVERALL and within each level of a
#' stratifying column (default the infant arm). Each call reuses the existing
#' within_visit_cor() engine, so r / p / BH / n semantics are identical; the
#' overall ("all infants") result is always retained.
within_visit_cor_strata <- function(w, assays, strata_col = "infant_arm",
                                     method = "spearman", min_n = 4L) {
  overall <- within_visit_cor(w, assays, method = method, min_n = min_n)
  by <- list()
  if (strata_col %in% names(w)) {
    lv <- levels(factor(w[[strata_col]]))
    for (g in lv) {
      wg <- w[!is.na(w[[strata_col]]) & w[[strata_col]] == g, , drop = FALSE]
      if (nrow(wg) >= min_n)
        by[[g]] <- within_visit_cor(wg, assays, method = method, min_n = min_n)
    }
  }
  list(overall = overall, by = by, strata_col = strata_col)
}

#' Flatten a within_visit_cor_strata() result to one tidy row per (pair, stratum),
#' with a `stratum` column ("all" plus each arm level). BH p-values are those
#' from within_visit_cor(), i.e. adjusted WITHIN each (timepoint x stratum) pair
#' set.
cor_strata_long <- function(sc, timepoint) {
  mk <- function(cc, stratum) {
    d <- cor_long(cc, timepoint)
    if (nrow(d)) d$stratum <- stratum
    d
  }
  out <- list(mk(sc$overall, "all"))
  for (g in names(sc$by)) out[[length(out) + 1]] <- mk(sc$by[[g]], g)
  dplyr::bind_rows(out)
}

#' Faceted scatter panel of every unique assay pair on one draw, points coloured
#' by `color_col` (the infant arm). Per-arm LM lines (coloured) and an overall
#' LM line (dashed black) make the stratified-vs-combined slopes visible; each
#' facet is annotated with the overall and per-arm Spearman r so the figure
#' corresponds directly to the correlation tables. The LM line is OLS (for
#' readability) while the annotated r is Spearman (matching the tables).
cor_scatter_panel <- function(w, assays, color_col = "infant_arm",
                              method = "spearman", title = NULL, ncol = NULL) {
  cols <- intersect(assays, names(w))
  if (length(cols) < 2) {
    cat("*Fewer than two assays available - no scatter panel.*\n\n")
    return(invisible())
  }
  has_color <- color_col %in% names(w) && any(!is.na(w[[color_col]]))

  pieces <- list(); k <- 0
  for (i in seq_along(cols)) for (j in seq_along(cols)) if (j > i) {
    k  <- k + 1
    xi <- w[[cols[i]]]; yj <- w[[cols[j]]]
    grp <- if (has_color) as.character(w[[color_col]]) else "all"
    ok <- is.finite(xi) & is.finite(yj)
    pieces[[k]] <- data.frame(
      pair = paste(cols[i], "vs", cols[j]),
      x = xi[ok], y = yj[ok], grp = grp[ok], stringsAsFactors = FALSE)
  }
  df <- dplyr::bind_rows(pieces)
  if (!nrow(df)) { cat("*No finite pairs to plot.*\n\n"); return(invisible()) }
  df$pair <- factor(df$pair, levels = unique(df$pair))
  if (has_color)
    df$grp <- factor(df$grp, levels = intersect(c("wP", "aP"), unique(df$grp)))

  # per-facet Spearman r annotations (overall + per arm)
  pair_levels <- levels(df$pair)
  lab_df <- data.frame(pair = factor(pair_levels, levels = pair_levels),
                       lab = NA_character_, stringsAsFactors = FALSE)
  for (pl in pair_levels) {
    sub <- df[df$pair == pl, , drop = FALSE]
    r_all <- suppressWarnings(stats::cor(sub$x, sub$y, method = method,
                                         use = "complete.obs"))
    txt <- sprintf("all r=%.2f", r_all)
    if (has_color) for (g in levels(df$grp)) {
      sg <- sub[sub$grp == g, , drop = FALSE]
      if (nrow(sg) >= 4) {
        rg <- suppressWarnings(stats::cor(sg$x, sg$y, method = method,
                                          use = "complete.obs"))
        txt <- paste0(txt, sprintf("\n%s r=%.2f", g, rg))
      }
    }
    lab_df$lab[lab_df$pair == pl] <- txt
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y))
  if (has_color) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = grp), alpha = 0.5, size = 1.1) +
      ggplot2::geom_smooth(ggplot2::aes(color = grp, group = grp),
                           method = "lm", formula = y ~ x, se = FALSE,
                           linewidth = 0.7) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                           color = "black", linetype = "dashed", linewidth = 0.7) +
      ggplot2::scale_color_manual(values = c(wP = "#2166ac", aP = "#d01c8b"),
                                  name = "Infant arm", drop = FALSE)
  } else {
    p <- p +
      ggplot2::geom_point(alpha = 0.5, size = 1.1, color = "#2166ac") +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                           color = "black")
  }
  p <- p +
    ggplot2::geom_text(data = lab_df,
                       ggplot2::aes(x = -Inf, y = Inf, label = lab),
                       hjust = -0.05, vjust = 1.15, size = 2.5,
                       inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ pair, scales = "free", ncol = ncol) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "top")
  print(p)
  invisible(df)
}


# ---- (2) single-assay responder subgroup detection -------------------------

#' Build a responder panel from ONE assay across `visits`, score it, and detect
#' subgroups on that single feature's feature@visit columns. Whole-bacterial
#' assays are gated by the availability map; binding assays are taken everywhere.
#' min_feat = 1 lets PAM run on the (1-3 column) single-assay matrix instead of
#' being forced to the tertile fallback.
single_feature_subgroup <- function(data_raw, visits, feature, avail_feats,
                                     whole_map = WHOLE_ASSAY_VISITS_DEFAULT,
                                     breadth_thresh = 0.6,
                                     has_cluster = requireNamespace("cluster",
                                                                    quietly = TRUE)) {
  rfeat <- resolve_feat(feature, avail_feats)
  if (is.na(rfeat))
    return(list(ok = FALSE, feature = feature,
                reason = sprintf("feature '%s' absent from panel", feature)))
  is_whole <- toupper(rfeat) %in% toupper(names(whole_map))
  binding  <- if (is_whole) character(0) else rfeat
  whole    <- if (is_whole) rfeat       else character(0)
  long <- build_responder_panel(data_raw, visits, binding, whole,
                                avail_feats, whole_map)
  if (!nrow(long))
    return(list(ok = FALSE, feature = rfeat,
                reason = "no data for feature at requested visits"))
  subject  <- add_tertile_bands(summarise_responder_subjects(long), breadth_thresh)
  wide     <- build_resp_wide(long)
  by_visit <- responder_score_by_visit(long)
  sg <- detect_subgroups(wide, subject, has_cluster = has_cluster, min_feat = 1L)
  list(ok = TRUE, feature = rfeat, long = long, subject = subject,
       wide = wide, by_visit = by_visit, subgroups = sg)
}


# ---- (3) linkage-class collapsing & per-class serological profiles ----------

#' Collapse the detailed subgroup labels ("High responders"/"High score",
#' "Low responders"/"Low score", "Mid"/"Mid1"...) to {High, Low, Mid}.
collapse_resp_level <- function(group) {
  g <- as.character(group)
  out <- ifelse(grepl("High", g), "High",
         ifelse(grepl("Low",  g), "Low", "Mid"))
  out
}

#' Per (class, feature@visit) standardised profile: mean z, sd, n. `resp_long`
#' is a (combined) responder long table with columns subject_accession, key, z;
#' `class_tbl` carries subject_accession + a class column.
class_profile_table <- function(resp_long, class_tbl,
                                id_col = "subject_accession",
                                class_col = "link_class") {
  ct <- class_tbl[, c(id_col, class_col)]
  ct <- ct[!duplicated(ct[[id_col]]), , drop = FALSE]
  d  <- dplyr::left_join(resp_long, ct, by = id_col)
  names(d)[names(d) == class_col] <- "class"
  d <- d[!is.na(d$class), , drop = FALSE]
  d %>%
    dplyr::group_by(class, key) %>%
    dplyr::summarize(mean_z = mean(z, na.rm = TRUE),
                     sd_z   = stats::sd(z, na.rm = TRUE),
                     n      = sum(is.finite(z)), .groups = "drop")
}

#' For each feature@visit key, contrast a focal class (default "Low_Low")
#' against all other subjects: mean z in each set, the difference, and a
#' Wilcoxon rank-sum p. Answers "which assays most distinguish the focal class?"
distinguishing_features <- function(resp_long, class_tbl, focal = "Low_Low",
                                    id_col = "subject_accession",
                                    class_col = "link_class") {
  ct <- class_tbl[, c(id_col, class_col)]
  ct <- ct[!duplicated(ct[[id_col]]), , drop = FALSE]
  d  <- dplyr::left_join(resp_long, ct, by = id_col)
  names(d)[names(d) == class_col] <- "class"
  d <- d[!is.na(d$class), , drop = FALSE]
  d$is_focal <- d$class == focal
  keys <- unique(d$key)
  rows <- lapply(keys, function(kk) {
    dk <- d[d$key == kk, , drop = FALSE]
    zf <- dk$z[dk$is_focal];  zo <- dk$z[!dk$is_focal]
    zf <- zf[is.finite(zf)];  zo <- zo[is.finite(zo)]
    p <- NA_real_
    if (length(zf) >= 2 && length(zo) >= 2)
      p <- tryCatch(stats::wilcox.test(zf, zo)$p.value, error = function(e) NA_real_)
    data.frame(key = kk,
               n_focal = length(zf), n_other = length(zo),
               mean_focal = if (length(zf)) mean(zf) else NA_real_,
               mean_other = if (length(zo)) mean(zo) else NA_real_,
               diff = (if (length(zf)) mean(zf) else NA_real_) -
                      (if (length(zo)) mean(zo) else NA_real_),
               p = p, stringsAsFactors = FALSE)
  })
  res <- dplyr::bind_rows(rows)
  if (nrow(res)) {
    res$p_BH <- stats::p.adjust(res$p, method = "BH")
    res <- res[order(res$diff), , drop = FALSE]
  }
  res
}


# =============================================================================
# ADDED FOR THE SECOND REFACTOR (Stage B2 figures + concordance + low/low map)
#   1. Per single-assay distribution figure: Low/High subgroup + maternal arm
#   2. Cross-assay concordance of the single-assay High/Low classifications
#   3. (used in Stage D) per-subject single-assay label table for the low/low map
# =============================================================================

# ---- (1) single-assay distribution figure ----------------------------------

#' Distribution of a single-assay responder score, with subjects coloured by
#' their Low/High subgroup and the maternal arm (TT vs TdaP) shown as the two
#' rows of the same response axis. Per-arm violins show the distribution shape;
#' dashed lines mark each subgroup's median score. Falls back to the score-band
#' labels if PAM was skipped for this assay.
single_feature_dist_plot <- function(sf) {
  if (!isTRUE(sf$ok)) return(invisible())
  sg <- sf$subgroups
  if (isTRUE(sg$ok)) {
    d <- sg$subgroup_tbl
    d$level <- collapse_resp_level(d$subgroup)
  } else {
    d <- sf$subject
    d$level <- collapse_resp_level(d$score_band)
  }
  d <- d[!is.na(d$universal_z) & !is.na(d$arm_name), , drop = FALSE]
  if (!nrow(d)) { cat("*No classified subjects to plot.*\n\n"); return(invisible()) }
  d$level    <- factor(d$level, levels = c("Low", "Mid", "High"))
  d$arm_name <- factor(d$arm_name)

  meds <- stats::aggregate(universal_z ~ level, data = d, FUN = stats::median)
  meds <- meds[is.finite(meds$universal_z), , drop = FALSE]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = universal_z, y = arm_name)) +
    ggplot2::geom_violin(ggplot2::aes(group = arm_name), orientation = "y",
                         fill = "grey92", color = "grey75",
                         width = 0.85, alpha = 0.5) +
    ggplot2::geom_vline(data = meds,
                        ggplot2::aes(xintercept = universal_z, color = level),
                        linetype = "dashed", linewidth = 0.6,
                        show.legend = FALSE) +
    ggplot2::geom_jitter(ggplot2::aes(color = level),
                         height = 0.22, width = 0, size = 1.7, alpha = 0.75) +
    ggplot2::scale_color_manual(values = c(Low = "#d01c8b", Mid = "grey60",
                                           High = "#2166ac"),
                                name = "Subgroup", drop = FALSE) +
    ggplot2::labs(
      title = sprintf("%s - single-assay response by subgroup and maternal arm",
                      sf$feature),
      subtitle = "rows = maternal arm; colour = Low/High subgroup; dashed = subgroup medians",
      x = "Single-assay responder score (mean within-timepoint z)",
      y = "Maternal arm") +
    ggplot2::theme_bw(base_size = 10)
  print(p)
  invisible(d)
}


# ---- (2) cross-assay concordance of single-assay subgroups ------------------

#' Long + wide per-subject Low/High label table across the single-assay
#' subgroup results. PAM labels ("High/Low responders") and the tertile-band
#' fallback ("High/Low score") are both collapsed to High/Low; Mid is dropped.
single_feature_label_table <- function(single_feature_results) {
  rows <- list(); k <- 0
  for (f in names(single_feature_results)) {
    sf <- single_feature_results[[f]]
    if (!isTRUE(sf$ok)) next
    sg <- sf$subgroups
    if (isTRUE(sg$ok)) {
      ids <- sg$subgroup_tbl$subject_accession
      lev <- collapse_resp_level(sg$subgroup_tbl$subgroup)
    } else {
      ids <- sf$subject$subject_accession
      lev <- collapse_resp_level(sf$subject$score_band)
    }
    keep <- !is.na(lev) & lev %in% c("High", "Low")
    if (!any(keep)) next
    k <- k + 1
    rows[[k]] <- data.frame(subject_accession = ids[keep],
                            assay = sf$feature, level = lev[keep],
                            stringsAsFactors = FALSE)
  }
  long <- dplyr::bind_rows(rows)
  wide <- if (nrow(long))
    tidyr::pivot_wider(long, id_cols = subject_accession,
                       names_from = assay, values_from = level)
  else dplyr::tibble(subject_accession = character())
  list(long = long, wide = wide)
}

#' Cohen's kappa and raw % agreement for two categorical label vectors
#' (pairwise-complete). Returns c(kappa, agree, n).
cohen_kappa <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- a[ok]; b <- b[ok]
  if (!length(a)) return(c(kappa = NA_real_, agree = NA_real_, n = 0))
  lv  <- union(unique(a), unique(b))
  tab <- table(factor(a, levels = lv), factor(b, levels = lv))
  n   <- sum(tab)
  po  <- sum(diag(tab)) / n
  pe  <- sum((rowSums(tab) / n) * (colSums(tab) / n))
  kap <- if (abs(1 - pe) < 1e-12) NA_real_ else (po - pe) / (1 - pe)
  c(kappa = kap, agree = po, n = n)
}

#' Pairwise concordance (n, % agreement, Cohen's kappa) for every pair of the
#' single-assay label columns in `labels_wide`.
pairwise_concordance <- function(labels_wide, assays) {
  assays <- intersect(assays, names(labels_wide))
  out <- list(); k <- 0
  for (i in seq_along(assays)) for (j in seq_along(assays)) if (j > i) {
    r <- cohen_kappa(labels_wide[[assays[i]]], labels_wide[[assays[j]]])
    k <- k + 1
    out[[k]] <- data.frame(assay_x = assays[i], assay_y = assays[j],
                           n = unname(r["n"]), pct_agree = unname(r["agree"]),
                           kappa = unname(r["kappa"]), stringsAsFactors = FALSE)
  }
  if (!length(out))
    return(data.frame(assay_x = character(), assay_y = character(),
                      n = integer(), pct_agree = numeric(), kappa = numeric()))
  dplyr::bind_rows(out)
}

#' Square symmetric kappa matrix (diagonal 1) from a pairwise_concordance() df.
kappa_matrix <- function(pc, assays) {
  M <- matrix(NA_real_, length(assays), length(assays),
              dimnames = list(assays, assays))
  diag(M) <- 1
  for (r in seq_len(nrow(pc))) {
    M[pc$assay_x[r], pc$assay_y[r]] <- pc$kappa[r]
    M[pc$assay_y[r], pc$assay_x[r]] <- pc$kappa[r]
  }
  M
}


# =============================================================================
# ADDED FOR STAGE E (pan-assay concordance classes + low/low vs always-low)
#   Reusable assay-level profile table and heatmap, used both for the
#   always-high / mixed / always-low classes and for the low/low-vs-always-low
#   contrast.
# =============================================================================

#' Per (group, assay) mean within-timepoint z, averaging over visits. `prof_long`
#' is the combined responder long table (columns subject_accession, key, z);
#' `group_tbl` maps subject_accession -> a grouping column `group_col`.
assay_profile_table <- function(prof_long, group_tbl, group_col,
                                id_col = "subject_accession") {
  gt <- group_tbl[, c(id_col, group_col)]
  gt <- gt[!duplicated(gt[[id_col]]), , drop = FALSE]
  d  <- dplyr::left_join(prof_long, gt, by = id_col)
  names(d)[names(d) == group_col] <- "grp"
  d <- d[!is.na(d$grp), , drop = FALSE]
  d$assay <- sub("@.*$", "", d$key)
  d %>%
    dplyr::group_by(grp, assay) %>%
    dplyr::summarize(mean_z = mean(z, na.rm = TRUE),
                     n = sum(is.finite(z)), .groups = "drop")
}

#' Tile heatmap of an assay_profile_table()-style frame (columns grp, assay,
#' mean_z). `grp_levels` optionally fixes the row order. Returns a ggplot.
profile_heatmap <- function(prof_assay, title, grp_levels = NULL,
                            assay_levels = NULL) {
  d <- prof_assay
  if (!is.null(grp_levels))   d$grp   <- factor(d$grp,   levels = grp_levels)
  if (!is.null(assay_levels)) d$assay <- factor(d$assay, levels = assay_levels)
  ggplot2::ggplot(d, ggplot2::aes(assay, grp, fill = mean_z)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mean_z)), size = 2.6) +
    ggplot2::scale_fill_gradient2(low = "#d01c8b", mid = "white", high = "#2166ac",
                                  midpoint = 0, name = "mean z") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


# =============================================================================
# TRANSITION-RESPONSE SUBGROUP HETEROGENEITY
#   Added for the per-antigen x per-transition heterogeneity output.
#   For one antigen and one within-subject transition (from -> to), cluster
#   subjects on their per-analyte change (delta = log10 later - earlier) to ask
#   whether the transition response splits into distinct subgroups rather than a
#   single population. PAM with silhouette-chosen k when 'cluster' is available
#   and enough subjects exist; otherwise a transparent tertile fallback on the
#   total-IgG delta. Reuses paired_change() from R/serology_helpers.R.
# =============================================================================

#' Subgroup detection for ONE antigen across ONE transition.
#' Returns list(ok, antigen, from, to, method, best_k, silhouette, n, assignments).
transition_subgroups <- function(data_raw, antigen, from_visit, to_visit,
                                  analytes = c("IgG", "IgG1", "IgG2", "IgG3", "IgG4"),
                                  k_range = 2:4, min_complete = 20L,
                                  has_cluster = requireNamespace("cluster", quietly = TRUE)) {
  feats <- paste(antigen, analytes, sep = "_")
  pc <- paired_change(data_raw, from_visit, to_visit)
  pc <- pc[as.character(pc$feature) %in% feats, , drop = FALSE]
  if (!nrow(pc))
    return(list(ok = FALSE, antigen = antigen,
                reason = "no paired data for this antigen/transition"))

  arm  <- unique(pc[, c("subject_accession", "arm_name")])
  wide <- tidyr::pivot_wider(
            pc[, c("subject_accession", "analyte", "delta")],
            names_from = analyte, values_from = delta,
            values_fn = function(x) mean(x, na.rm = TRUE))
  feat_cols <- setdiff(names(wide), "subject_accession")
  total_col <- if ("IgG" %in% feat_cols) "IgG" else feat_cols[1]
  M <- as.matrix(wide[, feat_cols, drop = FALSE])
  td <- M[, total_col]

  Z <- suppressWarnings(scale(M))
  Z[!is.finite(Z)] <- 0

  method <- NA_character_; best_k <- NA_integer_; sil <- NA_real_; cl <- NULL
  if (has_cluster && nrow(Z) >= min_complete) {
    sils <- vapply(k_range, function(k) {
      if (k >= nrow(Z)) return(NA_real_)
      tryCatch(cluster::pam(Z, k = k)$silinfo$avg.width,
               error = function(e) NA_real_)
    }, numeric(1))
    if (any(is.finite(sils))) {
      best_k <- k_range[which.max(replace(sils, !is.finite(sils), -Inf))]
      fit    <- cluster::pam(Z, k = best_k)
      cl     <- factor(fit$clustering)
      sil    <- max(sils, na.rm = TRUE)
      method <- sprintf("PAM k=%d", best_k)
    }
  }
  if (is.null(cl)) {                                   # tertile / median fallback
    cl <- tryCatch({
      qs <- stats::quantile(td, c(1/3, 2/3), na.rm = TRUE)
      br <- unique(c(-Inf, qs, Inf))
      if (length(br) < 3) stop("degenerate")
      cut(td, br, labels = as.character(seq_len(length(br) - 1)))
    }, error = function(e)
      factor(ifelse(td > stats::median(td, na.rm = TRUE), "2", "1")))
    method <- "tertiles (total-IgG delta)"
  }

  res <- data.frame(subject_accession = wide$subject_accession,
                    cluster = cl, total_delta = td, stringsAsFactors = FALSE)
  res <- dplyr::left_join(res, arm, by = "subject_accession")
  lab <- res %>% dplyr::group_by(cluster) %>%
    dplyr::summarize(m = mean(total_delta, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(m)
  k <- nrow(lab)
  lab$subgroup <- if (k == 2) c("Low", "High") else
    c("Low", paste0("Mid", seq_len(max(k - 2, 0))), "High")[seq_len(k)]
  res$subgroup <- factor(lab$subgroup[match(res$cluster, lab$cluster)],
                         levels = lab$subgroup)
  list(ok = TRUE, antigen = antigen, from = from_visit, to = to_visit,
       method = method, best_k = best_k, silhouette = sil,
       n = nrow(res), assignments = res)
}

#' Scan antigens x transitions. `transitions` is a list of named character
#' vectors c(from=, to=, label=). Returns list(summary, assignments).
transition_subgroup_scan <- function(data_raw, antigens, transitions,
                                     analytes = c("IgG", "IgG1", "IgG2", "IgG3", "IgG4"),
                                     ...) {
  rows <- list(); assigns <- list()
  for (tr in transitions) for (a in antigens) {
    sg <- transition_subgroups(data_raw, a, tr["from"], tr["to"], analytes, ...)
    if (!isTRUE(sg$ok)) next
    asn  <- sg$assignments
    bysg <- asn %>% dplyr::group_by(subgroup) %>%
      dplyr::summarize(n = dplyr::n(),
                       md = round(stats::median(total_delta, na.rm = TRUE), 3),
                       .groups = "drop")
    tab   <- table(asn$subgroup, asn$arm_name)
    p_arm <- tryCatch(
      if (all(dim(tab) >= 2))
        suppressWarnings(stats::fisher.test(tab, simulate.p.value = TRUE,
                                            B = 2000)$p.value) else NA_real_,
      error = function(e) NA_real_)
    rows[[length(rows) + 1]] <- data.frame(
      transition  = unname(tr["label"]),
      antigen     = a,
      method      = sg$method,
      n           = sg$n,
      silhouette  = ifelse(is.na(sg$silhouette), NA_real_, round(sg$silhouette, 3)),
      subgroups   = paste(sprintf("%s: n=%d (med \u0394%.2f)",
                                  bysg$subgroup, bysg$n, bysg$md), collapse = "; "),
      p_arm_assoc = signif(p_arm, 3),
      stringsAsFactors = FALSE)
    asn$transition <- unname(tr["label"]); asn$antigen <- a
    assigns[[length(assigns) + 1]] <- asn
  }
  list(summary     = dplyr::bind_rows(rows),
       assignments = dplyr::bind_rows(assigns))
}
