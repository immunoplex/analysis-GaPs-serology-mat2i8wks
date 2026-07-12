# =====================================================================
# R/infmon2_effector_helpers.R
# ---------------------------------------------------------------------
# 8-weeks (infant baseline) maternal-vaccination-effect component:
# whole-bacterium assays (SBA, WT_IgG, PTNA) + antibody effector
# functions (ADCD, ADCP, ADNP) tested for TdaP-vs-TT differences and
# rendered as a Figure-2C/D-style forest plot.
#
# Sourced AFTER config/endpoints.R and R/serology_helpers.R (it reuses
# ANTIGEN_ORDER, MATERNAL_ARMS, pvalue_to_label(), and the nt tracker).
#
# The statistical test is IDENTICAL to the response-bubble grid
# (04_response_bubbles.Rmd / serology_helpers::bubble_arm_effect):
#   * Wilcoxon rank-sum of TT vs TdaP on log_assay_value (exact = FALSE)
#   * multiplicity via Benjamini-Hochberg within each assay class (default)
# For the forest plot we additionally take the point estimate + interval
# that NATIVELY pairs with that test: the Hodges-Lehmann location shift
# and its Wilcoxon conf.int (wilcox.test(conf.int = TRUE)). The p-value
# from that call is exactly the bubble p-value. Orientation is TdaP - TT
# (positive = higher with maternal TdaP), matching how Saso et al. 2026
# Figure 2C/D sign the maternal-vaccination coefficient.
#
# NOTE on antigen/analyte: the whole-bacterium features are stored as
# antigen = "WHOLE", analyte = SBA / WT_IgG / PTNA, so the feature string
# WHOLE_WT_IgG must NOT be split on "_" (that would give analyte "WT").
# We therefore read the data's own antigen/analyte columns rather than
# re-deriving them from the feature string.
# =====================================================================


# ---------------------------------------------------------------------
# 1. FEATURE VOCABULARY
# ---------------------------------------------------------------------

# The three whole-bacterium assays and their display labels.
INFMON2_WHOLE_ASSAYS <- c("WHOLE_SBA", "WHOLE_WT_IgG", "WHOLE_PTNA")
INFMON2_WHOLE_LABELS <- c(WHOLE_SBA = "SBA", WHOLE_WT_IgG = "WT IgG",
                          WHOLE_PTNA = "PTNA")

# The full possible antigen list, and the antigens each effector actually
# has results for. ADCP / ADNP have NO FHA and NO FIM assays.
INFMON2_POSSIBLE_ANTIGENS <- c("PT", "FHA", "PRN", "FIM", "TT", "DT")
INFMON2_EFFECTOR_ANTIGENS <- list(
  ADCD = c("PT", "FHA", "PRN", "FIM", "TT", "DT"),
  ADCP = c("PT", "PRN", "TT", "DT"),
  ADNP = c("PT", "PRN", "TT", "DT")
)

# Assemble the full antigen_analyte feature list this component tests.
infmon2_effector_features <- function(
    whole     = INFMON2_WHOLE_ASSAYS,
    effectors = INFMON2_EFFECTOR_ANTIGENS) {
  eff <- unlist(lapply(names(effectors), function(fn)
    paste(effectors[[fn]], fn, sep = "_")), use.names = FALSE)
  c(whole, eff)
}

# ---------------------------------------------------------------------
# 2. DISPLAY CONSTANTS (colours, panel + row order)
# ---------------------------------------------------------------------
# One colour per whole-bacterium assay and one per effector function, so
# the legend differentiates the whole-bacterium assays AND each effector.
INFMON2_SERIES_COLORS <- c(
  "SBA"    = "#B2182B",   # whole-bacterium assays (warm)
  "WT IgG" = "#E8853A",
  "PTNA"   = "#F2C14E",
  "ADCD"   = "#2166AC",   # effector functions (cool)
  "ADCP"   = "#1B9E77",
  "ADNP"   = "#762A83"
)

# Vertical panel order (top-to-bottom) and within-panel row order.
INFMON2_CLASS_ORDER <- c("Whole bacterium assays", "ADCD", "ADCP", "ADNP")
INFMON2_ROW_ORDER   <- c("SBA", "WT IgG", "PTNA",              # whole
                         "PT", "FHA", "PRN", "FIM", "TT", "DT") # effector antigens


# ---------------------------------------------------------------------
# 3. THE TEST  (same test as the response bubbles) + HL point/interval
# ---------------------------------------------------------------------
# Returns one row per feature with: per-arm n, per-arm medians, the bubble
# effect (median(TdaP) - median(TT)) and direction, the Hodges-Lehmann shift
# and its Wilcoxon CI, the (identical) Wilcoxon p, the BH-adjusted p,
# and the class / colour-series / row-label used by the forest plot.
# Orientation is TdaP - TT: with TT as the reference level this is the
# effect of TdaP relative to TT (positive = higher with maternal TdaP).
infmon2_arm_effect <- function(data,
                               visit           = "InfMon2",  # data value; displayed as "8 weeks"
                               features        = infmon2_effector_features(),
                               infant_arm_sel  = NULL,       # NULL = pool aP + wP
                               p_adjust_method = "BH",       # Benjamini-Hochberg
                               by_class        = TRUE,       # correct within each assay class
                               conf_level      = 0.95) {

  ua <- unique(as.character(data$arm_name))
  if (!all(c("TT", "TdaP") %in% ua))
    stop("infmon2_arm_effect(): arm_name must contain both 'TT' and 'TdaP'.")

  d <- data %>%
    dplyr::filter(as.character(visit_name) == visit,
                  as.character(arm_name)   %in% c("TT", "TdaP"),
                  as.character(feature)    %in% features)
  if (!is.null(infant_arm_sel))
    d <- d[as.character(d$infant_arm) == infant_arm_sel, , drop = FALSE]
  d <- dplyr::distinct(d, subject_accession, feature, .keep_all = TRUE)

  feats <- intersect(features, unique(as.character(d$feature)))
  med   <- function(v) if (length(v)) stats::median(v) else NA_real_

  rows <- lapply(feats, function(ft) {
    sub <- d[as.character(d$feature) == ft, , drop = FALSE]
    x <- stats::na.omit(sub$log_assay_value[as.character(sub$arm_name) == "TdaP"]) # TdaP
    y <- stats::na.omit(sub$log_assay_value[as.character(sub$arm_name) == "TT"])   # TT
    n_TdaP <- length(x); n_TT <- length(y)

    est <- ci_lo <- ci_hi <- p <- NA_real_
    if (n_TT >= 2 && n_TdaP >= 2) {
      # first arg = TdaP so the shift is oriented TdaP - TT: with TT as the
      # reference level, this is the effect of TdaP relative to TT
      # (positive = higher with maternal TdaP), as in Saso et al. Fig 2C/D.
      wt <- tryCatch(
        suppressWarnings(stats::wilcox.test(x, y, conf.int = TRUE,
                                            conf.level = conf_level, exact = FALSE)),
        error = function(e) NULL)
      # if the CI estimator fails, still recover the (bubble) p-value
      if (is.null(wt))
        p <- tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value,
                      error = function(e) NA_real_)
      else {
        p     <- wt$p.value
        est   <- unname(wt$estimate)      # Hodges-Lehmann shift, TdaP - TT (log10)
        ci_lo <- wt$conf.int[1]
        ci_hi <- wt$conf.int[2]
      }
    }
    data.frame(
      feature     = ft,
      antigen     = as.character(sub$antigen[1]),
      analyte     = as.character(sub$analyte[1]),
      n_TT        = n_TT,
      n_TdaP      = n_TdaP,
      median_TT   = med(y),
      median_TdaP = med(x),
      effect_size = med(x) - med(y),      # median(TdaP) - median(TT), TT referent
      hl_shift    = est,                  # HL location shift, TdaP - TT (log10)
      ci_lo       = ci_lo,
      ci_hi       = ci_hi,
      p_value     = p,
      stringsAsFactors = FALSE)
  })
  stats_df <- do.call(rbind, rows)
  if (is.null(stats_df) || !nrow(stats_df)) return(stats_df)

  # direction (semantics as in bubble_arm_effect(); here effect = TdaP - TT)
  stats_df$direction <- dplyr::case_when(
    is.na(stats_df$effect_size) ~ "No effect",
    stats_df$effect_size > 0    ~ "Higher TdaP",
    stats_df$effect_size < 0    ~ "Lower TdaP",
    TRUE                        ~ "No effect")

  # classification / colour series / row label
  eff_names <- names(INFMON2_EFFECTOR_ANTIGENS)
  is_eff    <- stats_df$analyte %in% eff_names
  stats_df$class  <- ifelse(is_eff, stats_df$analyte, "Whole bacterium assays")
  stats_df$series <- ifelse(is_eff, stats_df$analyte,
                            unname(INFMON2_WHOLE_LABELS[stats_df$feature]))
  stats_df$row_label <- ifelse(is_eff, stats_df$antigen,
                               unname(INFMON2_WHOLE_LABELS[stats_df$feature]))

  # multiplicity — Benjamini-Hochberg WITHIN each assay class (Whole / ADCD /
  # ADCP / ADNP) by default, so each functional family is its own test family.
  # Set by_class = FALSE to correct across the whole panel instead.
  stats_df$p_adj <- NA_real_
  if (isTRUE(by_class)) {
    for (cl in unique(stats_df$class)) {
      idx <- which(stats_df$class == cl)
      stats_df$p_adj[idx] <- stats::p.adjust(stats_df$p_value[idx], method = p_adjust_method)
    }
  } else {
    stats_df$p_adj <- stats::p.adjust(stats_df$p_value, method = p_adjust_method)
  }
  stats_df$significant   <- !is.na(stats_df$p_adj) & stats_df$p_adj < 0.05
  stats_df$fold_TdaP_vs_TT <- 10^(stats_df$hl_shift)  # TdaP/TT ratio on the linear (MFI/titre) scale

  # ordering for tables + facets
  stats_df$class     <- factor(stats_df$class, levels = INFMON2_CLASS_ORDER)
  stats_df$series    <- factor(stats_df$series, levels = names(INFMON2_SERIES_COLORS))
  stats_df$row_label <- factor(stats_df$row_label,
                               levels = rev(INFMON2_ROW_ORDER))  # rev => first = top
  stats_df <- stats_df[order(stats_df$class,
                             match(as.character(stats_df$row_label),
                                   INFMON2_ROW_ORDER)), , drop = FALSE]
  attr(stats_df, "conf_level")      <- conf_level
  attr(stats_df, "p_adjust_method") <- p_adjust_method
  attr(stats_df, "by_class")        <- isTRUE(by_class)
  stats_df
}


# ---------------------------------------------------------------------
# 4. THE TABLE
# ---------------------------------------------------------------------
infmon2_effect_kable <- function(stats, caption) {
  padj_lab <- switch(attr(stats, "p_adjust_method") %||% "BH",
                     BH = "p (BH)", fdr = "p (BH)", bonferroni = "p (Bonf)", "p (adj)")
  tab <- stats %>%
    dplyr::arrange(class, match(as.character(row_label), INFMON2_ROW_ORDER)) %>%
    dplyr::transmute(
      Class    = as.character(class),
      Assay    = as.character(series),
      Antigen  = ifelse(class == "Whole bacterium assays", "\u2014", as.character(antigen)),
      n_TT     = n_TT,
      n_TdaP   = n_TdaP,
      med_TT   = round(median_TT,   3),
      med_TdaP = round(median_TdaP, 3),
      hl       = round(hl_shift,    3),
      ci       = ifelse(is.na(ci_lo), "\u2014", sprintf("[%.2f, %.2f]", ci_lo, ci_hi)),
      fold     = ifelse(is.na(fold_TdaP_vs_TT), "\u2014", sprintf("%.2f", fold_TdaP_vs_TT)),
      dir      = direction,
      p        = signif(p_value, 3),
      p_adj    = signif(p_adj,   3),
      sig      = pvalue_to_label(p_adj))
  knitr::kable(
    tab,
    col.names = c("Class", "Assay", "Antigen", "n TT", "n TdaP",
                  "Median TT", "Median TdaP", "HL shift (TdaP\u2212TT)",
                  paste0(round(100 * (attr(stats, "conf_level") %||% 0.95)), "% CI"),
                  "Fold (TdaP/TT)", "Direction", "p", padj_lab, ""),
    caption = caption,
    align = c("l", "l", "l", "r", "r", "r", "r", "r", "c", "r", "l", "r", "r", "c"))
}


# ---------------------------------------------------------------------
# 5. THE FOREST PLOT  (Saso et al. 2026 Figure 2C/D style)
# ---------------------------------------------------------------------
# Stacked panels (one per assay class) sharing the x-axis; each row is a
# point (Hodges-Lehmann shift) with a 95% CI line; dashed reference at 0;
# adjusted p-values printed at the right; colour differentiates the whole-
# bacterium assays and each effector function.
infmon2_forest_plot <- function(stats,
                                title  = "Maternal TdaP vs TT effect at 8 weeks (infant baseline)",
                                x_lab  = NULL,
                                colors = INFMON2_SERIES_COLORS,
                                show_p = c("adjusted", "raw")) {
  show_p <- match.arg(show_p)
  df <- stats[!is.na(stats$hl_shift), , drop = FALSE]
  if (!nrow(df)) stop("infmon2_forest_plot(): no estimable rows to plot.")

  conf_level <- attr(stats, "conf_level") %||% 0.95
  pcol       <- if (show_p == "adjusted") df$p_adj else df$p_value
  df$p_lab   <- ifelse(is.na(pcol), "",
                ifelse(pcol < 0.001, "p<0.001", paste0("p=", formatC(pcol, format = "g", digits = 2))))

  # deterministic x position for the right-hand p-value column
  rng  <- range(c(df$ci_lo, df$ci_hi, 0), na.rm = TRUE)
  span <- diff(rng); if (!is.finite(span) || span == 0) span <- 1
  p_x  <- rng[2] + 0.08 * span

  if (is.null(x_lab))
    x_lab <- expression("Hodges\u2013Lehmann shift, TdaP \u2212 TT  (log"[10] * " scale)")

  ggplot2::ggplot(df, ggplot2::aes(x = hl_shift, y = row_label, colour = series)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey45", linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(size = significant)) +
    ggplot2::geom_text(ggplot2::aes(label = p_lab), x = p_x, hjust = 0,
                       size = 2.9, colour = "grey25", show.legend = FALSE) +
    ggplot2::facet_grid(class ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_colour_manual(name = NULL, values = colors, drop = FALSE) +
    ggplot2::scale_size_manual(values = c(`FALSE` = 2, `TRUE` = 3.4), guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.30))) +
    ggplot2::labs(
      title = title, x = x_lab, y = NULL,
      caption = sprintf(paste0(
        "Wilcoxon rank-sum (TdaP vs TT) on log10 assay values \u2014 the same test as the response-bubble grid (TT referent).\n",
        "Point = Hodges\u2013Lehmann location shift (TdaP \u2212 TT); line = %d%% CI; %s p shown. Benjamini\u2013Hochberg within each assay class; ",
        "larger dot = adjusted p < 0.05.\nDashed line = no difference; positive = higher with maternal TdaP."),
        round(100 * conf_level), show_p)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.placement    = "outside",
      strip.background    = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text.y.left  = ggplot2::element_text(angle = 0, face = "bold"),
      plot.title         = ggplot2::element_text(face = "bold"),
      plot.caption       = ggplot2::element_text(hjust = 0, size = 8, colour = "grey30"),
      legend.position    = "bottom") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 3)))
}

# small null-coalescing helper (kept local so the file is self-contained)
`%||%` <- function(a, b) if (is.null(a)) b else a
