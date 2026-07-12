## =============================================================================
## R/C2_C2_interval_helpers.R  --  Phase 10 engine
## Vaccination -> delivery interval adjustment for the maternal responder analysis.
##
## SELF-CONTAINED by design. Depends only on: base R, stats, splines, and
## (optionally) glmnet / dplyr / ggplot2 / patchwork. It is sourced READ-ONLY
## beside the Phase 07/08 helpers and does NOT modify any of them. In particular
## it re-implements a forced-unpenalised elastic net rather than editing
## R/C3_responder_helpers.R::elastic_net_signature(), so the existing
## pipeline is untouched.
##
## Convention throughout: w = weeks between vaccination and delivery
##   w_weeks = vaccine_birth_interval_days / 7
## Antibody values are log10 (the `log_assay_value` field of the recoded long
## serology object, and the v_from / v_to fields of paired_change()).
## =============================================================================

suppressPackageStartupMessages({
  library(splines)
})

## ---- small internal utilities -----------------------------------------------

## Mann-Whitney AUC for a numeric score against a 0/1 outcome (self-contained so
## the helper does not depend on C1_C1_serology_helpers::auc_binary()).
.auc01 <- function(score, y01) {
  ok <- is.finite(score) & !is.na(y01)
  score <- score[ok]; y01 <- y01[ok]
  n1 <- sum(y01 == 1); n0 <- sum(y01 == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[y01 == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## Choose a natural-spline df that is safe for the sample at hand. ns(df=k) needs
## at least k+1 distinct predictor values and comfortably more rows than df; we
## step down rather than error out on thin strata.
.safe_df <- function(x, df) {
  x <- x[is.finite(x)]
  nu <- length(unique(x)); n <- length(x)
  cap <- max(1L, min(df, nu - 1L, floor(n / 8)))
  as.integer(min(df, cap))
}

## Detect the subject-id column of a serology feature matrix W.
.detect_id <- function(W, ref_ids = NULL) {
  if ("subject_accession" %in% names(W)) return("subject_accession")
  if (!is.null(ref_ids)) {
    hit <- names(W)[vapply(W, function(z) any(as.character(z) %in% as.character(ref_ids)),
                           logical(1))]
    if (length(hit)) return(hit[1])
  }
  nonnum <- names(W)[!vapply(W, is.numeric, logical(1))]
  if (length(nonnum)) return(nonnum[1])
  names(W)[1]
}

## ---- 0. interval construction ------------------------------------------------

## Attach w_weeks to a clinical frame. Keeps days as a secondary column.
add_w_weeks <- function(clin, days_col = "vaccine_birth_interval_days") {
  if (!days_col %in% names(clin))
    stop("add_w_weeks(): column '", days_col, "' not found in the clinical frame.")
  clin$w_days  <- suppressWarnings(as.numeric(clin[[days_col]]))
  clin$w_weeks <- clin$w_days / 7
  clin
}

## Build the per-subject interval table used everywhere downstream:
## subject_accession, arm_name (factor TT < TdaP), w_weeks (+ w_days).
## `arm_tbl` and `clin_w` are joined; rows missing arm or w_weeks are dropped
## only where noted by the caller.
make_interval_tbl <- function(arm_tbl, clin_w, maternal_arms = c("TT", "TdaP")) {
  keep <- c("subject_accession", "w_weeks", "w_days")
  keep <- intersect(keep, names(clin_w))
  it <- merge(arm_tbl, clin_w[, keep, drop = FALSE], by = "subject_accession",
              all.x = TRUE)
  it <- it[it$arm_name %in% maternal_arms, , drop = FALSE]
  it$arm_name <- factor(as.character(it$arm_name), levels = maternal_arms)
  it
}

## ---- 1. per-feature boost model: MatBirth ~ ns(w)*arm + PregEarly -------------

## Fit the interval model for ONE feature's paired frame.
## pdf must contain: v_from (PregEarly log10), v_to (MatBirth log10),
##                   w_weeks, arm_name (factor), subject_accession.
## Returns a list with the full and reduced fits, an anova for the whole w block
## and for the w x arm interaction, and within-arm partial R^2 for w.
fit_boost_spline <- function(pdf, spline_df = 3, min_n_arm = 12L) {
  d <- pdf[is.finite(pdf$v_from) & is.finite(pdf$v_to) & is.finite(pdf$w_weeks), ,
           drop = FALSE]
  d$arm_name <- droplevels(factor(d$arm_name))
  arms <- levels(d$arm_name)
  n_by_arm <- table(d$arm_name)
  out <- list(ok = FALSE, n = nrow(d), n_by_arm = n_by_arm, arms = arms,
              df_used = NA_integer_)
  if (nrow(d) < (min_n_arm + 4L) || length(arms) < 1L) return(out)

  df_use <- .safe_df(d$w_weeks, spline_df)
  out$df_used <- df_use
  ns_ok <- df_use >= 1L
  bx <- if (ns_ok) "ns(w_weeks, df = df_use)" else "w_weeks"

  # Full model: arm-specific spline + baseline covariate.
  f_full <- stats::as.formula(sprintf("v_to ~ %s * arm_name + v_from", bx))
  # Reduced (drop all w terms): baseline + arm only.
  f_red  <- stats::as.formula("v_to ~ arm_name + v_from")
  # Additive (w main effect, no interaction): to isolate the w x arm interaction.
  f_add  <- stats::as.formula(sprintf("v_to ~ %s + arm_name + v_from", bx))

  environment(f_full) <- environment(f_add) <- environment()
  m_full <- tryCatch(stats::lm(f_full, data = d), error = function(e) NULL)
  m_red  <- tryCatch(stats::lm(f_red,  data = d), error = function(e) NULL)
  m_add  <- tryCatch(stats::lm(f_add,  data = d), error = function(e) NULL)
  if (is.null(m_full) || is.null(m_red)) return(out)

  # LR/F test for the whole w block (full vs reduced) and for the interaction
  # (full vs additive), computed only where both arms are present.
  p_w   <- tryCatch(stats::anova(m_red, m_full)$`Pr(>F)`[2],  error = function(e) NA_real_)
  p_int <- if (!is.null(m_add) && length(arms) == 2L)
    tryCatch(stats::anova(m_add, m_full)$`Pr(>F)`[2], error = function(e) NA_real_)
    else NA_real_

  # Within-arm partial R^2 for w: fraction of residual variance (after baseline)
  # explained by adding the spline, computed separately in each arm.
  partial_r2 <- setNames(rep(NA_real_, length(arms)), arms)
  n_arm      <- setNames(as.integer(n_by_arm[arms]), arms)
  for (a in arms) {
    da <- d[d$arm_name == a, , drop = FALSE]
    if (nrow(da) < (min_n_arm + 2L)) next
    dfa <- .safe_df(da$w_weeks, spline_df)
    bxa <- if (dfa >= 1L) sprintf("ns(w_weeks, df = %d)", dfa) else "w_weeks"
    fa_full <- stats::as.formula(sprintf("v_to ~ %s + v_from", bxa))
    fa_base <- stats::as.formula("v_to ~ v_from")
    ma_full <- tryCatch(stats::lm(fa_full, data = da), error = function(e) NULL)
    ma_base <- tryCatch(stats::lm(fa_base, data = da), error = function(e) NULL)
    if (is.null(ma_full) || is.null(ma_base)) next
    rss_full <- sum(stats::residuals(ma_full)^2)
    rss_base <- sum(stats::residuals(ma_base)^2)
    partial_r2[a] <- if (rss_base > 0) (rss_base - rss_full) / rss_base else NA_real_
  }

  out$ok <- TRUE
  out$fit_full <- m_full; out$fit_red <- m_red; out$fit_add <- m_add
  out$p_w <- p_w; out$p_int <- p_int
  out$partial_r2 <- partial_r2; out$n_arm <- n_arm
  out$coef_vfrom <- unname(stats::coef(m_full)["v_from"])
  out$ref_from   <- stats::median(d$v_from, na.rm = TRUE)
  out$data <- d
  out
}

## Grid of fitted curves (at a reference PregEarly baseline) with 95% CI, per arm.
spline_curve_df <- function(fit, n_grid = 80L) {
  if (!isTRUE(fit$ok)) return(NULL)
  d <- fit$data
  ref_from <- fit$ref_from
  out <- lapply(levels(d$arm_name), function(a) {
    wa <- d$w_weeks[d$arm_name == a]
    if (length(wa) < 3L) return(NULL)
    g <- data.frame(
      w_weeks  = seq(min(wa), max(wa), length.out = n_grid),
      arm_name = factor(a, levels = levels(d$arm_name)),
      v_from   = ref_from)
    pr <- tryCatch(stats::predict(fit$fit_full, newdata = g, se.fit = TRUE),
                   error = function(e) NULL)
    if (is.null(pr)) return(NULL)
    g$fit <- pr$fit; g$lwr <- pr$fit - 1.96 * pr$se.fit; g$upr <- pr$fit + 1.96 * pr$se.fit
    g
  })
  do.call(rbind, out)
}

## Per-subject partial-residual points, adjusted to the reference baseline so
## they overlay the fitted curve honestly: adj = v_to - b_vfrom*(v_from - ref).
spline_points_df <- function(fit) {
  if (!isTRUE(fit$ok)) return(NULL)
  d <- fit$data
  b <- fit$coef_vfrom; if (!is.finite(b)) b <- 0
  d$adj <- d$v_to - b * (d$v_from - fit$ref_from)
  d[, c("subject_accession", "w_weeks", "arm_name", "adj")]
}

## Stage-1 headline table across features (one row per feature).
boost_spline_table <- function(paired, feats, spline_df = 3) {
  rows <- lapply(feats, function(ft) {
    pdf <- paired[paired$feature == ft, , drop = FALSE]
    if (!nrow(pdf)) return(NULL)
    fit <- fit_boost_spline(pdf, spline_df = spline_df)
    if (!isTRUE(fit$ok))
      return(data.frame(feature = ft, n = fit$n, df_used = fit$df_used,
                        n_TdaP = NA_integer_, n_TT = NA_integer_,
                        pR2_w_TdaP = NA_real_, pR2_w_TT = NA_real_,
                        p_w = NA_real_, p_interaction = NA_real_,
                        stringsAsFactors = FALSE))
    pr <- fit$partial_r2; na <- fit$n_arm
    data.frame(
      feature       = ft,
      n             = fit$n,
      df_used       = fit$df_used,
      n_TdaP        = if ("TdaP" %in% names(na)) na[["TdaP"]] else NA_integer_,
      n_TT          = if ("TT"   %in% names(na)) na[["TT"]]   else NA_integer_,
      pR2_w_TdaP    = if ("TdaP" %in% names(pr)) round(pr[["TdaP"]], 3) else NA_real_,
      pR2_w_TT      = if ("TT"   %in% names(pr)) round(pr[["TT"]], 3)   else NA_real_,
      p_w           = signif(fit$p_w, 3),
      p_interaction = signif(fit$p_int, 3),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## Supplemental multi-panel figure: fitted rise-and-decay curve (+95% CI) over
## partial-residual data points, one panel per feature. Returns a patchwork/ggplot
## object (or a single ggplot if patchwork is unavailable).
spline_fit_panel <- function(paired, feats, spline_df = 3,
                             arm_cols = c(TdaP = "#8A1538", TT = "#2C7A7B"),
                             ncol = 3L, base_size = 10) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  panels <- list()
  for (ft in feats) {
    pdf <- paired[paired$feature == ft, , drop = FALSE]
    fit <- if (nrow(pdf)) fit_boost_spline(pdf, spline_df = spline_df) else list(ok = FALSE)
    if (!isTRUE(fit$ok)) next
    cd <- spline_curve_df(fit); pt <- spline_points_df(fit)
    if (is.null(cd) || is.null(pt)) next
    g <- ggplot2::ggplot() +
      ggplot2::geom_point(data = pt,
        ggplot2::aes(w_weeks, adj, colour = arm_name), alpha = 0.35, size = 1.1) +
      ggplot2::geom_ribbon(data = cd,
        ggplot2::aes(w_weeks, ymin = lwr, ymax = upr, fill = arm_name),
        alpha = 0.15, colour = NA) +
      ggplot2::geom_line(data = cd,
        ggplot2::aes(w_weeks, fit, colour = arm_name), linewidth = 0.9) +
      ggplot2::scale_colour_manual(values = arm_cols, drop = FALSE, name = "Arm") +
      ggplot2::scale_fill_manual(values = arm_cols, drop = FALSE, guide = "none") +
      ggplot2::labs(title = ft, x = "weeks vaccination \u2192 delivery",
                    y = expression(log[10]~"(MatBirth, baseline-adj.)")) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = base_size))
    panels[[ft]] <- g
  }
  if (!length(panels)) return(NULL)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(panels, ncol = ncol, guides = "collect") &
      ggplot2::theme(legend.position = "bottom")
  } else {
    panels[[1]]
  }
}

## ---- 2. residualise the score-defining MatBirth features on ns(w)*arm --------

## Return a COPY of the long serology object in which, for each responder-defining
## feature at `to_visit` (MatBirth), log_assay_value has been replaced by
## residuals of  log_assay_value ~ ns(w_weeks, df) * arm_name  (mean added back to
## preserve scale). PregEarly rows and every other feature are left untouched.
## Only rows with non-missing w_weeks and arm are adjusted; others are returned
## unchanged (callers restrict to interval-available mothers for a clean compare).
residualise_matbirth <- function(data_raw, resp_feats, to_visit, interval_tbl,
                                  spline_df = 3,
                                  value_col = "log_assay_value",
                                  feature_col = "feature",
                                  visit_col = "visit_name",
                                  id_col = "subject_accession") {
  stopifnot(all(c(value_col, feature_col, visit_col, id_col) %in% names(data_raw)))
  it <- interval_tbl[, c("subject_accession", "w_weeks", "arm_name")]
  d <- data_raw
  for (ft in resp_feats) {
    idx <- which(as.character(d[[feature_col]]) == ft &
                 as.character(d[[visit_col]]) == to_visit)
    if (!length(idx)) next
    sub <- data.frame(subject_accession = d[[id_col]][idx],
                      y = d[[value_col]][idx], row = idx,
                      stringsAsFactors = FALSE)
    sub <- merge(sub, it, by = "subject_accession", all.x = TRUE)
    fit_rows <- with(sub, is.finite(y) & is.finite(w_weeks) & !is.na(arm_name))
    if (sum(fit_rows) < 16L) next
    fdat <- sub[fit_rows, , drop = FALSE]
    fdat$arm_name <- droplevels(factor(fdat$arm_name))
    dfa <- .safe_df(fdat$w_weeks, spline_df)
    bx  <- if (dfa >= 1L) sprintf("ns(w_weeks, df = %d)", dfa) else "w_weeks"
    frm <- if (nlevels(fdat$arm_name) >= 2L)
             stats::as.formula(sprintf("y ~ %s * arm_name", bx))
           else stats::as.formula(sprintf("y ~ %s", bx))
    m <- tryCatch(stats::lm(frm, data = fdat), error = function(e) NULL)
    if (is.null(m)) next
    # Remove ONLY the within-arm interval deviation, preserving each arm's mean
    # titre (the vaccine effect is biology, not an interval artifact, so it must
    # survive). Standardise every mother to her arm's mean interval: subtract
    # f(w_i, arm_i) - f(mean_w_in_arm, arm_i). NB: regressing out ns(w)*arm and
    # adding back the GLOBAL mean would instead delete the between-arm difference
    # and spuriously collapse the arm signal.
    fit_all <- stats::predict(m)
    ref <- fdat
    ref$w_weeks <- stats::ave(fdat$w_weeks, fdat$arm_name,
                              FUN = function(z) mean(z, na.rm = TRUE))
    fit_ref <- tryCatch(stats::predict(m, newdata = ref),
                        error = function(e) fit_all)
    adj <- fdat$y - (fit_all - fit_ref)
    d[[value_col]][fdat$row] <- adj    # write interval-standardised values back
  }
  d
}

## ---- 3. forced-unpenalised elastic net with the interval block ---------------

## Cross-validated logistic AUC via k-fold CV of a plain glm (used for the
## no-serology reference rows: arm-only, and arm + ns(w)*arm).
cv_auc_glm <- function(formula, data, positive = "High", nfolds = 10, seed = 1) {
  mf <- tryCatch(stats::model.frame(formula, data = data,
                                    na.action = stats::na.omit),
                 error = function(e) NULL)
  if (is.null(mf) || nrow(mf) < (nfolds + 5L)) return(NA_real_)
  y  <- mf[[1]]
  y01 <- as.numeric(as.character(y) == positive)
  set.seed(seed)
  fold <- sample(rep_len(seq_len(nfolds), nrow(mf)))
  pr <- rep(NA_real_, nrow(mf))
  for (k in seq_len(nfolds)) {
    tr <- fold != k; te <- fold == k
    fit <- tryCatch(stats::glm(formula, data = mf[tr, , drop = FALSE],
                               family = stats::binomial()),
                    error = function(e) NULL)
    if (is.null(fit)) next
    pr[te] <- tryCatch(stats::predict(fit, newdata = mf[te, , drop = FALSE],
                                      type = "response"),
                       error = function(e) NA_real_)
  }
  .auc01(pr, y01)
}

## Elastic net where the vaccine arm (always) and, optionally, the ns(w)*arm
## interval block are forced in UNPENALISED and the serology features compete for
## the remaining signal. Mirrors the contract of elastic_net_signature() but adds
## the interval block; implemented directly on glmnet so the existing helper is
## not touched.
##
## W          serology feature matrix (subject id column + numeric feature cols)
## mat_cls    data.frame(subject_accession, responder [Low/High], score)
## interval_tbl  data.frame(subject_accession, arm_name, w_weeks)
## force_interval  TRUE -> force arm + ns(w)*arm ; FALSE -> force arm only
en_forced <- function(W, mat_cls, interval_tbl, force_interval = TRUE,
                      spline_df = 3, alpha = 0.5, nfolds = 10, seed = 1,
                      positive = "High") {
  out <- list(ok = FALSE, cv_auc = NA_real_, n_serology_selected = NA_integer_,
              coefs = NULL, force_interval = force_interval)
  if (!requireNamespace("glmnet", quietly = TRUE)) return(out)

  idc <- .detect_id(W, ref_ids = mat_cls$subject_accession)
  feat_cols <- setdiff(names(W), idc)
  feat_cols <- feat_cols[vapply(W[feat_cols], is.numeric, logical(1))]
  Wm <- data.frame(subject_accession = as.character(W[[idc]]),
                   W[, feat_cols, drop = FALSE], check.names = FALSE,
                   stringsAsFactors = FALSE)

  cl <- data.frame(subject_accession = as.character(mat_cls$subject_accession),
                   responder = as.character(mat_cls$responder),
                   stringsAsFactors = FALSE)
  it <- data.frame(subject_accession = as.character(interval_tbl$subject_accession),
                   arm_name = as.character(interval_tbl$arm_name),
                   w_weeks = interval_tbl$w_weeks, stringsAsFactors = FALSE)

  dat <- merge(merge(cl, it, by = "subject_accession"), Wm, by = "subject_accession")
  dat <- dat[stats::complete.cases(dat[, c("responder", "arm_name", "w_weeks")]), ,
             drop = FALSE]
  # keep serology features that are fully observed on this aligned set
  fc <- feat_cols[colSums(is.na(dat[, feat_cols, drop = FALSE])) == 0]
  if (nrow(dat) < (nfolds + 10L) || !length(fc)) return(out)

  y01 <- as.numeric(dat$responder == positive)
  arm_tdap <- as.numeric(dat$arm_name == "TdaP")

  # forced (unpenalised) block
  if (force_interval) {
    dfa <- .safe_df(dat$w_weeks, spline_df)
    B <- if (dfa >= 1L) splines::ns(dat$w_weeks, df = dfa) else matrix(dat$w_weeks, ncol = 1)
    B <- as.matrix(B); colnames(B) <- paste0("w_ns", seq_len(ncol(B)))
    Bint <- B * arm_tdap; colnames(Bint) <- paste0(colnames(B), ":TdaP")
    forced <- cbind(armTdaP = arm_tdap, B, Bint)
  } else {
    forced <- cbind(armTdaP = arm_tdap)
  }
  Xser <- as.matrix(dat[, fc, drop = FALSE])
  X <- cbind(forced, Xser)
  pen <- c(rep(0, ncol(forced)), rep(1, ncol(Xser)))

  set.seed(seed)
  foldid <- sample(rep_len(seq_len(nfolds), nrow(X)))
  cvfit <- tryCatch(
    glmnet::cv.glmnet(X, y01, family = "binomial", alpha = alpha,
                      penalty.factor = pen, foldid = foldid,
                      type.measure = "auc", standardize = TRUE),
    error = function(e) NULL)
  if (is.null(cvfit)) return(out)

  i_min  <- which(cvfit$lambda == cvfit$lambda.min)
  cv_auc <- cvfit$cvm[i_min]
  co <- as.matrix(stats::coef(cvfit, s = "lambda.min"))
  co <- data.frame(term = rownames(co), coef = co[, 1], stringsAsFactors = FALSE)
  ser_terms <- co$term %in% fc
  n_sel <- sum(ser_terms & co$coef != 0)

  out$ok <- TRUE; out$cv_auc <- as.numeric(cv_auc)
  out$n_serology_selected <- as.integer(n_sel)
  out$coefs <- co[co$coef != 0 & co$term != "(Intercept)", , drop = FALSE]
  out$n <- nrow(X); out$n_serology_features <- length(fc)
  out
}
