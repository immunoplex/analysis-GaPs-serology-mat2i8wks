#' Within-TdaP vaccination-timing dose-response (Phase 4, Q4)
#' ---------------------------------------------------------------------------
#' analysis_plan_post_primary Q4: does the TIMING of maternal vaccination change
#' its impact, and is it detectable at all? Timing varies only within the TdaP
#' arm, so this is a within-arm continuous-exposure analysis with real power
#' limits — hence the detectability statement is mandatory, not optional.
#'
#'   fit_timing_lm()     linear timing slope (+ CI, p) AND a natural-spline fit
#'                       with a non-linearity test (Step 4.2), plus the minimum
#'                       detectable slope at the target power (Step 4.3)
#'   timing_modifier()   timing x baseline / timing x priming interaction
#'                       (Step 4.4 — does timing move the blunting, not the mean)
#'
#' The detectability logic (Step 4.3): MDES = (z_{1-a/2} + z_{power}) * SE(slope).
#' Read it as "the sample can detect a per-unit timing slope of >= MDES at the
#' target power". A non-significant slope with a SMALL MDES is an informative
#' null ("no timing effect larger than X"); a non-significant slope with a LARGE
#' MDES is uninformative and must be labelled so.
#'
#' Dependency-light: base lm + splines::ns (both base R).
#' ---------------------------------------------------------------------------

.r2 <- function(fit) { s <- summary(fit); if (is.null(s$r.squared)) NA_real_ else s$r.squared }

## ---- timing dose-response for one cell ------------------------------------
## d        : data frame with outcome `y`, timing column `timing`, and optional
##            covariates. timing should already be on an interpretable unit
##            (e.g. weeks) so the slope reads per-unit.
## Returns the linear slope with CI and p, the residual SD, the spline
## non-linearity p, and the minimum detectable slope at `power`.
fit_timing_lm <- function(d, y = "y", timing = "timing", covars = character(0),
                          spline_df = 3, min_n = 25, power = 0.8, alpha = 0.05) {
  vars <- c(y, timing, covars); vars <- vars[vars %in% names(d)]
  d <- d[stats::complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
  for (v in covars) if (v %in% names(d) && is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
  covars <- covars[vapply(covars, function(v)
    v %in% names(d) && (!is.factor(d[[v]]) || nlevels(d[[v]]) > 1), logical(1))]
  if (nrow(d) < min_n || stats::sd(d[[timing]], na.rm = TRUE) == 0) return(NULL)
  d$.y <- d[[y]]; d$.t <- d[[timing]]

  lin <- stats::lm(stats::reformulate(c(".t", covars), ".y"), data = d)
  cf  <- summary(lin)$coefficients
  if (!(".t" %in% rownames(cf))) return(NULL)
  b <- cf[".t", 1]; se <- cf[".t", 2]; p <- cf[".t", 4]
  zt <- stats::qnorm(1 - alpha / 2) + stats::qnorm(power)
  mdes <- zt * se

  nonlin_p <- NA_real_
  sp <- tryCatch(stats::lm(stats::reformulate(
           c(sprintf("splines::ns(.t, %d)", spline_df), covars), ".y"), data = d),
         error = function(e) NULL)
  if (!is.null(sp)) {
    av <- tryCatch(stats::anova(lin, sp), error = function(e) NULL)
    if (!is.null(av) && "Pr(>F)" %in% names(av)) nonlin_p <- av[["Pr(>F)"]][2]
  }

  ## detectability: `detected` = slope significant; the Rmd classifies a
  ## non-significant slope as an informative vs uninformative null by comparing
  ## `mdes` against a pre-specified meaningful slope.
  detected <- !is.na(p) && p < alpha
  data.frame(
    n = nrow(d), timing_sd = stats::sd(d[[timing]], na.rm = TRUE),
    slope = b, se = se, ci_lo = b - 1.96 * se, ci_hi = b + 1.96 * se,
    ci_halfwidth = 1.96 * se, p_value = p,
    resid_sd = stats::sigma(lin), mdes = mdes, nonlin_p = nonlin_p,
    detected = detected, stringsAsFactors = FALSE)
}

## ---- timing as an effect-modifier -----------------------------------------
## modifier may be continuous (e.g. baseline PT_IgG) or a factor (e.g. infant_arm).
timing_modifier <- function(d, y = "y", timing = "timing", modifier,
                            covars = character(0), min_n = 25) {
  vars <- c(y, timing, modifier, covars); vars <- vars[vars %in% names(d)]
  d <- d[stats::complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
  if (modifier %in% names(d) && is.factor(d[[modifier]]))
    d[[modifier]] <- droplevels(d[[modifier]])
  if (nrow(d) < min_n) return(NULL)
  if (is.factor(d[[modifier]]) && nlevels(d[[modifier]]) < 2) return(NULL)
  if (stats::sd(d[[timing]], na.rm = TRUE) == 0) return(NULL)
  d$.y <- d[[y]]; d$.t <- d[[timing]]; d$.m <- d[[modifier]]
  red  <- stats::lm(stats::reformulate(c(".t", ".m", covars), ".y"), data = d)
  full <- stats::lm(stats::reformulate(c(".t", ".m", ".t:.m", covars), ".y"), data = d)
  av <- tryCatch(stats::anova(red, full), error = function(e) NULL)
  p  <- if (!is.null(av) && "Pr(>F)" %in% names(av)) av[["Pr(>F)"]][2] else NA_real_
  data.frame(n = nrow(d), modifier = modifier,
             dR2_interaction = .r2(full) - .r2(red),
             p_interaction = p, stringsAsFactors = FALSE)
}
