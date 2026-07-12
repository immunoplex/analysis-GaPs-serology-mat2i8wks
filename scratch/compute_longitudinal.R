#' Longitudinal synthesis: mixed model on the 2->5->9 response (Phase 5)
#' ---------------------------------------------------------------------------
#' analysis_plan_post_primary Phase 5: one model per outcome that unifies the
#' arm, baseline, priming and visit structure across the post-primary visits.
#' Per the chosen design the OUTCOME is the change from the 2-month baseline
#' (resp = level(visit) - level(M02)), modelled at the 5- and 9-month visits
#' with a subject random intercept and ADDITIVE fixed effects:
#'
#'   resp ~ maternal_arm + infant_arm + baseline + visit + (1 | subject)
#'
#' so the maternal-arm coefficient is the adjusted arm effect on the response
#' (averaged over visits), the baseline coefficient is the longitudinal blunting
#' axis (negative = higher 2-month antibody -> smaller response), the visit term
#' is the 9-vs-5 durability shift, and the random intercept absorbs the
#' within-infant correlation between the two change points.
#'
#'   fit_lmm()    primary additive model: tidy fixed effects (+ Wald CI, p),
#'                ICC, marginal/conditional R^2, AIC, n_obs/n_subj
#'   lrt_term()   likelihood-ratio test for an added fixed term (e.g. the
#'                maternal x priming or maternal x baseline interaction)
#'
#' {lme4} is required; {lmerTest} is used for Satterthwaite p-values if present
#' (otherwise a large-sample normal approximation is used); {MuMIn} is used for
#' R^2 if present (otherwise R^2 is returned NA).
#' ---------------------------------------------------------------------------

.has_lme4    <- requireNamespace("lme4", quietly = TRUE)
.has_lmtest  <- requireNamespace("lmerTest", quietly = TRUE)
.has_mumin   <- requireNamespace("MuMIn", quietly = TRUE)
.lmer <- function(...) if (.has_lmtest) lmerTest::lmer(...) else lme4::lmer(...)

## ---- ICC and R^2 from a fitted merMod -------------------------------------
.lmm_fitstats <- function(fit, id) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  v_sub <- vc$vcov[vc$grp == id & is.na(vc$var1) == FALSE][1]
  if (is.na(v_sub)) v_sub <- vc$vcov[vc$grp == id][1]
  v_res <- vc$vcov[vc$grp == "Residual"][1]
  icc   <- v_sub / (v_sub + v_res)
  r2m <- r2c <- NA_real_
  if (.has_mumin) {
    rr <- tryCatch(MuMIn::r.squaredGLMM(fit), error = function(e) NULL)
    if (!is.null(rr)) { r2m <- unname(rr[1, "R2m"]); r2c <- unname(rr[1, "R2c"]) }
  }
  data.frame(icc = icc, r2_marginal = r2m, r2_conditional = r2c,
             aic = stats::AIC(fit))
}

## ---- primary additive mixed model -----------------------------------------
fit_lmm <- function(d, response = "resp",
                    fixed = c("maternal_arm","infant_arm","baseline","visit"),
                    id = "subject_accession", min_n = 40) {
  if (!.has_lme4) stop("lme4 is required for the Phase-5 mixed model.")
  vars <- c(response, fixed, id); vars <- vars[vars %in% names(d)]
  d <- d[stats::complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
  for (v in c("maternal_arm","infant_arm","visit"))
    if (v %in% names(d)) d[[v]] <- droplevels(factor(d[[v]]))
  fixed <- fixed[vapply(fixed, function(v)
    v %in% names(d) && (!is.factor(d[[v]]) || nlevels(d[[v]]) > 1), logical(1))]
  if (length(unique(d[[id]])) < min_n || !length(fixed)) return(NULL)
  d$.resp <- d[[response]]
  form <- stats::reformulate(c(fixed, sprintf("(1|%s)", id)), ".resp")
  fit  <- tryCatch(.lmer(form, data = d, REML = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  cf <- stats::coef(summary(fit))
  est <- cf[, "Estimate"]; se <- cf[, "Std. Error"]
  p  <- if ("Pr(>|t|)" %in% colnames(cf)) cf[, "Pr(>|t|)"]
        else 2 * stats::pnorm(-abs(cf[, "t value"]))
  tidy <- data.frame(
    term = rownames(cf), estimate = unname(est), se = unname(se),
    ci_lo = unname(est - 1.96 * se), ci_hi = unname(est + 1.96 * se),
    p_value = unname(p), stringsAsFactors = FALSE)
  fs <- .lmm_fitstats(fit, id)
  list(tidy = tidy,
       fit  = cbind(n_obs = nrow(d), n_subj = length(unique(d[[id]])), fs))
}

## ---- LRT for an added fixed-effect term (REML must be FALSE) ---------------
lrt_term <- function(d, response = "resp",
                     base_fixed = c("maternal_arm","infant_arm","baseline","visit"),
                     extra, id = "subject_accession", min_n = 40) {
  if (!.has_lme4) return(NULL)
  comp <- unique(unlist(strsplit(extra, ":")))
  vars <- c(response, base_fixed, comp, id); vars <- vars[vars %in% names(d)]
  d <- d[stats::complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
  for (v in c("maternal_arm","infant_arm","visit"))
    if (v %in% names(d)) d[[v]] <- droplevels(factor(d[[v]]))
  if (length(unique(d[[id]])) < min_n) return(NULL)
  d$.resp <- d[[response]]
  bf <- base_fixed[base_fixed %in% names(d)]
  f0 <- stats::reformulate(c(bf, sprintf("(1|%s)", id)), ".resp")
  f1 <- stats::reformulate(c(bf, extra, sprintf("(1|%s)", id)), ".resp")
  m0 <- tryCatch(lme4::lmer(f0, data = d, REML = FALSE), error = function(e) NULL)
  m1 <- tryCatch(lme4::lmer(f1, data = d, REML = FALSE), error = function(e) NULL)
  if (is.null(m0) || is.null(m1)) return(NULL)
  av <- tryCatch(stats::anova(m0, m1), error = function(e) NULL)
  if (is.null(av)) return(NULL)
  data.frame(term = extra, n_subj = length(unique(d[[id]])),
             dAIC = av$AIC[2] - av$AIC[1],
             chisq = av$Chisq[2], df = av$Df[2],
             p_value = av[["Pr(>Chisq)"]][2], stringsAsFactors = FALSE)
}
