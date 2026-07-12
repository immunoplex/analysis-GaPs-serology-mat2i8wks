#' Carry-forward prediction and mediation for post-primary function (Phase 3, Q1)
#' ---------------------------------------------------------------------------
#' analysis_plan_post_primary Q1 / Phase 3: does the 2-month subclass PROGRAMME
#' carry forward to predict the infant's post-primary Fc-receptor binding,
#' effector function and whole-bacterial responses, and how much of the
#' maternal-arm effect on function is MEDIATED by the 2-month composition?
#'
#'   fit_hier_blocks()       sequential block entry: arm -> +total IgG ->
#'                           +IgG1/composition -> +effector/FcR, with DR^2 and
#'                           nested-F increment tests (Step 1.4 decisive contrasts)
#'   arm_priming_interaction()  the maternal x infant interaction test
#'                           (differential carry-forward by priming type)
#'   mediate_simple()        product-of-coefficients mediation of the maternal-arm
#'                           effect through a 2-month mediator, subject-bootstrap CIs
#'                           (closes the "no mediation model" gap)
#'
#' All predictors are the PRIOR 2-month baseline, so the forward prediction is
#' non-circular by construction (the circularity lesson from summary_06).
#'
#' Dependency-light: base lm + anova. {mediation} is NOT required (manual
#' product-of-coefficients with a subject bootstrap is used instead).
#' ---------------------------------------------------------------------------

.r2  <- function(fit) { s <- summary(fit); if (is.null(s$r.squared)) NA_real_ else s$r.squared }

## ---- sequential hierarchical block entry ----------------------------------
## d        : data frame with the outcome `y`, the arm factors, and the named
##            baseline predictor columns referenced in `blocks`.
## arm_terms: character vector of model terms entered first (e.g.
##            c("maternal_arm","infant_arm","maternal_arm:infant_arm")).
## blocks   : NAMED list; each element a character vector of predictor columns
##            added cumulatively after arm_terms (e.g.
##            list(total="b_IGG", composition=c("b_IGG1","b_comp"),
##                 effector=c("b_FCGR2A","b_FCGR3B","b_ADCD"))).
## Returns one row per model ("arm", then each block) with R^2, incremental
## DR^2 over the previous model, and the nested-F increment p-value. All models
## are fit on the SAME complete-case rows so the nested F tests are valid.
fit_hier_blocks <- function(d, y, arm_terms, blocks, min_n = 30) {
  pred_vars <- unique(c(unlist(lapply(strsplit(arm_terms, ":"), identity)),
                        unlist(blocks)))
  pred_vars <- pred_vars[pred_vars %in% names(d)]
  cc <- stats::complete.cases(d[, c(y, pred_vars), drop = FALSE])
  d  <- d[cc, , drop = FALSE]
  for (v in c("maternal_arm","infant_arm"))
    if (v %in% names(d)) d[[v]] <- droplevels(factor(d[[v]]))
  if (nrow(d) < min_n) return(NULL)
  d$.y <- d[[y]]

  arm_terms <- arm_terms[vapply(arm_terms, function(t)
    all(strsplit(t, ":")[[1]] %in% names(d)) &&
    all(vapply(strsplit(t, ":")[[1]],
               function(v) !is.factor(d[[v]]) || nlevels(d[[v]]) > 1, logical(1))),
    logical(1))]

  rhs0 <- if (length(arm_terms)) arm_terms else "1"
  fits <- list(arm = stats::lm(stats::reformulate(rhs0, ".y"), data = d))
  rhs  <- rhs0
  for (bn in names(blocks)) {
    add <- blocks[[bn]][blocks[[bn]] %in% names(d)]
    add <- add[vapply(add, function(v) {
      x <- d[[v]]; sum(!is.na(x)) > 0 && stats::sd(x, na.rm = TRUE) > 0 }, logical(1))]
    if (!length(add)) { fits[[bn]] <- fits[[length(fits)]]; next }
    rhs <- c(rhs[rhs != "1"], add)
    fits[[bn]] <- stats::lm(stats::reformulate(rhs, ".y"), data = d)
  }

  out <- vector("list", length(fits)); nm <- names(fits)
  prev <- NULL
  for (i in seq_along(fits)) {
    r2 <- .r2(fits[[i]])
    dr2 <- if (i == 1) NA_real_ else r2 - .r2(fits[[i - 1]])
    pinc <- NA_real_
    if (i > 1) {
      av <- tryCatch(stats::anova(fits[[i - 1]], fits[[i]]), error = function(e) NULL)
      if (!is.null(av) && "Pr(>F)" %in% names(av)) pinc <- av[["Pr(>F)"]][2]
    }
    out[[i]] <- data.frame(model = nm[i], n = nrow(d), r2 = r2,
                           dR2 = dr2, p_increment = pinc,
                           stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

## ---- maternal x infant interaction (differential carry-forward) -----------
arm_priming_interaction <- function(d, y, predictors, min_n = 30) {
  pv <- unique(c("maternal_arm","infant_arm", predictors))
  pv <- pv[pv %in% names(d)]
  cc <- stats::complete.cases(d[, c(y, pv), drop = FALSE]); d <- d[cc, , drop = FALSE]
  d$maternal_arm <- droplevels(factor(d$maternal_arm))
  d$infant_arm   <- droplevels(factor(d$infant_arm))
  if (nrow(d) < min_n || nlevels(d$maternal_arm) < 2 || nlevels(d$infant_arm) < 2)
    return(NULL)
  d$.y <- d[[y]]
  base_terms <- c(predictors[predictors %in% names(d)], "maternal_arm", "infant_arm")
  red  <- stats::lm(stats::reformulate(base_terms, ".y"), data = d)
  full <- stats::lm(stats::reformulate(c(base_terms, "maternal_arm:infant_arm"), ".y"), data = d)
  av <- tryCatch(stats::anova(red, full), error = function(e) NULL)
  p  <- if (!is.null(av) && "Pr(>F)" %in% names(av)) av[["Pr(>F)"]][2] else NA_real_
  data.frame(n = nrow(d), dR2_interaction = .r2(full) - .r2(red),
             p_interaction = p, stringsAsFactors = FALSE)
}

## ---- mediation: arm -> 2-month mediator -> post-primary function ----------
## Product-of-coefficients with a subject bootstrap. exposure must be a 2-level
## factor with the reference first (TT) and the treated level (TdaP) second.
## covars: optional adjusters (e.g. "infant_arm" when pooled across priming).
mediate_simple <- function(d, y, mediator, exposure = "maternal_arm",
                           covars = character(0), id = "subject_accession",
                           min_n = 30, R = 1000, seed = 1L, conf = 0.95) {
  vars <- unique(c(y, mediator, exposure, covars, id))
  vars <- vars[vars %in% names(d)]
  cc <- stats::complete.cases(d[, setdiff(vars, id), drop = FALSE])
  d  <- d[cc, , drop = FALSE]
  d[[exposure]] <- droplevels(factor(d[[exposure]]))
  if (nrow(d) < min_n || nlevels(d[[exposure]]) < 2) return(NULL)
  exp_coef <- paste0(exposure, levels(d[[exposure]])[2])

  est <- function(dd) {
    am <- stats::lm(stats::reformulate(c(exposure, covars), mediator), data = dd)
    om <- stats::lm(stats::reformulate(c(mediator, exposure, covars), y), data = dd)
    tm <- stats::lm(stats::reformulate(c(exposure, covars), y), data = dd)
    cf_a <- stats::coef(am); cf_o <- stats::coef(om); cf_t <- stats::coef(tm)
    a <- if (exp_coef %in% names(cf_a)) cf_a[[exp_coef]] else NA_real_
    b <- if (mediator %in% names(cf_o)) cf_o[[mediator]] else NA_real_
    cprime <- if (exp_coef %in% names(cf_o)) cf_o[[exp_coef]] else NA_real_
    ctot   <- if (exp_coef %in% names(cf_t)) cf_t[[exp_coef]] else NA_real_
    ind <- a * b
    prop <- if (is.na(ctot) || abs(ctot) < 1e-8) NA_real_ else ind / ctot
    c(a = a, b = b, c_total = ctot, c_direct = cprime,
      indirect = ind, prop_mediated = prop)
  }
  point <- est(d)
  a2 <- (1 - conf) / 2
  ## R <= 0: instant dry-run (point estimates only, no CIs) to validate the grid.
  if (R <= 0)
    return(data.frame(
      n = nrow(d), a = point[["a"]], b = point[["b"]],
      c_total = point[["c_total"]], c_direct = point[["c_direct"]],
      indirect = point[["indirect"]], indirect_lo = NA_real_, indirect_hi = NA_real_,
      prop_mediated = point[["prop_mediated"]], prop_lo = NA_real_, prop_hi = NA_real_,
      c_total_lo = NA_real_, c_total_hi = NA_real_,
      boot_R_effective = 0L, mediator = mediator, stringsAsFactors = FALSE))
  ## fast subject bootstrap: precompute the row indices for each subject ONCE,
  ## then resample whole subjects by index (O(rows) per rep, vs the previous
  ## O(n^2) rbind-of-subsets that made this chunk run for hours).
  idx <- split(seq_len(nrow(d)), d[[id]]); ids <- names(idx)
  set.seed(seed)
  bs <- matrix(NA_real_, R, length(point), dimnames = list(NULL, names(point)))
  for (r in seq_len(R)) {
    rows <- unlist(idx[sample(ids, length(ids), replace = TRUE)], use.names = FALSE)
    db <- d[rows, , drop = FALSE]
    db[[exposure]] <- droplevels(factor(db[[exposure]]))
    if (nlevels(db[[exposure]]) < 2) next
    bs[r, ] <- tryCatch(est(db), error = function(e) rep(NA_real_, length(point)))
  }
  ci <- apply(bs, 2, stats::quantile, probs = c(a2, 1 - a2), na.rm = TRUE)
  data.frame(
    n = nrow(d),
    a = point[["a"]], b = point[["b"]],
    c_total = point[["c_total"]], c_direct = point[["c_direct"]],
    indirect = point[["indirect"]],
    indirect_lo = ci[1, "indirect"], indirect_hi = ci[2, "indirect"],
    prop_mediated = point[["prop_mediated"]],
    prop_lo = ci[1, "prop_mediated"], prop_hi = ci[2, "prop_mediated"],
    c_total_lo = ci[1, "c_total"], c_total_hi = ci[2, "c_total"],
    boot_R_effective = sum(stats::complete.cases(bs)),
    mediator = mediator, stringsAsFactors = FALSE)
}
