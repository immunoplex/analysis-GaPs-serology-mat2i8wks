#' Blunting models for the infant primary response (Phase 2, Q2 inferential)
#' ---------------------------------------------------------------------------
#' analysis_plan_post_primary Q2: "Is there a level of the 2-month assay above
#' which the infant's response to primary vaccination is lower than otherwise
#' expected (maternal-antibody interference / blunting)?"
#'
#' The GMR battery (Phase 0) located the threshold descriptively; Phase 1
#' quantified the maternal variance share. Phase 2 supplies the *inferential*
#' test of blunting, which needs an EXPECTED-RESPONSE reference:
#'
#'   build_response_frame()        2->5-month response per subject x feature
#'   fit_blunting_interaction()    baseline x arm interaction (the formal test)
#'   fit_changepoint_by_arm()      segmented + Davies on the response (per arm)
#'   observed_vs_expected()        high-baseline deficit vs the TT counterfactual
#'
#' "Blunting" = a NEGATIVE dependence of the response on the 2-month baseline
#' that is steeper when maternal antibody is present (TdaP). The interaction
#' coefficient baseline:maternal_armTdaP < 0 is the formal signature.
#'
#' No hard external dependency: {segmented} is used for the changepoint if
#' installed, otherwise that engine degrades to NA with a flag.
#' ---------------------------------------------------------------------------

.have_segmented <- requireNamespace("segmented", quietly = TRUE)

## ---- response frame: baseline (2mo), level (post), fold-rise --------------
## dxdata must carry: standstat, subject_accession, arm_name, feature,
## visit_name, log_assay_value (the working scale; residual for standardised).
build_response_frame <- function(dxdata,
                                 features,
                                 baseline_visit = "M02",
                                 post_visit      = "M05") {
  d <- dxdata[as.character(dxdata$feature) %in% features &
              as.character(dxdata$visit_name) %in% c(baseline_visit, post_visit),
              c("standstat","subject_accession","arm_name","feature",
                "visit_name","log_assay_value")]
  d$visit_name <- as.character(d$visit_name)
  id <- c("standstat","subject_accession","arm_name","feature")
  base <- d[d$visit_name == baseline_visit, c(id, "log_assay_value")]
  post <- d[d$visit_name == post_visit,     c(id, "log_assay_value")]
  names(base)[names(base) == "log_assay_value"] <- "baseline"
  names(post)[names(post) == "log_assay_value"] <- "level"
  ## one row per id (defensive against any residual duplicates)
  base <- base[!duplicated(base[, id]), , drop = FALSE]
  post <- post[!duplicated(post[, id]), , drop = FALSE]
  w <- merge(base, post, by = id)
  w <- w[stats::complete.cases(w[, c("baseline","level")]), , drop = FALSE]
  w$foldrise      <- w$level - w$baseline
  w$maternal_arm  <- factor(ifelse(w$arm_name %in% c("TdaP_wP","TdaP_aP"),"TdaP","TT"),
                            levels = c("TT","TdaP"))
  w$infant_arm    <- factor(ifelse(w$arm_name %in% c("TdaP_wP","TT_wP"),"wP","aP"),
                            levels = c("aP","wP"))
  w$antigen       <- sub("_.*$","", as.character(w$feature))
  rownames(w) <- NULL
  w
}

## ---- the formal blunting test: baseline x arm interaction -----------------
## d = one (standstat, feature, priming, visit) cell. response in {"foldrise",
## "level"}; "level" is the baseline-adjusted form (baseline enters as a term).
fit_blunting_interaction <- function(d, response = "foldrise", min_n = 30) {
  d <- d[stats::complete.cases(d[, c("baseline", response)]), , drop = FALSE]
  d$maternal_arm <- droplevels(factor(d$maternal_arm, levels = c("TT","TdaP")))
  if (nrow(d) < min_n || nlevels(d$maternal_arm) < 2) return(NULL)
  d$y <- d[[response]]
  fit <- stats::lm(y ~ baseline * maternal_arm, data = d)
  cf  <- summary(fit)$coefficients
  b_base <- if ("baseline" %in% rownames(cf)) cf["baseline", 1] else NA_real_
  irow   <- grep("^baseline:maternal_armTdaP$", rownames(cf))
  if (!length(irow)) return(NULL)
  b_int  <- cf[irow, 1]; se_int <- cf[irow, 2]; p_int <- cf[irow, 4]
  data.frame(
    n            = nrow(d),
    slope_TT     = b_base,                    # baseline-response slope, TT
    slope_TdaP   = b_base + b_int,            # baseline-response slope, TdaP
    slope_diff   = b_int,                     # TdaP - TT  (< 0 = extra blunting)
    diff_lo      = b_int - 1.96 * se_int,
    diff_hi      = b_int + 1.96 * se_int,
    p_interaction = p_int,
    response_type = response,
    stringsAsFactors = FALSE)
}

## ---- segmented (changepoint) on the response, per maternal arm ------------
fit_changepoint_by_arm <- function(d, response = "foldrise", min_n = 30) {
  out <- list()
  for (lvl in levels(droplevels(factor(d$maternal_arm)))) {
    dd <- d[d$maternal_arm == lvl &
            stats::complete.cases(d[, c("baseline", response)]), , drop = FALSE]
    if (nrow(dd) < min_n) next
    dd$y <- dd[[response]]
    lin <- stats::lm(y ~ baseline, data = dd)
    rec <- data.frame(maternal_arm = lvl, n = nrow(dd),
                      breakpoint = NA_real_, slope_below = NA_real_,
                      slope_above = NA_real_, davies_p = NA_real_,
                      converged = FALSE, stringsAsFactors = FALSE)
    if (.have_segmented) {
      dv <- tryCatch(segmented::davies.test(lin, ~ baseline),
                     error = function(e) NULL)
      if (!is.null(dv)) rec$davies_p <- unname(dv$p.value)
      sg <- tryCatch(segmented::segmented(lin, seg.Z = ~ baseline),
                     error = function(e) NULL)
      if (!is.null(sg) && !is.null(sg$psi)) {
        sl <- tryCatch(segmented::slope(sg)$baseline[, 1], error = function(e) NULL)
        if (!is.null(sl) && length(sl) >= 2) {
          rec$breakpoint  <- unname(sg$psi[1, "Est."])
          rec$slope_below <- unname(sl[1])
          rec$slope_above <- unname(sl[2])
          rec$converged   <- TRUE
        }
      }
    }
    out[[lvl]] <- rec
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

## ---- observed vs expected: the counterfactual blunting deficit ------------
## Fit the baseline->response relationship in the TT (low maternal-antibody)
## arm; use it to PREDICT each TdaP infant's expected response from its own
## baseline; report observed - expected. A negative deficit among high-baseline
## TdaP infants is "lower than expected". high_q: baseline quantile above which
## an infant counts as high-baseline (default top tertile of TdaP baselines).
observed_vs_expected <- function(d, response = "level", high_q = 2/3,
                                 min_n = 30, R = 1000, seed = 1L, conf = 0.95) {
  d <- d[stats::complete.cases(d[, c("baseline", response)]), , drop = FALSE]
  d$maternal_arm <- droplevels(factor(d$maternal_arm, levels = c("TT","TdaP")))
  tt   <- d[d$maternal_arm == "TT",   , drop = FALSE]
  tdap <- d[d$maternal_arm == "TdaP", , drop = FALSE]
  if (nrow(tt) < min_n || nrow(tdap) < min_n) return(NULL)
  cut <- stats::quantile(tdap$baseline, high_q, na.rm = TRUE)
  hi  <- tdap[tdap$baseline >= cut, , drop = FALSE]
  if (nrow(hi) < 5) return(NULL)

  deficit <- function(tt_df, hi_df) {
    fit <- stats::lm(stats::reformulate("baseline", response), data = tt_df)
    exp_hi <- stats::predict(fit, newdata = hi_df)
    mean(hi_df[[response]] - exp_hi)
  }
  point <- deficit(tt, hi)

  ## subject bootstrap (resample within TT and within high-baseline TdaP)
  set.seed(seed)
  bs <- replicate(R, {
    tb <- tt[sample(nrow(tt), replace = TRUE), , drop = FALSE]
    hb <- hi[sample(nrow(hi), replace = TRUE), , drop = FALSE]
    tryCatch(deficit(tb, hb), error = function(e) NA_real_)
  })
  a  <- (1 - conf) / 2
  ci <- stats::quantile(bs, c(a, 1 - a), na.rm = TRUE)
  data.frame(
    n_TT          = nrow(tt),
    n_TdaP_high   = nrow(hi),
    high_baseline_cut = unname(cut),
    deficit       = point,                   # observed - expected (log scale)
    deficit_lo    = unname(ci[1]),
    deficit_hi    = unname(ci[2]),
    fold_vs_expected = 10^point,             # multiplicative observed/expected
    boot_R_effective = sum(!is.na(bs)),
    response_type = response,
    stringsAsFactors = FALSE)
}
