#' Variance attribution for the maternal arm in post-primary outcome models
#' ---------------------------------------------------------------------------
#' Phase 1 of analysis_plan_post_primary (Q3): "How much of the variance in the
#' 5-month whole-bacterial assays is attributable to maternal TdaP vs TT?"
#'
#' GMRs are effect sizes (ratios) and do NOT answer a variance-share question;
#' this helper supplies the decomposition the GMR battery cannot. For a fitted
#' outcome ~ maternal_arm * infant_arm (+ baseline) model it returns, for each
#' model term:
#'   - partial eta^2   = SS_term / (SS_term + SS_resid)   (Type-II SS)
#'   - LMG / Shapley share of R^2 (average over predictor orderings) for the
#'     three main-effect groups {maternal_arm, infant_arm, baseline}, with the
#'     maternal x infant interaction reported as its incremental R^2 over the
#'     additive model (this side-steps the marginality problem that makes a
#'     naive LMG-with-interaction ill-defined).
#' A subject-level nonparametric bootstrap gives CIs on every share.
#'
#' No hard dependency on {relaimpo}: LMG is computed directly. {car} is used for
#' Type-II SS if present, otherwise a sequential-SS fallback is used.
#' ---------------------------------------------------------------------------

.have_car <- requireNamespace("car", quietly = TRUE)

## ---- R^2 of an lm fit (fraction of variance explained) --------------------
.r2 <- function(fit) {
  s <- summary(fit)
  if (is.null(s$r.squared) || is.na(s$r.squared)) return(NA_real_)
  s$r.squared
}

## ---- Build the design data for one (feature, visit, standstat) cell -------
## Expects columns: subject_accession, maternal_arm, infant_arm, response, and
## (optionally) baseline. Drops incomplete rows and unused factor levels.
prep_cell <- function(df, response = "response", use_baseline = TRUE) {
  keep <- c("subject_accession", "maternal_arm", "infant_arm", response,
            if (use_baseline) "baseline")
  d <- df[, intersect(keep, names(df)), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d$maternal_arm <- droplevels(factor(d$maternal_arm))
  d$infant_arm   <- droplevels(factor(d$infant_arm))
  d
}

## ---- Partial eta^2 (Type-II SS) for each term -----------------------------
partial_eta2 <- function(fit) {
  if (.have_car) {
    aa <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
    if (!is.null(aa)) {
      ss   <- aa[["Sum Sq"]]
      rn   <- rownames(aa)
      ssr  <- ss[rn == "Residuals"]
      keep <- rn != "Residuals"
      return(stats::setNames(ss[keep] / (ss[keep] + ssr), rn[keep]))
    }
  }
  ## fallback: sequential anova() SS (order-dependent; flagged in output)
  aa  <- stats::anova(fit)
  ss  <- aa[["Sum Sq"]]; rn <- rownames(aa)
  ssr <- ss[rn == "Residuals"]; keep <- rn != "Residuals"
  stats::setNames(ss[keep] / (ss[keep] + ssr), rn[keep])
}

## ---- LMG / Shapley shares over the additive main-effect groups ------------
## groups: character vector of predictor names that each form a "group".
## Returns named vector of R^2 contributions summing to the additive-model R^2.
lmg_shares <- function(data, response, groups) {
  groups <- groups[vapply(groups, function(g) g %in% names(data), logical(1))]
  k <- length(groups)
  if (k == 0) return(numeric(0))
  perms <- if (k == 1) matrix(1, 1, 1) else {
    do.call(rbind, combinat_perm(seq_len(k)))
  }
  contr <- stats::setNames(numeric(k), groups)
  r2_of <- function(gset) {
    if (!length(gset)) return(0)
    f <- stats::reformulate(groups[gset], response = response)
    .r2(stats::lm(f, data = data))
  }
  for (p in seq_len(nrow(perms))) {
    ord <- perms[p, ]
    prev <- 0; inset <- integer(0)
    for (pos in ord) {
      inset <- c(inset, pos)
      now <- r2_of(inset)
      contr[pos] <- contr[pos] + (now - prev)
      prev <- now
    }
  }
  contr / nrow(perms)
}

## small permutation generator (avoids a {combinat} dependency)
combinat_perm <- function(x) {
  if (length(x) == 1) return(list(x))
  out <- list()
  for (i in seq_along(x))
    for (p in combinat_perm(x[-i])) out[[length(out) + 1]] <- c(x[i], p)
  out
}

## ---- One full decomposition for a prepared cell ---------------------------
## Returns a one-row data.frame of shares (proportion-of-variance, R^2 units)
## and the same expressed as % of the full-model R^2.
arm_variance_decomp <- function(d, response = "response",
                                with_interaction = TRUE, use_baseline = TRUE) {
  if (nlevels(d$maternal_arm) < 2)
    return(NULL)                                   # no maternal contrast in cell
  mains <- c("maternal_arm",
             if (nlevels(d$infant_arm) > 1) "infant_arm",
             if (use_baseline && "baseline" %in% names(d)) "baseline")
  add_f  <- stats::reformulate(mains, response = response)
  add_fit <- stats::lm(add_f, data = d)
  r2_add  <- .r2(add_fit)

  full_fit <- add_fit; r2_full <- r2_add; int_inc <- 0
  if (with_interaction && "infant_arm" %in% mains) {
    full_f  <- stats::reformulate(c(mains, "maternal_arm:infant_arm"),
                                  response = response)
    full_fit <- stats::lm(full_f, data = d)
    r2_full  <- .r2(full_fit)
    int_inc  <- max(r2_full - r2_add, 0)
  }

  lmg <- lmg_shares(d, response, mains)            # sums to r2_add
  pe  <- partial_eta2(full_fit)

  share <- function(nm) if (nm %in% names(lmg)) unname(lmg[nm]) else NA_real_
  data.frame(
    n                 = nrow(d),
    r2_full           = r2_full,
    r2_additive       = r2_add,
    lmg_maternal      = share("maternal_arm"),
    lmg_infant        = share("infant_arm"),
    lmg_baseline      = share("baseline"),
    interaction_R2    = int_inc,
    maternal_footprint = sum(share("maternal_arm"), int_inc, na.rm = TRUE),
    peta2_maternal    = unname(if ("maternal_arm" %in% names(pe)) pe["maternal_arm"] else NA),
    peta2_interaction = unname(if ("maternal_arm:infant_arm" %in% names(pe)) pe["maternal_arm:infant_arm"] else NA),
    ## shares as % of full-model R^2 (handy for reporting)
    pct_maternal      = 100 * share("maternal_arm")      / r2_full,
    pct_infant        = 100 * share("infant_arm")        / r2_full,
    pct_baseline      = 100 * share("baseline")          / r2_full,
    pct_interaction   = 100 * int_inc                    / r2_full,
    car_used          = .have_car,
    stringsAsFactors  = FALSE
  )
}

## ---- Subject bootstrap of the shares --------------------------------------
## Resamples whole subjects (the clustering unit). Returns the point estimate
## plus percentile CIs for the maternal-arm LMG share, the maternal footprint
## (main + interaction), partial-eta^2 of the maternal main effect, and R^2_full.
boot_arm_variance <- function(d, response = "response", with_interaction = TRUE,
                              use_baseline = TRUE, R = 1000, seed = 1L,
                              conf = 0.95) {
  point <- arm_variance_decomp(d, response, with_interaction, use_baseline)
  if (is.null(point)) return(NULL)
  ids <- unique(d$subject_accession)
  set.seed(seed)
  grab <- c("lmg_maternal", "maternal_footprint", "peta2_maternal", "r2_full",
            "interaction_R2")
  mat <- matrix(NA_real_, nrow = R, ncol = length(grab),
                dimnames = list(NULL, grab))
  ## fast subject resample: precompute row indices per subject once, then index
  ## (the decomposition does not use subject_accession, so no rename is needed).
  idx <- split(seq_len(nrow(d)), d$subject_accession); ids <- names(idx)
  for (b in seq_len(R)) {
    rows <- unlist(idx[sample(ids, length(ids), replace = TRUE)], use.names = FALSE)
    db <- d[rows, , drop = FALSE]
    db$maternal_arm <- droplevels(factor(db$maternal_arm))
    db$infant_arm   <- droplevels(factor(db$infant_arm))
    est <- tryCatch(
      arm_variance_decomp(db, response, with_interaction, use_baseline),
      error = function(e) NULL)
    if (!is.null(est)) mat[b, ] <- as.numeric(est[1, grab])
  }
  a  <- (1 - conf) / 2
  ci <- apply(mat, 2, stats::quantile, probs = c(a, 1 - a), na.rm = TRUE)
  beff <- sum(stats::complete.cases(mat))
  out <- point
  for (g in grab) {
    out[[paste0(g, "_lo")]] <- ci[1, g]
    out[[paste0(g, "_hi")]] <- ci[2, g]
  }
  out$boot_R_effective <- beff
  out
}
