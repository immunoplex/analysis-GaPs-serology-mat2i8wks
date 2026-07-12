# =============================================================================
# maternal_responder_helpers.R   (Phase 08 branch)
#
# Engine for: who are the LOW vs HIGH maternal responders, and is there a
# clinical covariate or maternal-serology signature that predicts the class?
#
# The maternal Low/High classification is the SAME construct as Phase 07 Stage A
# (Blom within-(feature,visit) z of the maternal binding features across the two
# maternal visits -> mean responder score -> PAM k=2). This file reuses the
# Phase 07 engine in R/within_visit_helpers.R when it is on the search path, and
# otherwise reconstructs the identical score with a self-contained fallback.
#
# All functions take inputs as explicit arguments. Stratification is on the
# MATERNAL arm (TdaP vs TT) only; there is deliberately no infant-arm logic here.
# =============================================================================


# ---- 0. maternal responder classification ----------------------------------

.blom <- function(x) {
  r <- rank(x, na.last = "keep"); n <- sum(!is.na(x))
  stats::qnorm((r - 3/8) / (n + 1/4))
}

#' Self-contained reconstruction of the maternal responder score + Low/High,
#' used only if the Phase 07 engine functions are not on the search path.
.fallback_responder <- function(data_raw, mat_vis, feats,
                                has_cluster = requireNamespace("cluster", quietly = TRUE)) {
  d <- data_raw[data_raw$visit_name %in% mat_vis &
                  as.character(data_raw$feature) %in% feats,
                c("subject_accession", "feature", "visit_name", "log_assay_value")]
  if (!nrow(d))
    return(data.frame(subject_accession = character(), responder = character(),
                      score = numeric(), method = character(),
                      stringsAsFactors = FALSE))
  d$key <- paste(as.character(d$feature), as.character(d$visit_name), sep = "@")
  agg <- stats::aggregate(log_assay_value ~ subject_accession + key, d,
                          FUN = function(z) mean(z, na.rm = TRUE))
  agg$z <- stats::ave(agg$log_assay_value, agg$key, FUN = .blom)
  sc <- stats::aggregate(z ~ subject_accession, agg,
                         FUN = function(z) mean(z, na.rm = TRUE))
  names(sc)[2] <- "score"
  lab <- rep(NA_character_, nrow(sc))
  if (has_cluster && nrow(sc) >= 20) {
    km <- tryCatch(cluster::pam(matrix(sc$score, ncol = 1), k = 2),
                   error = function(e) NULL)
    if (!is.null(km)) {
      hi <- which.max(tapply(sc$score, km$clustering, mean))
      lab <- ifelse(km$clustering == hi, "High", "Low")
    }
  }
  if (all(is.na(lab)))
    lab <- ifelse(sc$score >= stats::median(sc$score, na.rm = TRUE), "High", "Low")
  data.frame(subject_accession = sc$subject_accession, responder = lab,
             score = sc$score, method = "fallback", stringsAsFactors = FALSE)
}

#' Build the maternal Low/High responder table (one row per mother):
#' subject_accession, responder (factor Low<High), score (responder z), method.
build_maternal_responder_class <- function(data_raw, mat_vis, resp_feats,
                                            breadth_thresh = 0.6,
                                            has_cluster = requireNamespace("cluster", quietly = TRUE)) {
  avail <- sort(unique(as.character(data_raw$feature)))
  rf <- if (exists("resolve_feats", mode = "function"))
    resolve_feats(resp_feats, avail) else intersect(resp_feats, avail)

  engine <- all(vapply(c("build_responder_panel", "summarise_responder_subjects",
                         "add_tertile_bands", "build_resp_wide", "detect_subgroups",
                         "collapse_resp_level"),
                       exists, logical(1), mode = "function"))
  if (engine) {
    wmap <- if (exists("WHOLE_ASSAY_VISITS_DEFAULT")) WHOLE_ASSAY_VISITS_DEFAULT else list()
    long <- build_responder_panel(data_raw, mat_vis, rf, character(0), avail, wmap)
    subj <- add_tertile_bands(summarise_responder_subjects(long), breadth_thresh)
    sg   <- detect_subgroups(build_resp_wide(long), subj, has_cluster = has_cluster)
    if (isTRUE(sg$ok)) {
      cls <- data.frame(subject_accession = sg$subgroup_tbl$subject_accession,
                        responder = collapse_resp_level(sg$subgroup_tbl$subgroup),
                        score = sg$subgroup_tbl$universal_z, method = "PAM",
                        stringsAsFactors = FALSE)
    } else {
      cls <- data.frame(subject_accession = subj$subject_accession,
                        responder = collapse_resp_level(subj$score_band),
                        score = subj$universal_z, method = "tertile",
                        stringsAsFactors = FALSE)
    }
  } else {
    cls <- .fallback_responder(data_raw, mat_vis, rf, has_cluster)
  }
  cls <- cls[cls$responder %in% c("Low", "High"), , drop = FALSE]
  cls$responder <- factor(cls$responder, levels = c("Low", "High"))
  cls
}

#' Maternal arm (TdaP / TT) per subject from the antibody long frame.
subject_maternal_arm <- function(data_raw, arm_col = "arm_name") {
  d <- unique(data_raw[!is.na(data_raw[[arm_col]]),
                       c("subject_accession", arm_col)])
  d <- d[!duplicated(d$subject_accession), , drop = FALSE]
  names(d)[2] <- "arm_name"
  d
}


# ---- 1. small stats utilities ----------------------------------------------

# pvalue_to_label lives in R/serology_helpers.R (the shared utility). It is
# defined here ONLY as a fallback for the rare case this engine is sourced
# without the shared helpers; when serology_helpers.R is present (the normal
# pipeline), its definition is used and this one does not clobber it.
if (!exists("pvalue_to_label", mode = "function"))
  pvalue_to_label <- function(p)
    ifelse(is.na(p), "\u2014",
    ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns"))))

make_tertile <- function(x) {
  qs <- stats::quantile(x, c(1/3, 2/3), na.rm = TRUE)
  cut(x, c(-Inf, qs[1], qs[2], Inf),
      labels = c("T1 (low)", "T2 (mid)", "T3 (high)"), include.lowest = TRUE)
}

#' AUC that a numeric score ranks the second factor level above the first
#' (Mann-Whitney; no pROC dependency). y is a 2-level factor.
auc_binary <- function(score, y) {
  ok <- is.finite(score) & !is.na(y); s <- score[ok]; yy <- droplevels(factor(y[ok]))
  if (nlevels(yy) != 2) return(NA_real_)
  g1 <- s[yy == levels(yy)[2]]; g0 <- s[yy == levels(yy)[1]]
  if (!length(g1) || !length(g0)) return(NA_real_)
  r <- rank(c(g1, g0)); n1 <- length(g1)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(g0))
}

#' ROC points for a numeric score against a 2-level factor y (positive = level 2).
roc_points <- function(score, y) {
  ok <- is.finite(score) & !is.na(y); s <- score[ok]; yy <- factor(y[ok])
  pos <- yy == levels(yy)[2]
  thr <- sort(unique(s), decreasing = TRUE)
  tpr <- fpr <- numeric(length(thr) + 1)
  P <- sum(pos); N <- sum(!pos)
  for (i in seq_along(thr)) {
    pred <- s >= thr[i]
    tpr[i + 1] <- if (P) sum(pred & pos) / P else NA
    fpr[i + 1] <- if (N) sum(pred & !pos) / N else NA
  }
  data.frame(fpr = fpr, tpr = tpr)
}

.tidy_fit <- function(m, logit = FALSE) {
  co <- summary(m)$coefficients
  d <- data.frame(Predictor = rownames(co),
                  est = co[, 1], SE = co[, 2], stat = co[, 3], p = co[, 4],
                  stringsAsFactors = FALSE, row.names = NULL)
  if (logit) { d$OR <- exp(d$est); d$OR_lo <- exp(d$est - 1.96 * d$SE)
               d$OR_hi <- exp(d$est + 1.96 * d$SE) }
  d$sig <- pvalue_to_label(d$p)
  d
}

.null_fit <- function(label, predictor) list(
  label = label, predictor = predictor, fit = NA, R2 = NA, adj_R2 = NA,
  auc = NA_real_, n_obs = NA, p_model = NA, tidy = NULL, model = NULL)


# ---- 2. univariate fitter (lm for the score, logistic for the class) --------

#' Fit one predictor. mode = "lm" (continuous responder score) or
#' "logit" (binary Low/High class). Returns a common result structure where
#' R2 is adj-R^2 (lm) or McFadden pseudo-R^2 (logit), and auc is filled for logit.
fit_one <- function(data, outcome_col, predictor_col, mode, label = "") {
  f <- stats::reformulate(predictor_col, response = outcome_col)
  if (mode == "lm") {
    m <- tryCatch(stats::lm(f, data = data, na.action = stats::na.omit),
                  error = function(e) NULL)
    if (is.null(m)) return(.null_fit(label, predictor_col))
    s <- summary(m)
    pmod <- tryCatch(stats::pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                               lower.tail = FALSE), error = function(e) NA)
    list(label = label, predictor = predictor_col, fit = "lm",
         R2 = unname(s$r.squared), adj_R2 = unname(s$adj.r.squared),
         auc = NA_real_, n_obs = stats::nobs(m), p_model = unname(pmod),
         tidy = .tidy_fit(m), model = m)
  } else {
    m <- tryCatch(stats::glm(f, data = data, family = stats::binomial(),
                             na.action = stats::na.omit), error = function(e) NULL)
    if (is.null(m)) return(.null_fit(label, predictor_col))
    m0 <- tryCatch(stats::update(m, . ~ 1), error = function(e) NULL)
    pmod <- tryCatch(stats::anova(m0, m, test = "LRT")[["Pr(>Chi)"]][2],
                     error = function(e) NA)
    mcf <- tryCatch(1 - as.numeric(stats::logLik(m)) / as.numeric(stats::logLik(m0)),
                    error = function(e) NA)
    au <- tryCatch(auc_binary(stats::predict(m, type = "response"),
                              stats::model.frame(m)[[1]]), error = function(e) NA)
    list(label = label, predictor = predictor_col, fit = "logit",
         R2 = mcf, adj_R2 = mcf, auc = au, n_obs = stats::nobs(m),
         p_model = pmod, tidy = .tidy_fit(m, logit = TRUE), model = m)
  }
}

run_multi <- function(data, outcome_col, predictors, mode, label = "Multivariate") {
  if (length(predictors) < 2) return(NULL)
  f <- stats::reformulate(predictors, response = outcome_col)
  if (mode == "lm") {
    m <- tryCatch(stats::lm(f, data = data, na.action = stats::na.omit),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    s <- summary(m)
    pmod <- tryCatch(stats::pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                               lower.tail = FALSE), error = function(e) NA)
    list(label = label, R2 = unname(s$r.squared), adj_R2 = unname(s$adj.r.squared),
         auc = NA_real_, n_obs = stats::nobs(m), p_model = unname(pmod),
         tidy = .tidy_fit(m), model = m)
  } else {
    m <- tryCatch(stats::glm(f, data = data, family = stats::binomial(),
                             na.action = stats::na.omit), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    m0 <- stats::update(m, . ~ 1)
    pmod <- tryCatch(stats::anova(m0, m, test = "LRT")[["Pr(>Chi)"]][2],
                     error = function(e) NA)
    mcf <- tryCatch(1 - as.numeric(stats::logLik(m)) / as.numeric(stats::logLik(m0)),
                    error = function(e) NA)
    au <- tryCatch(auc_binary(stats::predict(m, type = "response"),
                              stats::model.frame(m)[[1]]), error = function(e) NA)
    list(label = label, R2 = mcf, adj_R2 = mcf, auc = au, n_obs = stats::nobs(m),
         p_model = pmod, tidy = .tidy_fit(m, logit = TRUE), model = m)
  }
}

#' Univariate screen over a named covariate list, with continuous + tertile
#' forms for continuous variables and form selection by the fit's R2 field.
run_clinical_screen <- function(d, outcome_col, cov_list, cat_vars,
                                mode, alpha = 0.05) {
  cont_v <- setdiff(unname(unlist(cov_list)), cat_vars)
  rows <- list(); res <- list()
  for (lbl in names(cov_list)) {
    v <- cov_list[[lbl]]
    if (v %in% cat_vars) {
      r <- fit_one(d, outcome_col, v, mode, label = lbl)
      rows[[paste0(v, "_cat")]] <- data.frame(
        Label = lbl, Variable = v, Form = "categorical",
        n = ifelse(is.na(r$n_obs), NA_integer_, r$n_obs),
        R2 = r$R2, adj_R2 = r$adj_R2, auc = r$auc, p_model = r$p_model,
        sig = pvalue_to_label(r$p_model), stringsAsFactors = FALSE)
      res[[paste0(v, "_cat")]] <- r
    } else {
      vT <- paste0(v, "_T")
      rc <- fit_one(d, outcome_col, v,  mode, label = paste0(lbl, " (continuous)"))
      rt <- fit_one(d, outcome_col, vT, mode, label = paste0(lbl, " (tertile)"))
      for (tag in c("cont", "tert")) {
        r <- if (tag == "cont") rc else rt
        rows[[paste0(v, "_", tag)]] <- data.frame(
          Label = lbl, Variable = if (tag == "cont") v else vT,
          Form = if (tag == "cont") "continuous" else "tertile",
          n = ifelse(is.na(r$n_obs), NA_integer_, r$n_obs),
          R2 = r$R2, adj_R2 = r$adj_R2, auc = r$auc, p_model = r$p_model,
          sig = pvalue_to_label(r$p_model), stringsAsFactors = FALSE)
      }
      res[[paste0(v, "_cont")]] <- rc; res[[paste0(v, "_tert")]] <- rt
    }
  }
  univ_table <- do.call(rbind, rows); rownames(univ_table) <- NULL

  selected <- list(); notes <- list()
  for (lbl in names(cov_list)) {
    v <- cov_list[[lbl]]
    if (v %in% cat_vars) {
      r <- res[[paste0(v, "_cat")]]
      if (!is.na(r$p_model) && r$p_model < alpha) {
        selected[[lbl]] <- v; notes[[lbl]] <- "categorical \u2014 significant" }
    } else {
      rc <- res[[paste0(v, "_cont")]]; rt <- res[[paste0(v, "_tert")]]
      sc <- !is.na(rc$p_model) && rc$p_model < alpha
      st <- !is.na(rt$p_model) && rt$p_model < alpha
      if (sc && st) {
        if (!is.na(rc$adj_R2) && !is.na(rt$adj_R2) && rc$adj_R2 >= rt$adj_R2) {
          selected[[lbl]] <- v
          notes[[lbl]] <- sprintf("continuous (R2=%.3f) over tertile (R2=%.3f) \u2014 both sig",
                                  rc$adj_R2, rt$adj_R2)
        } else {
          selected[[lbl]] <- paste0(v, "_T")
          notes[[lbl]] <- sprintf("tertile (R2=%.3f) over continuous (R2=%.3f) \u2014 both sig",
                                  rt$adj_R2, rc$adj_R2)
        }
      } else if (sc) { selected[[lbl]] <- v;            notes[[lbl]] <- "continuous \u2014 only form significant"
      } else if (st) { selected[[lbl]] <- paste0(v,"_T"); notes[[lbl]] <- "tertile \u2014 only form significant" }
    }
  }
  sel <- unname(unlist(selected))
  multi <- if (length(sel) > 1) run_multi(d, outcome_col, sel, mode,
                                          "Multivariate (all significant predictors)") else NULL
  list(univ_table = univ_table, selected = selected, form_notes = notes,
       multi_result = multi, univ_res = res, mode = mode)
}


# ---- 3. maternal serology signature -----------------------------------------

#' Tag a feature name (antigen_analyte) into coarse classes for grouping.
classify_feature <- function(feat) {
  f <- toupper(feat)
  analyte <- dplyr::case_when(
    grepl("IGG[1-4]", f) ~ "IgG subclass",
    grepl("IGG", f)      ~ "IgG total",
    grepl("IGA|IGM", f)  ~ "other Ig",
    grepl("FCG|FCR|FCGR|FCR?N", f) ~ "Fc receptor",
    grepl("ADCD|ADCP|ADNP|ADNKA|NK", f) ~ "effector function",
    grepl("PTNA|SBA", f) ~ "functional",
    TRUE ~ "other")
  antigen <- sub("_.*$", "", feat)
  data.frame(feature = feat, antigen = antigen, analyte_class = analyte,
             stringsAsFactors = FALSE)
}

#' Wide maternal serology matrix (subject x feature@visit), Blom-standardised
#' within each feature@visit, excluding the responder-defining features so the
#' signature is not circular with the score it predicts.
build_maternal_serology_matrix <- function(data_raw, visits, exclude_feats = character(0),
                                            min_coverage = 0.6, standardize = TRUE) {
  d <- data_raw[data_raw$visit_name %in% visits,
                c("subject_accession", "feature", "visit_name", "log_assay_value")]
  d$feature <- as.character(d$feature)
  d <- d[!(d$feature %in% exclude_feats), , drop = FALSE]
  if (!nrow(d)) return(data.frame(subject_accession = character()))
  d$key <- paste(d$feature, as.character(d$visit_name), sep = "@")
  agg <- stats::aggregate(log_assay_value ~ subject_accession + key, d,
                          FUN = function(z) mean(z, na.rm = TRUE))
  W <- as.data.frame(tidyr::pivot_wider(agg, id_cols = subject_accession,
                                        names_from = key, values_from = log_assay_value))
  feat_cols <- setdiff(names(W), "subject_accession")
  cov <- vapply(feat_cols, function(c_) mean(!is.na(W[[c_]])), numeric(1))
  keep <- feat_cols[cov >= min_coverage]
  W <- W[, c("subject_accession", keep), drop = FALSE]
  if (standardize) for (c_ in keep) W[[c_]] <- .blom(W[[c_]])
  W
}

#' Median-impute a numeric matrix column-wise (for complete-case-free modelling).
impute_median <- function(X) {
  X <- as.matrix(X)
  for (j in seq_len(ncol(X))) {
    v <- X[, j]; m <- stats::median(v[is.finite(v)], na.rm = TRUE)
    v[!is.finite(v)] <- if (is.finite(m)) m else 0; X[, j] <- v
  }
  X
}

#' Univariate AUC ranking of every serology feature for the Low/High class.
serology_univariate_auc <- function(W, class_tbl, id = "subject_accession") {
  m <- merge(W, class_tbl[, c(id, "responder")], by = id)
  feats <- setdiff(names(W), id)
  rows <- lapply(feats, function(f) {
    a <- auc_binary(m[[f]], m$responder)
    p <- tryCatch(stats::wilcox.test(m[[f]] ~ m$responder)$p.value,
                  error = function(e) NA_real_)
    data.frame(feature = f, auc = a, auc_dist = abs(a - 0.5),
               p = p, n = sum(is.finite(m[[f]]) & !is.na(m$responder)),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$p_BH <- stats::p.adjust(out$p, method = "BH")
  out[order(-out$auc_dist), , drop = FALSE]
}

#' Elastic-net logistic with CV-AUC over a serology matrix (optionally forcing
#' the maternal arm in unpenalised). Returns selected coefficients and CV-AUC.
#' Requires glmnet; degrades gracefully if absent.
elastic_net_signature <- function(W, class_tbl, id = "subject_accession",
                                   arm_tbl = NULL, alpha = 0.5, nfolds = 10,
                                   seed = 1) {
  if (!requireNamespace("glmnet", quietly = TRUE))
    return(list(ok = FALSE, reason = "glmnet not installed"))
  m <- merge(W, class_tbl[, c(id, "responder")], by = id)
  force_arm <- !is.null(arm_tbl)
  if (force_arm) {
    m <- merge(m, arm_tbl, by = id)
    m$arm_bin <- as.numeric(factor(m$arm_name)) - 1
  }
  feats <- setdiff(names(W), id)
  X <- impute_median(m[, feats, drop = FALSE]); colnames(X) <- feats
  if (force_arm) X <- cbind(arm_TdaP = m$arm_bin, X)
  y <- m$responder
  pf <- rep(1, ncol(X)); if (force_arm) pf[colnames(X) == "arm_TdaP"] <- 0
  set.seed(seed)
  cv <- tryCatch(glmnet::cv.glmnet(X, y, family = "binomial", alpha = alpha,
                                   nfolds = nfolds, type.measure = "auc",
                                   penalty.factor = pf), error = function(e) NULL)
  if (is.null(cv)) return(list(ok = FALSE, reason = "cv.glmnet failed"))
  auc <- cv$cvm[which(cv$lambda == cv$lambda.min)]
  co <- as.matrix(stats::coef(cv, s = "lambda.min"))
  co <- co[co[, 1] != 0, , drop = FALSE]
  list(ok = TRUE, cv = cv, cv_auc = auc, lambda = cv$lambda.min,
       n = nrow(X), n_feat = ncol(X),
       coefs = data.frame(term = rownames(co), coef = co[, 1],
                          stringsAsFactors = FALSE),
       X = X, y = y)
}
