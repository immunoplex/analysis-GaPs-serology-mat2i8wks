#this code must not be run before connection.R

#the standardized_set does not return any values for features that correspond to the standardized analyte
standardized_set <- function (standardized_analyte = "IgG",
                             connection_data = NULL,
                             residual_file = NULL,
                             include_feature    = NULL,
                             include_feature_wo = NULL
                             ){

    ### load residuals scores.
    load(residual_file)
    ## select residuals for each include_feature that match the model to the selected model
    feature_residuals <- data_residv_total[data_residv_total$selmodel == data_residv_total$modeln &
                                             toupper(data_residv_total$feature) %in% toupper(include_feature),
                                           c("subject_accession", "feature", "visit_name", "modeln", "models", "resid", "selmodel")]
    # view(filter(feature_residuals, if_any(everything(), ~ is.na(.))))
    table(feature_residuals$feature)
    table(feature_residuals$visit_name)

    arm_dat <- dplyr::distinct(connection_data[ , c("subject_accession", "arm_name")])
    data_res_chk <- merge(feature_residuals[!is.null(feature_residuals$resid), ], arm_dat, by = "subject_accession", all.x = TRUE)

    table(data_res_chk$subject_accession, data_res_chk$arm_name)
    table(connection_data$subject_accession, connection_data$arm_name)

    connection_data <- connection_data[toupper(connection_data$feature) %in% toupper(include_feature_wo),]

    ## check for duplicates

    keys <- colnames(feature_residuals)[!grepl('resid',colnames(feature_residuals))]
    X <- as.data.table(feature_residuals)
    feature_residuals <- as.data.frame(X[,list(resid= mean(resid)),keys])
    rm(X)
    dups_resid <- {feature_residuals} %>%
      dplyr::group_by(subject_accession, visit_name, feature) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(n > 1L)
    # check for na that could prevent matrix operations
    tfeature_residuals <- pivot_wider(feature_residuals,
                                      id_cols = c("subject_accession","visit_name"),
                                      names_from = "feature",
                                      values_from = "resid"
                                      )
    table(tfeature_residuals$visit_name)
    table(na.omit(tfeature_residuals)$visit_name)

    # feature_residuals <<- data_residv_total[data_residv_total$selmodel == data_residv_total$modeln, ]

    #code for after model selection for residuals (need to select resid before this would work)
    connection_data$upper_feature <- toupper(connection_data$feature)
    feature_residuals$upper_feature <- toupper(feature_residuals$feature)
    feature_residuals <- feature_residuals[ , !names(feature_residuals) %in% c("feature")]
    data_resid_merged <- merge(x = connection_data,
                        y = feature_residuals,
                        by = c("subject_accession","visit_name","upper_feature"),
                        all.x = TRUE
                        )
    table(data_resid_merged$subject_accession, data_resid_merged$arm_name)
    ### remove rows for feature was the data standardized by
    regex_pattern <- paste0(standardized_analyte, "$")
    selected_rows <- grepl(regex_pattern, data_resid_merged$feature)
    data_resid <- data_resid_merged[!selected_rows,]
    table(data_resid$subject_accession, data_resid$arm_name)

    ##  move standardized value (resid) to log_assay_value
    data_resid$log_assay_value <- data_resid$resid
    data_resid$value_reported <- 10^data_resid$log_assay_value
    # data_resid$value_imputed<-data_resid$value_reported

    ## keep variables that match original connection_data
    # connection_data_resid<-data_resid[, c("subject_accession", "study_accession",
    #                                       "experiment_accession", "analyte_reported", "value_reported", "value_imputed", "visit_name", "arm_name", "antigen",
    #                                       "feature", "log_assay_value", "antigen_family", covariates)]
    #
    # revised to drop because above step means covariates need to have been specified
    # drop step keeps all connection_data variables from above and removes all feature_residual variables not needed
    connection_data_resid <- data_resid[!names(data_resid) %in% c("modeln", "models", "resid",
                                                                "weight", "residw", "absresidw",
                                                                "selmodel", "upper_feature")]

    return(connection_data_resid)
}

## ---------------------------------------------------------------------------
## C0 IgG standardisation / residual transform
## Shared engine for analysis/parts/C0_standardise.Rmd.
## Hardwired to d_set.RData / ebaa_extra2 (c_set.RData is retired).
## ---------------------------------------------------------------------------

## --- fit-quality helpers (defined once) ------------------------------------
.r2 <- function(x) {
  SSe <- sum((x$resid)^2)
  obs <- x$resid + x$fitted
  1 - SSe / sum((obs - mean(obs))^2)
}
.rsquare_inlier <- function(x) {
  x <- x[x$wclas >= 0.5, ]
  SSe <- sum((x$resid)^2)
  obs <- x$resid + x$fitted
  1 - SSe / sum((obs - mean(obs))^2)
}
.adj_r2 <- function(nobs, ncoef, rsq) 1 - ((1 - rsq) * (nobs - 1) / (nobs - ncoef - 1))

## --- main entry point ------------------------------------------------------
## Returns invisibly a list(residuals=, modelfit=) and writes both to data/.
standardise_arm_visit <- function(selected_infant_arm,
                                  timeperiods,
                                  antigenlist, igglist, fulligglist,
                                  wholeassay_map,
                                  residual_path, modelfit_path,
                                  data_file = here::here("data", "d_set.RData")) {

  ## Biological contract: whole-assay regressor antigens must exist in the vocab,
  ## otherwise SBA/WT-IgG/PTNA would be standardised against the wrong (or no) IgG.
  stopifnot(all(wholeassay_map$igg_antigen %in% antigenlist))
  stopifnot(all(wholeassay_map$igg_antigen %in% sub("_IgG$", "", igglist)))

  load(data_file)                       # restores ebaa_extra2
  data <- ebaa_extra2[ebaa_extra2$infant_arm == selected_infant_arm, ]

  ## Arm labelling (maternal TdaP vs TT, ref = TT)
  data$arm_name <- factor(data$maternal_arm)
  levels(data$arm_name) <- c("TdaP", "TT")
  data$arm_name <- relevel(data$arm_name, ref = "TT")

  ## Rebuild feature from antigen + analyte (canonical mixed-case)
  data$analyte <- as.character(data$analyte)
  data$antigen <- as.character(data$antigen)
  data$feature <- factor(paste(data$antigen, data$analyte, sep = "_"))

  ## --- whole-assay relabel: reassign antigen to the standardising antigen ---
  data$feature <- as.character(data$feature)
  wa_idx <- data$feature %in% wholeassay_map$source_feature
  if (sum(wa_idx) == 0) {
    warning("wholeassay_map: none of {",
            paste(wholeassay_map$source_feature, collapse = ", "),
            "} present for arm ", selected_infant_arm,
            " -- whole-bacterium assays skipped.")
  } else {
    wa <- wholeassay_map[match(data$feature[wa_idx], wholeassay_map$source_feature), ]
    data$antigen[wa_idx] <- wa$igg_antigen   # regress against this antigen's total IgG
    data$feature[wa_idx] <- wa$out_feature   # residual key
    message("whole-assay rows relabelled (", selected_infant_arm, "): ", sum(wa_idx))
  }
  data$feature <- factor(data$feature)
  data$antigen <- factor(data$antigen)

  data_e <- data[data$antigen %in% antigenlist & data$visit_name %in% timeperiods, ]
  data_e$visit_name <- droplevels(data_e$visit_name)
  data_e$antigen    <- droplevels(data_e$antigen)
  data_e$log_assay_value <- log10(data_e$value_reported + 0.0001)

  ## Merge per-timepoint total IgG (t0/t1/t2) as the regressor
  igg_at <- function(tp) {
    d <- data_e[data_e$feature %in% igglist & data_e$visit_name == tp,
                c("subject_accession","visit_name","antigen",
                  "value_reported","log_assay_value")]
    names(d)[names(d) == "log_assay_value"] <- "log_total_IgG_MFI"
    names(d)[names(d) == "value_reported"]  <- "total_IgG_MFI"
    d
  }
  data_igg_MFI <- do.call(rbind, lapply(timeperiods, igg_at))
  data_exna <- merge(data_e, data_igg_MFI,
                     by = c("subject_accession","visit_name","antigen"), all = TRUE)

  ## Drop IgG-vs-IgG self comparisons
  data_exna_s <- as.data.frame(data_exna[!(data_exna$feature %in% fulligglist), ])
  data_exna_s$feature <- droplevels(data_exna_s$feature)

  data_residv_total  <- data.frame()
  data_lmstatsv_total <- data.frame()

  ## ---- per antigen x feature model ladder + selection --------------------
  ## (linear/gam x unweighted/discrete/continuous; identical logic to the
  ##  original six files, factored into fit_feature() for one code path.)
  # for (antigen in as.character(unique(data_exna_s$antigen))) {
  #   ad <- data_exna_s[data_exna_s$antigen == antigen, ]
  #   for (feat in unique(ad$feature)) {
  #     fit <- .fit_feature(ad[ad$feature == feat, ], feat)   # helper below
  #     data_residv_total   <- rbind(data_residv_total,  fit$resid)
  #     data_lmstatsv_total <- rbind(data_lmstatsv_total, fit$stats)
  #   }
  # }
  for (antigen in as.character(unique(data_exna_s$antigen))) {
    if (render) { cat(sprintf("\n#### Antigen: %s {.tabset}\n\n", antigen)) }
    ad <- data_exna_s[data_exna_s$antigen == antigen, ]
    for (feat in unique(ad$feature)) {
      fit <- .fit_feature(ad[ad$feature == feat, ], feat)
      data_residv_total   <- rbind(data_residv_total,  fit$resid)
      data_lmstatsv_total <- rbind(data_lmstatsv_total, fit$stats)
      if (render) {
        cat(sprintf("\n##### Feature: %s\n\n", feat))
        .render_feature(fit$data_m1, fit$selmodel, feat)
      }
    }
  }

  save(data_residv_total,  file = residual_path)
  save(data_lmstatsv_total, file = modelfit_path)
  message("wrote ", basename(residual_path), " and ", basename(modelfit_path))
  invisible(list(residuals = residual_path, modelfit = modelfit_path))
}
