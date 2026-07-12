# =====================================================================
# R/connection_transform.R
# ---------------------------------------------------------------------
# IgG-standardization helper for the response-bubble analysis (Part 4).
# Project copy of the original connection_transform.R, with two changes:
#   1. visit_name in the loaded residual object is recoded to the canonical
#      names (recode_visits) so it joins to the recoded antibody data;
#   2. the include-feature lists are passed in (defaulting to the BUBBLE_*
#      vocabulary in config/endpoints.R) instead of relying on globals.
# Requires data.table and tidyr (loaded by the bubble part).
#
# The standardized_set() returns NO rows for the analyte the data were
# standardized by (e.g. total IgG); those are re-added by the caller.
# =====================================================================

standardized_set <- function(standardized_analyte = "IgG",
                             connection_data = NULL,
                             residual_file = NULL,
                             include_feature = BUBBLE_INCLUDE_FEATURE,
                             include_feature_wo = BUBBLE_INCLUDE_FEATURE_WO) {

  ## load residual scores (object: data_residv_total)
  load(residual_file)

  ## recode the residuals' visit_name to the canonical names so it matches
  ## the recoded antibody data (the single global recode rule)
  if (exists("recode_visits"))
    data_residv_total$visit_name <- recode_visits(data_residv_total$visit_name)

  ## select residuals matching the selected model and the included features
  feature_residuals <- data_residv_total[
    data_residv_total$selmodel == data_residv_total$modeln &
      toupper(data_residv_total$feature) %in% toupper(include_feature),
    c("subject_accession", "feature", "visit_name", "modeln", "models", "resid", "selmodel")]

  arm_dat <- dplyr::distinct(connection_data[, c("subject_accession", "arm_name")])
  data_res_chk <- merge(feature_residuals[!is.null(feature_residuals$resid), ],
                        arm_dat, by = "subject_accession", all.x = TRUE)

  connection_data <- connection_data[
    toupper(connection_data$feature) %in% toupper(include_feature_wo), ]

  ## collapse duplicate residuals (mean) using data.table
  keys <- colnames(feature_residuals)[!grepl("resid", colnames(feature_residuals))]
  X <- data.table::as.data.table(feature_residuals)
  feature_residuals <- as.data.frame(X[, list(resid = mean(resid)), keys])
  rm(X)

  ## join residuals onto the connection data on subject x visit x feature
  connection_data$upper_feature  <- toupper(connection_data$feature)
  feature_residuals$upper_feature <- toupper(feature_residuals$feature)
  feature_residuals <- feature_residuals[, !names(feature_residuals) %in% c("feature")]
  data_resid_merged <- merge(x = connection_data, y = feature_residuals,
                             by = c("subject_accession", "visit_name", "upper_feature"),
                             all.x = TRUE)

  ## drop the analyte the data were standardized by (e.g. total IgG)
  selected_rows <- grepl(paste0(standardized_analyte, "$"), data_resid_merged$feature)
  data_resid <- data_resid_merged[!selected_rows, ]

  ## move the standardized value (resid) into log_assay_value
  data_resid$log_assay_value <- data_resid$resid
  data_resid$value_reported  <- 10^data_resid$log_assay_value

  ## keep only the original connection_data variables
  data_resid[!names(data_resid) %in% c("modeln", "models", "resid", "weight",
                                       "residw", "absresidw", "selmodel", "upper_feature")]
}
