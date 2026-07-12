#this code must not be run before connection.R

#the standardized_set does not return any values for features that correspond to the standardized analyte
standardized_set <- function (standardized_analyte = "IgG",
                             connection_data = NULL,
                             residual_file = NULL
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
