### Loading functions from I-SPI

get_db_connection <- function() {
  dbConnect(RPostgres::Postgres(),
            dbname = Sys.getenv("db"),
            host = Sys.getenv("db_host"),
            port = Sys.getenv("db_port"),
            user = Sys.getenv("db_userid_x"),
            password = Sys.getenv("db_pwd_x"),
            options = "-c search_path=madi_results"
  )
}

check_plate <- function(selected_study){
  local_conn <- get_db_connection()
  on.exit(DBI::dbDisconnect(local_conn), add = TRUE)
  query_plates <- glue::glue_sql(
  "SELECT xmap_header_id, project_id, study_accession, experiment_accession,
    auth0_user, plate, nominal_sample_dilution
  	FROM madi_results.xmap_header
  	WHERE study_accession = {selected_study};", .con = local_conn)
  plates <- dbGetQuery(local_conn, query_plates)
  dbDisconnect(local_conn)
  plates <- plates[ , c("project_id", "study_accession", "experiment_accession", "auth0_user", "plate", "nominal_sample_dilution")]
  return(plates)
}

load_curves <- function(selected_study){
  local_conn <- get_db_connection()
  on.exit(DBI::dbDisconnect(local_conn), add = TRUE)
  query_curves <- glue::glue_sql(
    "SELECT xmap_header_id, xmap_header.project_id, xmap_header.study_accession, xmap_header.experiment_accession,
    auth0_user, xmap_header.plate, xmap_header.nominal_sample_dilution,
	curve_lookup.curve_id, curve_lookup.feature, curve_lookup.antigen, curve_lookup.source, curve_lookup.wavelength
  	FROM madi_results.xmap_header
	  INNER JOIN madi_results.curve_lookup ON
	  curve_lookup.project_id = xmap_header.project_id
	  AND curve_lookup.study_accession = xmap_header.study_accession
	  AND curve_lookup.experiment_accession = xmap_header.experiment_accession
	  AND curve_lookup.plate = xmap_header.plate
	  AND curve_lookup.nominal_sample_dilution = xmap_header.nominal_sample_dilution
  	WHERE xmap_header.study_accession = {selected_study};", .con = local_conn)
  curves <- dbGetQuery(local_conn, query_curves)
  dbDisconnect(local_conn)
  return(curves)
}


load_standards <- function(){
  local_conn <- get_db_connection()
  on.exit(DBI::dbDisconnect(local_conn), add = TRUE)
  query_standards <- glue::glue_sql(
    "SELECT curve_id, xmap_standard.study_accession, xmap_standard.experiment_accession, well,
stype, sampleid, xmap_standard.source, dilution, pctaggbeads, samplingerrors, xmap_standard.antigen,
antibody_mfi AS mfi, dilution * 100000 AS concentration, xmap_standard.feature, predicted_mfi, xmap_standard.project_id,
xmap_standard.nominal_sample_dilution, xmap_standard.plate, xmap_standard.wavelength
	FROM madi_results.xmap_standard
	  INNER JOIN madi_results.curve_lookup ON
	  curve_lookup.project_id = xmap_standard.project_id
	  AND curve_lookup.study_accession = xmap_standard.study_accession
	  AND curve_lookup.experiment_accession = xmap_standard.experiment_accession
	  AND curve_lookup.plate = xmap_standard.plate
	  AND curve_lookup.nominal_sample_dilution = xmap_standard.nominal_sample_dilution
  	WHERE xmap_standard.study_accession IN ('MADI_P3_GAPS', 'Gaps subclasses');", .con = local_conn)
  standards <- dbGetQuery(local_conn, query_standards)
  dbDisconnect(local_conn)
  return(standards)
}

pull_samples <- function(selected_study) {
  # conn <- get_db_connection()
  local_conn <- get_db_connection()
  on.exit(DBI::dbDisconnect(local_conn), add = TRUE)
  query_samples <- glue::glue_sql(
  "SELECT best_sample_se_all.project_id, best_sample_se_all.study_accession, best_sample_se_all.experiment_accession,
best_sample_se_all.plate, best_sample_se_all.nominal_sample_dilution, best_sample_se_all.feature,
best_sample_se_all.antigen, best_sample_se_all.source, best_sample_se_all.wavelength,
	best_sample_se_all.well, best_sample_se_all.sampleid, best_sample_se_all.patientid, best_sample_se_all.timeperiod,
	best_sample_se_all.dilution,
	best_sample_se_all.assay_response, best_sample_se_all.raw_predicted_concentration, best_sample_se_all.se_concentration,
	best_sample_se_all.final_predicted_concentration, best_sample_se_all.pcov,
	bayes_samples.mfi AS bayes_assay_response,
	bayes_samples.raw_predicted_concentration AS bayes_raw_predicted_concentration,
	bayes_samples.se_concentration AS bayes_se_concentration,
--	bayes_samples.raw_predicted_concentration * best_sample_se_all.dilution AS bayes_final_predicted_concentration,
	bayes_samples.pcov AS bayes_pcov
  	FROM madi_results.best_sample_se_all
	LEFT OUTER JOIN madi_results.bayes_samples ON
          bayes_samples.project_id = best_sample_se_all.project_id
	  AND bayes_samples.study_accession = best_sample_se_all.study_accession
	  AND bayes_samples.experiment_accession = best_sample_se_all.experiment_accession
	  AND bayes_samples.feature = best_sample_se_all.feature
	  AND bayes_samples.antigen = best_sample_se_all.antigen
	  AND bayes_samples.source = best_sample_se_all.source
	  AND bayes_samples.wavelength = best_sample_se_all.wavelength
	  AND bayes_samples.plate = best_sample_se_all.plate
	  AND bayes_samples.nominal_sample_dilution = best_sample_se_all.nominal_sample_dilution
	  AND bayes_samples.well = best_sample_se_all.well
  	WHERE best_sample_se_all.study_accession = {selected_study};", .con = local_conn)

#   select_query <- glue::glue_sql(
#       		"SELECT DISTINCT
#       		xmap_sample_id, xmap_sample.study_accession, xmap_sample.experiment_accession, xmap_sample.plate_id, xmap_sample.timeperiod,
# patientid, well, stype, sampleid, id_imi, agroup, pctaggbeads, samplingerrors, antigen, antibody_mfi AS MFI, antibody_au AS AU,
#         		dilution AS sample_dilution_factor, antibody_n, antibody_name, feature, gate_class, antibody_au_se, reference_dilution,
#         		gate_class_dil, norm_mfi,
#         		CASE
#         		  WHEN gate_class IN ('Between_Limits','Acceptable') THEN 'Acceptable'
#               WHEN gate_class IN ('Below_Lower_Limit','Too Diluted') THEN 'Too Diluted'
#       		    WHEN gate_class IN ('Above_Upper_Limit','Too Concentrated') THEN 'Too Concentrated'
#               WHEN gate_class IN ('Not Evaluated') OR gate_class IS NULL THEN 'Not Evaluated' END AS gclod,
#             CASE
#               WHEN gate_class_linear_region IN ('Between_Limits','Acceptable') THEN 'Acceptable'
#               WHEN gate_class_linear_region IN ('Below_Lower_Limit','Too Diluted') THEN 'Too Diluted'
#               WHEN gate_class_linear_region IN ('Above_Upper_Limit','Too Concentrated') THEN 'Too Concentrated'
#               WHEN gate_class_linear_region IN ('Not Evaluated') OR gate_class IS NULL THEN 'Not Evaluated' END AS gclin,
#             CASE
#               WHEN gate_class_loq IN ('Between_Limits','Acceptable') THEN 'Acceptable'
#               WHEN gate_class_loq IN ('Below_Lower_Limit','Too Diluted') THEN 'Too Diluted'
#               WHEN gate_class_loq IN ('Above_Upper_Limit','Too Concentrated') THEN 'Too Concentrated'
#               WHEN gate_class_loq IN ('Not Evaluated') OR gate_class IS NULL THEN 'Not Evaluated' END AS gcloq,
#             CASE WHEN antibody_n < lower_bc_threshold THEN 'LowBeadN' ELSE 'Acceptable' END AS lowbeadn,
#             CASE WHEN pctaggbeads > pct_agg_threshold THEN 'PctAggBeads' ELSE 'Acceptable' END AS highbeadagg
#       		FROM madi_results.xmap_sample
#           INNER JOIN (
#             SELECT study_accession, param_integer_value AS lower_bc_threshold
#             FROM madi_results.xmap_study_config
#     		    WHERE study_accession = {selected_study} AND param_user = {current_user} AND param_name = 'lower_bc_threshold'
#   		    ) AS bct ON bct.study_accession = xmap_sample.study_accession
#           INNER JOIN (
#             SELECT study_accession, param_integer_value AS pct_agg_threshold
#             FROM madi_results.xmap_study_config
#     		    WHERE study_accession = {selected_study} AND param_user = {current_user} AND param_name = 'pct_agg_threshold'
#   		    ) AS pab ON pab.study_accession = xmap_sample.study_accession
#   		    WHERE xmap_sample.study_accession = {selected_study};", .con = conn)



  sample_data <- dbGetQuery(local_conn, query_samples)

  # sample_data$plate_id <- str_trim(str_replace_all(sample_data$plate_id, "\\s", ""), side = "both")
  # sample_data$plate_id <- toupper(sample_data$plate_id)
  # sample_data$feature <- ifelse(is.na(sample_data$feature), sample_data$experiment_accession, sample_data$feature)

  # sample_data <- merge(sample_data,
  #                      plates[ , ! names(plates) %in% c("feature")],
  #                      by=c("plate_id","sample_dilution_factor"),
  #                      all.x = TRUE)
  sample_data <- distinct(sample_data)
  dbDisconnect(local_conn)
  return(sample_data)
}

split_well_column <- function(df, well_col = "well",
                              new_letters = "well_column",
                              new_numbers = "well_row") {
  if (!is.data.frame(df)) stop("df must be a data.frame or tibble")
  if (!well_col %in% names(df)) stop(sprintf("Column '%s' not found in df", well_col))

  # Ensure well column is character
  well_vals <- as.character(df[[well_col]])

  # Use regex to extract leading letters and trailing numbers
  # ^([A-Za-z]+) captures leading letters; ([0-9]+)$ captures trailing digits
  letters_part <- sub("^([A-Za-z]+).*", "\\1", well_vals)
  numbers_part <- sub("^.*?([0-9]+)$", "\\1", well_vals)

  # If extraction didn't match, set to NA
  letters_part[!grepl("^[A-Za-z]+", well_vals)] <- NA
  numbers_part[!grepl("[0-9]+$", well_vals)] <- NA

  # Convert numbers to integer (or NA where not present)
  numbers_part_num <- as.integer(numbers_part)

  # Add new columns to a copy of df
  df[[new_letters]] <- letters_part
  df[[new_numbers]] <- numbers_part_num

  return(df)
}
