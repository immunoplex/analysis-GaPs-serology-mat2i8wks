library(readxl)
library(dplyr)
library(tidyr)
library(here)
library(stringr)
library(purrr)
library(RPostgres)
library(glue)
library(janitor)
library(data.table)
library(DBI)
library(DT)
library(papeR)
library(foreach)
library(stringr)
library(ggplot2)

source(here::here("./code/db_load_functions.R"))
current_user <- "mscotzens@gmail.com"

standards_raw <- load_standards()
standards <- standards_raw


standards$set <- ifelse(standards$study_accession=='MADI_P3_GAPS',"inf_baa","mat_baa")
standards$source <- ifelse(standards$source %in% c("NIBSC06_140","Sando","NIBSC"),"NIBSC06_140", standards$source)

adcd_mat <- standards[standards$set == "mat_baa"
                      & standards$experiment_accession %in% c("ADCD")
                      & standards$dilution == 80, ]

standards$set_sel <- ifelse(standards$set == "inf_baa", 1, 0)
standards$set_sel <- ifelse(standards$set == "mat_baa"
                            & standards$experiment_accession %in% c("ADCD")
                            & standards$nominal_sample_dilution == 80, 1, standards$set_sel)
standards$set_sel <- ifelse(standards$set == "mat_baa"
                            & ! standards$experiment_accession %in% c("ADCD"), 1, standards$set_sel)


standards$feature <- ifelse(standards$experiment_accession == "IgG1","IgG1",standards$feature)
standards$feature <- ifelse(standards$experiment_accession == "IgG2","IgG2",standards$feature)
standards$feature <- ifelse(standards$experiment_accession == "IgG3","IgG3",standards$feature)
standards$feature <- ifelse(standards$experiment_accession == "IgG4","IgG4",standards$feature)
standards$feature <- ifelse(standards$experiment_accession == "ADCD","ADCD",standards$feature)
standards <- standards[standards$source != "V2",]
standards <- standards[standards$set_sel == 1 ,]

standards_sel <- standards[ ,
                          c("set","feature","antigen","plate","source","mfi","concentration", "dilution")]

summary(standards_sel)
table(standards_sel$set)
table(standards_sel$feature)
table(standards_sel$antigen)
table(standards_sel$plate)
table(standards_sel$source)
table(standards_sel$concentration)

save(standards_sel, file = here::here("./data/standards_sel.RData"))


### ---- Read sample data from bead arrays for infants ----
# selected_study <- 'MADI_P3_GAPS'
# plates <- check_plate(selected_study)
# curves <- load_curves(selected_study)
# src_inf_baa <- pull_samples(selected_study)
# src_inf_baa$sampleid <- as.numeric(src_inf_baa$sampleid)
# src_inf_baa$patientid <- as.numeric(src_inf_baa$patientid)
# table(src_inf_baa$plate, src_inf_baa$feature)
# table(src_inf_baa$plate, src_inf_baa$feature, src_inf_baa$timeperiod)
# inf_baa <- src_inf_baa
# inf_baa$visit_name <- ifelse(inf_baa$timeperiod %in% c("pre1st","pre2nd"), "M02", inf_baa$timeperiod)
# inf_baa$visit_name <- ifelse(inf_baa$visit_name %in% c("pre3rd","pre4th"), "M04", inf_baa$visit_name)
# inf_baa$visit_name <- ifelse(inf_baa$visit_name %in% c("post1st","post3rd"), "M05", inf_baa$visit_name)
# inf_baa$visit_name <- ifelse(inf_baa$visit_name %in% c("post5th", "post3rd5mo"), "M09", inf_baa$visit_name)
# inf_baa$visit_name <- factor(inf_baa$visit_name)
# inf_baa$subject_accession <- as.numeric(inf_baa$patientid)
# inf_baa$antigen <- sapply(strsplit(inf_baa$antigen, "[-_:]"), `[`, 1)
# inf_baa$sample_type <- "B"
# save(inf_baa, file = here("./data/inf_baa.RData"))

# table(inf_baa$visit_name)

# ---- Read sample data from bead arrays for mothers ----
selected_study <- 'Gaps subclasses'
mat_plates <- check_plate(selected_study)
mat_curves <- load_curves(selected_study)
src_mat_baa <- pull_samples(selected_study)
src_mat_baa$sampleid <- as.numeric(src_mat_baa$sampleid)
src_mat_baa$patientid <- as.numeric(src_mat_baa$patientid)
src_mat_baa <- src_mat_baa[src_mat_baa$experiment_accession !=
                             "IgG3_zensTest", ]
table(src_mat_baa$plate, src_mat_baa$feature)
table(src_mat_baa$timeperiod)
mat_baa <- src_mat_baa
mat_baa$visit_name <- ifelse(mat_baa$timeperiod %in% c("prevacc"), "P02", mat_baa$timeperiod)
mat_baa$visit_name <- ifelse(mat_baa$visit_name %in% c("delivery"), "P09", mat_baa$visit_name)
mat_baa$visit_name <- ifelse(mat_baa$visit_name %in% c("final"), "P18", mat_baa$visit_name)
mat_baa$visit_name <- ifelse(mat_baa$visit_name %in% c("cord"), "M00", mat_baa$visit_name)
mat_baa$visit_name <- factor(mat_baa$visit_name)
mat_baa$antigen <- sapply(strsplit(mat_baa$antigen, "[-_:]"), `[`, 1)
mat_baa$sample_type <- ifelse(mat_baa$timeperiod == "cord", "C", "M")


mat_baa$set_sel <- ifelse(mat_baa$experiment_accession %in% c("ADCD")
                            & mat_baa$nominal_sample_dilution == 80, 1, 0)
mat_baa$set_sel <- ifelse(!mat_baa$experiment_accession %in% "ADCD", 1, mat_baa$set_sel)

mat_baa <- mat_baa[mat_baa$set_sel  == 1, ]

save(mat_baa, file = here("./data/mat_baa.RData"))


load(file = here::here("./data/inf_baa.RData"))

inf_baa <- inf_baa[inf_baa$source %in% c("NIBSC06_140","Sando"),]
inf_baa$source <- ifelse(inf_baa$source %in% c("NIBSC06_140","Sando"),"NIBSC06_140", inf_baa$source)
table(inf_baa$source)


load(file = here::here("./data/plate_df.RData"))
plate_df$plate <- paste0("plate_",plate_df$plate_number)
plate_df$visit_name <- ifelse(plate_df$timeperiod=="cord","M00",plate_df$timeperiod)
plate_df$visit_name <- ifelse(plate_df$timeperiod=="maternal pre vacc","P02",plate_df$visit_name)
plate_df$visit_name <- ifelse(plate_df$timeperiod=="maternal at delivery","P09",plate_df$visit_name)
plate_df$visit_name <- ifelse(plate_df$timeperiod=="maternal final","P18",plate_df$visit_name)
names(plate_df)[names(plate_df) == "patientID"] <- "patientid"

plateidcorr <- plate_df[, c("plate","well","patientid","visit_name", "timeperiod")]
load(file = here::here("./data/mat_baa.RData"))
mat_baa <- mat_baa[ , c("project_id","study_accession","experiment_accession","plate",
                       "nominal_sample_dilution","feature","antigen","source","wavelength",
                       "well","sampleid","dilution","assay_response",
                       "raw_predicted_concentration","se_concentration",
                       "final_predicted_concentration","pcov","bayes_assay_response",
                       "bayes_raw_predicted_concentration","bayes_se_concentration",
                       "bayes_pcov","sample_type" )]

mat_baa <- inner_join(mat_baa,plateidcorr, by = c("plate","well"))
mat_baa$subject_accession <- as.numeric(mat_baa$patientid)

mat_baa$feature <- ifelse(mat_baa$experiment_accession == "IgG1","IgG1",mat_baa$feature)
mat_baa$feature <- ifelse(mat_baa$experiment_accession == "IgG2","IgG2",mat_baa$feature)
mat_baa$feature <- ifelse(mat_baa$experiment_accession == "IgG3","IgG3",mat_baa$feature)
mat_baa$feature <- ifelse(mat_baa$experiment_accession == "IgG4","IgG4",mat_baa$feature)
mat_baa$feature <- ifelse(mat_baa$experiment_accession == "ADCD","ADCD",mat_baa$feature)

mat_igg1_baa <- mat_baa[mat_baa$feature == 'IgG1', ]
mat_igg2_baa <- mat_baa[mat_baa$feature == 'IgG2', ]
mat_igg3_baa <- mat_baa[mat_baa$feature == 'IgG3', ]
mat_igg4_baa <- mat_baa[mat_baa$feature == 'IgG4', ]
mat_adcd_baa <- mat_baa[mat_baa$feature == 'ADCD', ]
table(mat_igg1_baa$plate, mat_igg1_baa$visit_name)
table(mat_igg2_baa$plate, mat_igg2_baa$visit_name)
table(mat_igg3_baa$plate, mat_igg3_baa$visit_name)
table(mat_igg4_baa$plate, mat_igg4_baa$visit_name)
table(mat_adcd_baa$plate, mat_adcd_baa$visit_name)
table(mat_igg1_baa$subject_accession, mat_igg1_baa$visit_name)
table(mat_igg2_baa$subject_accession, mat_igg2_baa$visit_name)
table(mat_igg3_baa$subject_accession, mat_igg3_baa$visit_name)
table(mat_igg4_baa$subject_accession, mat_igg4_baa$visit_name)
table(mat_adcd_baa$subject_accession, mat_adcd_baa$visit_name)

mat_baa$source = ifelse(mat_baa$source == "NIBSC","NIBSC06_140",mat_baa$source)
table(mat_baa$source)
mat_baa <- mat_baa[mat_baa$source %in% c("NIBSC06_140","Sando"),]
table(mat_baa$source)
mat_baa_x <- mat_baa
mat_baa_x$set <- "mat_baa"
mat_baa_x$mfi <- 10^mat_baa$assay_response

norm_lookup <- read.csv('normalization_params.csv')

# Normalize mat_baa samples onto the inf_baa scale:
#
# The join keys must match normalization_params.csv exactly. Two value-domain
# mismatches exist between the sample side and the params side and are handled here:
#   (1) source: samples use "Sando", params use "SD"  -> recode "Sando" -> "SD"
#   (2) antigen_base: params are lowercase (act, pt, ipv1, ...), sample antigen
#       casing comes straight from the DB -> force lowercase before joining.
# antigen_base also strips any trailing batch-code suffix (e.g. act_42 -> act),
# matching how standards_comparison.Rmd built the params.
mat_baa_normalised <- mat_baa_x %>%
  filter(set == 'mat_baa') %>%
  mutate(
    antigen_base = tolower(str_remove(antigen, '_\\d+$')),
    source_join  = ifelse(source == 'Sando', 'SD', source)
  ) %>%
  left_join(
    norm_lookup,
    by = c('feature' = 'feature',
           'antigen_base' = 'antigen_base',
           'plate' = 'plate',
           'source_join' = 'source')
  ) %>%
  mutate(
    mfi_norm = 10^(cross_intercept_mat2inf +
                   cross_slope_mat2inf * log10(pmax(mfi, 0.5)))
  )

# FIX: write the normalised log-MFI INTO mat_baa_normalised (the frame that flows
# downstream into mat_baa_s -> baa), NOT into mat_baa_x (which is discarded).
# Previously assay_response in mat_baa_normalised still held the ORIGINAL
# un-normalised value carried in from the left_join, so the normalisation never
# reached baa / value_reported and the analysis used un-normalised mat_baa values.
mat_baa_normalised$assay_response <- log10(mat_baa_normalised$mfi_norm)

# Safety check: rows that failed the norm_lookup join have NA mfi_norm and would
# silently become NA assay_response. Flag them so the join coverage is explicit.
if (any(is.na(mat_baa_normalised$mfi_norm))) {
  warning(sprintf("%d of %d mat_baa rows had no matching normalization_params row (NA mfi_norm).",
                  sum(is.na(mat_baa_normalised$mfi_norm)), nrow(mat_baa_normalised)))
}

mat_baa_s <- mat_baa_normalised[ , names(mat_baa)]

baa <- rbind(inf_baa,mat_baa_s)
names(baa)
table(baa$feature)
table(baa$visit_name)
table(baa$sample_type)
table(baa$antigen)
table(baa$source)
baa$feature <- ifelse(baa$feature == "Total_IgG", "IgG", baa$feature)

baa$analyte <- baa$feature
baa$feature <- paste(baa$antigen, baa$feature, sep = "_")

table(baa$feature, baa$visit_name)
table(baa$feature, baa$visit_name, baa$source)


### use for MFI responses
baa$value_reported <- 10^baa$assay_response ### calculated this way because the log10 value is always available and not dependent on any curve fitting of any standard curve process.
baa$value_reported <- ifelse(baa$value_reported < 0.001, 0.001, baa$value_reported)
baa$value_reported <- ifelse(baa$value_reported > 50000, 50000, baa$value_reported)
baa$log_assay_value <- log10(baa$value_reported + 1)
### end use for MFI responses
### use for concentration responses
# baa$value_reported <- ifelse(is.na(baa$bayes_raw_predicted_concentration),
#                              baa$final_predicted_concentration,
#                              baa$bayes_raw_predicted_concentration * baa$dilution)
# baa$value_reported <- ifelse(!is.na(baa$value_reported) & baa$value_reported < 0.00000001, 0.00000001, baa$value_reported)
# baa$value_reported <- ifelse(!is.na(baa$value_reported) & baa$value_reported > 1000000, 1000000, baa$value_reported)
# baa <- baa[!is.na(baa$assay_response),]
# baa$value_reported <- ifelse(is.na(baa$value_reported) & baa$assay_response > 2, 1000000, baa$value_reported)
# # baa$value_reported <- (10^baa$assay_response) - 1
# baa$log_assay_value <- log10(baa$value_reported + 0.00001)
### end use for concentration responses


baa$feature <- factor(toupper(baa$feature))
baa$antigen <- factor(toupper(baa$antigen))

table(baa$visit_name, baa$sample_type)
table(baa$antigen)
table(baa$feature)
unique(baa$subject_accession)
table(baa$feature,baa$source, baa$visit_name, baa$sample_type)
table(baa$analyte,baa$antigen,baa$visit_name)
table(baa$analyte,baa$antigen,baa$visit_name)

baa_select <- baa %>%
  dplyr::filter(analyte=="IgG1" & visit_name=="M00")

table(baa_select$subject_accession, baa_select$analyte)



### ---- load RIVM total IgG data ----
varnames <- c("Participant", "Time_point", "Participant_number", "SampleCode", "laborder",
              "sample_type", "PT", "FHA", "PRN", "FIM", "DT", "TT")
datastringpath <- here::here("./raw_files/Raw_MIA_for_R_msz.xlsx")
totIgGraw <- read_excel(paste0(datastringpath), sheet = "Sheet1", col_names = TRUE)
totIgGw <- totIgGraw[totIgGraw$sample_type=="B",]
table(totIgGw$Time_point)

tigg20 <- totIgGw[totIgGw$Time_point =="20", ]$Participant_number
tigg21 <- totIgGw[totIgGw$Time_point =="21", ]$Participant_number
tigg23 <- totIgGw[totIgGw$Time_point =="23", ]$Participant_number
extratigg21 <- !tigg21 %in% tigg20
extratigg23 <- !tigg23 %in% tigg20

tigg40 <- totIgGw[totIgGw$Time_point =="40", ]$Participant_number
tigg42 <- totIgGw[totIgGw$Time_point =="42", ]$Participant_number
tigg45 <- totIgGw[totIgGw$Time_point =="45", ]$Participant_number
tigg47 <- totIgGw[totIgGw$Time_point =="47", ]$Participant_number ### 4months plus 7 days
tigg50 <- totIgGw[totIgGw$Time_point =="50", ]$Participant_number
extratigg42 <- !tigg42 %in% tigg40
extratigg52 <- !tigg42 %in% tigg50 ### no records
extratigg45 <- !tigg45 %in% tigg40
extratigg55 <- !tigg45 %in% tigg50
extratigg47 <- !tigg47 %in% tigg40
extratigg57 <- !tigg47 %in% tigg50

# as.data.frame(tigg43[extratigg43])



totIgGw$visit_name <- ifelse(totIgGw$Time_point %in% c("20"), "M02", totIgGw$Time_point)
totIgGw$visit_name <- ifelse(totIgGw$Time_point %in% c("23") & totIgGw$Participant_number %in% tigg23[extratigg23], "M04", totIgGw$visit_name)
totIgGw$visit_name <- ifelse(totIgGw$Time_point %in% c("40"), "M04", totIgGw$visit_name)
# totIgGw$visit_name <- ifelse(totIgGw$Time_point %in% c("47") & totIgGw$Participant_number %in% tigg47[extratigg47], "mo4", totIgGw$visit_name)

totIgGw$visit_name <- ifelse(totIgGw$Time_point %in% c("50"), "M05", totIgGw$visit_name)
# totIgGw$visit_name <- ifelse(totIgGw$Time_point %in% c("47") & totIgGw$Participant_number %in% tigg47[extratigg57], "vaccinated", totIgGw$visit_name)
# totIgGw$visit_name <- ifelse(totIgGw$Time_point %in% c("45") & totIgGw$Participant_number %in% tigg45[extratigg55], "vaccinated", totIgGw$visit_name)
totIgGw$visit_name <- ifelse(totIgGw$Time_point %in% c("60"), "M09", totIgGw$visit_name)
totIgGw <- totIgGw[totIgGw$visit_name %in% c("M02","M04","M05", "M09"),]
table(totIgGw$visit_name,totIgGw$Time_point)

totIgGw$subject_accession <- as.numeric(totIgGw$Participant_number)
totIgGtw <- pivot_longer(totIgGw,
                       cols = c("PT","PRN","FHA","FIM","DT","TT"),
                       names_to = "antigen",
                       values_to = "value_reported")

totIgGtw$antigen <- ifelse(totIgGtw$antigen=="FIM","FIM",totIgGtw$antigen)
totIgGtw$antigen <- factor(totIgGtw$antigen)
totIgGtw$analyte <- factor("IgG")
totIgGtw$feature <- factor(paste(totIgGtw$antigen, totIgGtw$analyte, sep = "_"))
totIgGtw <- totIgGtw[ , c("subject_accession", "sample_type", "visit_name", "antigen", "analyte", "feature", "value_reported")]


totIgGm <- totIgGraw[totIgGraw$sample_type=="M",]
totIgGm$visit_name <- ifelse(totIgGm$Time_point=='00','P09',totIgGm$Time_point)
totIgGm$visit_name <- ifelse(totIgGm$visit_name=='02','P02',totIgGm$visit_name)
totIgGm$visit_name <- ifelse(totIgGm$visit_name=='50','P05',totIgGm$visit_name)
totIgGm$visit_name <- ifelse(totIgGm$visit_name=='60','P18',totIgGm$visit_name)
table(totIgGm$Time_point)
totIgGm$subject_accession <- as.numeric(totIgGm$Participant_number)
totIgGtm <- pivot_longer(totIgGm,
                       cols = c("PT","PRN","FHA","FIM","DT","TT"),
                       names_to = "antigen",
                       values_to = "value_reported")
totIgGtm$analyte <- factor("IgG")
totIgGtm$feature <- factor(paste(totIgGtm$antigen, totIgGtm$analyte, sep = "_"))
totIgGtm <- totIgGtm[ , c("subject_accession", "sample_type", "visit_name", "antigen", "analyte", "feature", "value_reported")]

totIgGc <- totIgGraw[totIgGraw$sample_type=="C" & totIgGraw$SampleCode != 'HC13200X',]
totIgGc$Time_point <- ifelse(totIgGc$Participant_number=="204",'00',totIgGc$Time_point)
# totIgGc$Time_point <- ifelse(totIgGc$SampleCode == 'HC13200X','00',totIgGc$Time_point)
# totIgGc$Participant <- ifelse(totIgGc$SampleCode == 'HC13200X','HC',totIgGc$Participant)
# totIgGc$Participant_number <- ifelse(totIgGc$SampleCode == 'HC13200X','132',totIgGc$Participant_number)
totIgGc$Participant <- ifelse(totIgGc$SampleCode == 'HC64700X','HC', totIgGc$Participant)
totIgGc$visit_name <- 'M00'
table(totIgGc$Time_point)
totIgGc$subject_accession <- as.numeric(totIgGc$Participant_number)
totIgGtc <- pivot_longer(totIgGc,
                        cols = c("PT","PRN","FHA","FIM","DT","TT"),
                        names_to = "antigen",
                        values_to = "value_reported")
totIgGtc$analyte <- factor("IgG")
totIgGtc$feature <- factor(paste(totIgGtc$antigen, totIgGtc$analyte, sep = "_"))
totIgGtc <- totIgGtc[ , c("subject_accession", "sample_type", "visit_name", "antigen", "analyte", "feature", "value_reported")]
totIgG <- rbind(totIgGtw, totIgGtc, totIgGtm)
totIgG$log_assay_value <- log10(totIgG$value_reported+0.001)
save(totIgG, file = here::here("./data/totalIgG.RData"))


load(here::here("./data/totalIgG.RData"))
table(totIgG$visit_name, totIgG$sample_type)
table(totIgG$antigen)
table(totIgG$feature)

unique(totIgG$subject_accession)

# Find common values between the two vectors
intersect(baa$subject_accession, totIgG$subject_accession)

# Or to see which values are in baa but NOT in totIgG
setdiff(baa$subject_accession, totIgG$subject_accession)

# Or to see which values are in totIgG but NOT in baa
setdiff(totIgG$subject_accession, baa$subject_accession)

baa$assay <- "mbaa"
totIgG$assay <- "elisa"
xbaa <- baa[baa$subject_accession %in% intersect(baa$subject_accession, totIgG$subject_accession), c("subject_accession","sample_type","assay","visit_name","antigen","analyte",
                                                                                                     "feature","value_reported","log_assay_value")]
xtotIgG <- totIgG[totIgG$subject_accession %in% intersect(baa$subject_accession, totIgG$subject_accession),]
ebaa <- rbind(xbaa,
              xtotIgG)
ebaa$feature <- toupper(ebaa$feature)
ebaa$assay_priority <- case_when(
  tolower(ebaa$assay) == "elisa" ~ 1,
  tolower(ebaa$assay) == "mbaa"  ~ 2,
  TRUE                      ~ 3
)

summary(ebaa)

table(ebaa$feature, ebaa$visit_name)

# key <- c("subject_accession","visit_name","sample_type","antigen","feature","analyte","value_reported", "log_assay_value")
key <- c("subject_accession","visit_name","sample_type","antigen","feature","analyte")


nrow(ebaa)
ebaa %>% distinct(across(all_of(key))) %>% nrow()   # this is your ceiling after dedup

ebaa %>%
  filter(visit_name %in% c("M00","P02","P09","P18")) %>%
  group_by(across(all_of(key))) %>%
  dplyr::summarise(n_in_key = n(), .groups = "drop") %>%
  dplyr::mutate(feat = sub(".*_", "", feature)) %>%     # IGG1 / IGG2 / ...
  dplyr::count(feat, n_in_key) %>%
  dplyr::arrange(feat, n_in_key)

ebaa %>%
  filter(feature == "PT_IGG1", visit_name == "M00") %>%
  dplyr::group_by(across(all_of(key))) %>%
  dplyr::summarise(across(everything(), ~ n_distinct(.)), .groups = "drop") %>%
  dplyr::summarise(across(-any_of(key), max)) %>%
  dplyr::glimpse()

ex_key <- ebaa %>%
  filter(feature == "PT_IGG1", visit_name == "M00") %>%
  add_count(across(all_of(key)), name = "n") %>%
  dplyr::filter(n == 16) %>% slice(1) %>%
  dplyr::select(all_of(key))

ebaa %>%
  semi_join(ex_key, by = key) %>%
  dplyr::arrange(value_reported) %>%
  as.data.frame()

ebaa %>%
  dplyr::filter(grepl("_IGG1$", feature)) %>%
  dplyr::count(assay, assay_priority)

# (1) Does subject 38 explode only for IgG1, or across all subclasses?
ebaa %>%
  dplyr::filter(subject_accession == 38, visit_name == "M00", antigen == "PT") %>%
  dplyr::count(feature, analyte)

# (2) Rows vs unique subjects per isotype — the ID-collapse signature
ebaa %>%
  dplyr::mutate(iso = sub(".*_", "", feature)) %>%
  dplyr::filter(visit_name == "M00", antigen == "PT") %>%
  group_by(iso) %>%
  dplyr::summarise(n_rows = n(), n_subjects = n_distinct(subject_accession))

ebaa_clean <- ebaa %>%
  # Sort so elisa comes first within each group
  arrange(subject_accession, visit_name, sample_type, antigen, feature, analyte, assay_priority) %>%
  # Keep first record per natural key (which will be elisa if it exists)
  distinct(subject_accession, visit_name, sample_type, antigen, feature, analyte,
           .keep_all = TRUE) %>%
  # Drop the helper column
  dplyr::select(-assay_priority)

table(ebaa_clean$feature, ebaa_clean$visit_name)

# ### ---- create the other features data from bead array spreadsheet ----
# load(file = here::here("./data/frombaa.RData"))
# otherfeat <- frombaa
# table(otherfeat$feature,otherfeat$antigen,otherfeat$timeperiod)
# # otherfeat$analyte <- case_when(
# #   otherfeat$analyte == 'ADCP_127' ~ "ADCP_90",
# #   otherfeat$analyte == 'ADCP_630' ~ "ADCP_90",
# #   otherfeat$analyte == 'ADNP_127' ~ "ADNP_90",
# #   otherfeat$analyte == 'ADNP_630' ~ "ADNP_90",
# #   TRUE ~ otherfeat$analyte)
# otherfeat$subject_accession = as.numeric(otherfeat$patientid)
# # ADCDC time periods
# # post 3rd dose (+1 mo) post 3rd dose (+5 mo) pre 1st dose (2 mo) pre 3rd dose (4 mo)
# # post1st                  3652                     0                   0                   0
# # post5th                     0                  3542                   0                   0
# # pre2nd                      0                     0                3454                   0
# # pre4th                      0                     0                   0                2310
# otherfeat$timeperiod <- ifelse(otherfeat$timeperiod %in% c("pre1st","pre2nd"), "prevaccinated", otherfeat$timeperiod)
# otherfeat$timeperiod <- ifelse(otherfeat$timeperiod %in% c("pre3rd","pre4th"), "mo4", otherfeat$timeperiod)
# otherfeat$timeperiod <- ifelse(otherfeat$timeperiod %in% c("post1st","post3rd"), "vaccinated", otherfeat$timeperiod)
# otherfeat$timeperiod <- ifelse(otherfeat$timeperiod %in% c("post5th", "post3rd5mo"), "vaccinated9", otherfeat$timeperiod)
# otherfeat$visit_name <- factor(otherfeat$timeperiod)
# otherfeat$antigen <- factor(toupper(otherfeat$antigen))
# otherfeat$analyte <- factor(otherfeat$feature)
# otherfeat$feature <- factor(paste(otherfeat$antigen, otherfeat$analyte, sep = "_"))
#
# table(otherfeat$analyte,otherfeat$antigen,otherfeat$visit_name)
#
# otherfeat$agroup <- ifelse(otherfeat$agroup=="Tdap_aP","TdaP_aP",otherfeat$agroup)
# otherfeat$maternal_arm <- factor(ifelse(otherfeat$agroup %in% c("TdaP_aP", "TdaP_wP"), "TdaP","TT"))
# otherfeat$infant_arm <- factor(ifelse(otherfeat$agroup %in% c("TdaP_aP", "TT_aP"), "aP","wP"))
# otherfeat$arm_name <- factor(paste(otherfeat$maternal_arm, otherfeat$infant_arm, sep = "_"))
# otherfeat$value_reported = as.numeric(otherfeat$mfi)
# otherfeat <- otherfeat[ , c("subject_accession","visit_name","maternal_arm","infant_arm", "arm_name", "antigen", "analyte", "feature", "value_reported")]
# save(otherfeat, file = here("./data/otherfeat.RData"))
# table(otherfeat$analyte,otherfeat$antigen,otherfeat$visit_name)
#
# ### add other antigens to totalIgG
# add_antig <- otherfeat[otherfeat$antigen %in% c("ACT","IPV1","IPV2","IPV3") & otherfeat$analyte == "Total_IgG",names(totIgG)]
# add_antig$analyte <- "IgG"
# add_antig$feature <- paste(add_antig$antigen,add_antig$analyte, sep = "_")
# table (add_antig$feature)
# table(add_antig$visit_name)
#
# names(add_antig)


### ---- combine totalIgG from RIVM and other features from bead array and join the maternal and infant arm designations ----
# load(file = here("./data/totalIgG.RData"))
# load(file = here("./data/otherfeat.RData"))

visit_list <- unique(ebaa_clean$visit_name)
# antigen_list <- unique(totIgG$antigen)
antigen_list <- c("PT","PRN","FHA","FIM","DT","TT","ACT","IPV1","IPV2","IPV3")
bindfeat <- ebaa_clean[ebaa_clean$visit_name %in% visit_list & ebaa_clean$antigen %in% antigen_list, names(ebaa_clean)]

datastringpath <- here::here("./raw_files/AT_GaPs data base_msz.xlsx")
subjfeat <- read_excel(paste0(datastringpath), sheet = "Sheet1", col_names = TRUE)
# subjfeat <- distinct(otherfeat[otherfeat$visit_name %in% visit_list & otherfeat$antigen %in% antigen_list, c("subject_accession", "arm_name", "infant_arm", "maternal_arm")])

# data <- rbind(totIgG,bindfeat)
# data <- rbind(data,add_antig)
data <- merge(bindfeat, subjfeat[, c("subject_accession","arm_name","maternal_arm","infant_arm")],
              by = "subject_accession", all.x = TRUE)
# data$visit_name <- factor(data$visit_name)
data$arm_name <- factor(data$arm_name)
data$maternal_arm <- factor(data$maternal_arm)
data$infant_arm <- factor(data$infant_arm)

data$visit_name <- droplevels(data$visit_name)
data$antigen <- droplevels(data$antigen)
data$feature <- droplevels(data$feature)

table(data$analyte,data$antigen,data$visit_name)

### ---- create the clinical_assessments dataset ----
# datastringpath <- here("./raw_files/GaPs_clinical_data_per_protocol.xlsx")
# raw_asses <- read_excel(paste0(datastringpath), sheet = "GaPs_clinical_data_per_protocol", col_names = TRUE)
# save(raw_asses, file = here("./data/raw_assess.RData"))
load(file = here::here("./data/raw_assess.RData"))
clin_assess <- raw_asses
clin_assess$subject_accession <- as.numeric(clin_assess$Periscope_stem)
clin_assess$maternal_age <- as.numeric(clin_assess$Age)
clin_assess$gestational_age_vaccination <- as.numeric(clin_assess$Gestation_M02)
clin_assess$gestational_age_birth <- as.numeric(clin_assess$Gestation_DEL)
clin_assess$vaccine_birth_interval_days <- as.numeric(clin_assess$vaccine_birth_interval_days)
clin_assess$maternal_bmi <- as.numeric(clin_assess$Maternal_BMI)
clin_assess$maternal_Hb <- as.numeric(clin_assess$Hb)
clin_assess$infant_sex <- factor(clin_assess$Sex)
clin_assess$parity_recode <- ifelse(clin_assess$Parity %in% c("2","3","4"),"2+",
                                    ifelse(clin_assess$Parity == "1", "1",
                                           ifelse(clin_assess$Parity == "0", "0", "NA")))
table(clin_assess$Parity, clin_assess$parity_recode)
clin_assess$parity <- factor(clin_assess$parity_recode)

clin_assess$delivery_mode <- factor(ifelse(clin_assess$Delivery_mode %in% c("Planned CS", "Emergency CS"), "Caesarian Section", "Vaginal Delivery"))
clin_assess$mat_id <- clin_assess$Maternal_ID
clin_assess$inf_id <- clin_assess$GaPs_ID
clin_assess$birth_weight <- as.numeric(clin_assess$Weight_DEL)
clin_assess <- clin_assess[ , c("subject_accession", "mat_id", "inf_id", "maternal_age", "gestational_age_vaccination", "gestational_age_birth",
                                "vaccine_birth_interval_days", "maternal_bmi", "infant_sex", "parity", "delivery_mode", "birth_weight",
                                "maternal_Hb"
                                )]
save(clin_assess, file = here("./data/clin_assess.RData"))

### ---- merge assay data to clinical assessment: a_set ----
load(file = here::here("./data/clin_assess.RData"))

data <- merge (data, clin_assess, by = "subject_accession", all.x = TRUE)

duplicates <- data |>
  group_by_all() |>
  filter(n() > 1) |>
  ungroup()

dist_data <- distinct(data)

keys <- colnames(dist_data)[!grepl('value_reported',colnames(dist_data))]
library(plyr)
dedup_data <-ddply(dist_data,keys,summarize, value=mean(value_reported))
names(dedup_data)[names(dedup_data) == "value"] <- "value_reported"


duplicate_counts <- dedup_data |>
  add_count(subject_accession, visit_name, feature) |>
  filter(n > 1) |>
  distinct()

data <- dedup_data

table(data$analyte,data$antigen,data$visit_name)


save(data, file = paste0(getwd(),"/data/a_set.RData"))



### ---- add Bmem PTNA SBA WT_IgG: b_set ----
load(file = paste0(getwd(),"/data/a_set.RData"))
datastringpath <- here::here("./raw_files/GaPs_immunogenicity_data_per_protocol_ALL_20_08_24_plusBrem_MSZ.xlsx")
extrabca <- readxl::read_excel(paste0(datastringpath), sheet = "GaPs_immunogenicity_data_per_pr", col_names = TRUE)
names(extrabca)[names(extrabca) == "Periscope_stem"] <- "subject_accession"
extrabca$visit_name <- dplyr::case_when (
  extrabca$Time_point == "V2" ~ "M02",
  extrabca$Time_point == "V5" ~ "M05",
  extrabca$Time_point == "V6" ~ "M09",
  TRUE ~ "OTHER"
)
table(extrabca$visit_name, extrabca$Time_point)

rextrabca <- extrabca[extrabca$Antigen %in% c("PTNA","SBA", "WT_IgG"),c("subject_accession","visit_name","Antigen","MIA_IgG")]
names(rextrabca)[names(rextrabca) == "Antigen"] <- "analyte"
names(rextrabca)[names(rextrabca) == "MIA_IgG"] <- "value_reported"
rextrabca$antigen <- "PT"
rextrabca$value_reported <- as.numeric(rextrabca$value_reported)
table(rextrabca$visit_name, rextrabca$analyte)
rextrabca$log_assay_value <- log10(rextrabca$value_reported)
rextrabca$antigen <- "WHOLE"
rextrabca$assay <- "WT"
rextrabca$sample_type <- "B"
rextrabca$feature <- factor(paste(rextrabca$antigen, rextrabca$analyte, sep = "_"))

table(rextrabca$visit_name, rextrabca$feature)

# rsbextraprn <- rextrabca[rextrabca$analyte %in% c("SBA", "WT_IgG"),]
# rsbextraprn$antigen <- "PRN"
# rsbextrafha <- rextrabca[rextrabca$analyte %in% c("SBA", "WT_IgG"),]
# rsbextrafha$antigen <- "FHA"
# rsbextrafim <- rextrabca[rextrabca$analyte %in% c("SBA", "WT_IgG"),]
# rsbextrafim$antigen <- "FIM"
# rsbextratt <- rextrabca[rextrabca$analyte %in% c("SBA", "WT_IgG"),]
# rsbextratt$antigen <- "TT"
# rsbextradt <- rextrabca[rextrabca$analyte %in% c("SBA", "WT_IgG"),]
# rsbextradt$antigen <- "DT"

sextrabca <- extrabca[extrabca$Bmem_IgG != 'NA',c("subject_accession","visit_name","Antigen","Bmem_IgG")]
names(sextrabca)[names(sextrabca) == "Antigen"] <- "antigen"
names(sextrabca)[names(sextrabca) == "Bmem_IgG"] <- "value_reported"
sextrabca$analyte <- "Bmem"
sextrabca$value_reported <- as.numeric(sextrabca$value_reported)
sextrabca$log_assay_value <- log10(sextrabca$value_reported)
sextrabca$assay <- "WT"
sextrabca$sample_type <- "B"
sextrabca$feature <- factor(paste(sextrabca$antigen, sextrabca$analyte, sep = "_"))



# id_data <- dplyr::distinct(data[ ,c("subject_accession","visit_name",
#                     "arm_name","maternal_arm","infant_arm","mat_id","inf_id","maternal_age","gestational_age_vaccination",
#                     "gestational_age_birth","vaccine_birth_interval_days","maternal_bmi","infant_sex",
#                     "parity","delivery_mode","birth_weight","maternal_Hb")])
id_data <- dplyr::distinct(data[ ,c("subject_accession","visit_name","arm_name","maternal_arm","infant_arm")])


datastringpath <- here::here("./raw_files/GaPs_serum_deposition_msz.xlsx")
wtIgGIgA <- read_excel(paste0(datastringpath), sheet = "Sheet1", col_names = TRUE)
wtIgGIgA$analyte <- paste(wtIgGIgA$wt_mt, wtIgGIgA$antibody, sep = "_")
wtIgGIgA$antigen <- "WHOLE"
wtIgGIgA <- wtIgGIgA %>%
  mutate(sample_type = str_sub(Sample_ID, 2, 2))
wtIgGIgA$subject_accession <- as.numeric(wtIgGIgA$ID)
wtIgGIgA$value_reported <- as.numeric(wtIgGIgA$Conc)
wtIgGIgA$log_assay_value <- log10(wtIgGIgA$value_reported)
wtIgGIgA$feature <- wtIgGIgA$analyte
wtIgGIgA$assay <- "WT"
wtIgGIgA$visit_name <- ifelse(wtIgGIgA$timepoint == "00", "M00", wtIgGIgA$timepoint)
wtIgGIgA$visit_name <- ifelse(wtIgGIgA$timepoint == "20", "M02", wtIgGIgA$visit_name)
wtIgGIgA$visit_name <- ifelse(wtIgGIgA$timepoint == "50", "M05", wtIgGIgA$visit_name)
wtIgGIgA$visit_name <- ifelse(wtIgGIgA$timepoint == "60", "M09", wtIgGIgA$visit_name)
cordWTIgG <- wtIgGIgA[wtIgGIgA$timepoint %in% c("00","20","50", "60") &
                       !is.na(wtIgGIgA$value_reported) &
                       wtIgGIgA$analyte == "WT_IgG" & wtIgGIgA$visit_name == "M00"
                       , c(
  "subject_accession","sample_type","visit_name","assay","antigen","analyte","feature","value_reported","log_assay_value"
)]


table(cordWTIgG$visit_name)
table(cordWTIgG$sample_type)
table(cordWTIgG$feature)
table(cordWTIgG$feature,cordWTIgG$visit_name, cordWTIgG$sample_type)


textrabca <- rbind(rextrabca, sextrabca, cordWTIgG)
# textrabca <- rbind(textrabca, rsbextraprn)
# textrabca <- rbind(textrabca, rsbextrafha)
# textrabca <- rbind(textrabca, rsbextrafim)
# textrabca <- rbind(textrabca, rsbextratt)
# textrabca <- rbind(textrabca, rsbextradt)
table(textrabca$analyte,textrabca$antigen,textrabca$visit_name)
textrabca$feature <- paste(textrabca$antigen,textrabca$analyte, sep = "_")
table(textrabca$feature,textrabca$visit_name)

# extrabca_f <- merge(clin_assess, textrabca, by = c("subject_accession"),all.y = TRUE)
# extrabca_f <- merge(id_data, extrabca_f, by = c("subject_accession","visit_name"),all.y = TRUE)
#
# table(extrabca_f$feature)
# table(extrabca_f$feature,extrabca_f$visit_name)


ebaa_extra <- rbind(ebaa_clean, textrabca)

ebaa_extra <- inner_join(ebaa_extra, subjfeat[, c("subject_accession","arm_name","maternal_arm","infant_arm")],
              by = "subject_accession")

# clin_assess <- distinct(clin_assess)
# ebaa_extra <- left_join(ebaa_extra, clin_assess,
#                          by = "subject_accession")

table(ebaa_extra$infant_arm)
save(ebaa_extra, file = here::here("./data/c_set.RData"))




load(file = here::here("./data/c_set.RData"))



### ---- cytokines ----
datastringpath <- here::here("./raw_files/20251215_IMI-Periscope2022_results PBLSGI_copie 19122025_msz.xlsx")
cytok <- readxl::read_excel(paste0(datastringpath), sheet = "GaPs_cytokines", col_names = TRUE)
# names(extrabca)[names(extrabca) == "Periscope_stem"] <- "subject_accession"
cytok$visit_name <- dplyr::case_when (
  cytok$timepoint == "V02" ~ "prevaccinated",
  cytok$timepoint == "V05" ~ "vaccinated",
  cytok$timepoint == "V06" ~ "vaccinated9",
  TRUE ~ "OTHER"
)
table(cytok$visit_name, cytok$timepoint)

names(cytok)

cytok_long <- pivot_longer(data = cytok, cols = c("G-CSF","GM-CSF","IFN-g","IL-1b","IL-2","IL-4","IL-5","IL-6","IL-7","IL-8","IL-10","IL12 p70","IL-13","IL-17","MCP-1","MIP-1b","TNFa"), names_to = "analyte", values_to = "value_reported")
cytok_long <- cytock_long[ , c("subject_accession", "timepoint", "visit_name","sample_type","analyte","value_reported")]
names(cytok_long)
cytok_long$sample_type <- ifelse(cytok_long$sample_type=="Mock", "MOCK", cytok_long$sample_type)


cytok_w <- pivot_wider(data = cytok_long, id_cols = c("subject_accession", "timepoint", "visit_name","analyte"), names_from = "sample_type", values_from = "value_reported")

cytok_w$calc_conc_PT = as.numeric(cytok_w$PT) - as.numeric(cytok_w$MOCK)

cytok_w$calc_conc_Bppp = as.numeric(cytok_w$Bppp) - as.numeric(cytok_w$MOCK)

cytok_w$calc_conc = as.numeric(cytok_w$PT) + as.numeric(cytok_w$Bppp) - as.numeric(cytok_w$MOCK)

cytok_w <- cytok_w[ , c("subject_accession","timepoint","analyte","calc_conc_PT","calc_conc_Bppp", "calc_conc")]
names(cytok_w)

cytok_ww <- pivot_wider(data = cytok_w, id_cols = c("subject_accession","analyte"), names_from = "timepoint", values_from = "calc_conc")
cytok_ww$V02 <- ifelse(cytok_ww$V02 < 0 , 0.0000001, cytok_ww$V02)
cytok_ww$V05 <- ifelse(cytok_ww$V05 < 0 , 0.0000001, cytok_ww$V05)
cytok_ww$cyto_ratio <- cytok_ww$V05 / cytok_ww$V02
cytok_ww <- cytok_ww[, c("subject_accession","analyte","cyto_ratio", "V02", "V05")]
names(cytok_ww)
summary(cytok_ww)
### ---- v02 quintiles ----

add_v02_quintiles <- function(df, na.rm = TRUE) {
  df %>%
    group_by(analyte) %>%
    mutate(
      V02_quintile = case_when(
        is.na(V02) ~ NA_integer_,
        is.infinite(V02) ~ NA_integer_,
        TRUE ~ ntile(V02[is.finite(V02)], 5)[match(
          row_number(),
          which(is.finite(V02))
        )]
      )
    ) %>%
    group_by(analyte, V02_quintile) %>%
    mutate(
      # V02_quintile_label = case_when(
      #   is.na(V02_quintile) ~ NA_character_,
      #   TRUE ~ paste0(
      #     round(min(V02, na.rm = TRUE), 2),
      #     "-",
      #     round(max(V02, na.rm = TRUE), 2)
      #   )
      # )

      V02_quintile_label = case_when(
        is.na(V02_quintile) ~ NA_character_,
        TRUE ~ sprintf("%.3f-%.3f",
                       min(V02, na.rm = TRUE),
                       max(V02, na.rm = TRUE))
      )
    ) %>%
    ungroup()
}

cytok_wwq <- as.data.frame(add_v02_quintiles(cytok_ww, na.rm=FALSE))
summary(cytok_wwq)


plot_cyto_by_v02_quintile <- function(df, ncol = 3, point_alpha = 0.3,
                                      show_points = TRUE, log_scale = FALSE) {

  # Create a label lookup within each analyte
  df_plot <- df %>%
    filter(!is.na(V02_quintile), !is.infinite(cyto_ratio)) %>%
    group_by(analyte, V02_quintile) %>%
    mutate(y_label = first(V02_quintile_label)) %>%
    ungroup() %>%
    mutate(y_label = factor(y_label, levels = unique(y_label[order(analyte, V02_quintile)])))

  # Reorder y_label within each analyte
  df_plot <- df_plot %>%
    arrange(analyte, V02_quintile) %>%
    mutate(
      y_label_ordered = factor(
        V02_quintile_label,
        levels = unique(V02_quintile_label)
      )
    )

  p <- ggplot(df_plot, aes(x = cyto_ratio, y = reorder(V02_quintile_label, V02_quintile))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
    geom_boxplot(orientation = "y", outlier.shape = NA, fill = "lightgray") +
    facet_wrap(~ analyte, scales = "free_y", ncol = ncol) +
    labs(
      x = "Cytokine Ratio (V05/V02)",
      y = "V02 Quintile"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightblue"),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 7),
      panel.grid.minor = element_blank()
    )

  if (show_points) {
    p <- p + geom_jitter(height = 0.2, alpha = point_alpha, size = 1, color = "steelblue")
  }

  if (log_scale) {
    p <- p + scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
      labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")
    )
  }

  return(p)
}


plot_cyto_by_v02_quintile(cytok_wwq, ncol = 4, point_alpha = 0.5, log_scale = TRUE)


unique(data_extra$feature)




