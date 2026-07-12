# ---
# title: "PTNA associations with IgG total and IgG subclasses binding to pertussis antigens"
# subtitle: "Stratified by infant arm (aP vs wP) — parallel to SBA analysis"
# author: "Scot Zens"
# ---
library(data.table)           # for the function in the connection_transform.R code
library(tidyverse)            # dplyr, ggplot2, forcats, tidyr, stringr, tibble, magrittr pipe
library(here)                 # here::here used in file paths
library(reshape2)             # reshape2::pivot_wider
library(RColorBrewer)         # brewer.pal color palettes
library(corrplot)
library(lme4)
library(NetworkComparisonTest)
library(networktools)
library(bootnet)
library(qgraph)
library(broom)
library(relaimpo)
library(tidyr)

source(here::here("./code/colors.R"), local = TRUE)
iarm.apwps.color <- c('aP' = "dodgerblue2", 'wP' = "forestgreen")
marm.apwps.color <- c('TdaP' = "skyblue",   'TT'  = "salmon")

color.Welch.fill <- brewer.pal(3, "YlGn")[2]

covariate.colors <- c(
  "arm_name"                    = "tomato",
  "birth_weight"                = "lightblue",
  "infant_sex"                  = "forestgreen",
  "delivery_mode"               = "#9467bd",
  "parity"                      = "#8c564b",
  "gravidity"                   = "firebrick",
  "gestational_age_vaccination" = "darkorange",
  "vaccine_birth_interval_days" = "steelblue",
  "gestational_age_birth"       = "mediumseagreen",
  "maternal_age"                = "salmon",
  "maternal_bmi"                = "#FEE08B",
  "maternal_Hb"                 = "#C51B7D"
)

antigen.colors <- c(
  "zmixture" = "firebrick",
  "zPC1"     = "salmon",
  "zPC2"     = "#E04C5C",
  "DT"       = "#FEE08B",
  "FHA"      = "#23AECE",
  "PRN"      = "steelblue",
  "PT"       = "#AF8DC3",
  "TT"       = "mediumseagreen"
)

include_feature <- c("DT_IgG","DT_IgG1","DT_IgG2","DT_IgG3","DT_IgG4","DT_ADCD","DT_ADNP","DT_ADCP","DT_FcgR2a","DT_FcgR3b",
      "TT_IgG","TT_IgG1","TT_IgG2","TT_IgG3","TT_IgG4","TT_ADCD","TT_ADNP","TT_ADCP","TT_FcgR2a","TT_FcgR3b","TT_Bmem",
      "PT_IgG","PT_IgG1","PT_IgG2","PT_IgG3","PT_IgG4","PT_ADCD","PT_ADNP","PT_ADCP","PT_FcgR2a","PT_FcgR3b","PT_Bmem",
      "FHA_IgG","FHA_IgG1","FHA_IgG2","FHA_IgG3","FHA_IgG4","FHA_ADCD","FHA_ADNP","FHA_ADCP","FHA_FcgR2a","FHA_FcgR3b","FHA_Bmem",
      "PRN_IgG","PRN_IgG1","PRN_IgG2","PRN_IgG3","PRN_IgG4","PRN_ADCD","PRN_ADNP","PRN_ADCP","PRN_FcgR2a","PRN_FcgR3b",
      "WHOLE_WT_IgG", "WHOLE_PTNA", "WHOLE_SBA")


### select covariates, set time points for difference and fold-change, name variables
covariates <- c("birth_weight", "infant_sex", "delivery_mode", "parity", "gravidity",
                "gestational_age_vaccination", "vaccine_birth_interval_days",
                "gestational_age_birth", "maternal_age", "log_igg_t0")

covariates_named_list <- list(
  "Arm Name"                                              = "arm_name",
  "Birth Weight"                                          = "birth_weight",
  "Infant Sex"                                            = "infant_sex",
  "Delivery Mode"                                         = "delivery_mode",
  "Parity"                                                = "parity",
  # "Gravidity"                                           = "gravidity",
  "Gestational Age at Vaccine"                            = "gestational_age_vaccination",
  "Infant Interval Between Vaccination and Delivery Days" = "vaccine_birth_interval_days",
  "Gestational Age at Delivery"                           = "gestational_age_birth",
  "Maternal Age at Vaccination"                           = "maternal_age",
  "Maternal BMI at Vaccination"                           = "maternal_bmi",
  "Maternal hemoglobin"                                   = "maternal_Hb"
  # "Baseline IgG"     = "log_igg_t0",
  # "Baseline DT IgG"  = "DT_IgGt0",
  # "Baseline FHA IgG" = "FHA_IgGt0",
  # "Baseline PRN IgG" = "PRN_IgGt0",
  # "Baseline PT IgG"  = "PT_IgGt0",
  # "Baseline TT IgG"  = "TT_IgGt0"
)

covar_list <- c("arm_name","infant_sex", "delivery_mode", "parity", "gravidity")
selected_covar <- 1

covar_cont_list <- c("birth_weight",
                     "gestational_age_vaccination",
                     "vaccine_birth_interval_days",
                     "gestational_age_birth",
                     "maternal_age", "log_igg_t0")
selected_covar_cont <- 1

covar_cat_list <- c("arm_name","infant_sex", "delivery_mode", "parity", "gravidity")
selected_covar_cat <- 1

igglist <- c("DT_IgG", "FHA_IgG", "PRN_IgG", "PT_IgG", "TT_IgG")

# t0 <- "M00"
t0 <- "M02"
t1 <- "M05"
t2 <- "M09"
visit_list <- c(t0, t1, t2)
# visit_list <- c(t0,t1,t2,t3)

visit_label <- function(visit = "vaccinated") {
  case_when(
    visit == "M02" ~ "month 2",
    # visit == "M00" ~ "birth",
    visit == "M05" ~ "month 5",
    visit == "M09" ~ "month 9",
    TRUE ~ NA_character_
  )
}

visit_number <- function(visit = "vaccinated") {
  case_when(
    visit == "M02" ~ 2,
    # visit == "M00" ~ 0,
    visit == "M05" ~ 5,
    visit == "M09" ~ 9
  )
}

# Helper function to convert p-values to significance labels
pvalue_to_label <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")
}

antigen_list   <- c("PT", "FHA", "PRN", "DT", "TT", "WHOLE")
analytes_order <- c("IgG","IgG3","IgG1","IgG2","IgG4","FcgR2a","FcgR3b","ADCD","ADCP","ADNP","WT_IgG","PTNA","SBA")
arms    <- c("TdaP", "TT")
refarms <- c("TT",   "TdaP")

# ============================================================
# PTNA Network Analysis — Stratified by Infant Arm
# Parallel to SBA analysis in response_matpm_SBA_model_v1.R
# — outcome switched to WHOLE_PTNA; all predictors are
# pertussis-antigen IgG total and IgG subclasses
# (PT, FHA, PRN × IgG1–4)
# ============================================================

# ---- 0. Setup: infant arm loop scaffold --------------------
load(file = here::here("./data/c_set.RData"))
data <- ebaa_extra
data$arm_name <- data$maternal_arm
data$arm_name <- fct_recode(data$arm_name, "TT" = "TT", "TdaP" = "TdaP")
data$feature  <- factor(paste(data$antigen, data$analyte, sep = "_"))

table(data$maternal_arm)
table(data$arm_name)
table(data$infant_arm)
table(data$visit_name)
table(data$feature)

data <- data[data$feature %in% include_feature,
             c("visit_name","subject_accession","antigen","arm_name","analyte",
               "infant_arm","feature","value_reported","log_assay_value")]
data$feature <- droplevels(data$feature)
table(data$arm_name)
table(data$visit_name)
table(data$feature)

ddata <- distinct(
  data[data$feature %in% include_feature,
       c("infant_arm","visit_name","subject_accession","antigen","arm_name",
         "analyte","feature","value_reported","log_assay_value")],
  subject_accession, visit_name, feature, .keep_all = TRUE
)

ddata$arm_name <- factor(ddata$arm_name)
table(ddata$arm_name)
ddata$arm_name <- relevel(ddata$arm_name, ref = as.character(refarms[1]))
table(ddata$arm_name)
data <- ddata
data <- data[data$arm_name   %in% arms,       ]
data <- data[data$visit_name %in% visit_list, ]

data$arm_name   <- droplevels(data$arm_name)
data$visit_name <- droplevels(data$visit_name)

data$study_accession <- ifelse(data$arm_name %in% arms, "SDY2818", NA)
table(data$feature)

data$visit_number <- visit_number(data$visit_name)
table(data$feature, data$visit_number)

infant_arms  <- c("aP", "wP")
iarm_colors  <- c("aP" = "dodgerblue2", "wP" = "forestgreen")

# Storage lists for cross-arm comparison at the end
networks_by_iarm    <- list()
data_wide_by_iarm   <- list()
model_full_by_iarm  <- list()
relimp_full_by_iarm <- list()
bridge_df_by_iarm   <- list()
ratio_model_by_iarm <- list()

# ---- 1. Feature selection ----------------------------------
# Outcome:    WHOLE_PTNA
# Predictors: PT, FHA, PRN total IgG and subclasses IgG1–IgG4

ptna_features <- c("WHOLE_PTNA",
                   "PT_IgG",  "PT_IgG1",  "PT_IgG2",  "PT_IgG3",  "PT_IgG4",
                   "PRN_IgG", "PRN_IgG1", "PRN_IgG2", "PRN_IgG3", "PRN_IgG4",
                   "FHA_IgG", "FHA_IgG1", "FHA_IgG2", "FHA_IgG3", "FHA_IgG4")

# Network nodes: PTNA + 12 pertussis-antigen subclasses
net_vars <- c("WHOLE_PTNA",
              "PT_IgG1",  "PT_IgG2",  "PT_IgG3",  "PT_IgG4",
              "FHA_IgG1", "FHA_IgG2", "FHA_IgG3", "FHA_IgG4",
              "PRN_IgG1", "PRN_IgG2", "PRN_IgG3", "PRN_IgG4")

node_names <- net_vars
groups_idx <- list(
  "PTNA" = which(node_names == "WHOLE_PTNA"),
  "PT"   = which(node_names %in% c("PT_IgG1",  "PT_IgG2",  "PT_IgG3",  "PT_IgG4")),
  "FHA"  = which(node_names %in% c("FHA_IgG1", "FHA_IgG2", "FHA_IgG3", "FHA_IgG4")),
  "PRN"  = which(node_names %in% c("PRN_IgG1", "PRN_IgG2", "PRN_IgG3", "PRN_IgG4"))
)
group_colors <- c("gold", "tomato", "steelblue", "seagreen")

communities <- c(1, rep(2,4), rep(3,4), rep(4,4))   # PTNA=1, PT=2, FHA=3, PRN=4

# ============================================================
# ---- 2. Main loop over infant arms -------------------------
# ============================================================

for (iarm in infant_arms) {

  cat("\n\n============================================================\n")
  cat(paste0("  INFANT ARM: ", iarm, "\n"))
  cat("============================================================\n\n")

  # ---- 2a. Subset and pivot to wide -----------------------

  data0259_iarm <- data[data$feature %in% ptna_features &
                          data$infant_arm == iarm, ]
  data0259_iarm$feature <- droplevels(data0259_iarm$feature)

  cat("--- Feature × Visit counts ---\n")
  print(table(data0259_iarm$feature, data0259_iarm$visit_number))

  data_wide_iarm <- data0259_iarm %>%
    dplyr::select(subject_accession, visit_name, arm_name,
                  infant_arm, feature, log_assay_value) %>%
    pivot_wider(
      id_cols     = c(subject_accession, visit_name, arm_name, infant_arm),
      names_from  = feature,
      values_from = log_assay_value
    )

  cat("--- Missingness ---\n")
  print(colSums(is.na(data_wide_iarm)))

  # ---- 2b. Pro-inflammatory / tolerogenic ratios ----------

  data_wide_iarm <- data_wide_iarm %>%
    mutate(
      PT_proinflam  = 10^PT_IgG1  + 10^PT_IgG3,
      PT_tolerogen  = 10^PT_IgG2  + 10^PT_IgG4,
      FHA_proinflam = 10^FHA_IgG1 + 10^FHA_IgG3,
      FHA_tolerogen = 10^FHA_IgG2 + 10^FHA_IgG4,
      PRN_proinflam = 10^PRN_IgG1 + 10^PRN_IgG3,
      PRN_tolerogen = 10^PRN_IgG2 + 10^PRN_IgG4,
      PT_ratio      = log10(PT_proinflam  / PT_tolerogen),
      FHA_ratio     = log10(FHA_proinflam / FHA_tolerogen),
      PRN_ratio     = log10(PRN_proinflam / PRN_tolerogen)
    )

  data_wide_iarm$visit_number <- visit_number(data_wide_iarm$visit_name)

  cat("--- Ratio summary ---\n")
  print(summary(data_wide_iarm[, c("PT_ratio","FHA_ratio","PRN_ratio","WHOLE_PTNA")]))

  # Store wide data for later cross-arm comparisons
  data_wide_by_iarm[[iarm]] <- data_wide_iarm

  # ---- 2c. Antigen-level relative importance --------------

  cat("\n--- Antigen-level model (total IgG → PTNA) ---\n")
  model_antigen_i <- lm(WHOLE_PTNA ~ PT_IgG + FHA_IgG + PRN_IgG,
                         data = data_wide_iarm)
  print(summary(model_antigen_i))
  print(calc.relimp(model_antigen_i, type = "lmg"))

  # ---- 2d. Subclass decomposition per antigen -------------

  cat("\n--- PT subclass model (→ PTNA) ---\n")
  model_PT_i <- lm(WHOLE_PTNA ~ PT_IgG1  + PT_IgG2  + PT_IgG3  + PT_IgG4,
                    data = data_wide_iarm)
  print(calc.relimp(model_PT_i, type = "lmg"))

  cat("\n--- FHA subclass model (→ PTNA) ---\n")
  model_FHA_i <- lm(WHOLE_PTNA ~ FHA_IgG1 + FHA_IgG2 + FHA_IgG3 + FHA_IgG4,
                     data = data_wide_iarm)
  print(calc.relimp(model_FHA_i, type = "lmg"))

  cat("\n--- PRN subclass model (→ PTNA) ---\n")
  model_PRN_i <- lm(WHOLE_PTNA ~ PRN_IgG1 + PRN_IgG2 + PRN_IgG3 + PRN_IgG4,
                     data = data_wide_iarm)
  print(calc.relimp(model_PRN_i, type = "lmg"))

  # ---- 2e. Full 12-subclass model -------------------------

  cat("\n--- Full 12-subclass model (→ PTNA) ---\n")
  model_full_i <- lm(
    WHOLE_PTNA ~ PT_IgG1  + PT_IgG2  + PT_IgG3  + PT_IgG4  +
                 FHA_IgG1 + FHA_IgG2 + FHA_IgG3 + FHA_IgG4 +
                 PRN_IgG1 + PRN_IgG2 + PRN_IgG3 + PRN_IgG4,
    data = data_wide_iarm)
  print(summary(model_full_i))

  relimp_i <- calc.relimp(model_full_i, type = "lmg")
  print(relimp_i)

  model_full_by_iarm[[iarm]]  <- model_full_i
  relimp_full_by_iarm[[iarm]] <- relimp_i

  # ---- 2f. Pro-inflammatory/tolerogenic ratio models ------

  cat("\n--- Ratio model (all timepoints combined, → PTNA) ---\n")
  model_ratio_i <- lm(WHOLE_PTNA ~ PT_ratio + FHA_ratio + PRN_ratio,
                       data = data_wide_iarm)
  print(summary(model_ratio_i))
  ratio_model_by_iarm[[iarm]] <- model_ratio_i

  cat("\n--- Ratio model by timepoint ---\n")
  ratio_by_visit <- data_wide_iarm %>%
    group_by(visit_number) %>%
    do(tidy(lm(WHOLE_PTNA ~ PT_ratio + FHA_ratio + PRN_ratio, data = .))) %>%
    dplyr::filter(term != "(Intercept)")
  print(ratio_by_visit, n = 30)

  # ---- 2g. Antigen × subclass interaction (mixed model) ---

  data_long_sub_i <- data_wide_iarm %>%
    dplyr::select(subject_accession, visit_number, arm_name,
                  WHOLE_PTNA, matches("IgG[1-4]$")) %>%
    pivot_longer(cols          = matches("IgG[1-4]$"),
                 names_to      = c("antigen", "subclass"),
                 names_pattern = "(PT|FHA|PRN)_(IgG[1-4])",
                 values_to     = "subclass_value") %>%
    drop_na()

  cat("\n--- Antigen × Subclass × Value interaction (lmer, → PTNA) ---\n")
  model_int_i <- lmer(
    WHOLE_PTNA ~ antigen * subclass * subclass_value + (1 | subject_accession),
    data = data_long_sub_i)
  print(anova(model_int_i))

  # ---- 2h. GGM network at Visit 5 -------------------------

  data_net_i <- data_wide_iarm %>%
    filter(visit_number == 5) %>%
    dplyr::select(all_of(net_vars)) %>%
    drop_na()

  cat(paste0("\n--- GGM network (Visit 5, n = ", nrow(data_net_i), ") ---\n"))

  if (nrow(data_net_i) >= 50) {
    network_i <- estimateNetwork(data_net_i, default = "EBICglasso",
                                  threshold = TRUE)
    networks_by_iarm[[iarm]] <- network_i

    plot(network_i,
         groups       = groups_idx,
         layout       = "spring",
         labels       = node_names,
         color        = group_colors,
         title        = paste0("Visit 5 — ", iarm, " infants (PTNA network)"),
         legend       = TRUE,
         legend.cex   = 0.4,
         vsize        = 8,
         border.width = 1.5,
         edge.labels  = FALSE)

    centralityPlot(network_i, include = c("Strength", "ExpectedInfluence"))

    # ---- 2i. Bridge centrality --------------------------

    bridge_i <- bridge(network_i$graph, communities = communities)

    bridge_df_i <- data.frame(
      Variable        = names(bridge_i$`Bridge Strength`),
      Bridge_Strength = bridge_i$`Bridge Strength`,
      Bridge_EI       = bridge_i$`Bridge Expected Influence (1-step)`,
      Community       = communities,
      infant_arm      = iarm,
      row.names       = NULL
    ) %>% arrange(desc(Bridge_Strength))

    bridge_df_by_iarm[[iarm]] <- bridge_df_i

    cat("\n--- Bridge centrality ---\n")
    print(bridge_df_i)

    # ---- 2j. Corrplot heatmap ---------------------------

    adj_i <- network_i$graph
    rownames(adj_i) <- colnames(adj_i) <- net_vars
    corrplot(adj_i, method = "color", type = "lower",
             tl.col = "black", tl.cex = 0.7,
             col    = colorRampPalette(c("blue","white","red"))(200),
             title  = paste0("Partial Correlation — Visit 5 — ", iarm, " (PTNA network)"),
             mar    = c(0,0,2,0), order = "hclust", addrect = 4)

    # ---- 2k. Bootstrap stability (optional, slow) -------
    # Uncomment if you want stability estimates per arm
    # boot_case_i <- bootnet(network_i, nBoots = 1000, type = "case",
    #                        statistics = c("strength","expectedInfluence"))
    # cat(paste0("\n--- CS coefficient (", iarm, ") ---\n"))
    # print(corStability(boot_case_i))

  } else {
    cat(paste0("Insufficient n for network estimation in ", iarm,
               " at visit 5 (n = ", nrow(data_net_i), ")\n"))
  }

  # ---- 2l. Networks at all four timepoints ----------------

  par(mfrow = c(2, 2))
  for (v in c(0, 2, 5, 9)) {
    data_tmp_i <- data_wide_iarm %>%
      filter(visit_number == v) %>%
      dplyr::select(all_of(net_vars)) %>%
      drop_na()

    if (nrow(data_tmp_i) >= 30) {
      net_tmp_i <- estimateNetwork(data_tmp_i, default = "EBICglasso",
                                    threshold = TRUE)
      plot(net_tmp_i, groups = groups_idx, color = group_colors,
           layout = "spring",
           title  = paste0("Visit ", v, " — ", iarm, " (PTNA)"),
           vsize  = 8, legend = FALSE)
    } else {
      plot.new()
      title(paste0("Visit ", v, " — ", iarm, " (n<30, skipped)"))
    }
  }
  par(mfrow = c(1, 1))

  # ---- 2m. NCT: timepoint comparison within infant arm ----

  data_v2_i <- data_wide_iarm %>%
    filter(visit_number == 2) %>%
    dplyr::select(all_of(net_vars)) %>% drop_na()
  data_v9_i <- data_wide_iarm %>%
    filter(visit_number == 9) %>%
    dplyr::select(all_of(net_vars)) %>% drop_na()

  cat(paste0("\n--- NCT: Visit 2 vs Visit 9 (", iarm, ") ---\n"))
  cat(paste0("n V2 = ", nrow(data_v2_i), "  n V9 = ", nrow(data_v9_i), "\n"))

  if (nrow(data_v2_i) >= 30 & nrow(data_v9_i) >= 30) {
    nct_time_i <- NCT(data_v2_i, data_v9_i, it = 1000,
                      test.edges      = TRUE,
                      test.centrality = TRUE,
                      centrality      = c("strength","expectedInfluence"))
    print(summary(nct_time_i))
  }

  # ---- 2n. NCT: maternal arm comparison within infant arm -

  data_TdaP_i <- data_wide_iarm %>%
    filter(visit_number == 5, arm_name == "TdaP") %>%
    dplyr::select(all_of(net_vars)) %>% drop_na()
  data_TT_i <- data_wide_iarm %>%
    filter(visit_number == 5, arm_name == "TT") %>%
    dplyr::select(all_of(net_vars)) %>% drop_na()

  cat(paste0("\n--- NCT: TdaP vs TT at Visit 5 (", iarm, ") ---\n"))
  cat(paste0("n TdaP = ", nrow(data_TdaP_i), "  n TT = ", nrow(data_TT_i), "\n"))

  if (nrow(data_TdaP_i) >= 30 & nrow(data_TT_i) >= 30) {
    nct_arms_i <- NCT(data_TdaP_i, data_TT_i, it = 1000,
                      test.edges      = TRUE,
                      test.centrality = TRUE,
                      centrality      = c("strength","expectedInfluence"))
    print(summary(nct_arms_i))
  }

} # end infant_arm loop


# ============================================================
# ---- 3. Cross-arm comparisons (aP vs wP) ------------------
# ============================================================

# ---- 3a. NCT: aP vs wP at Visit 5 -----------------------

cat("\n\n============================================================\n")
cat("  CROSS-ARM: aP vs wP network comparison at Visit 5 (PTNA)\n")
cat("============================================================\n\n")

data_net_aP <- data_wide_by_iarm[["aP"]] %>%
  filter(visit_number == 5) %>%
  dplyr::select(all_of(net_vars)) %>% drop_na()

data_net_wP <- data_wide_by_iarm[["wP"]] %>%
  filter(visit_number == 5) %>%
  dplyr::select(all_of(net_vars)) %>% drop_na()

cat(paste0("n aP = ", nrow(data_net_aP), "  n wP = ", nrow(data_net_wP), "\n"))

nct_iarm <- NCT(data_net_aP, data_net_wP, it = 1000,
                test.edges      = TRUE,
                test.centrality = TRUE,
                centrality      = c("strength","expectedInfluence"))
summary(nct_iarm)

# ---- 3b. Side-by-side network plots ----------------------

par(mfrow = c(1, 2))
for (iarm in infant_arms) {
  if (!is.null(networks_by_iarm[[iarm]])) {
    plot(networks_by_iarm[[iarm]],
         groups     = groups_idx, color = group_colors,
         layout     = "spring",   labels = node_names,
         title      = paste0("Visit 5 — ", iarm, " (PTNA)"),
         legend     = (iarm == "wP"),   # legend on right panel only
         legend.cex = 0.4, vsize = 8)
  }
}
par(mfrow = c(1, 1))

# ---- 3c. Bridge centrality comparison table --------------

cat("\n--- Bridge centrality comparison: aP vs wP (PTNA network) ---\n")
bridge_combined <- bind_rows(bridge_df_by_iarm) %>%
  dplyr::select(Variable, infant_arm, Bridge_Strength, Bridge_EI) %>%
  pivot_wider(names_from  = infant_arm,
              values_from = c(Bridge_Strength, Bridge_EI),
              names_sep   = "_") %>%
  mutate(
    delta_BS = Bridge_Strength_aP - Bridge_Strength_wP,
    delta_EI = Bridge_EI_aP       - Bridge_EI_wP
  ) %>%
  arrange(desc(abs(delta_BS)))
print(bridge_combined)

# ---- 3d. Relative importance comparison table ------------

cat("\n--- Relative importance comparison: aP vs wP (PTNA) ---\n")
relimp_aP <- data.frame(
  feature = names(relimp_full_by_iarm[["aP"]]@lmg),
  lmg_aP  = relimp_full_by_iarm[["aP"]]@lmg
)
relimp_wP <- data.frame(
  feature = names(relimp_full_by_iarm[["wP"]]@lmg),
  lmg_wP  = relimp_full_by_iarm[["wP"]]@lmg
)
relimp_compare <- merge(relimp_aP, relimp_wP, by = "feature") %>%
  mutate(delta_lmg = lmg_aP - lmg_wP) %>%
  arrange(desc(abs(delta_lmg)))
print(relimp_compare)

# ---- 3e. Ratio model coefficient comparison plot ---------

ratio_coef <- bind_rows(
  tidy(ratio_model_by_iarm[["aP"]]) %>%
    filter(term != "(Intercept)") %>% mutate(infant_arm = "aP"),
  tidy(ratio_model_by_iarm[["wP"]]) %>%
    filter(term != "(Intercept)") %>% mutate(infant_arm = "wP")
)

ggplot(ratio_coef, aes(x = term, y = estimate, fill = infant_arm,
                        ymin = estimate - 1.96 * std.error,
                        ymax = estimate + 1.96 * std.error)) +
  geom_col(position = position_dodge(0.6), width = 0.5) +
  geom_errorbar(position = position_dodge(0.6), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = iarm_colors) +
  labs(x     = "Pro-inflammatory/Tolerogenic Ratio",
       y     = "Regression Coefficient (→ PTNA)",
       title = "Ratio model: aP vs wP infants (outcome = PTNA)",
       fill  = "Infant arm") +
  theme_bw()

# ---- 3f. Interaction plot: antigen × subclass slopes -----
#          separately for aP and wP

data_long_combined <- bind_rows(
  data_wide_by_iarm[["aP"]] %>%
    dplyr::select(subject_accession, visit_number, arm_name,
                  WHOLE_PTNA, matches("IgG[1-4]$")) %>%
    mutate(infant_arm = "aP"),
  data_wide_by_iarm[["wP"]] %>%
    dplyr::select(subject_accession, visit_number, arm_name,
                  WHOLE_PTNA, matches("IgG[1-4]$")) %>%
    mutate(infant_arm = "wP")
) %>%
  pivot_longer(cols          = matches("IgG[1-4]$"),
               names_to      = c("antigen", "subclass"),
               names_pattern = "(PT|FHA|PRN)_(IgG[1-4])",
               values_to     = "subclass_value") %>%
  drop_na()

ggplot(data_long_combined,
       aes(x = subclass_value, y = WHOLE_PTNA,
           color = subclass, linetype = infant_arm)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  facet_grid(antigen ~ infant_arm) +
  scale_color_brewer(palette = "Set1") +
  labs(x        = "Log10 Subclass-Specific IgG",
       y        = "Log10 PTNA (Pertussis Toxin Neutralization Activity)",
       color    = "Subclass",
       linetype = "Infant arm",
       title    = "Antigen × Subclass → PTNA: aP vs wP infants") +
  theme_bw(base_size = 11)
