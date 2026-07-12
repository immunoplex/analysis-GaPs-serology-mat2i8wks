# ---
# title:    "Predicting Infant M02 PTNA from Maternal and Cord IgG"
# subtitle: "Six predictor blocks from maternal vaccination through cord transfer to M02"
# author:   "Scot Zens"
# note:     Set outcome_feature / outcome_label in Section 1 to switch
#           between WHOLE_PTNA (this script), WHOLE_SBA, or WHOLE_WT_IgG
# ---
#
# PREDICTOR BLOCKS
# ----------------
# Block 1 — Quantity:       Maternal total IgG + subclasses at P02 and P09
# Block 2 — Class switch:   Maternal pro-/tolerogenic ratio change P02 → P09
# Block 3 — Quantity:       IgG level change P09 (maternal) → M00 (cord)  [transfer]
# Block 4 — Class switch:   Ratio change P09 (maternal) → M00 (cord)       [transfer]
# Block 5 — Quantity:       IgG level change M00 (cord)  → M02 (infant)   [decay]
# Block 6 — Class switch:   Ratio change M00 (cord)      → M02 (infant)   [decay]
#
# OUTCOME:  WHOLE_PTNA at M02 (infant, before first vaccination)
# STRATA:   infant_arm (aP vs wP)
# ============================================================

library(data.table)
library(tidyverse)
library(here)
library(RColorBrewer)
library(broom)
library(relaimpo)

# ============================================================
# SECTION 1 — Parameters (only section that needs editing to
#             switch to SBA or WT_IgG)
# ============================================================

outcome_feature <- "WHOLE_PTNA"       # WHOLE_PTNA | WHOLE_SBA | WHOLE_WT_IgG
outcome_label   <- "PTNA"             # used in cat() and plot titles
outcome_y_label <- "Log10 PTNA (Pertussis Toxin Neutralization Activity)"

# Antigen-level total IgG predictor names
antigen_totals  <- c("PT_IgG",  "FHA_IgG",  "PRN_IgG")

# Subclass predictor names (without visit suffix)
pt_subclasses   <- c("PT_IgG1",  "PT_IgG2",  "PT_IgG3",  "PT_IgG4")
fha_subclasses  <- c("FHA_IgG1", "FHA_IgG2", "FHA_IgG3", "FHA_IgG4")
prn_subclasses  <- c("PRN_IgG1", "PRN_IgG2", "PRN_IgG3", "PRN_IgG4")
all_subclasses  <- c(pt_subclasses, fha_subclasses, prn_subclasses)

predictor_features <- c(antigen_totals, all_subclasses)

# Visits required
PRED_VISITS    <- c("P02", "P09", "M00", "M02")
OUTCOME_VISIT  <- "M02"

# Maternal arms and infant arms
arms        <- c("TdaP", "TT")
infant_arms <- c("aP", "wP")
iarm_colors <- c("aP" = "dodgerblue2", "wP" = "forestgreen")

# ============================================================
# SECTION 2 — Helper functions
# ============================================================

# ---- 2a. p-value labels ------------------------------------
pvalue_to_label <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

# ---- 2b. Ratio calculator (log10 pro-inflam / tolerogenic) --
ratio_col <- function(igG1, igG2, igG3, igG4) {
  log10((10^igG1 + 10^igG3) / (10^igG2 + 10^igG4))
}

# ---- 2c. Core model runner ---------------------------------
# Returns a list(model, relimp) and prints summaries.
run_lm <- function(outcome_col, predictors, data, label = "") {
  cat(paste0("\n--- ", label, " ---\n"))
  f <- reformulate(predictors, response = outcome_col)
  m <- lm(f, data = data, na.action = na.omit)
  print(summary(m))
  ri <- tryCatch(
    calc.relimp(m, type = "lmg"),
    error = function(e) {
      cat("  [relimp skipped:", conditionMessage(e), "]\n"); NULL
    }
  )
  if (!is.null(ri)) print(ri)
  invisible(list(model = m, relimp = ri, label = label,
                 R2     = summary(m)$r.squared,
                 adj_R2 = summary(m)$adj.r.squared,
                 n_obs  = nobs(m)))
}

# ---- 2d. Collect R² from a run_lm result into a data frame row --
r2_row <- function(res, infant_arm, block, model_name) {
  data.frame(
    infant_arm = infant_arm,
    block      = block,
    model      = model_name,
    n_obs      = res$n_obs,
    R2         = round(res$R2,     4),
    adj_R2     = round(res$adj_R2, 4),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# SECTION 3 — Load data and pivot to wide (one row per subject)
# ============================================================

load(file = here::here("./data/c_set.RData"))
data_raw           <- ebaa_extra
data_raw$arm_name  <- data_raw$maternal_arm
data_raw$arm_name  <- fct_recode(data_raw$arm_name, "TT" = "TT", "TdaP" = "TdaP")
data_raw$feature   <- factor(paste(data_raw$antigen, data_raw$analyte, sep = "_"))

# Diagnostic overview
cat("=== Full dataset visit counts ===\n")
print(table(data_raw$visit_name))

cat(paste0("\n=== ", outcome_feature, " counts by visit ===\n"))
print(table(data_raw$visit_name[data_raw$feature == outcome_feature]))

# Filter to relevant features and visits
all_features_needed <- c(predictor_features, outcome_feature)

data_long <- data_raw %>%
  dplyr::filter(
    feature    %in% all_features_needed,
    visit_name %in% PRED_VISITS,
    arm_name   %in% arms
  ) %>%
  dplyr::select(subject_accession, visit_name, arm_name, infant_arm,
                feature, log_assay_value) %>%
  distinct(subject_accession, visit_name, feature, .keep_all = TRUE)

cat("\n=== Feature × Visit counts (after filtering) ===\n")
print(table(data_long$feature, data_long$visit_name))

# Pivot wide: columns = {feature}_{visit}
# e.g. PT_IgG_P02, PT_IgG1_P09, WHOLE_PTNA_M02
data_wide <- data_long %>%
  pivot_wider(
    id_cols     = c(subject_accession, arm_name, infant_arm),
    names_from  = c(feature, visit_name),
    values_from = log_assay_value,
    names_sep   = "_"
  )

outcome_col <- paste0(outcome_feature, "_", OUTCOME_VISIT)   # "WHOLE_PTNA_M02"

cat("\n=== Wide data dimensions ===\n")
cat(nrow(data_wide), "subjects ×", ncol(data_wide), "columns\n")
cat("Outcome column:", outcome_col, "—",
    sum(!is.na(data_wide[[outcome_col]])), "non-missing\n")

# ============================================================
# SECTION 4 — Compute derived quantities
# ============================================================

data_wide <- data_wide %>%
  mutate(

    # ---- Ratios at each visit --------------------------------

    # P02  (maternal, pre-vaccination)
    PT_ratio_P02  = ratio_col(PT_IgG1_P02,  PT_IgG2_P02,  PT_IgG3_P02,  PT_IgG4_P02),
    FHA_ratio_P02 = ratio_col(FHA_IgG1_P02, FHA_IgG2_P02, FHA_IgG3_P02, FHA_IgG4_P02),
    PRN_ratio_P02 = ratio_col(PRN_IgG1_P02, PRN_IgG2_P02, PRN_IgG3_P02, PRN_IgG4_P02),

    # P09  (maternal, at birth)
    PT_ratio_P09  = ratio_col(PT_IgG1_P09,  PT_IgG2_P09,  PT_IgG3_P09,  PT_IgG4_P09),
    FHA_ratio_P09 = ratio_col(FHA_IgG1_P09, FHA_IgG2_P09, FHA_IgG3_P09, FHA_IgG4_P09),
    PRN_ratio_P09 = ratio_col(PRN_IgG1_P09, PRN_IgG2_P09, PRN_IgG3_P09, PRN_IgG4_P09),

    # M00  (cordblood at birth)
    PT_ratio_M00  = ratio_col(PT_IgG1_M00,  PT_IgG2_M00,  PT_IgG3_M00,  PT_IgG4_M00),
    FHA_ratio_M00 = ratio_col(FHA_IgG1_M00, FHA_IgG2_M00, FHA_IgG3_M00, FHA_IgG4_M00),
    PRN_ratio_M00 = ratio_col(PRN_IgG1_M00, PRN_IgG2_M00, PRN_IgG3_M00, PRN_IgG4_M00),

    # M02  (infant at 2 months — needed for Block 6)
    PT_ratio_M02  = ratio_col(PT_IgG1_M02,  PT_IgG2_M02,  PT_IgG3_M02,  PT_IgG4_M02),
    FHA_ratio_M02 = ratio_col(FHA_IgG1_M02, FHA_IgG2_M02, FHA_IgG3_M02, FHA_IgG4_M02),
    PRN_ratio_M02 = ratio_col(PRN_IgG1_M02, PRN_IgG2_M02, PRN_IgG3_M02, PRN_IgG4_M02),

    # ---- Block 2: Δ ratio P02 → P09  (class switching during pregnancy) ----
    delta_PT_ratio_P02_P09  = PT_ratio_P09  - PT_ratio_P02,
    delta_FHA_ratio_P02_P09 = FHA_ratio_P09 - FHA_ratio_P02,
    delta_PRN_ratio_P02_P09 = PRN_ratio_P09 - PRN_ratio_P02,

    # ---- Block 3: Δ quantity P09 → M00  (maternal-to-cord transfer) --------
    delta_PT_IgG_P09_M00   = PT_IgG_M00   - PT_IgG_P09,
    delta_PT_IgG1_P09_M00  = PT_IgG1_M00  - PT_IgG1_P09,
    delta_PT_IgG2_P09_M00  = PT_IgG2_M00  - PT_IgG2_P09,
    delta_PT_IgG3_P09_M00  = PT_IgG3_M00  - PT_IgG3_P09,
    delta_PT_IgG4_P09_M00  = PT_IgG4_M00  - PT_IgG4_P09,
    delta_FHA_IgG_P09_M00  = FHA_IgG_M00  - FHA_IgG_P09,
    delta_FHA_IgG1_P09_M00 = FHA_IgG1_M00 - FHA_IgG1_P09,
    delta_FHA_IgG2_P09_M00 = FHA_IgG2_M00 - FHA_IgG2_P09,
    delta_FHA_IgG3_P09_M00 = FHA_IgG3_M00 - FHA_IgG3_P09,
    delta_FHA_IgG4_P09_M00 = FHA_IgG4_M00 - FHA_IgG4_P09,
    delta_PRN_IgG_P09_M00  = PRN_IgG_M00  - PRN_IgG_P09,
    delta_PRN_IgG1_P09_M00 = PRN_IgG1_M00 - PRN_IgG1_P09,
    delta_PRN_IgG2_P09_M00 = PRN_IgG2_M00 - PRN_IgG2_P09,
    delta_PRN_IgG3_P09_M00 = PRN_IgG3_M00 - PRN_IgG3_P09,
    delta_PRN_IgG4_P09_M00 = PRN_IgG4_M00 - PRN_IgG4_P09,

    # ---- Block 4: Δ ratio P09 → M00  (subclass-selective transfer) ---------
    delta_PT_ratio_P09_M00  = PT_ratio_M00  - PT_ratio_P09,
    delta_FHA_ratio_P09_M00 = FHA_ratio_M00 - FHA_ratio_P09,
    delta_PRN_ratio_P09_M00 = PRN_ratio_M00 - PRN_ratio_P09,

    # ---- Block 5: Δ quantity M00 → M02  (cord to 2 months, decay) ----------
    delta_PT_IgG_M00_M02   = PT_IgG_M02   - PT_IgG_M00,
    delta_PT_IgG1_M00_M02  = PT_IgG1_M02  - PT_IgG1_M00,
    delta_PT_IgG2_M00_M02  = PT_IgG2_M02  - PT_IgG2_M00,
    delta_PT_IgG3_M00_M02  = PT_IgG3_M02  - PT_IgG3_M00,
    delta_PT_IgG4_M00_M02  = PT_IgG4_M02  - PT_IgG4_M00,
    delta_FHA_IgG_M00_M02  = FHA_IgG_M02  - FHA_IgG_M00,
    delta_FHA_IgG1_M00_M02 = FHA_IgG1_M02 - FHA_IgG1_M00,
    delta_FHA_IgG2_M00_M02 = FHA_IgG2_M02 - FHA_IgG2_M00,
    delta_FHA_IgG3_M00_M02 = FHA_IgG3_M02 - FHA_IgG3_M00,
    delta_FHA_IgG4_M00_M02 = FHA_IgG4_M02 - FHA_IgG4_M00,
    delta_PRN_IgG_M00_M02  = PRN_IgG_M02  - PRN_IgG_M00,
    delta_PRN_IgG1_M00_M02 = PRN_IgG1_M02 - PRN_IgG1_M00,
    delta_PRN_IgG2_M00_M02 = PRN_IgG2_M02 - PRN_IgG2_M00,
    delta_PRN_IgG3_M00_M02 = PRN_IgG3_M02 - PRN_IgG3_M00,
    delta_PRN_IgG4_M00_M02 = PRN_IgG4_M02 - PRN_IgG4_M00,

    # ---- Block 6: Δ ratio M00 → M02  (class-selective decay) ---------------
    delta_PT_ratio_M00_M02  = PT_ratio_M02  - PT_ratio_M00,
    delta_FHA_ratio_M00_M02 = FHA_ratio_M02 - FHA_ratio_M00,
    delta_PRN_ratio_M00_M02 = PRN_ratio_M02 - PRN_ratio_M00
  )

cat("\n=== Derived variable missingness summary ===\n")
derived_cols <- c(
  "PT_ratio_P02", "PT_ratio_P09", "PT_ratio_M00", "PT_ratio_M02",
  "delta_PT_ratio_P02_P09", "delta_PT_ratio_P09_M00", "delta_PT_ratio_M00_M02",
  "delta_PT_IgG_P09_M00",  "delta_PT_IgG_M00_M02",
  outcome_col
)
print(colSums(is.na(data_wide[, intersect(derived_cols, names(data_wide))])))

# ============================================================
# SECTION 5 — Storage for cross-arm comparison
# ============================================================

r2_all     <- list()      # collects r2_row() data frames
ratio_coef <- list()      # tidy() from each ratio model, by block × arm

# Per-arm model storage (for cross-arm relimp comparison)
relimp_b1_full_by_iarm <- list()  # Block 1 full-12 at P09
relimp_b3_full_by_iarm <- list()  # Block 3 full-12 transfer
relimp_b5_full_by_iarm <- list()  # Block 5 full-12 decay
ratio_models_by_block_arm <- list()  # ratio model objects keyed "B2_aP" etc.

# ============================================================
# SECTION 6 — Main loop over infant arms
# ============================================================

for (iarm in infant_arms) {

  cat("\n\n============================================================\n")
  cat(paste0("  INFANT ARM: ", iarm, "\n"))
  cat("============================================================\n\n")

  d <- data_wide %>% dplyr::filter(infant_arm == iarm)

  cat("n subjects =", nrow(d), "| outcome non-missing =",
      sum(!is.na(d[[outcome_col]])), "\n")

  # ----------------------------------------------------------
  # BLOCK 1 — Maternal IgG quantity at P02 and P09
  # Biological question: Do maternal antibody levels before or
  # at the time of delivery predict infant PTNA at M02?
  # ----------------------------------------------------------

  cat("\n\n===== BLOCK 1: Maternal IgG Quantity (P02 and P09) =====\n")
  cat("Missingness in Block 1 predictors:\n")
  b1_cols <- c(paste0(antigen_totals, "_P02"), paste0(antigen_totals, "_P09"),
               paste0(all_subclasses,  "_P02"), paste0(all_subclasses,  "_P09"))
  b1_cols_present <- intersect(b1_cols, names(d))
  print(colSums(is.na(d[, c(b1_cols_present, outcome_col)])))

  # --- 1a. Antigen totals at P02 ---
  r1a <- run_lm(outcome_col,
                paste0(antigen_totals, "_P02"), d,
                paste0("Block 1a — Antigen totals at P02 → ", outcome_label, "_M02"))
  r2_all[[paste0("B1a_", iarm)]] <- r2_row(r1a, iarm, "B1_quantity", "antigen_P02")

  # --- 1b. Antigen totals at P09 ---
  r1b <- run_lm(outcome_col,
                paste0(antigen_totals, "_P09"), d,
                paste0("Block 1b — Antigen totals at P09 → ", outcome_label, "_M02"))
  r2_all[[paste0("B1b_", iarm)]] <- r2_row(r1b, iarm, "B1_quantity", "antigen_P09")

  # --- 1c. Antigen totals: P02 + P09 combined ---
  r1c <- run_lm(outcome_col,
                c(paste0(antigen_totals, "_P02"), paste0(antigen_totals, "_P09")), d,
                paste0("Block 1c — Antigen totals P02+P09 combined → ", outcome_label, "_M02"))
  r2_all[[paste0("B1c_", iarm)]] <- r2_row(r1c, iarm, "B1_quantity", "antigen_P02+P09")

  # --- 1d–1f. PT subclasses at P02 and P09 ---
  r1d <- run_lm(outcome_col, paste0(pt_subclasses, "_P02"), d,
                "Block 1d — PT subclasses at P02")
  r2_all[[paste0("B1d_", iarm)]] <- r2_row(r1d, iarm, "B1_quantity", "PT_subclass_P02")

  r1e <- run_lm(outcome_col, paste0(pt_subclasses, "_P09"), d,
                "Block 1e — PT subclasses at P09")
  r2_all[[paste0("B1e_", iarm)]] <- r2_row(r1e, iarm, "B1_quantity", "PT_subclass_P09")

  # --- 1g–1h. FHA subclasses at P02 and P09 ---
  r1g <- run_lm(outcome_col, paste0(fha_subclasses, "_P02"), d,
                "Block 1g — FHA subclasses at P02")
  r2_all[[paste0("B1g_", iarm)]] <- r2_row(r1g, iarm, "B1_quantity", "FHA_subclass_P02")

  r1h <- run_lm(outcome_col, paste0(fha_subclasses, "_P09"), d,
                "Block 1h — FHA subclasses at P09")
  r2_all[[paste0("B1h_", iarm)]] <- r2_row(r1h, iarm, "B1_quantity", "FHA_subclass_P09")

  # --- 1i–1j. PRN subclasses at P02 and P09 ---
  r1i <- run_lm(outcome_col, paste0(prn_subclasses, "_P02"), d,
                "Block 1i — PRN subclasses at P02")
  r2_all[[paste0("B1i_", iarm)]] <- r2_row(r1i, iarm, "B1_quantity", "PRN_subclass_P02")

  r1j <- run_lm(outcome_col, paste0(prn_subclasses, "_P09"), d,
                "Block 1j — PRN subclasses at P09")
  r2_all[[paste0("B1j_", iarm)]] <- r2_row(r1j, iarm, "B1_quantity", "PRN_subclass_P09")

  # --- 1k. Full 12-subclass model at P02 ---
  r1k <- run_lm(outcome_col, paste0(all_subclasses, "_P02"), d,
                "Block 1k — Full 12-subclass model at P02")
  r2_all[[paste0("B1k_", iarm)]] <- r2_row(r1k, iarm, "B1_quantity", "full12_P02")

  # --- 1l. Full 12-subclass model at P09 ---
  r1l <- run_lm(outcome_col, paste0(all_subclasses, "_P09"), d,
                "Block 1l — Full 12-subclass model at P09")
  r2_all[[paste0("B1l_", iarm)]] <- r2_row(r1l, iarm, "B1_quantity", "full12_P09")
  relimp_b1_full_by_iarm[[iarm]] <- r1l$relimp   # store P09 full model relimp

  # --- 1m. Full 24-predictor model: all subclasses at both P02 and P09 ---
  r1m <- run_lm(outcome_col,
                c(paste0(all_subclasses, "_P02"), paste0(all_subclasses, "_P09")), d,
                "Block 1m — Full 24-predictor model (P02 + P09 subclasses)")
  r2_all[[paste0("B1m_", iarm)]] <- r2_row(r1m, iarm, "B1_quantity", "full24_P02+P09")

  # ----------------------------------------------------------
  # BLOCK 2 — Maternal IgG ratio change P02 → P09
  # Biological question: Does the SHIFT in pro-/tolerogenic
  # subclass balance during pregnancy predict infant PTNA?
  # Positive delta_ratio = shift toward pro-inflammatory
  # ----------------------------------------------------------

  cat("\n\n===== BLOCK 2: Ratio Change P02 → P09 (Class Switching in Pregnancy) =====\n")
  cat("Missingness in Block 2 predictors:\n")
  b2_cols <- c("delta_PT_ratio_P02_P09", "delta_FHA_ratio_P02_P09",
               "delta_PRN_ratio_P02_P09", outcome_col)
  print(colSums(is.na(d[, intersect(b2_cols, names(d))])))

  cat("\n--- Block 2 ratio delta summaries ---\n")
  print(summary(d[, intersect(
    c("delta_PT_ratio_P02_P09","delta_FHA_ratio_P02_P09","delta_PRN_ratio_P02_P09"),
    names(d))]))

  r2b <- run_lm(outcome_col,
                c("delta_PT_ratio_P02_P09", "delta_FHA_ratio_P02_P09",
                  "delta_PRN_ratio_P02_P09"),
                d, paste0("Block 2 — Δ ratio P02→P09 → ", outcome_label, "_M02"))
  r2_all[[paste0("B2_", iarm)]] <- r2_row(r2b, iarm, "B2_ratio_switch", "delta_ratio_P02_P09")
  ratio_models_by_block_arm[[paste0("B2_", iarm)]] <- r2b$model

  cat("\n--- Block 2: by-antigen ratio delta models ---\n")
  for (ant in c("PT","FHA","PRN")) {
    dc <- paste0("delta_", ant, "_ratio_P02_P09")
    if (dc %in% names(d)) {
      cat(paste0("  ", ant, ": β = "))
      m_tmp <- lm(reformulate(dc, response = outcome_col), data = d, na.action = na.omit)
      s_tmp <- summary(m_tmp)
      b_tmp <- coef(s_tmp)
      if (nrow(b_tmp) >= 2) {
        cat(round(b_tmp[2,1], 4), " p =", round(b_tmp[2,4], 4),
            pvalue_to_label(b_tmp[2,4]), "\n")
      }
    }
  }

  # ----------------------------------------------------------
  # BLOCK 3 — IgG quantity change P09 (maternal) → M00 (cord)
  # Biological question: Does the amount of antibody transferred
  # across the placenta predict infant PTNA at M02?
  # Positive delta = cord higher than maternal (efficient transfer)
  # ----------------------------------------------------------

  cat("\n\n===== BLOCK 3: Quantity Change P09 → M00 (Placental Transfer) =====\n")
  b3_antigen  <- paste0("delta_", antigen_totals, "_P09_M00")
  b3_pt       <- paste0("delta_", pt_subclasses,  "_P09_M00")
  b3_fha      <- paste0("delta_", fha_subclasses, "_P09_M00")
  b3_prn      <- paste0("delta_", prn_subclasses, "_P09_M00")
  b3_full     <- c(b3_pt, b3_fha, b3_prn)

  cat("Missingness in Block 3 predictors:\n")
  print(colSums(is.na(d[, intersect(c(b3_full, outcome_col), names(d))])))

  # --- 3a. Antigen totals ---
  r3a <- run_lm(outcome_col, b3_antigen, d,
                "Block 3a — Δ antigen totals P09→M00 → PTNA_M02")
  r2_all[[paste0("B3a_", iarm)]] <- r2_row(r3a, iarm, "B3_transfer_qty", "antigen_delta")

  # --- 3b–3d. Per-antigen subclass models ---
  r3b <- run_lm(outcome_col, b3_pt,  d, "Block 3b — Δ PT subclasses P09→M00")
  r3c <- run_lm(outcome_col, b3_fha, d, "Block 3c — Δ FHA subclasses P09→M00")
  r3d <- run_lm(outcome_col, b3_prn, d, "Block 3d — Δ PRN subclasses P09→M00")
  r2_all[[paste0("B3b_", iarm)]] <- r2_row(r3b, iarm, "B3_transfer_qty", "PT_sub_delta")
  r2_all[[paste0("B3c_", iarm)]] <- r2_row(r3c, iarm, "B3_transfer_qty", "FHA_sub_delta")
  r2_all[[paste0("B3d_", iarm)]] <- r2_row(r3d, iarm, "B3_transfer_qty", "PRN_sub_delta")

  # --- 3e. Full 12-subclass delta model ---
  r3e <- run_lm(outcome_col, b3_full, d,
                "Block 3e — Full 12-subclass Δ P09→M00")
  r2_all[[paste0("B3e_", iarm)]] <- r2_row(r3e, iarm, "B3_transfer_qty", "full12_delta")
  relimp_b3_full_by_iarm[[iarm]] <- r3e$relimp

  # ----------------------------------------------------------
  # BLOCK 4 — Ratio change P09 (maternal) → M00 (cord)
  # Biological question: Is there subclass-selective placental
  # transfer that predicts infant PTNA at M02?
  # ----------------------------------------------------------

  cat("\n\n===== BLOCK 4: Ratio Change P09 → M00 (Class-Selective Transfer) =====\n")
  b4_cols <- c("delta_PT_ratio_P09_M00", "delta_FHA_ratio_P09_M00",
               "delta_PRN_ratio_P09_M00")
  cat("Missingness in Block 4 predictors:\n")
  print(colSums(is.na(d[, intersect(c(b4_cols, outcome_col), names(d))])))

  cat("\n--- Block 4 ratio delta summaries ---\n")
  print(summary(d[, intersect(b4_cols, names(d))]))

  r4 <- run_lm(outcome_col, b4_cols, d,
               paste0("Block 4 — Δ ratio P09→M00 (transfer) → ", outcome_label, "_M02"))
  r2_all[[paste0("B4_", iarm)]] <- r2_row(r4, iarm, "B4_transfer_ratio", "delta_ratio_P09_M00")
  ratio_models_by_block_arm[[paste0("B4_", iarm)]] <- r4$model

  cat("\n--- Block 4: by-antigen ratio delta models ---\n")
  for (ant in c("PT","FHA","PRN")) {
    dc <- paste0("delta_", ant, "_ratio_P09_M00")
    if (dc %in% names(d)) {
      m_tmp <- lm(reformulate(dc, response = outcome_col), data = d, na.action = na.omit)
      s_tmp <- summary(m_tmp); b_tmp <- coef(s_tmp)
      if (nrow(b_tmp) >= 2)
        cat(paste0("  ", ant, ": β=", round(b_tmp[2,1],4),
                   "  p=", round(b_tmp[2,4],4),
                   " ", pvalue_to_label(b_tmp[2,4]), "\n"))
    }
  }

  # ----------------------------------------------------------
  # BLOCK 5 — IgG quantity change M00 (cord) → M02 (infant)
  # Biological question: Do the decay kinetics of maternal
  # antibodies from birth to 2 months predict infant PTNA?
  # Negative deltas expected (antibody waning after birth)
  # ----------------------------------------------------------

  cat("\n\n===== BLOCK 5: Quantity Change M00 → M02 (Cord to 2 Months, Decay) =====\n")
  b5_antigen <- paste0("delta_", antigen_totals, "_M00_M02")
  b5_pt      <- paste0("delta_", pt_subclasses,  "_M00_M02")
  b5_fha     <- paste0("delta_", fha_subclasses, "_M00_M02")
  b5_prn     <- paste0("delta_", prn_subclasses, "_M00_M02")
  b5_full    <- c(b5_pt, b5_fha, b5_prn)

  cat("Missingness in Block 5 predictors:\n")
  print(colSums(is.na(d[, intersect(c(b5_full, outcome_col), names(d))])))

  cat("\n--- Decay magnitude check (medians of delta columns) ---\n")
  b5_present <- intersect(b5_antigen, names(d))
  if (length(b5_present) > 0)
    print(round(apply(d[, b5_present], 2,
                      function(x) median(x, na.rm = TRUE)), 3))

  # --- 5a. Antigen totals ---
  r5a <- run_lm(outcome_col, b5_antigen, d,
                "Block 5a — Δ antigen totals M00→M02 → PTNA_M02")
  r2_all[[paste0("B5a_", iarm)]] <- r2_row(r5a, iarm, "B5_decay_qty", "antigen_delta")

  # --- 5b–5d. Per-antigen subclass models ---
  r5b <- run_lm(outcome_col, b5_pt,  d, "Block 5b — Δ PT subclasses M00→M02")
  r5c <- run_lm(outcome_col, b5_fha, d, "Block 5c — Δ FHA subclasses M00→M02")
  r5d <- run_lm(outcome_col, b5_prn, d, "Block 5d — Δ PRN subclasses M00→M02")
  r2_all[[paste0("B5b_", iarm)]] <- r2_row(r5b, iarm, "B5_decay_qty", "PT_sub_delta")
  r2_all[[paste0("B5c_", iarm)]] <- r2_row(r5c, iarm, "B5_decay_qty", "FHA_sub_delta")
  r2_all[[paste0("B5d_", iarm)]] <- r2_row(r5d, iarm, "B5_decay_qty", "PRN_sub_delta")

  # --- 5e. Full 12-subclass delta model ---
  r5e <- run_lm(outcome_col, b5_full, d,
                "Block 5e — Full 12-subclass Δ M00→M02")
  r2_all[[paste0("B5e_", iarm)]] <- r2_row(r5e, iarm, "B5_decay_qty", "full12_delta")
  relimp_b5_full_by_iarm[[iarm]] <- r5e$relimp

  # ----------------------------------------------------------
  # BLOCK 6 — Ratio change M00 (cord) → M02 (infant)
  # Biological question: Do subclass-specific decay rates from
  # cord blood to 2 months predict infant PTNA?
  # ----------------------------------------------------------

  cat("\n\n===== BLOCK 6: Ratio Change M00 → M02 (Class-Selective Decay) =====\n")
  b6_cols <- c("delta_PT_ratio_M00_M02", "delta_FHA_ratio_M00_M02",
               "delta_PRN_ratio_M00_M02")
  cat("Missingness in Block 6 predictors:\n")
  print(colSums(is.na(d[, intersect(c(b6_cols, outcome_col), names(d))])))

  cat("\n--- Block 6 ratio delta summaries ---\n")
  print(summary(d[, intersect(b6_cols, names(d))]))

  r6 <- run_lm(outcome_col, b6_cols, d,
               paste0("Block 6 — Δ ratio M00→M02 (decay) → ", outcome_label, "_M02"))
  r2_all[[paste0("B6_", iarm)]] <- r2_row(r6, iarm, "B6_decay_ratio", "delta_ratio_M00_M02")
  ratio_models_by_block_arm[[paste0("B6_", iarm)]] <- r6$model

  cat("\n--- Block 6: by-antigen ratio delta models ---\n")
  for (ant in c("PT","FHA","PRN")) {
    dc <- paste0("delta_", ant, "_ratio_M00_M02")
    if (dc %in% names(d)) {
      m_tmp <- lm(reformulate(dc, response = outcome_col), data = d, na.action = na.omit)
      s_tmp <- summary(m_tmp); b_tmp <- coef(s_tmp)
      if (nrow(b_tmp) >= 2)
        cat(paste0("  ", ant, ": β=", round(b_tmp[2,1],4),
                   "  p=", round(b_tmp[2,4],4),
                   " ", pvalue_to_label(b_tmp[2,4]), "\n"))
    }
  }

} # ---- end infant_arm loop -----------------------------------


# ============================================================
# SECTION 7 — Cross-arm comparisons
# ============================================================

cat("\n\n============================================================\n")
cat(paste0("  CROSS-ARM COMPARISONS — Outcome: ", outcome_label, " at M02\n"))
cat("============================================================\n\n")

# ---- 7a. R² summary table across all blocks and models ------

cat("--- R² summary: all blocks × both arms ---\n")
r2_table <- bind_rows(r2_all) %>%
  arrange(block, model, infant_arm) %>%
  dplyr::select(infant_arm, block, model, n_obs, R2, adj_R2)
print(as.data.frame(r2_table))

# Pivot wider for side-by-side aP vs wP comparison
r2_wide <- r2_table %>%
  pivot_wider(names_from  = infant_arm,
              values_from = c(R2, adj_R2, n_obs),
              names_sep   = "_") %>%
  mutate(delta_R2_aP_wP = R2_aP - R2_wP) %>%
  arrange(block, model)
cat("\n--- R² comparison (aP vs wP) ---\n")
print(as.data.frame(r2_wide))

# ---- 7b. Ratio-model coefficient comparison plots -----------
#      One plot per ratio block (B2, B4, B6), side-by-side aP vs wP

for (blk in c("B2", "B4", "B6")) {

  block_labels <- c(
    B2 = "Block 2: Δ Ratio P02→P09 (class switching in pregnancy)",
    B4 = "Block 4: Δ Ratio P09→M00 (class-selective transfer)",
    B6 = "Block 6: Δ Ratio M00→M02 (class-selective decay)"
  )

  coef_rows <- list()
  for (iarm in infant_arms) {
    key <- paste0(blk, "_", iarm)
    if (!is.null(ratio_models_by_block_arm[[key]])) {
      coef_rows[[iarm]] <- tidy(ratio_models_by_block_arm[[key]]) %>%
        dplyr::filter(term != "(Intercept)") %>%
        mutate(infant_arm = iarm)
    }
  }

  if (length(coef_rows) > 0) {
    ratio_df <- bind_rows(coef_rows)
    p <- ggplot(ratio_df,
                aes(x = term, y = estimate, fill = infant_arm,
                    ymin = estimate - 1.96 * std.error,
                    ymax = estimate + 1.96 * std.error)) +
      geom_col(position = position_dodge(0.6), width = 0.5) +
      geom_errorbar(position = position_dodge(0.6), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = iarm_colors) +
      labs(
        x     = "Pro-inflammatory/Tolerogenic Ratio",
        y     = paste0("Regression Coefficient (→ ", outcome_label, "_M02)"),
        title = paste0(block_labels[blk], "\n(outcome = ", outcome_label, ")"),
        fill  = "Infant arm"
      ) +
      theme_bw(base_size = 11)
    print(p)
  }
}

# ---- 7c. Relative importance comparison — full-12-subclass models ---

cat("\n--- Relative importance comparison: aP vs wP (Block 1, P09 full model) ---\n")
if (!is.null(relimp_b1_full_by_iarm[["aP"]]) &&
    !is.null(relimp_b1_full_by_iarm[["wP"]])) {
  ri_compare <- merge(
    data.frame(feature = names(relimp_b1_full_by_iarm[["aP"]]@lmg),
               lmg_aP  = relimp_b1_full_by_iarm[["aP"]]@lmg),
    data.frame(feature = names(relimp_b1_full_by_iarm[["wP"]]@lmg),
               lmg_wP  = relimp_b1_full_by_iarm[["wP"]]@lmg),
    by = "feature"
  ) %>%
    mutate(delta_lmg = lmg_aP - lmg_wP) %>%
    arrange(desc(abs(delta_lmg)))
  print(ri_compare)
}

cat("\n--- Relative importance comparison: aP vs wP (Block 3 transfer full model) ---\n")
if (!is.null(relimp_b3_full_by_iarm[["aP"]]) &&
    !is.null(relimp_b3_full_by_iarm[["wP"]])) {
  ri_compare <- merge(
    data.frame(feature = names(relimp_b3_full_by_iarm[["aP"]]@lmg),
               lmg_aP  = relimp_b3_full_by_iarm[["aP"]]@lmg),
    data.frame(feature = names(relimp_b3_full_by_iarm[["wP"]]@lmg),
               lmg_wP  = relimp_b3_full_by_iarm[["wP"]]@lmg),
    by = "feature"
  ) %>%
    mutate(delta_lmg = lmg_aP - lmg_wP) %>%
    arrange(desc(abs(delta_lmg)))
  print(ri_compare)
}

cat("\n--- Relative importance comparison: aP vs wP (Block 5 decay full model) ---\n")
if (!is.null(relimp_b5_full_by_iarm[["aP"]]) &&
    !is.null(relimp_b5_full_by_iarm[["wP"]])) {
  ri_compare <- merge(
    data.frame(feature = names(relimp_b5_full_by_iarm[["aP"]]@lmg),
               lmg_aP  = relimp_b5_full_by_iarm[["aP"]]@lmg),
    data.frame(feature = names(relimp_b5_full_by_iarm[["wP"]]@lmg),
               lmg_wP  = relimp_b5_full_by_iarm[["wP"]]@lmg),
    by = "feature"
  ) %>%
    mutate(delta_lmg = lmg_aP - lmg_wP) %>%
    arrange(desc(abs(delta_lmg)))
  print(ri_compare)
}

# ---- 7d. Interaction plot: subclass slopes by arm (Block 1, P09) ---

cat("\n--- Block 1 (P09) subclass interaction plot: aP vs wP ---\n")

data_long_b1 <- bind_rows(
  data_wide %>%
    dplyr::filter(infant_arm == "aP") %>%
    dplyr::select(all_of(outcome_col),
                  matches("IgG[1-4]_P09$")) %>%
    mutate(infant_arm = "aP"),
  data_wide %>%
    dplyr::filter(infant_arm == "wP") %>%
    dplyr::select(all_of(outcome_col),
                  matches("IgG[1-4]_P09$")) %>%
    mutate(infant_arm = "wP")
) %>%
  pivot_longer(
    cols          = matches("IgG[1-4]_P09$"),
    names_to      = c("antigen", "subclass", "visit"),
    names_pattern = "(PT|FHA|PRN)_(IgG[1-4])_(P09)",
    values_to     = "subclass_value"
  ) %>%
  drop_na()

if (nrow(data_long_b1) > 0) {
  p_b1 <- ggplot(data_long_b1,
                 aes(x = subclass_value, y = .data[[outcome_col]],
                     color = subclass, linetype = infant_arm)) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    facet_grid(antigen ~ infant_arm) +
    scale_color_brewer(palette = "Set1") +
    labs(
      x        = "Log10 Maternal Subclass IgG at P09",
      y        = outcome_y_label,
      color    = "Subclass",
      linetype = "Infant arm",
      title    = paste0("Block 1 (P09): Antigen × Subclass → ",
                        outcome_label, "_M02: aP vs wP")
    ) +
    theme_bw(base_size = 11)
  print(p_b1)
}

# ---- 7e. Interaction plot: Block 5 (decay deltas) ----------

data_long_b5 <- bind_rows(
  data_wide %>%
    dplyr::filter(infant_arm == "aP") %>%
    dplyr::select(all_of(outcome_col),
                  matches("^delta_(PT|FHA|PRN)_IgG[1-4]_M00_M02$")) %>%
    mutate(infant_arm = "aP"),
  data_wide %>%
    dplyr::filter(infant_arm == "wP") %>%
    dplyr::select(all_of(outcome_col),
                  matches("^delta_(PT|FHA|PRN)_IgG[1-4]_M00_M02$")) %>%
    mutate(infant_arm = "wP")
) %>%
  pivot_longer(
    cols          = matches("^delta_(PT|FHA|PRN)_IgG[1-4]_M00_M02$"),
    names_to      = c("antigen", "subclass"),
    names_pattern = "^delta_(PT|FHA|PRN)_(IgG[1-4])_M00_M02$",
    values_to     = "delta_value"
  ) %>%
  drop_na()

if (nrow(data_long_b5) > 0) {
  p_b5 <- ggplot(data_long_b5,
                 aes(x = delta_value, y = .data[[outcome_col]],
                     color = subclass, linetype = infant_arm)) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    facet_grid(antigen ~ infant_arm) +
    scale_color_brewer(palette = "Set1") +
    labs(
      x        = "Δ Log10 Subclass IgG (M02 − M00, decay = negative)",
      y        = outcome_y_label,
      color    = "Subclass",
      linetype = "Infant arm",
      title    = paste0("Block 5 (decay): Antigen × Subclass → ",
                        outcome_label, "_M02: aP vs wP")
    ) +
    theme_bw(base_size = 11)
  print(p_b5)
}

cat("\n\n============================================================\n")
cat(paste0("  DONE — Outcome: ", outcome_label,
           " at M02 — adapt outcome_feature in Section 1 for SBA/WT_IgG\n"))
cat("============================================================\n")
