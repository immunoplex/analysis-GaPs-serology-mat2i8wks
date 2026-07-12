# =====================================================================
# config/endpoints.R   (REFACTORED 2026-07 — cord-predictor + ratio-stage)
# ---------------------------------------------------------------------
# THE single control surface for the whole project.
#
# To run the analysis for a different endpoint you change ONE thing:
# the `endpoint` parameter in the master .Rmd YAML header
# (params$endpoint = "SBA" | "WT_IgG" | "PTNA"). Everything downstream
# — outcome column, plot titles, which visits are in play, the predictor
# lists — is derived from the ENDPOINTS list below.
#
# Nothing in this file depends on the data being loaded; it is pure
# configuration and is sourced before anything else.
#
# ---------------------------------------------------------------------
# WHAT CHANGED IN THIS REFACTOR (all additive; nothing removed):
#   §5  Every ENDPOINTS entry gains four OPTIONAL fields consumed by the
#       two focused SBA analyses (see run_SBA_cord_and_ratios.Rmd):
#         cord_predictor        raw feature name used as the cord-blood
#                               temporal predictor of the InfMon2 outcome
#                               (Result #1). NULL disables the cord block.
#         cord_predictor_label  display string for that predictor.
#         chain_signal_visit    the maternal visit whose measurement of
#                               `cord_predictor` is the upstream "Block-1
#                               maternal chain signal" that the cord block
#                               tests mediation of (default "MatBirth").
#         ratio_stages          the upstream visits whose ABSOLUTE
#                               pro-/tolerogenic subclass ratios are
#                               regressed on the outcome in one place
#                               (Result #2; default = the three non-infant
#                               chain visits).
#       SBA turns the cord block ON with cord_predictor = "WHOLE_WT_IgG"
#       (the whole-cell binding IgG measured in cord blood). This is the
#       one-line change that flips the previously "Not applicable for SBA"
#       cord-blood predictor slot into an active analysis.
#   §6  get_endpoint_config() now surfaces those fields with safe
#       defaults, plus the derived cord/outcome column names, so the parts
#       never have to reconstruct a column string by hand.
#
# NOTE: this file is a faithful superset of the previous endpoints.R — the
# feature vocabulary, arms/colours, clinical covariates, bubble config and
# casing guard are unchanged, so it is a drop-in replacement.
#
# ---------------------------------------------------------------------
# !! CASING ASSUMPTION — VERIFY, DO NOT ASSUME !!
# This file names features in mixed case (PT_IgG, PT_IgG1, FcgR2a), which
# is what the WORKING pipeline (endpoints.R + C1_C1_serology_helpers.R) uses and
# what the rendered analyses produced. The chain/concurrent paths match
# feature names CASE-SENSITIVELY against levels(data$feature); a mismatch
# silently drops the feature. Call verify_predictor_casing(data_raw)
# immediately after load_serology_data() in each driver .Rmd. It ERRORS
# (not warns) if any feature is absent. Do not knit until it passes.
# =====================================================================


# ---------------------------------------------------------------------
# 1. GLOBAL VISIT RECODE MAP
# ---------------------------------------------------------------------
#   old   ->  new
#   P02   ->  PregEarly    (maternal, early pregnancy / pre-vaccination)
#   P09   ->  MatBirth     (maternal, at delivery)
#   M00   ->  CordBlood    (cord blood at birth)
#   M02   ->  InfMon2      (infant, 2 months — before first vaccination)
#   M04   ->  InfMon4      (infant, 4 months)
#   M05   ->  InfMon5      (infant, 5 months)
#   M09   ->  InfMon9      (infant, 9 months)
VISIT_RECODE <- c(
  P02 = "PregEarly",
  P09 = "MatBirth",
  M00 = "CordBlood",
  M02 = "InfMon2",
  M04 = "InfMon4",
  M05 = "InfMon5",
  M09 = "InfMon9"
)

# Canonical display order of the (new) visit names.
VISIT_ORDER <- c("PregEarly", "MatBirth", "CordBlood",
                 "InfMon2", "InfMon4", "InfMon5", "InfMon9")

# Numeric index per visit (used by the concurrent network part only).
VISIT_NUMBER <- c(
  PregEarly = -9, MatBirth = -1, CordBlood = 0,
  InfMon2   =  2, InfMon4   =  4, InfMon5   = 5, InfMon9 = 9
)


# ---------------------------------------------------------------------
# 2. PREDICTOR FEATURE VOCABULARY (shared by all endpoints)
# ---------------------------------------------------------------------
ANTIGEN_TOTALS <- c("PT_IgG", "FHA_IgG", "PRN_IgG",
                    "DT_IgG", "TT_IgG", "FIM_IgG")

PT_SUBCLASSES  <- c("PT_IgG1",  "PT_IgG2",  "PT_IgG3",  "PT_IgG4")
FHA_SUBCLASSES <- c("FHA_IgG1", "FHA_IgG2", "FHA_IgG3", "FHA_IgG4")
PRN_SUBCLASSES <- c("PRN_IgG1", "PRN_IgG2", "PRN_IgG3", "PRN_IgG4")
DT_SUBCLASSES  <- c("DT_IgG1",  "DT_IgG2",  "DT_IgG3",  "DT_IgG4")
TT_SUBCLASSES  <- c("TT_IgG1",  "TT_IgG2",  "TT_IgG3",  "TT_IgG4")
FIM_SUBCLASSES <- c("FIM_IgG1", "FIM_IgG2", "FIM_IgG3", "FIM_IgG4")

# Chain-safe subclass set (PT/FHA/PRN + DT/TT). NOT FIM — FIM subclasses
# are absent at maternal/cord visits and would wipe chain models by
# complete-case deletion. This drives the ratio construction in
# add_chain_derived(), which is what Result #2 reads.
ALL_SUBCLASSES <- c(PT_SUBCLASSES, FHA_SUBCLASSES, PRN_SUBCLASSES,
                    DT_SUBCLASSES, TT_SUBCLASSES)

# Network-node subclass set for the concurrent (infant-visit) GGM only.
NETWORK_SUBCLASSES <- c(PT_SUBCLASSES, FHA_SUBCLASSES, PRN_SUBCLASSES,
                        FIM_SUBCLASSES)

# Core pertussis set — drives the ratio loop in add_chain_derived().
ANTIGENS       <- c("PT", "FHA", "PRN")

# Canonical antigen display order (single source of truth).
ANTIGEN_ORDER  <- c("PT", "FHA", "PRN", "FIM", "TT", "DT")

# Selection vocabulary pulled into the wide data.
PREDICTOR_FEATURES <- c(ANTIGEN_TOTALS, ALL_SUBCLASSES, FIM_SUBCLASSES)


# ---------------------------------------------------------------------
# 3. STRATIFICATION ARMS + COLOURS (shared)
# ---------------------------------------------------------------------
MATERNAL_ARMS <- c("TdaP", "TT")     # TdaP = pertussis booster; TT = comparator
INFANT_ARMS   <- c("aP", "wP")       # acellular vs whole-cell infant priming

IARM_COLORS <- c(aP = "dodgerblue2", wP = "forestgreen")
MARM_COLORS <- c(TdaP = "skyblue",   TT = "salmon")


# ---------------------------------------------------------------------
# 4. CLINICAL COVARIATES (kept for compatibility with the clinical part)
# ---------------------------------------------------------------------
COVARIATES_NAMED <- list(
  "Birth Weight"                             = "birth_weight",
  "Infant Sex"                               = "infant_sex",
  "Delivery Mode"                            = "delivery_mode",
  "Parity"                                   = "parity",
  "Gestational Age at Vaccine"               = "gestational_age_vaccination",
  "Interval: Vaccination to Delivery (days)" = "vaccine_birth_interval_days",
  "Gestational Age at Delivery"              = "gestational_age_birth",
  "Maternal Age at Vaccination"              = "maternal_age",
  "Maternal BMI at Vaccination"              = "maternal_bmi",
  "Maternal Haemoglobin"                     = "maternal_Hb"
)
COVARIATES_CATEGORICAL <- c("infant_sex", "delivery_mode", "parity")


# ---------------------------------------------------------------------
# 5. PER-ENDPOINT CONFIGURATION
# ---------------------------------------------------------------------
# Each endpoint specifies only what differs between SBA / WT_IgG / PTNA.
# The four cord_/ratio_ fields are OPTIONAL and consumed only by the two
# focused analyses; omit them and get_endpoint_config() supplies defaults.
ENDPOINTS <- list(

  SBA = list(
    outcome_feature      = "WHOLE_SBA",
    outcome_label        = "SBA",
    outcome_y_label      = "Log10 SBA",
    outcome_long         = "serum bactericidal activity",
    chain_visits         = c("PregEarly", "MatBirth", "CordBlood", "InfMon2"),
    chain_outcome_visit  = "InfMon2",
    concurrent_visits    = c("InfMon2", "InfMon5", "InfMon9"),

    # --- cord-blood temporal predictor (Result #1) — NOW ACTIVE FOR SBA ---
    cord_predictor       = "WHOLE_WT_IgG",       # whole-cell IgG in cord blood
    cord_predictor_label = "Whole-cell binding IgG",
    cord_visit           = "CordBlood",
    chain_signal_visit   = "MatBirth",           # upstream Block-1 maternal signal
    # --- consolidated absolute subclass ratios (Result #2) ---
    ratio_stages         = c("PregEarly", "MatBirth", "CordBlood")
  ),

  WT_IgG = list(
    outcome_feature      = "WHOLE_WT_IgG",
    outcome_label        = "WT_IgG",
    outcome_y_label      = "Log10 whole-cell IgG",
    outcome_long         = "whole-cell-pertussis-binding IgG",
    chain_visits         = c("PregEarly", "MatBirth", "CordBlood", "InfMon2"),
    chain_outcome_visit  = "InfMon2",
    concurrent_visits    = c("CordBlood", "InfMon2", "InfMon5", "InfMon9"),

    cord_predictor       = "WHOLE_WT_IgG",       # its own cord level
    cord_predictor_label = "Whole-cell binding IgG",
    cord_visit           = "CordBlood",
    chain_signal_visit   = "MatBirth",
    ratio_stages         = c("PregEarly", "MatBirth", "CordBlood")
  ),

  PTNA = list(
    outcome_feature      = "WHOLE_PTNA",
    outcome_label        = "PTNA",
    outcome_y_label      = "Log10 PTNA",
    outcome_long         = "pertussis-toxin neutralization",
    chain_visits         = c("PregEarly", "MatBirth", "CordBlood", "InfMon2"),
    chain_outcome_visit  = "InfMon2",
    concurrent_visits    = c("InfMon2", "InfMon5", "InfMon9"),

    # cord_predictor left NULL -> cord block prints "Not applicable"
    cord_predictor       = NULL,
    cord_predictor_label = NULL,
    cord_visit           = "CordBlood",
    chain_signal_visit   = "MatBirth",
    ratio_stages         = c("PregEarly", "MatBirth", "CordBlood")
  )
)


# ---------------------------------------------------------------------
# 6. ACCESSOR
# ---------------------------------------------------------------------
# Returns the config list for one endpoint, augmented with derived
# convenience fields the parts rely on. New: safe defaults + derived
# column names for the cord-predictor and ratio-stage analyses.
get_endpoint_config <- function(endpoint) {
  if (!endpoint %in% names(ENDPOINTS)) {
    stop("Unknown endpoint '", endpoint, "'. Known: ",
         paste(names(ENDPOINTS), collapse = ", "))
  }
  cfg <- ENDPOINTS[[endpoint]]

  # small NULL-coalescing helper (avoid depending on rlang::`%||%`)
  `%or%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

  # Derived fields shared by all parts
  cfg$endpoint        <- endpoint
  cfg$outcome_col     <- paste0(cfg$outcome_feature, "_", cfg$chain_outcome_visit)
  cfg$antigen_totals  <- ANTIGEN_TOTALS
  cfg$all_subclasses  <- ALL_SUBCLASSES
  cfg$network_subclasses <- NETWORK_SUBCLASSES
  cfg$antigens        <- ANTIGENS
  cfg$maternal_arms   <- MATERNAL_ARMS
  cfg$infant_arms     <- INFANT_ARMS
  cfg$iarm_colors     <- IARM_COLORS
  cfg$marm_colors     <- MARM_COLORS

  # ---- cord-predictor block (Result #1) ----
  cfg$cord_visit          <- cfg$cord_visit          %or% "CordBlood"
  cfg$chain_signal_visit  <- cfg$chain_signal_visit  %or% "MatBirth"
  cfg$infant_visit        <- cfg$chain_outcome_visit                   # "InfMon2"
  # NULL cord_predictor stays NULL -> the cord part self-disables.
  if (!is.null(cfg$cord_predictor) && length(cfg$cord_predictor) &&
      nzchar(cfg$cord_predictor)) {
    cfg$cord_col   <- paste0(cfg$cord_predictor, "_", cfg$cord_visit)          # WHOLE_WT_IgG_CordBlood
    cfg$signal_col <- paste0(cfg$cord_predictor, "_", cfg$chain_signal_visit)  # WHOLE_WT_IgG_MatBirth
    cfg$cord_predictor_label <- cfg$cord_predictor_label %or% cfg$cord_predictor
  } else {
    cfg$cord_predictor <- NULL
    cfg$cord_col <- NULL; cfg$signal_col <- NULL
  }

  # ---- ratio-stage block (Result #2) ----
  # default: every chain visit except the infant outcome visit
  cfg$ratio_stages <- cfg$ratio_stages %or%
    setdiff(cfg$chain_visits, cfg$chain_outcome_visit)

  cfg
}


# ---------------------------------------------------------------------
# 7. PART 4 CONFIG — MATERNAL-VACCINATION-EFFECT "BUBBLE" ANALYSIS
# ---------------------------------------------------------------------
# (Unchanged. Not used by the two focused SBA analyses, but retained so
#  this file remains a drop-in replacement for the full pipeline.)
BUBBLE_ANTIGENS <- ANTIGEN_ORDER
BUBBLE_ANALYTES <- c("IgG", "IgG3", "IgG1", "IgG2", "IgG4",
                     "FcgR2a", "FcgR3b", "ADCD", "ADCP", "ADNP")
BUBBLE_VISITS   <- c("InfMon2", "InfMon5", "InfMon9")
BUBBLE_IGGLIST  <- c("DT_IgG", "FHA_IgG", "PRN_IgG", "PT_IgG", "TT_IgG",
                     "FIM_IgG")

BUBBLE_INCLUDE_FEATURE <- c(
  "DT_IgG","DT_IgG1","DT_IgG2","DT_IgG3","DT_IgG4","DT_ADCD","DT_ADNP","DT_ADCP","DT_FcgR2a","DT_FcgR3b",
  "TT_IgG","TT_IgG1","TT_IgG2","TT_IgG3","TT_IgG4","TT_ADCD","TT_ADNP","TT_ADCP","TT_FcgR2a","TT_FcgR3b","TT_Bmem",
  "PT_IgG","PT_IgG1","PT_IgG2","PT_IgG3","PT_IgG4","PT_ADCD","PT_ADNP","PT_ADCP","PT_FcgR2a","PT_FcgR3b","PT_Bmem",
  "FHA_IgG","FHA_IgG1","FHA_IgG2","FHA_IgG3","FHA_IgG4","FHA_ADCD","FHA_ADNP","FHA_ADCP","FHA_FcgR2a","FHA_FcgR3b","FHA_Bmem",
  "PRN_IgG","PRN_IgG1","PRN_IgG2","PRN_IgG3","PRN_IgG4","PRN_ADCD","PRN_ADNP","PRN_ADCP","PRN_FcgR2a","PRN_FcgR3b",
  "FIM_IgG","FIM_IgG1","FIM_IgG2","FIM_IgG3","FIM_IgG4","FIM_ADCD","FIM_FcgR2a","FIM_FcgR3b")

BUBBLE_INCLUDE_FEATURE_WO <- c(
  "DT_IgG1","DT_IgG2","DT_IgG3","DT_IgG4","DT_ADCD","DT_ADNP","DT_ADCP","DT_FcgR2a","DT_FcgR3b",
  "TT_IgG1","TT_IgG2","TT_IgG3","TT_IgG4","TT_ADCD","TT_ADNP","TT_ADCP","TT_FcgR2a","TT_FcgR3b","TT_Bmem",
  "PT_IgG1","PT_IgG2","PT_IgG3","PT_IgG4","PT_ADCD","PT_ADNP","PT_ADCP","PT_FcgR2a","PT_FcgR3b","PT_Bmem",
  "FHA_IgG1","FHA_IgG2","FHA_IgG3","FHA_IgG4","FHA_ADCD","FHA_ADNP","FHA_ADCP","FHA_FcgR2a","FHA_FcgR3b","FHA_Bmem",
  "PRN_IgG1","PRN_IgG2","PRN_IgG3","PRN_IgG4","PRN_ADCD","PRN_ADNP","PRN_ADCP","PRN_FcgR2a","PRN_FcgR3b",
  "FIM_IgG1","FIM_IgG2","FIM_IgG3","FIM_IgG4","FIM_ADCD","FIM_FcgR2a","FIM_FcgR3b")

BUBBLE_EFFECT_COLORS <- c("Lower TdaP" = "salmon", "Higher TdaP" = "skyblue", "No effect" = "darkgrey")

BUBBLE_RESIDUAL_AP <- "data/igg_standard_residuals_ap_matpm_prevacvacc_k.RData"
BUBBLE_RESIDUAL_WP <- "data/igg_standard_residuals_wp_matpm_prevacvacc_k.RData"


# ---------------------------------------------------------------------
# 8. CASING VERIFICATION  (unchanged)
# ---------------------------------------------------------------------
verify_predictor_casing <- function(data_long, feature_col = "feature",
                                    stop_on_fail = TRUE) {
  present <- as.character(unique(data_long[[feature_col]]))

  needed <- unique(c(
    ANTIGEN_TOTALS,
    ALL_SUBCLASSES, FIM_SUBCLASSES,
    DT_SUBCLASSES, TT_SUBCLASSES
  ))

  missing <- setdiff(needed, present)
  found   <- intersect(needed, present)

  case_clashes <- character(0)
  if (length(missing)) {
    up <- toupper(present)
    for (m in missing) {
      alt <- present[up == toupper(m)]
      if (length(alt)) case_clashes <- c(case_clashes, sprintf("%s -> data has %s", m, alt[1]))
    }
  }

  msg <- sprintf("verify_predictor_casing(): %d/%d requested features present.",
                 length(found), length(needed))
  if (length(case_clashes)) {
    msg <- paste0(msg,
      "\n  CASE MISMATCH (config casing != data casing) for ",
      length(case_clashes), " feature(s):\n    ",
      paste(utils::head(case_clashes, 12), collapse = "\n    "),
      if (length(case_clashes) > 12) "\n    ..." else "",
      "\n  FIX: make the config constants match the data casing (or recase",
      "\n       `analyte` once inside load_serology_data()), then re-source.")
  }
  truly_absent <- setdiff(missing, sub(" ->.*$", "", case_clashes))
  if (length(truly_absent)) {
    msg <- paste0(msg, "\n  ABSENT in this dataset (any case): ",
                  paste(utils::head(truly_absent, 12), collapse = ", "),
                  if (length(truly_absent) > 12) " ..." else "")
  }

  if (length(missing) && stop_on_fail) stop(msg, call. = FALSE)
  message(msg, if (!length(missing)) "  All present — casing OK." else "")
  invisible(list(found = found, missing = missing, case_clashes = case_clashes))
}
