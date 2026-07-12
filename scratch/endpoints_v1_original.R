# =====================================================================
# config/endpoints.R
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
# =====================================================================
#
# =====================================================================
# PATCH 2026-06 — add cross-reactive toxoid (DT, TT) + FIM features to
#                 answer the toxoid-driver and FIM-after-TdaP questions.
# ---------------------------------------------------------------------
# WHAT CHANGED (all additive; nothing removed):
#   §2  ANTIGEN_TOTALS gains DT_IgG, TT_IgG, FIM_IgG  -> flows into BOTH
#       the chain Block-1 antigen models (run_chain_blocks uses
#       cfg$antigen_totals) and the concurrent antigen-total relimp
#       (run_concurrent_models uses cfg$antigen_totals). All three totals
#       exist at every visit, so this is safe across the whole chain.
#   §2  NEW DT_SUBCLASSES / TT_SUBCLASSES / FIM_SUBCLASSES, and a NEW
#       NETWORK_SUBCLASSES vector for the concurrent GGM nodes.
#       ALL_SUBCLASSES is deliberately LEFT UNCHANGED (PT/FHA/PRN only):
#       it feeds the maternal-chain full-12 models, and FIM subclasses do
#       NOT exist at PregEarly/MatBirth/CordBlood — adding them there would
#       silently wipe those models via complete-case deletion. FIM nodes
#       therefore enter only the infant-visit network (see Part 3 edit).
#   §7  BUBBLE_* gains FIM (IgG, 4 subclasses, ADCD, FcgR2a, FcgR3b only —
#       FIM has no ADCP/ADNP/Bmem).
#
# ---------------------------------------------------------------------
# !! CASING ASSUMPTION — VERIFY, DO NOT ASSUME !!
# This file names features in mixed case (PT_IgG, PT_IgG1, FcgR2a), which
# is what the WORKING pipeline (endpoints.R + serology_helpers.R) uses and
# what the rendered analyses produced. BUT the raw `table(data$feature,...)`
# dump and config/endpoints_additions.R are UPPERCASE (PT_IGG1, FCGR2A).
# The chain/concurrent paths match feature names CASE-SENSITIVELY against
# levels(data$feature); a mismatch silently drops the feature (the bubble
# path is immune — it toupper()s both sides). So the new constants must
# match the post-load casing EXACTLY.
#
# Call verify_predictor_casing(data_raw) immediately after
# load_serology_data() in each driver .Rmd. It ERRORS (not warns) if any
# new feature is absent, and tells you if only the OTHER casing is present
# so the fix is obvious. Do not knit until it passes.
# =====================================================================


# ---------------------------------------------------------------------
# 1. GLOBAL VISIT RECODE MAP
# ---------------------------------------------------------------------
# The raw data uses the old P0x / M0x visit codes. The recipients of the
# analysis found these confusing, so every visit label is recoded ONCE,
# immediately after load (see R/serology_helpers.R::load_serology_data),
# and the new names then flow automatically into every column name,
# derived delta, regex, caption and figure label downstream.
#
#   old   ->  new
#   P02   ->  PregEarly    (maternal, early pregnancy / pre-vaccination)
#   P09   ->  MatBirth     (maternal, at delivery)
#   M00   ->  CordBlood    (cord blood at birth)
#   M02   ->  InfMon2      (infant, 2 months — before first vaccination)
#   M04   ->  InfMon4      (infant, 4 months)
#   M05   ->  InfMon5      (infant, 5 months)
#   M09   ->  InfMon9      (infant, 9 months)
#
# Order matters only for display; named lookup is used for recoding.
VISIT_RECODE <- c(
  P02 = "PregEarly",
  P09 = "MatBirth",
  M00 = "CordBlood",
  M02 = "InfMon2",
  M04 = "InfMon4",
  M05 = "InfMon5",
  M09 = "InfMon9"
)

# Canonical display order of the (new) visit names, used to keep tables
# and factor levels in chronological order.
VISIT_ORDER <- c("PregEarly", "MatBirth", "CordBlood",
                 "InfMon2", "InfMon4", "InfMon5", "InfMon9")

# A small numeric index per infant visit, used only by the concurrent
# network part for `visit_number`-style filtering. Maternal/cord visits
# are given negative/zero sentinels so they never collide with infant
# month numbers.
VISIT_NUMBER <- c(
  PregEarly = -9, MatBirth = -1, CordBlood = 0,
  InfMon2   =  2, InfMon4   =  4, InfMon5   = 5, InfMon9 = 9
)


# ---------------------------------------------------------------------
# 2. PREDICTOR FEATURE VOCABULARY (shared by all endpoints)
# ---------------------------------------------------------------------
# Antigen totals: PT/FHA/PRN are the pertussis drivers; DT/TT are the
# cross-reactive toxoids; FIM is fimbriae. All six TOTALS exist at every
# visit (maternal/cord and infant), so all six are safe in the chain
# Block-1 antigen models and the concurrent antigen-total relimp.
ANTIGEN_TOTALS <- c("PT_IgG", "FHA_IgG", "PRN_IgG",
                    "DT_IgG", "TT_IgG", "FIM_IgG")   # PATCH: +DT/TT/FIM totals

PT_SUBCLASSES  <- c("PT_IgG1",  "PT_IgG2",  "PT_IgG3",  "PT_IgG4")
FHA_SUBCLASSES <- c("FHA_IgG1", "FHA_IgG2", "FHA_IgG3", "FHA_IgG4")
PRN_SUBCLASSES <- c("PRN_IgG1", "PRN_IgG2", "PRN_IgG3", "PRN_IgG4")
DT_SUBCLASSES  <- c("DT_IgG1",  "DT_IgG2",  "DT_IgG3",  "DT_IgG4")    # PATCH (maternal+infant)
TT_SUBCLASSES  <- c("TT_IgG1",  "TT_IgG2",  "TT_IgG3",  "TT_IgG4")    # PATCH (maternal+infant)
FIM_SUBCLASSES <- c("FIM_IgG1", "FIM_IgG2", "FIM_IgG3", "FIM_IgG4")   # PATCH (INFANT VISITS ONLY)

# Chain-safe subclass set: PT/FHA/PRN only. Used by the maternal-chain
# full-12 models AND the concurrent full-12 models via cfg$all_subclasses.
# DO NOT add FIM here — FIM subclasses are absent at maternal/cord visits
# and would wipe the chain full-12 models by complete-case deletion.
ALL_SUBCLASSES <- c(PT_SUBCLASSES, FHA_SUBCLASSES, PRN_SUBCLASSES,DT_SUBCLASSES,TT_SUBCLASSES)

# Network-node subclass set for the concurrent (infant-visit) GGM only.
# FIM is estimable at InfMon2/5/9, so it can be a network node here. Add
# DT_SUBCLASSES / TT_SUBCLASSES to this vector too if you want them as
# nodes (they exist at infant visits); keep the order so the positional
# community vector in Part 3 stays aligned.
NETWORK_SUBCLASSES <- c(PT_SUBCLASSES, FHA_SUBCLASSES, PRN_SUBCLASSES,
                        FIM_SUBCLASSES)              # PATCH: +FIM nodes

ANTIGENS       <- c("PT", "FHA", "PRN")   # core pertussis set (unchanged on purpose;
                                          # this drives add_concurrent_ratios' loop)

# ---------------------------------------------------------------------
# CANONICAL ANTIGEN DISPLAY ORDER  (single source of truth)
# ---------------------------------------------------------------------
# Every figure / table that lays antigens out in rows or columns reads
# this one vector, so the order is changed in exactly one place. The
# order is: the three TdaP pertussis antigens (PT, FHA, PRN), then FIM,
# then the cross-reactive toxoids (TT, DT). Consumers: results_B_to_F.Rmd
# (Figure 3 pathway), 06_prediction_models.Rmd (ANTI6), BUBBLE_ANTIGENS
# below, and config/endpoints_additions.R (Parts 4/5 ANTIGENS shadow).
ANTIGEN_ORDER  <- c("PT", "FHA", "PRN", "FIM", "TT", "DT")

# Selection vocabulary pulled into the wide data. Includes FIM_SUBCLASSES
# so the concurrent network has the FIM columns available; the chain only
# *models* ALL_SUBCLASSES, so the extra FIM columns sit unused there.
PREDICTOR_FEATURES <- c(ANTIGEN_TOTALS, ALL_SUBCLASSES, FIM_SUBCLASSES)  # PATCH


# ---------------------------------------------------------------------
# 3. STRATIFICATION ARMS + COLOURS (shared)
# ---------------------------------------------------------------------
MATERNAL_ARMS <- c("TdaP", "TT")     # TdaP = pertussis booster; TT = comparator
INFANT_ARMS   <- c("aP", "wP")       # acellular vs whole-cell infant priming

IARM_COLORS <- c(aP = "dodgerblue2", wP = "forestgreen")
MARM_COLORS <- c(TdaP = "skyblue",   TT = "salmon")


# ---------------------------------------------------------------------
# 4. CLINICAL COVARIATES (shared by the clinical-covariate part)
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
# Each endpoint specifies only what differs between SBA / WT_IgG / PTNA:
#   outcome_feature   : the raw feature name in the data
#   outcome_label     : short label used in titles / prose
#   outcome_y_label   : axis label for plots
#   chain_visits      : ordered visits for the maternal->infant chain part
#                       (Block 1..6 framework). Identical across endpoints
#                       because the chain always predicts the infant at
#                       2 months from the upstream maternal/cord chain.
#   chain_outcome_visit : the infant visit the chain predicts.
#   concurrent_visits : the infant visits used by the concurrent network /
#                       mixed-model part. THIS differs: WT_IgG adds CordBlood.
#
# NOTE on chain_outcome_visit: the maternal chain is fundamentally a
# "predict the infant at 2 months from the maternal supply" analysis, so
# the outcome visit is InfMon2 for every endpoint. The InfMon5 / InfMon9
# (and CordBlood, for WT_IgG) outcomes live in the concurrent part.
ENDPOINTS <- list(

  SBA = list(
    outcome_feature      = "WHOLE_SBA",
    outcome_label        = "SBA",
    outcome_y_label      = "Log10 SBA",
    outcome_long         = "serum bactericidal activity",
    chain_visits         = c("PregEarly", "MatBirth", "CordBlood", "InfMon2"),
    chain_outcome_visit  = "InfMon2",
    concurrent_visits    = c("InfMon2", "InfMon5", "InfMon9")
  ),

  WT_IgG = list(
    outcome_feature      = "WHOLE_WT_IgG",
    outcome_label        = "WT_IgG",
    outcome_y_label      = "Log10 whole-cell IgG",
    outcome_long         = "whole-cell-pertussis-binding IgG",
    chain_visits         = c("PregEarly", "MatBirth", "CordBlood", "InfMon2"),
    chain_outcome_visit  = "InfMon2",
    concurrent_visits    = c("CordBlood", "InfMon2", "InfMon5", "InfMon9")
  ),

  PTNA = list(
    outcome_feature      = "WHOLE_PTNA",
    outcome_label        = "PTNA",
    outcome_y_label      = "Log10 PTNA",
    outcome_long         = "pertussis-toxin neutralization",
    chain_visits         = c("PregEarly", "MatBirth", "CordBlood", "InfMon2"),
    chain_outcome_visit  = "InfMon2",
    concurrent_visits    = c("InfMon2", "InfMon5", "InfMon9")
  )
)


# ---------------------------------------------------------------------
# 6. ACCESSOR
# ---------------------------------------------------------------------
# Returns the config list for one endpoint, augmented with a few derived
# convenience fields the parts rely on.
get_endpoint_config <- function(endpoint) {
  if (!endpoint %in% names(ENDPOINTS)) {
    stop("Unknown endpoint '", endpoint, "'. Known: ",
         paste(names(ENDPOINTS), collapse = ", "))
  }
  cfg <- ENDPOINTS[[endpoint]]

  # Derived fields shared by all parts
  cfg$endpoint        <- endpoint
  cfg$outcome_col     <- paste0(cfg$outcome_feature, "_", cfg$chain_outcome_visit)
  cfg$antigen_totals  <- ANTIGEN_TOTALS
  cfg$all_subclasses  <- ALL_SUBCLASSES
  cfg$network_subclasses <- NETWORK_SUBCLASSES   # PATCH: concurrent GGM nodes (incl. FIM)
  cfg$antigens        <- ANTIGENS
  cfg$maternal_arms   <- MATERNAL_ARMS
  cfg$infant_arms     <- INFANT_ARMS
  cfg$iarm_colors     <- IARM_COLORS
  cfg$marm_colors     <- MARM_COLORS
  cfg
}


# ---------------------------------------------------------------------
# 7. PART 4 CONFIG — MATERNAL-VACCINATION-EFFECT "BUBBLE" ANALYSIS
# ---------------------------------------------------------------------
# This part is NOT endpoint-parameterized. It surveys the maternal
# vaccination effect (TdaP vs TT, referent TT) across a broad antigen x
# analyte grid, in untransformed and IgG-standardized form, at the three
# infant visits. It reuses VISIT_RECODE, the arms, colours, the N tracker,
# pvalue_to_label, and mapnote from the shared library.

# Broader vocabulary than the endpoint analyses: adds DT/TT antigens and
# functional analytes (ADCD/ADCP/ADNP, FcgR2a/FcgR3b) alongside subclasses.
# PATCH: +FIM. NOTE FIM has IgG, IgG1-4, ADCD, FcgR2a, FcgR3b only — no
# ADCP/ADNP/Bmem; those FIM x analyte cells will simply be empty.
BUBBLE_ANTIGENS <- ANTIGEN_ORDER   # canonical order (PT, FHA, PRN, FIM, TT, DT)
BUBBLE_ANALYTES <- c("IgG", "IgG3", "IgG1", "IgG2", "IgG4",
                     "FcgR2a", "FcgR3b", "ADCD", "ADCP", "ADNP")
BUBBLE_VISITS   <- c("InfMon2", "InfMon5", "InfMon9")   # recoded t0 / t1 / t2
BUBBLE_IGGLIST  <- c("DT_IgG", "FHA_IgG", "PRN_IgG", "PT_IgG", "TT_IgG",
                     "FIM_IgG")                          # PATCH: +FIM_IgG

# Features included in the grid (antigen_analyte). Carried verbatim from the
# original analysis (note DT/PRN have no Bmem; Bmem is excluded from the grid).
BUBBLE_INCLUDE_FEATURE <- c(
  "DT_IgG","DT_IgG1","DT_IgG2","DT_IgG3","DT_IgG4","DT_ADCD","DT_ADNP","DT_ADCP","DT_FcgR2a","DT_FcgR3b",
  "TT_IgG","TT_IgG1","TT_IgG2","TT_IgG3","TT_IgG4","TT_ADCD","TT_ADNP","TT_ADCP","TT_FcgR2a","TT_FcgR3b","TT_Bmem",
  "PT_IgG","PT_IgG1","PT_IgG2","PT_IgG3","PT_IgG4","PT_ADCD","PT_ADNP","PT_ADCP","PT_FcgR2a","PT_FcgR3b","PT_Bmem",
  "FHA_IgG","FHA_IgG1","FHA_IgG2","FHA_IgG3","FHA_IgG4","FHA_ADCD","FHA_ADNP","FHA_ADCP","FHA_FcgR2a","FHA_FcgR3b","FHA_Bmem",
  "PRN_IgG","PRN_IgG1","PRN_IgG2","PRN_IgG3","PRN_IgG4","PRN_ADCD","PRN_ADNP","PRN_ADCP","PRN_FcgR2a","PRN_FcgR3b",
  "FIM_IgG","FIM_IgG1","FIM_IgG2","FIM_IgG3","FIM_IgG4","FIM_ADCD","FIM_FcgR2a","FIM_FcgR3b")   # PATCH: +FIM (no ADCP/ADNP/Bmem)

# Same list without the total-IgG features (used by the standardization step).
BUBBLE_INCLUDE_FEATURE_WO <- c(
  "DT_IgG1","DT_IgG2","DT_IgG3","DT_IgG4","DT_ADCD","DT_ADNP","DT_ADCP","DT_FcgR2a","DT_FcgR3b",
  "TT_IgG1","TT_IgG2","TT_IgG3","TT_IgG4","TT_ADCD","TT_ADNP","TT_ADCP","TT_FcgR2a","TT_FcgR3b","TT_Bmem",
  "PT_IgG1","PT_IgG2","PT_IgG3","PT_IgG4","PT_ADCD","PT_ADNP","PT_ADCP","PT_FcgR2a","PT_FcgR3b","PT_Bmem",
  "FHA_IgG1","FHA_IgG2","FHA_IgG3","FHA_IgG4","FHA_ADCD","FHA_ADNP","FHA_ADCP","FHA_FcgR2a","FHA_FcgR3b","FHA_Bmem",
  "PRN_IgG1","PRN_IgG2","PRN_IgG3","PRN_IgG4","PRN_ADCD","PRN_ADNP","PRN_ADCP","PRN_FcgR2a","PRN_FcgR3b",
  "FIM_IgG1","FIM_IgG2","FIM_IgG3","FIM_IgG4","FIM_ADCD","FIM_FcgR2a","FIM_FcgR3b")   # PATCH: +FIM (no IgG total, no ADCP/ADNP/Bmem)

# Bubble effect-direction colours (sign convention: effect = median(TT) - median(TdaP)).
BUBBLE_EFFECT_COLORS <- c("Lower TdaP" = "salmon", "Higher TdaP" = "skyblue", "No effect" = "darkgrey")

# IgG-standardized residual files (resolved at runtime via find_proj_file).
# PATCH NOTE: the IgG-standardized FIM bubbles only populate if these
# residual objects were built WITH FIM. If they predate the FIM addition,
# the untransformed FIM bubbles will render but the standardized FIM cells
# will be empty until the residuals are regenerated with FIM included.
BUBBLE_RESIDUAL_AP <- "data/igg_standard_residuals_ap_matpm_prevacvacc_k.RData"
BUBBLE_RESIDUAL_WP <- "data/igg_standard_residuals_wp_matpm_prevacvacc_k.RData"


# ---------------------------------------------------------------------
# 8. CASING VERIFICATION  (PATCH)
# ---------------------------------------------------------------------
# Explicit, abortive guard against the silent case-mismatch described in
# the patch header. Call ONCE per driver immediately after
# load_serology_data(), e.g.:
#
#     data_raw <- load_serology_data(find_proj_file("data/c_set.RData"))
#     verify_predictor_casing(data_raw)        # errors loudly on mismatch
#
# It checks the feature strings this config will actually request against
# the feature levels present in the loaded data, CASE-SENSITIVELY (which is
# how the chain/concurrent paths match). For any miss, it reports whether a
# different-case variant exists, so a casing problem is named, not guessed.
verify_predictor_casing <- function(data_long, feature_col = "feature",
                                    stop_on_fail = TRUE) {
  present <- as.character(unique(data_long[[feature_col]]))

  # Features this config will request from the maternal/cord + infant data.
  # (Bubble matching is case-insensitive, so it is not gated here; FcR/effector
  #  and FIM subclasses are infant-visit only, but the *strings* must still
  #  exist among the data's feature levels.)
  needed <- unique(c(
    ANTIGEN_TOTALS,                                   # incl. DT/TT/FIM totals
    ALL_SUBCLASSES, FIM_SUBCLASSES,                   # chain + network subclasses
    DT_SUBCLASSES, TT_SUBCLASSES                      # available; harmless to check
  ))

  missing <- setdiff(needed, present)
  found   <- intersect(needed, present)

  # For each miss, see whether ONLY a different-case form is present.
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
