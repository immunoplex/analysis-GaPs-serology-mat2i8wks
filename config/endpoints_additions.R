# =============================================================================
# ADDITIONS TO config/endpoints.R  (v2 — corrected for actual data vocabulary)
#
# Data confirmed long-format: one row per subject × visit × feature.
# Key columns: subject_accession | visit_name | feature | log_assay_value
#              maternal_arm | infant_arm
#
# Paste SECTION A constants near the top of config/endpoints.R (alongside
# VISIT_RECODE). Merge SECTION B fields into each endpoint list inside
# get_endpoint_config(). SECTION C holds a pivot helper used by Parts 4 & 5.
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# SECTION A  —  Feature-group constants (correct data vocabulary)
# ─────────────────────────────────────────────────────────────────────────────

## Core antigen panel ----------------------------------------------------------
## Use the canonical display order from config/endpoints.R when it has been
## sourced first (it always is, by both Part 4/5 and results_B_to_F). This
## removes the previous order disagreement (was PT,PRN,FHA,DT,TT,FIM here vs
## PT,FHA,PRN,FIM,TT,DT in endpoints.R). The 6-antigen ANTIGENS used by
## Parts 4/5 intentionally shadows the core-3 ANTIGENS that endpoints.R uses
## for add_concurrent_ratios(); they are only both in scope inside Parts 4/5.
if (!exists("ANTIGEN_ORDER"))
  ANTIGEN_ORDER <- c("PT", "FHA", "PRN", "FIM", "TT", "DT")
ANTIGENS <- ANTIGEN_ORDER

## IgG subclasses (uppercase in data, no hyphen) --------------------------------
SUBCLASSES            <- c("IGG1", "IGG2", "IGG3", "IGG4")
ANTIGEN_SUBCLASS_FEATS <- as.vector(outer(ANTIGENS, SUBCLASSES, paste, sep = "_"))
# 24 features: PT_IGG1 ... FIM_IGG4

## FcG/Fc-binding — only FCGR2A and FCGR3B are in this dataset ----------------
FCG_RECEPTORS <- c("FCGR2A", "FCGR3B")
FCG_FEATS     <- as.vector(outer(ANTIGENS, FCG_RECEPTORS, paste, sep = "_"))
# 12 features: PT_FCGR2A, PT_FCGR3B, ... FIM_FCGR3B

## Effector functions ----------------------------------------------------------
# ADCD  → available for ALL 6 antigens
# ADCP, ADNP → available for PT, PRN, DT, TT ONLY (FHA and FIM are ADCD-only)
ANTIGENS_FULL_EFF  <- c("PT", "PRN", "DT", "TT")   # have ADCD + ADCP + ADNP
ANTIGENS_ADCD_ONLY <- c("FHA", "FIM")               # ADCD only

ADCD_FEATS         <- paste0(ANTIGENS, "_ADCD")                              # 6
ADCP_ADNP_FEATS    <- as.vector(outer(ANTIGENS_FULL_EFF, c("ADCP","ADNP"),
                                      paste, sep = "_"))                      # 8
EFFECTOR_FEATS     <- c(ADCD_FEATS, ADCP_ADNP_FEATS)                        # 14 total


# ─────────────────────────────────────────────────────────────────────────────
# SECTION B  —  Per-endpoint additions to get_endpoint_config()
#
# Add these five fields to each endpoint's list inside get_endpoint_config().
# Confirm the recoded INFMon2 visit label matches your VISIT_RECODE mapping
# (check: paste(paste0(names(VISIT_RECODE),"→",VISIT_RECODE), collapse=", "))
# ─────────────────────────────────────────────────────────────────────────────

## SBA -------------------------------------------------------------------------
#   outcome_feature        = "WHOLE_SBA",
#   antigen_subclass_feats = ANTIGEN_SUBCLASS_FEATS,
#   fcg_feats              = FCG_FEATS,
#   effector_feats         = EFFECTOR_FEATS,
#   infmon2_visit          = "INFMon2",   # recoded visit label — verify against VISIT_RECODE
#   cord_predictor         = NULL,
#   cord_visit             = NULL,
#   cord_predictor_label   = NULL

## WT_IgG  (cord blood Block 10 is active for this endpoint only) -------------
#   outcome_feature        = "WHOLE_WT_IgG",
#   antigen_subclass_feats = ANTIGEN_SUBCLASS_FEATS,
#   fcg_feats              = FCG_FEATS,
#   effector_feats         = EFFECTOR_FEATS,
#   infmon2_visit          = "INFMon2",   # recoded visit label — verify
#   cord_predictor         = "WHOLE_WT_IgG",  # same feature, different visit
#   cord_visit             = "Delivery",  # ← VERIFY: recoded name for the cord visit
#   cord_predictor_label   = "Cord blood WT IgG"

## PTNA ------------------------------------------------------------------------
#   outcome_feature        = "WHOLE_PTNA",
#   antigen_subclass_feats = ANTIGEN_SUBCLASS_FEATS,
#   fcg_feats              = FCG_FEATS,
#   effector_feats         = EFFECTOR_FEATS,
#   infmon2_visit          = "INFMon2",   # recoded visit label — verify
#   cord_predictor         = NULL,
#   cord_visit             = NULL,
#   cord_predictor_label   = NULL


# ─────────────────────────────────────────────────────────────────────────────
# SECTION C  —  Pivot helper
#
# Call make_wide() wherever a child part needs a participant × feature matrix.
# data_long  : the long-format data object (data_raw after load_serology_data)
# visit_label: recoded visit name to filter on (e.g. cfg$infmon2_visit)
# features   : character vector of feature values to keep (NULL = keep all)
# ─────────────────────────────────────────────────────────────────────────────

make_wide <- function(data_long,
                      visit_label,
                      features = NULL,
                      visit_col   = "visit_name",
                      feature_col = "feature",
                      value_col   = "log_assay_value",
                      id_cols     = c("subject_accession",
                                      "maternal_arm", "infant_arm")) {
  stopifnot(requireNamespace("tidyr", quietly = TRUE),
            requireNamespace("dplyr", quietly = TRUE))

  df <- dplyr::filter(data_long, .data[[visit_col]] == visit_label)

  if (!is.null(features))
    df <- dplyr::filter(df, .data[[feature_col]] %in% features)

  tidyr::pivot_wider(
    df,
    id_cols     = dplyr::all_of(id_cols),
    names_from  = dplyr::all_of(feature_col),
    values_from = dplyr::all_of(value_col),
    values_fn   = mean    # averages any technical duplicates
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION D  —  Column-name verification (run interactively)
# ─────────────────────────────────────────────────────────────────────────────

verify_new_feature_columns <- function(data_long,
                                       feature_col = "feature",
                                       verbose = TRUE) {
  present <- unique(data_long[[feature_col]])
  check   <- function(feats, label) {
    found   <- intersect(feats, present)
    missing <- setdiff(feats, present)
    if (verbose)
      message(label, ": ", length(found), " / ", length(feats), " found",
              if (length(missing))
                paste0("\n  MISSING: ",
                       paste(missing[seq_len(min(10L, length(missing)))],
                             collapse = ", "),
                       if (length(missing) > 10) " ...")
              else " — all present")
    invisible(list(found = found, missing = missing))
  }
  list(
    antigen_subclass = check(ANTIGEN_SUBCLASS_FEATS, "Antigen × IgG subclass"),
    fcg              = check(FCG_FEATS,              "FcG binding"),
    effector         = check(EFFECTOR_FEATS,         "Effector functions")
  )
}

## ---------------------------------------------------------------------------
## C0 data-harmonisation / standardisation constants
## Sourced by analysis/C0_data_harmonisation.Rmd and its child fragments.
## ANTIGENS <- ANTIGEN_ORDER is assumed already defined in endpoints.R.
## ---------------------------------------------------------------------------

## Antigen / IgG vocabularies, keyed by window.
## Antigen strings match ebaa_extra2$antigen EXACTLY (verified against d_set.RData).
VOCAB <- list(
  ## Infant window (M02/M05/M09), both arms: 9 antigens, full analyte panel
  full = list(
    antigenlist = c("DT","FHA","FIM","PRN","PT","TT","IPV1","IPV2","IPV3"),
    igglist     = c("DT_IgG","FHA_IgG","FIM_IgG","PRN_IgG","PT_IgG","TT_IgG",
                    "IPV1_IgG","IPV2_IgG","IPV3_IgG"),
    fulligglist = c("PT_IgG","PRN_IgG","FHA_IgG","FIM_IgG","TT_IgG","DT_IgG",
                    "IPV1_IgG","IPV2_IgG","IPV3_IgG")
  ),
  ## Maternal (P02/P09/M00) + cord/transfer (M00/M02/M05) windows: core 5 antigens
  core5 = list(
    antigenlist = c("DT","FHA","PRN","PT","TT"),
    igglist     = c("DT_IgG","FHA_IgG","PRN_IgG","PT_IgG","TT_IgG"),
    fulligglist = c("PT_IgG","PRN_IgG","FHA_IgG","TT_IgG","DT_IgG")
  )
)

## Whole-bacterium functional assays are standardised against a SINGLE, biologically
## chosen antigen's total IgG -- NOT their own antigen:
##   SBA    ~ PRN total IgG
##   WT-IgG ~ PRN total IgG
##   PTNA   ~ PT  total IgG
## The engine reassigns each WHOLE_* row's `antigen` to igg_antigen before the
## per-antigen fit loop, so it is regressed on (and displayed under) that antigen.
wholeassay_map <- tibble::tribble(
  ~source_feature, ~igg_antigen, ~out_feature,
  "WHOLE_SBA",     "PRN",        "PRN_SBA",
  "WHOLE_WT_IgG",  "PRN",        "PRN_WT",
  "WHOLE_PTNA",    "PT",         "PT_PTNA"
)

## Driver table: one row per arm x window. `vocab` selects a VOCAB entry.
## `suffix` builds the output filenames; residual/modelfit names are held
## IDENTICAL to what C1/C4/C5/C7 already expect (note wP-inf stays `_k`, not `_j`).
std_config <- tibble::tribble(
  ~arm, ~window, ~t0,   ~t1,   ~t2,   ~suffix,  ~vocab,
  "aP", "inf",   "M02", "M05", "M09", "k",      "full",
  "wP", "inf",   "M02", "M05", "M09", "k",      "full",
  "aP", "mat",   "P02", "P09", "M00", "k_mat",  "core5",
  "wP", "mat",   "P02", "P09", "M00", "k_mat",  "core5",
  "aP", "trn",   "M00", "M02", "M05", "k_trn",  "core5",
  "wP", "trn",   "M00", "M02", "M05", "k_trn",  "core5"
)

## All C0 intermediate outputs live in data/ (git-ignored).
c0_out <- function(...) here::here("data", ...)

# Usage:
#   source("config/endpoints_additions.R")
#   verify_new_feature_columns(ebaa_extra2)   # or data_raw
