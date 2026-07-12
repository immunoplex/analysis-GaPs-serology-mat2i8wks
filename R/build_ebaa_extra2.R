## =====================================================================
## build_ebaa_extra2()
## ---------------------------------------------------------------------
## Returns a data.frame that contains EVERY row of `ebaa_extra` (the object
## stored in c_set.RData) plus appended ADCD rows taken from `mat_baa_norm`.
##
## Rows are pulled from mat_baa_norm where ALL of the following hold:
##     experiment_accession == "ADCD"
##     feature              == "ADCD"
##     source               == "NIBSC"                 (one row per key)
##     visit_name %in% c("P18","P09","M00","P02","P05")
##
## The appended rows are translated into the ebaa_extra schema following the
## same conventions used in format_raw_data.R:
##
##   subject_accession <- as.numeric(patientid)                # line 140
##   sample_type       <- sample_type                          # already "M"/"C" (line 102)
##   assay             <- "mbaa"                                # bead-array (line 377)
##   visit_name        <- visit_name                           # already P02/P09/P18/M00 (lines 96-99)
##   antigen           <- toupper(antigen)                      # "Pentamer" -> "PENTAMER" (line 253)
##   analyte           <- "ADCD"                                # the feature carried the isotype (line 226)
##   feature           <- toupper(paste(antigen,"ADCD",sep="_"))# e.g. "PT_ADCD" (lines 227,252)
##   value_reported    <- mfi_inf_raw                           # (your chosen value column)
##   log_assay_value   <- log10(value_reported)
##   arm_name / maternal_arm / infant_arm <- joined from ebaa_extra by subject_accession
##
## NOTE on de-duplication: mat_baa_norm holds TWO normalised rows per
## (subject, visit, antigen) — one for source "NIBSC" and one for "SD".
## Restricting to source == "NIBSC" yields exactly one row per key, matching
## the one-row-per-key structure of the ADCD data already in ebaa_extra.
## =====================================================================

build_ebaa_extra2 <- function(ebaa_extra,
                              mat_baa_norm,
                              source    = "NIBSC",
                              value_col = "mfi_inf_raw",
                              visits    = c("P18", "P09", "M00", "P02", "P05")) {

  ## ---- 0. sanity checks -------------------------------------------------
  target_cols <- c("subject_accession", "sample_type", "assay", "visit_name",
                   "antigen", "analyte", "feature", "value_reported",
                   "log_assay_value", "arm_name", "maternal_arm", "infant_arm")
  stopifnot(all(target_cols %in% names(ebaa_extra)))

  need <- c("experiment_accession", "feature", "source", "visit_name",
            "patientid", "antigen", "sample_type", value_col)
  missing_cols <- setdiff(need, names(mat_baa_norm))
  if (length(missing_cols))
    stop("mat_baa_norm is missing required column(s): ",
         paste(missing_cols, collapse = ", "))

  ## ---- 1. filter mat_baa_norm ------------------------------------------
  ## as.character() guards against visit_name / source being factors.
  keep <- mat_baa_norm$experiment_accession == "ADCD" &
          mat_baa_norm$feature              == "ADCD" &
          as.character(mat_baa_norm$source) == source &
          as.character(mat_baa_norm$visit_name) %in% visits

  src <- mat_baa_norm[keep, , drop = FALSE]
  if (nrow(src) == 0L)
    warning("No rows in mat_baa_norm matched the ADCD / source / visit filter.")

  ## ---- 2. translate into the ebaa_extra schema -------------------------
  antigen_up <- toupper(sub("_[0-9]+$", "", as.character(src$antigen)))  # strip trailing _<n>, upper-case
  value_rep  <- as.numeric(src[[value_col]])

  if (any(!is.na(value_rep) & value_rep <= 0))
    warning("Some ", value_col,
            " values are <= 0; log10() will yield -Inf/NaN for those rows.")

  add <- data.frame(
    subject_accession = as.numeric(src$patientid),
    sample_type       = as.character(src$sample_type),
    assay             = "mbaa",
    visit_name        = as.character(src$visit_name),
    antigen           = antigen_up,
    analyte           = "ADCD",
    feature           = paste(antigen_up, "ADCD", sep = "_"),
    value_reported    = value_rep,
    log_assay_value   = log10(value_rep),
    stringsAsFactors  = FALSE
  )

  ## ---- 3. attach arm designations from ebaa_extra ----------------------
  ## ebaa_extra already carries a clean subject -> arm mapping; reuse it so the
  ## function stays self-contained (no external subject spreadsheet needed).
  arm_lookup <- unique(ebaa_extra[, c("subject_accession", "arm_name",
                                      "maternal_arm", "infant_arm")])
  dup <- arm_lookup$subject_accession[duplicated(arm_lookup$subject_accession)]
  if (length(dup))
    warning("Multiple arm assignments found for subject_accession(s): ",
            paste(sort(unique(dup)), collapse = ", "))

  add <- merge(add, arm_lookup, by = "subject_accession", all.x = TRUE, sort = FALSE)

  n_no_arm <- sum(is.na(add$arm_name))
  if (n_no_arm)
    warning(n_no_arm, " appended ADCD row(s) had no matching arm assignment ",
            "in ebaa_extra (arm columns are NA).")

  ## ---- 4. match column order & classes to ebaa_extra -------------------
  add <- add[, target_cols, drop = FALSE]

  ## Keep factor columns as factors, aligned to ebaa_extra's levels. If the
  ## appended data introduces a level not already present, extend the levels
  ## (with a warning) rather than silently turning the value into NA.
  align_factor <- function(new_vals, template) {
    if (!is.factor(template)) return(new_vals)
    lv  <- levels(template)
    new <- setdiff(unique(new_vals[!is.na(new_vals)]), lv)
    if (length(new)) {
      warning("New level(s) added while appending: ",
              paste(new, collapse = ", "))
      lv <- c(lv, new)
    }
    factor(new_vals, levels = lv)
  }
  add$visit_name <- align_factor(add$visit_name, ebaa_extra$visit_name)
  add$antigen    <- align_factor(add$antigen,    ebaa_extra$antigen)

  ## Coerce any remaining columns to the class ebaa_extra uses (character vs
  ## numeric) so rbind() does not complain.
  for (nm in target_cols) {
    if (is.factor(ebaa_extra[[nm]])) next
    if (is.numeric(ebaa_extra[[nm]]))    add[[nm]] <- as.numeric(add[[nm]])
    else if (is.character(ebaa_extra[[nm]])) add[[nm]] <- as.character(add[[nm]])
  }

  ## ---- 5. combine ------------------------------------------------------
  ebaa_extra2 <- rbind(ebaa_extra[, target_cols, drop = FALSE], add)
  rownames(ebaa_extra2) <- NULL

  message(sprintf("ebaa_extra: %d rows  +  appended ADCD (%s / %s): %d rows  =  ebaa_extra2: %d rows",
                  nrow(ebaa_extra), source, value_col, nrow(add), nrow(ebaa_extra2)))
  ebaa_extra2
}


## =====================================================================
## Example usage
## =====================================================================
load(here::here("./data/c_set.RData"))        # provides `ebaa_extra`
load(here::here("./data/mat_baa_norm.RData")) # provides `mat_baa_norm`

ebaa_extra2 <- build_ebaa_extra2(ebaa_extra, mat_baa_norm)

# Quick checks:
table(ebaa_extra2$visit_name)
table(ebaa_extra2$feature[ebaa_extra2$analyte == "ADCD"],
      ebaa_extra2$visit_name[ebaa_extra2$analyte == "ADCD"])
summary(ebaa_extra2)
save(ebaa_extra2, file = here::here("./data/d_set.RData"))
