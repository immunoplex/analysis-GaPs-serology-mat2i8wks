# diagnose_censoring.R
# ===========================================================================
# Read-only diagnostic on the normalized ADCD output (mat_baa_norm).
#
# Question it answers: is the below-floor (left-censored) clamping DIFFERENTIAL
# across timepoints/antigens? If censoring rates are similar across the groups
# you plan to compare, "filter to oor_flag == 'ok' and compare" is defensible.
# If they differ markedly, a naive comparison of the in-range values (or of the
# clamped values) will confound detectability with biology, and a left-censored
# method is needed instead.
#
# Also shows HOW FAR below the infant floor the censored samples sit, which
# tells you whether the floored points are marginally low (censored-analysis
# worth it) or effectively zero (treat as absent).
#
# Run after mat_baa_norm exists. READ-ONLY. Copy the whole output back.
# ===========================================================================

suppressMessages({library(dplyr); library(tidyr)})
# Defeat papeR masking regardless of load order
summarise <- dplyr::summarise; summarize <- dplyr::summarize
mutate <- dplyr::mutate; filter <- dplyr::filter
count <- dplyr::count; arrange <- dplyr::arrange; group_by <- dplyr::group_by

rule <- function(t) cat("\n", strrep("=",70), "\n", t, "\n", strrep("=",70), "\n", sep="")
stopifnot(exists("mat_baa_norm"))

A <- mat_baa_norm %>% filter(feature == "ADCD")

# A single, ordered timepoint label. Prefer 'timeperiod'; fall back to sample_type.
tp_col <- if ("timeperiod" %in% names(A)) "timeperiod" else "sample_type"
A$tp <- A[[tp_col]]
cat("Using timepoint column:", tp_col, "\n")
cat("Timepoint levels:", paste(sort(unique(A$tp)), collapse = ", "), "\n")
cat("Total ADCD rows:", nrow(A), "\n")

# ---------------------------------------------------------------------------
# 1. Censoring rate by TIMEPOINT (the confounding check)
# ---------------------------------------------------------------------------
rule("1. CENSORING RATE BY TIMEPOINT")
by_tp <- A %>%
  group_by(tp) %>%
  summarise(
    n            = n(),
    ok           = sum(oor_flag == "ok"),
    below        = sum(oor_flag == "below_clamped"),
    above        = sum(oor_flag == "above_clamped"),
    no_transfer  = sum(oor_flag == "no_transfer"),
    pct_below    = round(100 * below / n, 1),
    pct_ok       = round(100 * ok    / n, 1),
    .groups = "drop"
  ) %>% arrange(desc(pct_below))
print(as.data.frame(by_tp), row.names = FALSE)
cat("\n  -> If pct_below varies a lot across the timepoints you will COMPARE,\n")
cat("     censoring is differential and a naive comparison is biased.\n")
cat("     Spread in pct_below:",
    round(diff(range(by_tp$pct_below)), 1), "percentage points.\n")

# ---------------------------------------------------------------------------
# 2. Censoring rate by ANTIGEN x TIMEPOINT (finer view)
# ---------------------------------------------------------------------------
rule("2. % BELOW-CLAMPED BY ANTIGEN x TIMEPOINT")
at <- A %>%
  group_by(antigen, tp) %>%
  summarise(n = n(), pct_below = round(100 * mean(oor_flag == "below_clamped"), 0),
            .groups = "drop") %>%
  tidyr::pivot_wider(names_from = tp, values_from = pct_below, values_fill = NA,
                     id_cols = antigen)
print(as.data.frame(at), row.names = FALSE)
cat("\n  (cells are % of that antigen x timepoint that are below-floor censored)\n")

# ---------------------------------------------------------------------------
# 3. HOW FAR below the floor do censored samples sit?
#    Compare each below-clamped sample's own MFI to the infant floor it was
#    clamped to. gap_log10 = log10(inf_floor) - log10(sample_mfi) >= 0.
#    Small gap -> just under LLOQ (censored methods meaningful);
#    large gap -> effectively absent.
# ---------------------------------------------------------------------------
rule("3. DEPTH BELOW FLOOR (log10 units) FOR CENSORED SAMPLES")
bc <- A %>%
  filter(oor_flag == "below_clamped",
         is.finite(mfi), is.finite(inf_mfi_min), mfi > 0) %>%
  mutate(gap_log10 = log10(inf_mfi_min) - log10(mfi))
if (nrow(bc) == 0) {
  cat("  no below-clamped samples.\n")
} else {
  cat(sprintf("  n below-clamped with usable mfi: %d\n", nrow(bc)))
  cat("  gap_log10 (how many log10 units below the infant floor):\n")
  print(summary(bc$gap_log10))
  cat("\n  gap distribution (log10 bands):\n")
  br <- cut(bc$gap_log10, breaks = c(-Inf, 0.1, 0.3, 0.5, 1, Inf),
            labels = c("<0.1 (~just under)", "0.1-0.3", "0.3-0.5", "0.5-1", ">1 (~10x under)"))
  print(as.data.frame(table(band = br)), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# 4. In-range (ok) signal spread by timepoint -- is there dynamic range left?
#    If the 'ok' values are themselves nearly flat, even the usable fraction
#    carries little information.
# ---------------------------------------------------------------------------
rule("4. IN-RANGE (ok) infant-scale MFI SPREAD BY TIMEPOINT")
okr <- A %>%
  filter(oor_flag == "ok", is.finite(mfi_inf_scale), mfi_inf_scale > 0) %>%
  group_by(tp) %>%
  summarise(n = n(),
            median_log10 = round(median(log10(mfi_inf_scale)), 2),
            iqr_log10    = round(IQR(log10(mfi_inf_scale)), 2),
            min_log10    = round(min(log10(mfi_inf_scale)), 2),
            max_log10    = round(max(log10(mfi_inf_scale)), 2),
            .groups = "drop")
print(as.data.frame(okr), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Paired-subject view: for within-subject comparisons (e.g. cord vs
#    maternal delivery, or prevacc vs delivery), how many subjects have a
#    NON-censored value at BOTH timepoints? Censoring at either end drops the pair.
# ---------------------------------------------------------------------------
rule("5. PAIR RETENTION FOR WITHIN-SUBJECT COMPARISONS")
subj_col <- intersect(c("patientid","subject_accession","sampleid"), names(A))[1]
if (is.na(subj_col)) {
  cat("  no subject id column found; skipping pair retention.\n")
} else {
  cat("  subject id column:", subj_col, "\n")
  A$subj <- A[[subj_col]]
  # per subject x antigen, which timepoints are 'ok'
  pr <- A %>%
    group_by(subj, antigen) %>%
    summarise(
      tps_ok = paste(sort(unique(tp[oor_flag == "ok"])), collapse = "+"),
      n_tp_ok = length(unique(tp[oor_flag == "ok"])),
      .groups = "drop")
  cat("\n  distribution of # timepoints with an in-range ADCD value, per subject x antigen:\n")
  print(as.data.frame(table(n_timepoints_ok = pr$n_tp_ok)), row.names = FALSE)
  cat("\n  (subject x antigen combos with >=2 in-range timepoints can support a\n")
  cat("   within-subject change; those with <2 cannot without censored handling.)\n")
}

rule("END - copy everything above back")
