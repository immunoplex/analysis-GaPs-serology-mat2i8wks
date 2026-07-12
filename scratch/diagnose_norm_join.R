# diagnose_norm_join.R
# ---------------------------------------------------------------------------
# Diagnoses why 9,063 / 49,323 mat_baa rows failed the normalization_params
# left_join (NA mfi_norm). Run this AFTER the line that creates
# mat_baa_normalised (i.e. it expects mat_baa_normalised, mat_baa, and
# norm_lookup to already be in the session).
#
# It prints a series of labelled blocks. Copy the WHOLE console output back.
# Nothing is modified; this is read-only.
# ---------------------------------------------------------------------------

suppressMessages({library(dplyr)})

rule <- function(txt) cat("\n", strrep("=", 70), "\n", txt, "\n",
                          strrep("=", 70), "\n", sep = "")

# --- 0. Objects present? ----------------------------------------------------
rule("0. OBJECTS IN SESSION")
for (o in c("mat_baa_normalised", "mat_baa", "norm_lookup")) {
  cat(sprintf("  %-20s exists=%s\n", o, exists(o)))
}
stopifnot(exists("mat_baa_normalised"), exists("norm_lookup"))

# --- 1. Overall join coverage ----------------------------------------------
rule("1. OVERALL JOIN COVERAGE")
n_tot <- nrow(mat_baa_normalised)
n_na  <- sum(is.na(mat_baa_normalised$mfi_norm))
cat(sprintf("  total rows       : %d\n", n_tot))
cat(sprintf("  NA mfi_norm      : %d (%.1f%%)\n", n_na, 100 * n_na / n_tot))
cat(sprintf("  matched          : %d (%.1f%%)\n", n_tot - n_na,
            100 * (n_tot - n_na) / n_tot))

# --- 2. Which KEY COLUMN is responsible: per-column setdiff -----------------
# Sample-side key values with NO counterpart in norm_lookup = guaranteed misses
rule("2. SAMPLE-SIDE KEY VALUES WITH NO MATCH IN norm_lookup (setdiff)")

show_setdiff <- function(sample_vals, param_vals, label) {
  miss <- sort(setdiff(unique(as.character(sample_vals)),
                       unique(as.character(param_vals))))
  cat(sprintf("\n  [%s] %d sample value(s) absent from params:\n", label, length(miss)))
  if (length(miss)) cat("    ", paste(miss, collapse = ", "), "\n")
  else              cat("     (none - all sample values exist in params)\n")
}

# feature
show_setdiff(mat_baa_normalised$feature, norm_lookup$feature, "feature")
# antigen_base  (param side column is 'antigen_base')
if ("antigen_base" %in% names(norm_lookup)) {
  show_setdiff(mat_baa_normalised$antigen_base, norm_lookup$antigen_base, "antigen_base")
} else {
  cat("\n  [antigen_base] WARNING: norm_lookup has no 'antigen_base' column! names:\n    ",
      paste(names(norm_lookup), collapse = ", "), "\n")
}
# plate
show_setdiff(mat_baa_normalised$plate, norm_lookup$plate, "plate")
# source_join  (param side column is 'source')
show_setdiff(mat_baa_normalised$source_join, norm_lookup$source, "source_join vs params$source")

# --- 3. Value inventories (so mismatches in casing/format are visible) ------
rule("3. VALUE INVENTORIES (sample side vs params side)")

cat("\n  -- feature --\n")
cat("  sample:", paste(sort(unique(as.character(mat_baa_normalised$feature))), collapse = ", "), "\n")
cat("  params:", paste(sort(unique(as.character(norm_lookup$feature))),        collapse = ", "), "\n")

cat("\n  -- source (sample source_join vs params source) --\n")
cat("  sample:", paste(sort(unique(as.character(mat_baa_normalised$source_join))), collapse = ", "), "\n")
cat("  params:", paste(sort(unique(as.character(norm_lookup$source))),             collapse = ", "), "\n")

cat("\n  -- plate --\n")
cat("  sample:", paste(sort(unique(as.character(mat_baa_normalised$plate))), collapse = ", "), "\n")
cat("  params:", paste(sort(unique(as.character(norm_lookup$plate))),        collapse = ", "), "\n")

cat("\n  -- antigen_base (sample) --\n")
cat("  sample:", paste(sort(unique(as.character(mat_baa_normalised$antigen_base))), collapse = ", "), "\n")
if ("antigen_base" %in% names(norm_lookup))
  cat("  params:", paste(sort(unique(as.character(norm_lookup$antigen_base))), collapse = ", "), "\n")

# --- 4. Where are the failures concentrated? (per-key crosstabs) ------------
rule("4. FAILURES BY KEY (n and n_failed per level)")

fail_by <- function(col) {
  if (!col %in% names(mat_baa_normalised)) {
    cat(sprintf("\n  [%s] not present\n", col)); return(invisible())
  }
  tab <- mat_baa_normalised %>%
    mutate(failed = is.na(mfi_norm)) %>%
    group_by(.data[[col]]) %>%
    summarise(n = n(), n_failed = sum(failed), .groups = "drop") %>%
    mutate(pct_failed = round(100 * n_failed / n, 1)) %>%
    arrange(desc(n_failed))
  cat(sprintf("\n  -- failures by %s --\n", col))
  print(as.data.frame(tab), row.names = FALSE)
}
fail_by("feature")
fail_by("source_join")
fail_by("plate")

# --- 5. Distinct FULL failing key combinations ------------------------------
rule("5. DISTINCT FAILING KEY COMBINATIONS (up to 60 shown)")
unmatched <- mat_baa_normalised %>%
  filter(is.na(mfi_norm)) %>%
  distinct(feature, antigen_base, plate, source_join) %>%
  arrange(feature, antigen_base, plate, source_join)
cat(sprintf("  %d distinct failing (feature, antigen_base, plate, source_join) combos\n",
            nrow(unmatched)))
print(head(as.data.frame(unmatched), 60), row.names = FALSE)

# --- 6. Upstream source filter check (the suspected SD-drop bug) ------------
rule("6. UPSTREAM SOURCE FILTER CHECK")
cat("  norm_lookup$source levels :",
    paste(sort(unique(as.character(norm_lookup$source))), collapse = ", "), "\n")
cat("  mat_baa$source levels     :",
    if (exists("mat_baa")) paste(sort(unique(as.character(mat_baa$source))), collapse = ", ") else "mat_baa not present", "\n")
cat("  NOTE: if params expect 'SD' but mat_baa$source has no SD/Sando rows,\n")
cat("        the SD rows were dropped upstream by the %in% c('NIBSC06_140','Sando') filter.\n")

# --- 7. norm_lookup structure ----------------------------------------------
rule("7. norm_lookup STRUCTURE")
cat("  dim:", paste(dim(norm_lookup), collapse = " x "), "\n")
cat("  names:", paste(names(norm_lookup), collapse = ", "), "\n")
cat("  duplicated join keys in params (feature,antigen_base,plate,source):\n")
if (all(c("feature","antigen_base","plate","source") %in% names(norm_lookup))) {
  dup <- norm_lookup %>% count(feature, antigen_base, plate, source) %>% filter(n > 1)
  cat(sprintf("    %d duplicated key combos (>1 param row -> row duplication on join)\n", nrow(dup)))
  if (nrow(dup)) print(head(as.data.frame(dup), 20), row.names = FALSE)
} else {
  cat("    (cannot check - expected key columns missing from norm_lookup)\n")
}

rule("END OF DIAGNOSTIC - copy everything above back")
