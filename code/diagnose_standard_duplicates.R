# diagnose_standard_duplicates.R
# ---------------------------------------------------------------------------
# WHY: a plate should carry exactly ONE standard reading per
#   (set, plate, antigen, feature, source, dilution).
# But standards_sel shows MULTIPLE distinct MFI values at some such keys
# (e.g. plate_1 / PT / ADCD / NIBSC06_140 / dilution 10 -> 841 AND 1758).
# This script finds every duplicated key, quantifies how far apart the
# duplicate MFIs are, and looks for a hidden grouping that explains them
# (e.g. two curve runs, two wavelengths, two plates collapsed, a wells axis
# that got dropped). It is READ-ONLY.
#
# Run after standards_sel is loaded. Copy the whole console output back.
# ---------------------------------------------------------------------------

suppressMessages({library(dplyr); library(tidyr)})

# papeR (loaded elsewhere in the session) masks dplyr::summarise/summarize and
# breaks grouped summaries silently. Force the dplyr versions here regardless of
# load order. All summarise() calls below are already dplyr::summarise().
summarise <- dplyr::summarise
summarize <- dplyr::summarize
mutate    <- dplyr::mutate
filter    <- dplyr::filter
count     <- dplyr::count

rule <- function(t) cat("\n", strrep("=",70), "\n", t, "\n", strrep("=",70), "\n", sep="")

stopifnot(exists("standards_sel"))
S <- standards_sel

rule("0. STRUCTURE")
cat("  dim:", paste(dim(S), collapse=" x "), "\n")
cat("  names:", paste(names(S), collapse=", "), "\n")
cat("\n  class per column:\n")
for (nm in names(S)) cat(sprintf("    %-14s %s\n", nm, class(S[[nm]])[1]))

# The intended unique key for a standard point:
key_cols <- c("set","plate","antigen","feature","source","dilution")
key_cols <- key_cols[key_cols %in% names(S)]
cat("\n  key used:", paste(key_cols, collapse=" + "), "\n")

# --- 1. How many rows per key, and how many DISTINCT mfi per key ------------
rule("1. ROWS AND DISTINCT MFI PER KEY")
per_key <- S %>%
  group_by(across(all_of(key_cols))) %>%
  dplyr::summarise(n_rows = n(),
            n_distinct_mfi = n_distinct(mfi),
            .groups = "drop")

cat("  distribution of n_rows per key:\n")
print(as.data.frame(table(n_rows = per_key$n_rows)), row.names = FALSE)

cat("\n  distribution of n_distinct_mfi per key (THIS is the real question):\n")
print(as.data.frame(table(n_distinct_mfi = per_key$n_distinct_mfi)), row.names = FALSE)

cat(sprintf("\n  keys with >1 distinct mfi: %d of %d (%.1f%%)\n",
            sum(per_key$n_distinct_mfi > 1), nrow(per_key),
            100*mean(per_key$n_distinct_mfi > 1)))

# --- 2. Are the extra rows EXACT duplicates or DIFFERENT values? ------------
# If n_rows>1 but n_distinct_mfi==1 -> harmless exact repeats (e.g. long format).
# If n_distinct_mfi>1 -> genuinely conflicting readings at one standard point.
rule("2. EXACT-REPEAT vs CONFLICTING")
per_key <- per_key %>%
  mutate(kind = case_when(
    n_rows == 1                      ~ "single",
    n_distinct_mfi == 1              ~ "exact_repeat",
    TRUE                             ~ "conflicting"))
print(as.data.frame(table(kind = per_key$kind)), row.names = FALSE)

# --- 3. Characterise the CONFLICTING keys: spread of their mfi values -------
rule("3. SPREAD WITHIN CONFLICTING KEYS")
conf_keys <- per_key %>% filter(kind == "conflicting") %>% select(all_of(key_cols))
if (nrow(conf_keys) == 0) {
  cat("  none - no conflicting keys.\n")
} else {
  conf <- S %>% inner_join(conf_keys, by = key_cols) %>%
    group_by(across(all_of(key_cols))) %>%
    dplyr::summarise(n = n(),
              n_distinct_mfi = n_distinct(mfi),
              mfi_min = min(mfi), mfi_max = max(mfi),
              ratio_max_min = round(max(mfi)/pmax(min(mfi),1e-9), 2),
              .groups = "drop")
  cat(sprintf("  %d conflicting keys. Ratio (max/min mfi) summary:\n", nrow(conf)))
  print(summary(conf$ratio_max_min))
  cat("\n  first 25 conflicting keys:\n")
  print(head(as.data.frame(conf), 25), row.names = FALSE)
}

# --- 4. Look for a HIDDEN grouping that separates the duplicates ------------
# For each candidate column NOT in the key, test whether adding it makes the
# key unique (i.e. that column explains the duplication).
rule("4. DOES ANY OTHER COLUMN EXPLAIN THE DUPLICATES?")
other_cols <- setdiff(names(S), c(key_cols, "mfi", "concentration"))
if (length(other_cols) == 0) {
  cat("  no other columns available to explain duplicates.\n")
  cat("  -> duplication cannot be resolved from within standards_sel alone.\n")
} else {
  base_conf <- sum(per_key$kind == "conflicting")
  for (oc in other_cols) {
    k2 <- c(key_cols, oc)
    pk2 <- S %>% group_by(across(all_of(k2))) %>%
      dplyr::summarise(nd = n_distinct(mfi), .groups="drop")
    still_conf <- sum(pk2$nd > 1)
    cat(sprintf("  + %-20s : conflicting keys %d -> %d %s\n",
                oc, base_conf, still_conf,
                if (still_conf == 0) "***EXPLAINS IT***" else ""))
  }
  cat("\n  (A column that drives conflicting keys to 0 is the missing axis:\n")
  cat("   the standard was measured more than once per dilution, distinguished\n")
  cat("   by that column - it should be part of the key, or collapsed over.)\n")
}

# --- 5. Concentration vs dilution within a conflicting key ------------------
# Check whether the two mfi values also differ in 'concentration', which would
# mean the pair are actually different assigned potencies, not repeats.
rule("5. concentration BEHAVIOUR WITHIN CONFLICTING KEYS")
if (exists("conf_keys") && nrow(conf_keys) > 0) {
  ex <- S %>% inner_join(conf_keys[1,,drop=FALSE], by = key_cols) %>%
    arrange(mfi) %>%
    select(all_of(key_cols), all_of(intersect(c("concentration"), names(S))), mfi)
  cat("  one conflicting key, all its rows:\n")
  print(head(as.data.frame(ex), 60), row.names = FALSE)
}

# --- 6. Is the duplication tied to a specific set/feature? ------------------
rule("6. CONFLICTING KEYS BY set x feature")
if (nrow(conf_keys) > 0) {
  bt <- per_key %>% filter(kind=="conflicting") %>% count(set, feature)
  print(as.data.frame(bt), row.names = FALSE)
}

rule("END - copy everything above back")
