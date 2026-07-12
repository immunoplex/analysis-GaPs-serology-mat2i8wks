# normalize_mat_to_inf_mfi.R
# ===========================================================================
# Harmonise maternal/cord test-sample MFI onto the INFANT panel's MFI scale,
# per plate x antigen x feature x source, using the paired standard curves.
#
# METHOD (agreed):
#   * Signal is MFI = 10^assay_response (both sides; standards use their own mfi).
#   * For each plate x antigen(case-folded) x feature x source, pair the
#     maternal-set and infant-set standard points at matched dilution, giving
#     (MFI_mat, MFI_inf) couples. Both raw curves hook with dilution, but the
#     PAIRED relationship MFI_inf ~ f(MFI_mat) is near-monotone because the hook
#     cancels (same reference dilutions on both panels).
#   * Fit a MONOTONE-increasing transfer function  log10(MFI_inf) ~ f(log10(MFI_mat))
#     and predict, then back-transform to NATURAL MFI.
#       - >= 6 distinct paired points -> monotone smoothing spline
#                                        (scam mpi; base monoH.FC if scam absent)
#       - 3-5 distinct paired points  -> monotone LINEAR INTERPOLATION (no smoothing)
#       - < 3 distinct paired points  -> SKIP key, flag its test samples
#   * Apply f to each test sample's MFI (maternal panel) -> infant-scale MFI.
#   * CLAMP (Q.C): test MFI below the infant standard's per-key min -> infant min;
#     above the infant standard's per-key max -> infant max. (Clamp in MFI space
#     using the infant curve's own min/max MFI for that plate/antigen/feature/source.)
#
# KEY HARMONISATION:
#   * antigen: fold case so maternal 'PT' pairs with infant 'pt'  (antigen_key = tolower)
#   * source : recode sample 'NIBSC' -> 'NIBSC06_140' to match standards
#   * shared dilutions only (10-2560); maternal-only 5120 has no infant partner
#
# INPUTS in session : standards_sel, mat_baa   (inf_baa handled later, same fn)
# OUTPUT            : mat_baa_norm  (mat_baa + mfi, mfi_inf_scale, assay_response_inf,
#                                    transfer_method, oor_flag)
# ===========================================================================

suppressMessages({library(dplyr); library(tidyr)})

# Defeat papeR masking regardless of load order
summarise <- dplyr::summarise; summarize <- dplyr::summarize
mutate <- dplyr::mutate; filter <- dplyr::filter
count <- dplyr::count; arrange <- dplyr::arrange; group_by <- dplyr::group_by

SHARED_DIL_MAX <- 2560          # infant ladder tops out here; drop maternal 5120
MIN_SPLINE_PTS <- 6             # >=6 -> monotone spline
MIN_INTERP_PTS <- 3             # 3-5 -> monotone linear interpolation; <3 -> skip
HAVE_SCAM <- requireNamespace("scam", quietly = TRUE)
if (!HAVE_SCAM)
  message("NOTE: package 'scam' not installed; the >=6-point case will use a ",
          "monotone INTERPOLATING spline (splinefun monoH.FC) instead of a ",
          "monotone smoothing spline. Install 'scam' for smoothing.")

# ---------------------------------------------------------------------------
# 1. Build the standards long table, collapse exact repeats, split by set.
#    (Diagnostics confirmed exactly one MFI per key; distinct() is loss-free.)
# ---------------------------------------------------------------------------
std <- standards_sel %>%
  transmute(
    set,
    feature,
    antigen_key = tolower(antigen),
    plate,
    source,
    dilution,
    mfi
  ) %>%
  filter(is.finite(mfi), dilution <= SHARED_DIL_MAX) %>%
  distinct()                                   # collapse identical long-format repeats

std_mat <- std %>% filter(set == "mat_baa") %>%
  transmute(feature, antigen_key, plate, source, dilution, mfi_mat = mfi)
std_inf <- std %>% filter(set == "inf_baa") %>%
  transmute(feature, antigen_key, plate, source, dilution, mfi_inf = mfi)

# Paired couples at matched dilution, per plate x antigen x feature x source
pairs <- inner_join(std_mat, std_inf,
                    by = c("feature", "antigen_key", "plate", "source", "dilution")) %>%
  filter(is.finite(mfi_mat), is.finite(mfi_inf), mfi_mat > 0, mfi_inf > 0)

# Infant per-key MFI range for clamping (from the FULL infant ladder within shared dil)
inf_range <- std_inf %>%
  group_by(feature, antigen_key, plate, source) %>%
  summarise(inf_mfi_min = min(mfi_inf, na.rm = TRUE),
            inf_mfi_max = max(mfi_inf, na.rm = TRUE),
            .groups = "drop")

# ---------------------------------------------------------------------------
# 2. Fit one monotone transfer function per key. Returns a closure:
#    predict_inf(mfi_mat_vec) -> mfi_inf_vec (natural MFI), plus method label.
#    Fit is on log10-log10; prediction back-transformed with 10^.
# ---------------------------------------------------------------------------
build_transfer <- function(df) {
  # df: rows for ONE key, columns mfi_mat, mfi_inf (natural, >0)
  d <- df %>%
    transmute(x = log10(mfi_mat), y = log10(mfi_inf)) %>%
    filter(is.finite(x), is.finite(y)) %>%
    arrange(x)
  # collapse duplicate x (same maternal MFI at >1 dilution): average y on log scale.
  # (Rare; keeps x strictly increasing for the interpolators/spline.)
  d <- d %>% group_by(x) %>% summarise(y = mean(y), .groups = "drop") %>% arrange(x)
  npts <- nrow(d)

  if (npts < MIN_INTERP_PTS) {
    return(list(method = "skip", n = npts, fn = function(xin) rep(NA_real_, length(xin))))
  }

  if (npts >= MIN_SPLINE_PTS) {
    if (HAVE_SCAM) {
      fit <- tryCatch(
        scam::scam(y ~ s(x, k = min(10, npts - 1), bs = "mpi"), data = d),
        error = function(e) NULL)
      if (!is.null(fit)) {
        rng <- range(d$x)
        fn <- function(xin) {
          xc <- pmin(pmax(xin, rng[1]), rng[2])           # hold flat outside fit range
          as.numeric(predict(fit, newdata = data.frame(x = xc)))
        }
        return(list(method = "scam_mpi_spline", n = npts, fn = fn))
      }
    }
    # scam unavailable or failed -> monotone interpolating spline (base R)
    sf <- stats::splinefun(d$x, d$y, method = "monoH.FC")
    rng <- range(d$x)
    fn <- function(xin) { xc <- pmin(pmax(xin, rng[1]), rng[2]); sf(xc) }
    return(list(method = "monoH.FC_spline", n = npts, fn = fn))
  }

  # 3-5 points -> monotone LINEAR interpolation between the couples (no smoothing)
  rng <- range(d$x)
  fn <- function(xin) {
    xc <- pmin(pmax(xin, rng[1]), rng[2])
    stats::approx(d$x, d$y, xout = xc, method = "linear", rule = 2)$y
  }
  list(method = "monotone_linear_interp", n = npts, fn = fn)
}

# Build all transfer functions, keyed
key_cols <- c("feature", "antigen_key", "plate", "source")
pair_keys <- pairs %>% distinct(across(all_of(key_cols)))
transfer <- vector("list", nrow(pair_keys))
method_tbl <- pair_keys %>% mutate(transfer_method = NA_character_, n_pairs = NA_integer_)

for (i in seq_len(nrow(pair_keys))) {
  k <- pair_keys[i, ]
  di <- pairs %>% semi_join(k, by = key_cols)
  tf <- build_transfer(di)
  transfer[[i]] <- c(as.list(k), list(tf = tf))
  method_tbl$transfer_method[i] <- tf$method
  method_tbl$n_pairs[i]         <- tf$n
}
key_index <- pair_keys %>% mutate(.row = row_number())

# ---------------------------------------------------------------------------
# 3. Apply to the maternal/cord test samples.
#    Test MFI = 10^assay_response ; maternal-panel antigens are uppercase, so
#    antigen_key = tolower(antigen); recode source NIBSC -> NIBSC06_140.
# ---------------------------------------------------------------------------
samp <- mat_baa %>%
  mutate(
    mfi         = 10^assay_response,
    antigen_key = tolower(antigen),
    source_key  = ifelse(source == "NIBSC", "NIBSC06_140", source)
  )

# attach transfer + clamp range by key
samp <- samp %>%
  left_join(key_index, by = c("feature", "antigen_key", "plate", "source_key" = "source")) %>%
  left_join(inf_range, by = c("feature", "antigen_key", "plate", "source_key" = "source"))

apply_key <- function(row_idx, mfi_vec) {
  out <- rep(NA_real_, length(mfi_vec))
  for (ri in unique(row_idx[!is.na(row_idx)])) {
    sel <- which(row_idx == ri)
    tf  <- transfer[[ri]]$tf
    out[sel] <- 10^ tf$fn(log10(pmax(mfi_vec[sel], .Machine$double.eps)))
  }
  out
}

samp$mfi_inf_raw <- apply_key(samp$.row, samp$mfi)

# ---------------------------------------------------------------------------
# 4. Clamp to infant per-key MFI range (Q.C) and set flags.
#    oor_flag: "below" / "above" (clamped), "no_transfer" (key skipped/absent),
#    "ok" otherwise.
# ---------------------------------------------------------------------------
samp <- samp %>%
  mutate(
    transfer_method = method_tbl$transfer_method[match(.row, seq_len(nrow(method_tbl)))],
    below = is.finite(mfi) & is.finite(inf_mfi_min) & mfi < inf_mfi_min,
    above = is.finite(mfi) & is.finite(inf_mfi_max) & mfi > inf_mfi_max,
    mfi_inf_scale = dplyr::case_when(
      is.na(.row) | is.na(mfi_inf_raw) ~ NA_real_,           # no usable transfer
      below                            ~ inf_mfi_min,        # clamp low
      above                            ~ inf_mfi_max,        # clamp high
      TRUE                             ~ mfi_inf_raw
    ),
    oor_flag = dplyr::case_when(
      is.na(.row)                      ~ "no_transfer",
      is.na(mfi_inf_raw)               ~ "no_transfer",
      below                            ~ "below_clamped",
      above                            ~ "above_clamped",
      TRUE                             ~ "ok"
    ),
    assay_response_inf = ifelse(is.finite(mfi_inf_scale) & mfi_inf_scale > 0,
                                log10(mfi_inf_scale), NA_real_)
  )

mat_baa_norm <- samp %>% select(-.row, -below, -above)

# ---------------------------------------------------------------------------
# 5. Coverage report
# ---------------------------------------------------------------------------
rule <- function(t) cat("\n", strrep("=",66), "\n", t, "\n", strrep("=",66), "\n", sep="")
rule("TRANSFER-FUNCTION METHODS PER KEY")
print(as.data.frame(count(method_tbl, transfer_method)), row.names = FALSE)
cat("\n  keys with <3 pairs (skipped):",
    sum(method_tbl$transfer_method == "skip"), "of", nrow(method_tbl), "\n")

rule("TEST-SAMPLE OUTCOME (all features)")
print(as.data.frame(count(mat_baa_norm, feature, oor_flag)), row.names = FALSE)

rule("ADCD ONLY: outcome by sample_type")
print(as.data.frame(
  mat_baa_norm %>% filter(feature == "ADCD") %>% count(sample_type, oor_flag)),
  row.names = FALSE)

rule("DONE")
cat("  Output object: mat_baa_norm\n")
cat("  New columns: mfi (10^assay_response), mfi_inf_scale (infant-scale MFI, clamped),\n")
cat("               assay_response_inf (log10 of infant-scale MFI),\n")
cat("               transfer_method, oor_flag.\n")
