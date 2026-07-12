# build_interval_days.R
#
# Computes the exact number of days between the cord-blood sample (C00)
# and the infant 8-week sample (B20) for each subject in the GaPs database.
#
# Date column problem & fix
# --------------------------
# 'Date of receipt' (col R) contains a mix of:
#   (a) Excel date serials stored correctly as date cells -> readxl "guess"
#       sees the earliest string cells (row ~107) and downgrades the whole
#       column to text, returning serials as character strings "43509" etc.
#   (b) DD/MM/YYYY strings entered as plain text (correct dates, wrong format)
#   (c) DD/MM/YYYY strings with a stray dot: "07/.07/2021" -> "07/07/2021"
#   (d) Unrecoverable typos: "22/072020", "20/10/020" -> NA with a message
#
# Because readxl will always classify this column as text (due to the early
# string cells), we accept that and use col_types = "text" throughout.
# A single parse_receipt_date() function handles all four cases explicitly,
# converting serial strings via the Excel epoch (1899-12-30).
#
# Output: INTERVAL_DAYS_BY_SUBJECT  (subject_accession, interval_days)
#         Ready for the half-life chunk in cord_to_infmon2_analysis_refactored.Rmd,
#         which detects this object and uses per-subject intervals automatically.
# ---------------------------------------------------------------------------

library(readxl)
library(dplyr)

XLSX_PATH  <- here::here("./raw_files/AT_GaPs_database.xlsx")
SHEET_NAME <- "Final GAP MIA list"

# ---------------------------------------------------------------------------
# Date parser — handles all observed formats in col R.
# Input:  character vector (as returned by readxl col_types = "text")
# Output: Date vector (NA where unrecoverable)
# ---------------------------------------------------------------------------

EXCEL_EPOCH <- as.Date("1899-12-30")   # Excel's date serial origin

parse_receipt_date <- function(x) {

  out <- as.Date(rep(NA, length(x)))

  # (a) Excel serial strings: 4–6 digit integers ("43509", "43591", ...)
  is_serial <- grepl("^\\d{4,6}$", x, perl = TRUE)
  if (any(is_serial, na.rm = TRUE)) {
    out[is_serial] <- EXCEL_EPOCH + as.integer(x[is_serial])
  }

  # (b+c) DD/MM/YYYY strings, with optional stray dot before the month
  #       "14/08/2019", "07/.07/2021" -> strip the dot then parse
  is_dmy <- grepl("^\\d{2}/\\.?\\d{2}/\\d{4}$", x, perl = TRUE)
  if (any(is_dmy, na.rm = TRUE)) {
    cleaned        <- gsub("/\\.", "/", x[is_dmy], perl = TRUE)
    out[is_dmy]    <- as.Date(cleaned, format = "%d/%m/%Y")
  }

  # (d) Anything still NA that was non-missing: unrecoverable typo — warn once
  unrecoverable <- !is.na(x) & is.na(out) & !is_serial & !is_dmy
  if (any(unrecoverable, na.rm = TRUE)) {
    bad <- unique(x[unrecoverable])
    warning(sprintf(
      "parse_receipt_date: %d value(s) could not be parsed and were set to NA: %s",
      sum(unrecoverable),
      paste(shQuote(bad), collapse = ", ")
    ), call. = FALSE)
  }

  out
}

# ---------------------------------------------------------------------------
# Read — col_types = "text" throughout (readxl would guess text anyway
# because of the early plain-text cells, so this is explicit not lossy).
# Select columns by name so the script is robust to column count changes.
# ---------------------------------------------------------------------------

raw <- read_excel(XLSX_PATH, sheet = SHEET_NAME, col_types = "text") |>
  select(
    subject_accession = `sample ID`,
    timepoint         = `time point`,
    date_raw          = `Date of receipt`
  ) |>
  filter(timepoint %in% c("C00", "B20")) |>
  mutate(date = parse_receipt_date(date_raw))

# ---------------------------------------------------------------------------
# Split, join, compute interval.
# ---------------------------------------------------------------------------

cord <- raw |>
  filter(timepoint == "C00") |>
  select(subject_accession, cord_date = date)

inf8w <- raw |>
  filter(timepoint == "B20") |>
  select(subject_accession, b20_date = date)

INTERVAL_DAYS_BY_SUBJECT <- cord |>
  inner_join(inf8w, by = "subject_accession") |>
  mutate(interval_days = as.numeric(b20_date - cord_date)) |>
  select(subject_accession, interval_days)

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

n_cord   <- nrow(cord)
n_b20    <- nrow(inf8w)
n_paired <- nrow(INTERVAL_DAYS_BY_SUBJECT)
n_na     <- sum(is.na(INTERVAL_DAYS_BY_SUBJECT$interval_days))
n_neg    <- sum(INTERVAL_DAYS_BY_SUBJECT$interval_days <= 0, na.rm = TRUE)

message("=== build_interval_days.R diagnostics ===")
message(sprintf("C00 rows:               %d", n_cord))
message(sprintf("B20 rows:               %d", n_b20))
message(sprintf("Paired subjects:        %d", n_paired))
message(sprintf("C00 only (no B20):      %d", n_cord - n_paired))
message(sprintf("B20 only (no C00):      %d", n_b20  - n_paired))
message(sprintf("NA interval_days:       %d  (date missing or unparseable)", n_na))
message(sprintf("Non-positive intervals: %d  (data check)", n_neg))
message("")
message("Interval summary (days):")
print(summary(INTERVAL_DAYS_BY_SUBJECT$interval_days))

flag <- INTERVAL_DAYS_BY_SUBJECT |>
  filter(is.na(interval_days) | interval_days <= 0 | interval_days > 120)
if (nrow(flag) > 0) {
  message("\nRows flagged for review (NA, <= 0, or > 120 days):")
  print(as.data.frame(flag))
} else {
  message("\nNo rows flagged for review.")
}

med_days <- median(INTERVAL_DAYS_BY_SUBJECT$interval_days, na.rm = TRUE)
message(sprintf(
  "\nINTERVAL_DAYS_BY_SUBJECT ready: %d subjects, median %.1f days (%.1f weeks).",
  n_paired, med_days, med_days / 7
))
