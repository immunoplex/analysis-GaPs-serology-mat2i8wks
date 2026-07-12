# =====================================================================
# R/find_proj_file.R
# ---------------------------------------------------------------------
# SINGLE source of truth for the self-locating project-file resolver.
# Previously this ~15-line function was copy-pasted into every driver
# .Rmd; it now lives here and each driver sources it via a small stub
# (see the "shared path resolver" block at the top of each file).
#
# Finds config/, R/, parts/, data/ whether they sit beside the calling
# .Rmd, under an analysis/ subfolder, or at the project root — and works
# even when rmarkdown renders from a temp directory. Returns the first
# candidate that exists; errors (or returns NULL when optional=TRUE) if
# nothing matches.
# =====================================================================

find_proj_file <- function(rel, optional = FALSE) {
  cand <- character(0)
  if (requireNamespace("here", quietly = TRUE))
    cand <- c(cand, here::here(rel), here::here("analysis", rel))
  rmd_dir <- tryCatch(dirname(knitr::current_input(dir = TRUE)),
                      error = function(e) NA_character_)
  if (!is.na(rmd_dir)) cand <- c(cand, file.path(rmd_dir, rel),
                                 file.path(dirname(rmd_dir), rel))
  cand <- unique(c(cand, file.path(getwd(), rel), rel))
  hit  <- cand[file.exists(cand)]
  if (length(hit)) return(normalizePath(hit[1], winslash = "/", mustWork = FALSE))
  if (optional) return(NULL)
  stop("Could not locate '", rel, "'.\nTried:\n  ", paste(cand, collapse = "\n  "),
       "\nPlace config/, R/, parts/ (and data/) beside this file, under analysis/, ",
       "or at the project root.", call. = FALSE)
}
