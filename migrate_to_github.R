#!/usr/bin/env Rscript

## Migrate the contents of a local directory into an existing GitHub repo.
##
## The remote already contains an AGPL LICENSE (and possibly a README).
## This script PRESERVES that existing content and layers your local files
## on top of it, rather than overwriting the remote history.
##
## Run from R:     source("migrate_to_github.R")
## or from shell:  Rscript migrate_to_github.R
##
## NOTE ON AUTH: reading the remote (public repo) needs no credentials, but
## the final `git push` does. Make sure git is already set up to authenticate
## to GitHub on this machine (a credential helper + Personal Access Token with
## "Contents: write", or an SSH key). See notes at the bottom of this file.

## ----------------------------- CONFIG -----------------------------
repo_url   <- "https://github.com/immunoplex/analysis-GaPs-serology-mat2i8wks.git"
local_dir  <- "/gaps_system_serology"
commit_msg <- "Migrate gaps_system_serology contents"
## ------------------------------------------------------------------

## Run a git command inside local_dir; echo output; stop on non-zero exit.
git <- function(...) {
  args <- c("-C", local_dir, ...)
  out  <- suppressWarnings(system2("git", args, stdout = TRUE, stderr = TRUE))
  st   <- attr(out, "status")
  if (length(out)) cat(out, sep = "\n")
  if (!is.null(st) && st != 0L)
    stop(sprintf("git %s  (exit %s)", paste(..., collapse = " "), st), call. = FALSE)
  invisible(out)
}

## Run a git command and return TRUE/FALSE for success (never stops).
git_ok <- function(...) {
  args <- c("-C", local_dir, ...)
  st <- suppressWarnings(system2("git", args, stdout = FALSE, stderr = FALSE))
  identical(st, 0L)
}

## 0. Sanity checks -------------------------------------------------
if (nchar(Sys.which("git")) == 0L)
  stop("git is not installed or not on PATH.", call. = FALSE)
if (!dir.exists(local_dir))
  stop(sprintf("Local directory does not exist: %s", local_dir), call. = FALSE)

## A commit needs an author identity.
name  <- suppressWarnings(system2("git", c("config", "--get", "user.name"),
                                  stdout = TRUE, stderr = FALSE))
email <- suppressWarnings(system2("git", c("config", "--get", "user.email"),
                                  stdout = TRUE, stderr = FALSE))
if (length(name) == 0L || length(email) == 0L)
  stop(paste0(
    "git identity is not configured. Set it once, then re-run:\n",
    '  git config --global user.name  "Your Name"\n',
    '  git config --global user.email "you@example.com"'), call. = FALSE)

## 1. Detect the remote's default branch ---------------------------
symref <- suppressWarnings(
  system2("git", c("ls-remote", "--symref", repo_url, "HEAD"),
          stdout = TRUE, stderr = TRUE))
branch <- sub("^ref: refs/heads/(\\S+)\\s+HEAD.*$", "\\1",
              grep("^ref:", symref, value = TRUE))
if (length(branch) != 1L || branch == "") branch <- "main"
message(sprintf("Remote default branch: %s", branch))

## 2. Write .gitignore (keeps the data/ directory out of the repo) -
gitignore <- file.path(local_dir, ".gitignore")
writeLines(c(
  "# Keep the local data/ directory out of version control.",
  "# Anchored to the repository root (/gaps_system_serology).",
  "/data/"
), gitignore)
message("Wrote ", gitignore)

## 3. Initialise the repo in place (idempotent) --------------------
if (!dir.exists(file.path(local_dir, ".git"))) git("init")

## 4. Stage & commit local content (data/ is ignored) -------------
git("add", "-A")
if (!git_ok("diff", "--cached", "--quiet")) {
  git("commit", "-m", commit_msg)
} else {
  message("Nothing new to commit.")
}

## 5. Align the local branch name with the remote's default --------
if (git_ok("rev-parse", "--verify", "HEAD")) git("branch", "-M", branch)

## 6. Point 'origin' at the remote (add or update) -----------------
if (git_ok("remote", "get-url", "origin")) {
  git("remote", "set-url", "origin", repo_url)
} else {
  git("remote", "add", "origin", repo_url)
}

## 7. Merge the existing remote history so the AGPL LICENSE is kept -
has_remote <- length(suppressWarnings(
  system2("git", c("ls-remote", "--heads", repo_url, branch),
          stdout = TRUE, stderr = FALSE))) > 0L
if (has_remote) {
  git("fetch", "origin", branch)
  ## -X ours => on a file collision (e.g. both have README.md) keep YOUR
  ## local version; remote-only files such as LICENSE are still merged in.
  git("pull", "--allow-unrelated-histories", "--no-edit", "-X", "ours",
      "origin", branch)
} else {
  message("Remote branch has no commits yet; skipping merge.")
}

## 8. Push ----------------------------------------------------------
git("push", "-u", "origin", branch)
message("\nDone. Pushed to ", repo_url, " (", branch, ").")

## ------------------------------------------------------------------
## If `git push` fails with an auth error:
##   - HTTPS: configure a credential helper and use a Personal Access Token
##     (fine-grained token with Contents: read/write on this repo).
##   - or switch to SSH:
##       git -C /gaps_system_serology remote set-url origin \
##         git@github.com:immunoplex/analysis-GaPs-serology-mat2i8wks.git
##     then re-run this script (or just `git push -u origin <branch>`).
##
## If the merge stops on a conflict you'd rather resolve by hand, edit the
## marked files, then:  git add <files>;  git commit;  git push -u origin <branch>
## ------------------------------------------------------------------
