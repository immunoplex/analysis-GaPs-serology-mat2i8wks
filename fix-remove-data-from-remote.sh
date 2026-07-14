#!/usr/bin/env bash
#
# fix-remove-data-from-remote.sh
#
# Goal: keep the CLEAN remote history (data removed), reapply only your
# legitimate non-data work, and ensure /data holds nothing but README.md.
#
# Run this from inside your repo:
#   cd ~/Documents/R-git/gaps_system_serology
#   bash fix-remove-data-from-remote.sh
#
# It PAUSES before destructive actions. Read each prompt.

set -euo pipefail

# ----------------------------------------------------------------------
# Files from your old commit to KEEP (edit this list as needed).
# Deliberately EXCLUDES everything under data/ except README.md,
# and excludes docs/site_libs/ and .Rproj.user/.
# ----------------------------------------------------------------------
KEEP_FILES=(
  "analysis/C0_data_harmonisation.Rmd"
  "analysis/parts/C0_standardise.Rmd"
  "analysis/parts/C0_standards_comparison.Rmd"
  "R/C0_connection_transform.R"
  "README.md"
  "config/endpoints_additions.R"
  ".gitignore"
  "data/README.md"
  docs/C0_data_harmonisation.html
    docs/site_libs/crosstalk-1.2.2/css/crosstalk.min.css
    docs/site_libs/crosstalk-1.2.2/js/crosstalk.js
    docs/site_libs/crosstalk-1.2.2/js/crosstalk.js.map
    docs/site_libs/crosstalk-1.2.2/js/crosstalk.min.js
    docs/site_libs/crosstalk-1.2.2/js/crosstalk.min.js.map
    docs/site_libs/crosstalk-1.2.2/scss/crosstalk.scss
    docs/site_libs/datatables-binding-0.34.0/datatables.js
    docs/site_libs/datatables-css-0.0.0/datatables-crosstalk.css
    docs/site_libs/dt-core-1.13.6/css/jquery.dataTables.extra.css
    docs/site_libs/dt-core-1.13.6/css/jquery.dataTables.min.css
    docs/site_libs/dt-core-1.13.6/js/jquery.dataTables.min.js
    docs/site_libs/htmltools-fill-0.5.9/fill.css
    docs/site_libs/htmlwidgets-1.6.4/htmlwidgets.js
    docs/site_libs/nouislider-7.0.10/jquery.nouislider.min.css
    docs/site_libs/nouislider-7.0.10/jquery.nouislider.min.js
    docs/site_libs/selectize-0.12.0/selectize.bootstrap3.css
    docs/site_libs/selectize-0.12.0/selectize.min.js
)

BACKUP_BRANCH="backup-before-reset"
COMMIT_MSG="revised C0 (analysis + config, no data)"

pause() {
  echo
  read -r -p ">>> $1  (press Enter to continue, Ctrl-C to abort) "
  echo
}

# ----------------------------------------------------------------------
echo "=== 0. Sanity: confirm we are in a git repo on branch main ==="
git rev-parse --is-inside-work-tree >/dev/null
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
echo "Current branch: $CURRENT_BRANCH"
if [ "$CURRENT_BRANCH" != "main" ]; then
  echo "WARNING: you are not on 'main'. Switch to main before continuing."
  exit 1
fi

# ----------------------------------------------------------------------
echo "=== 1. Back up current local state to '$BACKUP_BRANCH' ==="
if git show-ref --verify --quiet "refs/heads/$BACKUP_BRANCH"; then
  echo "Backup branch '$BACKUP_BRANCH' already exists — leaving it as-is."
else
  git branch "$BACKUP_BRANCH"
  echo "Created backup branch '$BACKUP_BRANCH' pointing at current HEAD."
fi
git log --oneline -3 "$BACKUP_BRANCH"

# ----------------------------------------------------------------------
echo "=== 2. Fetch the clean remote history ==="
git fetch origin
echo
echo "Remote main (the history you WANT to keep):"
git log --oneline -10 origin/main
echo
echo "Local main (current):"
git log --oneline -10 main

pause "Step 3 will HARD RESET local main to origin/main. Your work is safe on '$BACKUP_BRANCH'."

# ----------------------------------------------------------------------
echo "=== 3. Reset local main to match clean remote ==="
git reset --hard origin/main
echo "Local main now matches origin/main."

# ----------------------------------------------------------------------
echo "=== 4. Reapply ONLY the kept (non-data) files from the backup ==="
for f in "${KEEP_FILES[@]}"; do
  if git cat-file -e "$BACKUP_BRANCH:$f" 2>/dev/null; then
    git checkout "$BACKUP_BRANCH" -- "$f"
    echo "  restored: $f"
  else
    echo "  SKIP (not in backup): $f"
  fi
done

# ----------------------------------------------------------------------
echo "=== 5. Stage changes and VERIFY no data files (except README) ==="
git add -A

echo
echo "Staged files under data/ (should be ONLY data/README.md, or nothing):"
DATA_STAGED=$(git diff --cached --name-only | grep '^data/' || true)
echo "${DATA_STAGED:-<none>}"

# Abort if any data file other than README.md is staged
BAD=$(echo "${DATA_STAGED}" | grep -v '^data/README.md$' | grep -v '^$' || true)
if [ -n "$BAD" ]; then
  echo
  echo "ERROR: Unexpected data files are staged:"
  echo "$BAD"
  echo "Unstage them before committing. Aborting."
  exit 1
fi

echo
echo "Full staged file list:"
git diff --cached --name-only

pause "Step 6 will COMMIT the above. Review the list."

# ----------------------------------------------------------------------
echo "=== 6. Commit ==="
if git diff --cached --quiet; then
  echo "Nothing to commit — working tree already matches. Skipping commit."
else
  git commit -m "$COMMIT_MSG"
fi

# ----------------------------------------------------------------------
echo "=== 7. Push (fast-forward, NO --force needed) ==="
pause "Step 7 will push to origin/main."
git push

echo
echo "=== DONE ==="
echo "If everything looks correct on GitHub, remove the backup branch with:"
echo "    git branch -D $BACKUP_BRANCH"
