#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# setup_repo.sh — build the GaPs systems-serology repository skeleton and move
# the (flat) downloaded files into their correct locations.
#
# Usage:
#   1. Put every downloaded file (the 22 refactored files + README.md +
#      CHANGELOG.md + your unchanged originals + data/figures/docx) into a
#      single folder, say  ~/Downloads/gaps_drop/
#   2. Run:   bash setup_repo.sh  ~/Downloads/gaps_drop  ~/gaps-systems-serology
#
# Files not present in the drop folder are simply skipped (a notice is printed),
# so you can run it repeatedly as you gather the unchanged originals.
# ---------------------------------------------------------------------------
set -euo pipefail

SRC="${1:?usage: setup_repo.sh <source_dir> <target_repo_dir>}"
DST="${2:?usage: setup_repo.sh <source_dir> <target_repo_dir>}"

mkdir -p "$DST"/{config,R,parts,data,figures,docs,output}
: > "$DST/.here"   # here::here() sentinel

place() {  # place <dest_subdir> <file...>
  local sub="$1"; shift
  for f in "$@"; do
    if [ -f "$SRC/$f" ]; then
      cp "$SRC/$f" "$DST/$sub/$f"
    else
      echo "  (skip, not found: $f)"
    fi
  done
}

echo "config/";  place config endpoints.R endpoints_additions.R
echo "R/";       place R find_proj_file.R serology_helpers.R within_visit_helpers.R \
                        maternal_responder_helpers.R connection_transform.R
echo "parts/";   place parts 01_maternal_chain.Rmd 02_clinical_covariates.Rmd \
                        03_concurrent_networks.Rmd 04_response_bubbles.Rmd \
                        04_antigen_subclass_predictors.Rmd 05_effector_functions.Rmd
echo "data/";    place data c_set.RData clin_assess.RData \
                        igg_standard_residuals_ap_matpm_prevacvacc_k.RData \
                        igg_standard_residuals_wp_matpm_prevacvacc_k.RData
echo "figures/"; place figures Fig_coverage_matrix.png Fig_coverage_matrix.pdf \
                        Fig_samples_flow.png Fig_samples_flow.pdf
echo "docs/";    place docs GaPs_methods_and_results_section_baseline_prediction_chapter.docx \
                        statistical_methods_section-PartA_maternal_pre_vacc_to_infant_baseline.docx

echo "root (drivers + README/CHANGELOG)";
place . README.md CHANGELOG.md \
        SBA_analysis.Rmd WT_IgG_analysis.Rmd PTNA_analysis.Rmd \
        06_prediction_models.Rmd \
        07_within_visit_correlations_and_responder_subgroups.Rmd \
        08_maternal_responder_predictors.Rmd \
        matprevacc_to_matbirth_analysis.Rmd matbirth_to_cordblood_analysis.Rmd \
        cord_to_infmon2_analysis_refactored.Rmd interaction_and_deming.Rmd \
        responder_in_infant_model.Rmd standards_comparison.Rmd \
        results_B_to_F.Rmd response_bubbles_analysis.Rmd \
        response_bubbles_matpm_v4.Rmd transition_subgroup_heterogeneity.Rmd

echo
echo "Done. Repository skeleton at: $DST"
echo "Next:  cd \"$DST\"  &&  git init  &&  git add -A  &&  git commit -m 'Initial refactored codebase'"
