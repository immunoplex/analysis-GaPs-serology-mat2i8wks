# GaPs systems-serology — repository layout

This is the canonical directory structure for the refactored codebase. The
self-locating resolver `R/find_proj_file.R` finds `config/`, `R/`, `parts/`, and
`data/` whether they sit beside the calling `.Rmd`, under an `analysis/` subfolder,
or at the project root — so the **flat repo-root layout below is the simplest and
recommended one**. Knit/run from the repository root.

```
gaps-systems-serology/
│
├── .here                         # empty sentinel so here::here() anchors to the root
├── README.md                     # this file
├── CHANGELOG.md                  # what changed in the refactor (★)
│
├── config/                       # sourced first by every driver
│   ├── endpoints.R               # (●) ANTIGEN_ORDER, vocab, ENDPOINTS, BUBBLE_*, verify_predictor_casing
│   └── endpoints_additions.R     # (●) Parts 4/5 feature constants (ANTIGENS ← ANTIGEN_ORDER)
│
├── R/                            # sourced helpers (NOT knitted)
│   ├── find_proj_file.R          # (★ NEW) single canonical path resolver
│   ├── serology_helpers.R        # (●) loader+casing guard, demographics_by_arm, bubble/paired/chain helpers
│   ├── within_visit_helpers.R    # (●) within-visit engine + transition_subgroups()/_scan() (★ added)
│   ├── maternal_responder_helpers.R  # (●) responder classifier (pvalue_to_label de-duplicated)
│   └── connection_transform.R    # (○ unchanged) IgG-standardisation for the bubble residuals
│
├── parts/                        # CHILD documents — included via child=; do NOT knit directly
│   ├── 01_maternal_chain.Rmd            # (○)
│   ├── 02_clinical_covariates.Rmd       # (○)
│   ├── 03_concurrent_networks.Rmd       # (○)
│   ├── 04_response_bubbles.Rmd          # (○ reverted to original)
│   ├── 04_antigen_subclass_predictors.Rmd  # (○)  ← latent casing bug flagged in CHANGELOG
│   └── 05_effector_functions.Rmd        # (●) + inert ADCD scaffold (eval=FALSE)
│
├── data/                         # inputs (○ unchanged)
│   ├── c_set.RData
│   ├── clin_assess.RData
│   ├── igg_standard_residuals_ap_matpm_prevacvacc_k.RData
│   └── igg_standard_residuals_wp_matpm_prevacvacc_k.RData
│
├── figures/                      # static study-design assets (○)
│   ├── Fig_coverage_matrix.png   /  Fig_coverage_matrix.pdf
│   └── Fig_samples_flow.png      /  Fig_samples_flow.pdf
│
├── docs/                         # manuscript text — NOT edited; regenerate from analysis
│   ├── GaPs_methods_and_results_section_baseline_prediction_chapter.docx
│   └── statistical_methods_section-PartA_maternal_pre_vacc_to_infant_baseline.docx
│
├── output/                       # created at knit time (e.g. v4's M02 effect-table .docx)
│
└── (repo root) ── TOP-LEVEL DRIVERS — these are the files you knit ───────────────
    │
    │   Endpoint-parameterised pipeline (each includes parts/01–05):
    ├── SBA_analysis.Rmd                  # (●)
    ├── WT_IgG_analysis.Rmd               # (●)
    ├── PTNA_analysis.Rmd                 # (●)
    │
    │   Standalone analyses:
    ├── 06_prediction_models.Rmd          # (●) + Family 4 (IgG vs IgG1 vs IgG+IgG1)
    ├── 07_within_visit_correlations_and_responder_subgroups.Rmd   # (●) responder = PT/FHA/PRN
    ├── 08_maternal_responder_predictors.Rmd                       # (●) responder = PT/FHA/PRN
    ├── matprevacc_to_matbirth_analysis.Rmd      # (●)  production transition
    ├── matbirth_to_cordblood_analysis.Rmd       # (●)  transfer transition
    ├── cord_to_infmon2_analysis_refactored.Rmd  # (●)  decay transition
    ├── interaction_and_deming.Rmd        # (●) + IgG vs IgG1 slope tables
    ├── responder_in_infant_model.Rmd     # (●) responder = PT/FHA/PRN
    ├── standards_comparison.Rmd          # (○ unchanged)
    ├── results_B_to_F.Rmd                # (●) Table 1 demographics + expanded Figure 3 (+per-arm, +WT_IgG)
    ├── response_bubbles_analysis.Rmd     # (●) wraps parts/04_response_bubbles
    ├── response_bubbles_matpm_v4.Rmd     # (★ NEW/refactored) manuscript bubble Figs 1–2 (SBA removed, PTNA kept, 60°)
    └── transition_subgroup_heterogeneity.Rmd     # (★ NEW) per-antigen × per-transition subgroups
```

Legend: ★ new · ● changed in the refactor · ○ unchanged (use your original copy).

---

## Where each downloaded file goes

The files were delivered flat (no folders). Place them as follows:

| Downloaded file | → Repository path |
|---|---|
| `endpoints.R`, `endpoints_additions.R` | `config/` |
| `find_proj_file.R`, `serology_helpers.R`, `within_visit_helpers.R`, `maternal_responder_helpers.R`, `connection_transform.R` | `R/` |
| `01_maternal_chain.Rmd`, `02_clinical_covariates.Rmd`, `03_concurrent_networks.Rmd`, `04_response_bubbles.Rmd`, `04_antigen_subclass_predictors.Rmd`, `05_effector_functions.Rmd` | `parts/` |
| `*.RData` (4 files) | `data/` |
| `Fig_*.png` / `Fig_*.pdf` | `figures/` |
| `*.docx` (2 files) | `docs/` |
| every other `*.Rmd` (SBA/WT_IgG/PTNA, 06/07/08, transitions, interaction, responder, standards, results_B_to_F, response_bubbles_analysis, response_bubbles_matpm_v4, transition_subgroup_heterogeneity) | repo root |
| `CHANGELOG.md`, `README.md` | repo root |

**The one thing to watch:** `01`–`05` live in `parts/`, but `06`/`07`/`08` live at the
**root** — they are standalone drivers, not children.

Unchanged files (○) can come from your existing originals; the changed (●) and new
(★) ones are the 22 deliverables from this refactor.

---

## How `find_proj_file()` resolves paths

For a request like `find_proj_file("config/endpoints.R")` it tries, in order:
`here::here("config/endpoints.R")`, `here::here("analysis/config/endpoints.R")`, the
calling `.Rmd`'s own directory and its parent, then `getwd()`, then the bare relative
path — returning the first that exists. So as long as `config/`, `R/`, `parts/`,
`data/` sit at the repo root (or under `analysis/`), it resolves correctly from any of
those working directories. Add an empty `.here` file at the root so `here::here()`
anchors reliably (or open the folder as an RStudio Project).

## Render order / entry points

You don't knit `parts/` or `R/` directly. The documents you render are the root
drivers. A sensible order:

1. `results_B_to_F.Rmd` — study population (Table 1) + Figure 3 pathway.
2. `response_bubbles_matpm_v4.Rmd` — manuscript bubble Figures 1–2.
3. `SBA_analysis.Rmd`, `WT_IgG_analysis.Rmd`, `PTNA_analysis.Rmd` — full per-endpoint reports (pull in parts/01–05).
4. `06_prediction_models.Rmd`, `07_…`, `08_…` — prediction / responder analyses.
5. `matprevacc_to_matbirth`, `matbirth_to_cordblood`, `cord_to_infmon2`, `interaction_and_deming`, `responder_in_infant_model`, `standards_comparison`, `response_bubbles_analysis`, `transition_subgroup_heterogeneity` — supporting analyses.

Suggested R packages: tidyverse, broom, relaimpo, knitr, here, cluster (for #4),
glmnet/ppcor (Parts 4 importance / partial cor), data.table + flextable + officer
(v4 standardised grid + Lancet tables), deming (optional cross-check), RColorBrewer.

## Notes carried over from the refactor (see CHANGELOG.md)

- Nothing is knit-verified (no R in the build environment) — knit once before relying on output.
- `parts/04_antigen_subclass_predictors.Rmd` has a pre-existing uppercase-casing bug
  (Blocks 7–9 silently empty) flagged but not fixed.
- The two `docs/*.docx` still say "four defining features" (now three) and use old
  figure numbering — regenerate them from the updated analysis.
