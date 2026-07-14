# GaPs Systems Serology — Maternal-to-Infant Antibody Pipeline

> **Pertussis immunisation in pregnancy shapes the functional antibody
> landscape before primary infant pertussis vaccination**
> Anja Saso, Michael Scot Zens et al.

Analysis code for the systems-serology sub-study of the
[Gambian Pertussis Study (GaPs)](https://doi.org/10.1016/s1473-3099(25)00072-6),
a randomised 2×2 maternal–infant pertussis vaccine trial (n = 312
mother–infant dyads). The pipeline traces antigen-specific IgG quantity,
subclass composition and Fc-effector function from maternal vaccination
through placental transfer to the pre-primary infant baseline at 8 weeks.

📊 **[View the analysis site](https://immunoplex.github.io/analysis-GaPs-serology-mat2i8wks/)**

---

## Table of contents

- [Analysis components](#analysis-components)
- [Quick start](#quick-start)
- [Execution order](#execution-order)
- [Data requirements](#data-requirements)
  - [Summary of required data objects](#summary-of-required-data-objects)
- [Key design decisions](#key-design-decisions)
- [R environment](#r-environment)
- [Repository layout](#repository-layout)

---

## Analysis components

Each component maps to a manuscript section and one or more analysis files:

| Code | Rmd file(s) | Manuscript | What it does |
|------|-------------|------------|--------------|
| **C0** | `C0_data_harmonisation.Rmd` | Methods | Data harmonisation & IgG standardisation. Cross-batch MFI normalisation (standard-curve back-interpolation) and total-IgG standardisation of functional/subclass assays. Produces every `igg_standard_residuals_*` / `igg_standard_model_fit_*` input for C1, C4, C5, C7. Run once. |
| **C1** | `C1_arm_contrast.Rmd` `C1_arm_contrast_forest.Rmd` | Results A | Maternal arm contrast at the 8-week baseline — bubble plots (Fig 1B/C) and forest plot (Fig 1A). Runs in both untransformed and IgG-standardised representations. |
| **C2** | `C2_antibody_chain.Rmd` | Results B | Within-subject antibody changes across three transitions: maternal production (PregEarly→MatBirth), placental transfer (MatBirth→CordBlood), and postnatal decay (CordBlood→8 weeks). Fig 2, Table 1, half-life comparison with Oguti et al. |
| **C3** | `C3_responder_phenotype.Rmd` | Results C | Predictors of maternal High/Low responder status. Elastic-net logistic regression (CV-AUC), pre-vaccination model, interval-robustness check. Fig 3, Table 2. |
| **C4** | `C4_forward_prediction_PTNA.Rmd` `C4_forward_prediction_SBA.Rmd` `C4_forward_prediction_WTIgG.Rmd` `C4_forward_prediction.Rmd` | Results D–E | Forward prediction of infant function from maternal antibody at delivery. Six-block framework (quantity and class-switching at each chain stage). Arm×antibody interactions and Deming EIV regression. Tables 3, 4a/4b. |
| **C5** | `C5_concurrent_models.Rmd` | Results F | Concurrent (same-visit) infant models at 8 weeks. Hierarchical model building: arm → clinical → total IgG → IgG1 → effectors. Forward-vs-concurrent contrast. Table 5, Table 6a/6b. |
| **C6** | `C6_responder_in_model.Rmd` | Results J | Tests whether the maternal responder phenotype adds predictive information beyond vaccine arm and matched IgG1. Incremental CV-R², Shapley partition. Table S5. |
| **C7** | `C7_effector_pathway.Rmd` `C7_prn_pathway_figure4.Rmd` | Results H | Concurrent effector-pathway dissection: IgG1/IgG3 → FcγRIIa/FcγRIIIb → ADCD → SBA/WT-IgG. Graphical-lasso networks, serial mediation, layered LMG decomposition, NCT arm comparison. Fig 4, Tables S6–S7. |
| **C8** | `C8_subclass_balance.Rmd` | Results I | Compositional IgG subclass balance (ILR) and ADCD across the three chain stages. Mediation walk: balance→ADCD→SBA at maternal, cord and 8-week nodes. Fig 5, Table S8. |
| **CK** | `CK_block_ladder.Rmd` | Results K | Cross-validated block-ladder decomposition. Six blocks entered in temporal order; Shapley/LMG antigen attribution within the maternal-production block; commonality analysis. Fig 6. |

---

## Quick start

```r
# 1. Open the project — always via the .Rproj file to anchor here::here()
#    Double-click gaps_system_serology.Rproj in File Explorer

# 2. Build all pages
library(workflowr)
wflow_build()

# 3. To force a full rebuild of all pages
wflow_build(republish = TRUE)
```

---

## Execution order

```text
C0 (data prep, run once)
  ↓
C1 → C2 → C3 → C4 (PTNA / SBA / WTIgG / interaction, can run in parallel)
                ↓
               C5 → C6 → C7 → C8 → CK
```

All components are independent once data files are in place.
`wflow_build()` will only rebuild files that have changed since the last build.

---

## Data requirements

Data are **not** included in this repository. See [`data/README.md`](data/README.md)
for the full file list and per-component dependency table.

### Summary of required data objects

| File | Contents | Used by |
|------|----------|---------|
| `d_set.RData` | **Canonical dataset.** Full bead-array long-format object (`ebaa_extra2`). All timepoints, all antigens, all analytes. | C1–C8, CK |
| `clin_assess.RData` | Clinical and demographic variables per subject | C3, C6 |
| `igg_standard_residuals_ap_matpm_prevacvacc_k.RData` | aP-arm IgG-standardised residuals (infant visit) | C1, C5, C7 |
| `igg_standard_residuals_wp_matpm_prevacvacc_k.RData` | wP-arm IgG-standardised residuals (infant visit) | C1, C5, C7 |
| `igg_standard_residuals_*_mat.RData` | Maternal-visit IgG-standardised residuals (aP and wP arms) | C4 |
| `igg_standard_residuals_*_trn.RData` | Cord-blood IgG-standardised residuals (aP and wP arms) | C4 |
| `igg_standard_model_fit_*_mat.RData` | IgG standardisation model fits (for applying to new data) | C0 |
| `INTERVAL_DAYS_BY_SUBJECT.RData` | Exact cord-blood to 8-week interval per subject (days) | C2 |
| `S1_data.csv` | Supplementary Table S1 summary: paired within-subject changes across transitions | C7 |
| `figure3_abc_data.csv` | Pre-computed data for Figure 3 panels A/B/C | C3 |

---

## Key design decisions

| Decision | Rationale |
|----------|-----------|
| **workflowr** for pipeline management | Each HTML page records the git commit, R session, working directory and seed at render time — full reproducibility audit trail |
| **`here::here()`** for all file paths | Path-agnostic; works from any working directory as long as the `.Rproj` file is used to open the session |
| **`d_set.RData` / `ebaa_extra2`** as canonical data | Includes additional ADCD timepoints (P02, P09, M00, P18) not present in the earlier `c_set.RData` |
| **C-prefix naming** for all analysis files | Maps one-to-one to the eight prespecified analytical components in the manuscript (Supplementary Table S.10) |
| **`calc.relimp()` gated** to `compute_relimp = TRUE` | Without gating, Shapley/LMG decomposition runs on ~88 models per file (2^12 permutations each), causing >1 hour build time. Gated to the three full 12-subclass models that feed the comparison tables |
| **`dplyr::` explicit namespacing** in analysis chunks | Prevents silent masking of `select()`/`filter()` by `fs`, `MASS`, or other packages loaded in the knit environment |
| **`analysis/*_files/` gitignored** | knitr writes figure files beside the Rmd as a build cache; workflowr copies them to `docs/` for GitHub Pages. Only `docs/` figures are committed |

---

## R environment

- **R version:** 4.5.1 (2025-06-13, "Great Square Root")
- **Key packages:** `workflowr`, `tidyverse`, `here`, `broom`, `relaimpo`,
  `knitr`, `RColorBrewer`, `data.table`, `ppcor`, `glmnet`, `igraph`,
  `qgraph`, `NetworkComparisonTest`, `mediation`

---

## Repository layout

```text
gaps_system_serology/
│
├── _workflowr.yml                    # knit_root_dir: "."  seed: 1
├── gaps_system_serology.Rproj        # always open R via this file
├── README.md                         # this file
├── CHANGELOG.md                      # refactor history
├── LICENSE
│
├── analysis/                         # PRIMARY RMD FILES — knit via wflow_build()
│   ├── index.Rmd                     # Pipeline overview (GitHub Pages home)
│   │
│   ├── C0_standards_comparison.Rmd   # C0 · cross-batch MFI harmonisation
│   ├── C0_standardise.Rmd            # C0 · looped IgG-standardisation engine
│   │
│   ├── C1_arm_contrast.Rmd           # Results A · arm contrast bubble plots (Fig 1B/C)
│   ├── C1_arm_contrast_forest.Rmd    # Results A · forest plot (Fig 1A)
│   │
│   ├── C2_antibody_chain.Rmd         # Results B · production/transfer/decay (Fig 2, Table 1)
│   │
│   ├── C3_responder_phenotype.Rmd    # Results C · maternal responder predictors (Fig 3, Table 2)
│   │
│   ├── C4_forward_prediction.Rmd     # Results E · arm interaction & Deming EIV (Tables 4a/4b)
│   ├── C4_forward_prediction_PTNA.Rmd# Results D · maternal→infant PTNA (six-block framework)
│   ├── C4_forward_prediction_SBA.Rmd # Results D · maternal→infant SBA (six-block framework)
│   ├── C4_forward_prediction_WTIgG.Rmd# Results D · maternal→infant WT IgG binding
│   │
│   ├── C5_concurrent_models.Rmd      # Results F · concurrent infant models (Table 5)
│   │
│   ├── C6_responder_in_model.Rmd     # Results J · responder phenotype as predictor (Table S5)
│   │
│   ├── C7_effector_pathway.Rmd       # Results H · antibody chain transition summary (Fig 2)
│   ├── C7_prn_pathway_figure4.Rmd    # Results H · PRN effector pathway (Fig 4, Tables S6–S7)
│   │
│   ├── C8_subclass_balance.Rmd       # Results I · subclass balance & complement (Fig 5, Table S8)
│   │
│   ├── CK_block_ladder.Rmd           # Results K · CV block-ladder decomposition (Fig 6)
│   │
│   └── parts/                        # Child Rmd fragments — sourced via child=, NOT knitted directly
│       ├── C1_response_bubbles.Rmd          # Bubble plot engine for C1
│       ├── C2_chain_main.Rmd                # Chain transition engine for C2
│       ├── C3_clinical_covariates.Rmd       # Clinical covariate models for C3
│       ├── C4_antigen_subclass_predictors.Rmd # Subclass predictor models for C4
│       ├── C7_step00_networks.Rmd           # Effector pathway step 0 for C7
│       ├── C7_step01_pathway_by_arm.Rmd     # Effector pathway step 1 for C7
│       ├── C7_step02_link_regressions.Rmd   # Effector pathway step 2 for C7
│       ├── C7_step03_mediation.Rmd          # Effector pathway step 3 for C7
│       ├── C7_step04_variance_decomp.Rmd    # Effector pathway step 4 for C7
│       ├── C7_step05_fcgr_collinearity.Rmd  # Effector pathway step 5 for C7
│       ├── C7_step06_nct.Rmd                # Effector pathway step 6 for C7
│       ├── C7_step07_fcgr3b_unique.Rmd      # Effector pathway step 7 for C7
│       ├── C7_step08_agg_complement.Rmd     # Effector pathway step 8 for C7
│       ├── C7_step09_arm_interactions.Rmd   # Effector pathway step 9 for C7
│       ├── C7_step10_fcgr3b_vs_adcd.Rmd      # Effector pathway step 10 for C7
│       ├── C7_step11_effector_activation.Rmd # Effector pathway step 11 for C7
│       ├── C7_step12_effector_to_sba.Rmd    # Effector pathway step 12 for C7
│       ├── C7_step13_effector_aggregate.Rmd # Effector pathway step 13 for C7
│       └── C7_step14_effector_functions.Rmd # Effector pathway step 14 for C7
│
├── R/                                # Helper functions sourced at knit time (NOT knitted directly)
│   ├── shared_utils.R                # find_proj_file(), %||% operator
│   ├── C0_connection_transform.R     # IgG standardisation / residual transform
│   ├── C1_serology_helpers.R         # Arm contrast, bubble plots, paired-change helpers
│   ├── C2_interval_helpers.R         # Half-life and decay-rate helpers (Oguti method)
│   ├── C3_responder_helpers.R        # Responder classifier, elastic-net wrapper
│   ├── C4_within_visit_helpers.R     # Within-visit and forward-prediction modelling engine
│   ├── C7_effector_helpers.R         # Graphical-lasso network, bridge-centrality helpers
│   ├── C8_subclass_balance_helpers.R # ILR subclass balance and mediation helpers
│   └── colors.R                      # Colour palettes shared across components
│
├── code/                             # Build-once data preparation scripts
│   │                                 # (NOT sourced during knit — run manually to rebuild data)
│   ├── C0_build_interval_days.R      # Compute cord-to-8wk interval per subject
│   ├── C0_get_8week_interval_days.R  # Extract and format 8-week interval days
│   ├── C0_normalise_mat_to_inf_mfi.R # Cross-batch MFI normalisation (mat/cord → infant scale)
│   └── C0_db_load_functions.R        # Database load helpers (MFI data extraction)
│
├── config/
│   ├── endpoints.R                   # ANTIGEN_ORDER, ENDPOINTS, feature vocabulary
│   └── endpoints_additions.R         # Parts 4/5 feature constants (ANTIGENS ← ANTIGEN_ORDER)
│
├── data/                             # NOT committed — see data/README.md
│   └── README.md                     # Full data schema with per-component cross-reference
│
├── figures/                          # Static publication-quality figures
│   ├── C0_coverage_matrix.png/.pdf   # Figure S: assay coverage matrix
│   ├── C0_samples_flow.png/.pdf      # Figure S: sample flow diagram
│   ├── C3_interval_robustness.png/.pdf # Figure S1: interval-robustness plot
│   ├── C7_prn_pathway_figure4.png/.pdf # Figure 4: PRN effector pathway
│   └── CK_block_ladder.png/.pdf      # Figure 6: CV block-ladder decomposition
│
├── docs/                             # GitHub Pages output — auto-generated by wflow_build()
│   └── (do not edit manually)
│
└── scratch/                          # Archived pre-refactor files (146 files, not used)
```
