# GaPs systems-serology refactor — changelog

**Important:** none of this was knit-verified — there is no R in the environment
(no `Rscript`, network disabled), so every edit is static. Brace/paren balance and
call-site wiring were checked by hand; please knit once to confirm before relying on
output. Edits preserve the existing project layout (`config/`, `R/`, `parts/`, `data/`).

The data-casing question is resolved (see "Data note" at the bottom): the actual
features in `c_set.RData` are mixed-case (`PT_IgG`, `TT_IgG1`), the loader rebuilds
them to match the config, and TT/DT antibodies are present and correctly loaded —
the "flat TT" result is real data, not an artifact.

---

## Step 1 — Consolidation pass (no intended change to output)

- **`R/find_proj_file.R` (new):** single canonical definition of the self-locating
  path resolver. The ~15-line function previously copy-pasted into 13 driver `.Rmd`
  files is now sourced from here via a small guarded stub at the top of each driver
  (`if (!exists("find_proj_file", mode="function")) { ... source(...) }`). Drivers
  touched: `06`, `07`, `08`, `SBA_/WT_IgG_/PTNA_analysis`, `matprevacc_to_matbirth`,
  `matbirth_to_cordblood`, `cord_to_infmon2_analysis_refactored`,
  `interaction_and_deming`, `responder_in_infant_model`, `response_bubbles_analysis`,
  `results_B_to_F`.
- **Canonical antigen order (`config/endpoints.R`):** added
  `ANTIGEN_ORDER <- c("PT","FHA","PRN","FIM","TT","DT")` as the single source of
  truth for antigen layout. `BUBBLE_ANTIGENS` now points at it; `ANTI6` in
  `results_B_to_F.Rmd` and `06_prediction_models.Rmd` now point at it.
- **`ANTIGENS` collision reconciled (`config/endpoints_additions.R`):** the 6-antigen
  vector there was in a different order (`PT,PRN,FHA,DT,TT,FIM`) than `endpoints.R`.
  It now derives from `ANTIGEN_ORDER` (guarded), removing the order disagreement
  while keeping the intentional 6-antigen shadow used by Parts 4/5. The core-3
  `ANTIGENS` in `endpoints.R` (used by `add_concurrent_ratios`) is unchanged.
- **`pvalue_to_label` de-duplicated (`R/maternal_responder_helpers.R`):** the second
  copy is now `if (!exists("pvalue_to_label", mode="function"))`, so the shared
  definition in `R/serology_helpers.R` (sourced first) wins and the duplicate only
  acts as a standalone fallback.
- **Mandatory casing guard (`R/serology_helpers.R`):** `load_serology_data()` now
  calls `verify_predictor_casing()` on every load (guarded by `exists()`), so any
  future config/data casing drift errors loudly instead of silently dropping features.

## Step 2 — Change #5: responder defined by PT/FHA/PRN binding

- `responder_features` param changed to `c("PT_IgG","FHA_IgG","PRN_IgG")` in
  `07`, `08`, and `responder_in_infant_model.Rmd` (was PT/PRN total + PT/PRN IgG1).
- Hardcoded "four defining features" prose in `08` made dynamic via
  `` `r length(params$responder_features)` `` (4 sites), so it tracks the param.

## Step 3 — Change #1: Figure 3 expanded (`results_B_to_F.Rmd`)

- Antigen rows reordered to PT, FHA, PRN, FIM, TT, DT (via `ANTIGEN_ORDER`).
- Figure expanded from total IgG + IgG1 to **total IgG + IgG1–IgG4** (five series,
  offset within each antigen band; filled circle + open square/diamond/triangles).
- A **whole-cell IgG (`WT_IgG`) decay row** is pinned to the bottom (cord→InfMon2
  only; `WHOLE_WT_IgG` has no subclasses and is not measured pre-cord).
- Refactored the transition builders into `build_transition(tr, dat)` /
  `build_all_tr(dat)` / `build_fig_df(dat)` / `plot_pathway(fig_df, subtitle)` so the
  same code renders pooled and per-arm figures.
- **Two new per-arm figures** added: Figure 3 — Tdap-IPV mothers, Figure 3 — TT
  mothers (the pooled figure is retained). These resolve the TT/DT question
  directly: both arms contain tetanus toxoid (TT rises in both, subject to a
  pre-existing-immunity ceiling); diphtheria toxoid is Tdap-IPV-only (DT rises mainly
  there). FIM subclasses are unavailable pre-InfMon2, so the FIM row carries total
  IgG only in production/transfer — handled automatically by data availability.

## Step 4 — Change #2: bubble plots (`response_bubbles_matpm_v4.Rmd`, refactored)

The published bubble Figures 1–2 are produced by `response_bubbles_matpm_v4.Rmd`.
That file is now refactored into the project and is the single home for change #2.

- **`response_bubbles_matpm_v4.Rmd` (new/refactored):** the original v4 engines
  (`compute_arm_effect`, `plot_arm_effect`, `lancet_effect_table`) are preserved, but
  the file now uses the shared bootstrap (`find_proj_file` stub → `config/endpoints.R`,
  `R/serology_helpers.R`, `R/connection_transform.R`), loads via `load_serology_data()`
  (recoded visits + the mandatory casing guard), takes residual files through
  `find_proj_file`, and uses the recoded infant visits (`InfMon2/InfMon5/InfMon9`) and
  `BUBBLE_IGGLIST`. The old `./code/colors.R` dependency is dropped; flextable/officer
  are guarded so the doc still knits without them.
- **SBA removed, PTNA kept:** `functional_remap` now contains only
  `WHOLE_PTNA → PT_PTNA` (the `WHOLE_SBA → PRN_SBA` row is gone); `PRN_SBA` is dropped
  from `include_feature`/`include_feature_wo`; and `SBA` is removed from every analyte
  order and from the table labels. PTNA continues to display under PT.
- **60° rotation:** `plot_arm_effect()`'s theme rotates the x-axis (analyte) labels 60°.
- **Antigen order** is preserved from the published figure (`PT, FHA, PRN, FIM, DT, TT`);
  a comment notes it can be switched to the canonical `ANTIGEN_ORDER` for consistency.

NOTE: my earlier stand-in implementation of change #2 (a synthetic PTNA injection in
`parts/04_response_bubbles.Rmd`, plus `BUBBLE_WHOLECELL` and a `PTNA` entry in
`BUBBLE_ANALYTES`, and a rotation in `bubble_arm_effect()`) was made before the v4
source was available and has been **reverted** — `parts/04_response_bubbles.Rmd` is
back to its original, and those config/helper additions are removed — so there is a
single, authoritative implementation in the v4 file.

## Step 5 — Changes #3, #4, #6

- **#3 Demographic Table 1 (`R/serology_helpers.R` + `results_B_to_F.Rmd`):** new
  `demographics_by_arm()` helper (dependency-light; median [IQR] for continuous, n (%)
  for categorical) and a new "A. Study population" section rendering Table 1 **overall
  + Tdap-IPV + TT**. Arm membership is joined from the serology data (the clinical
  table has no arm column). Variables: maternal age/BMI/Hb, gestational ages,
  vaccine–birth interval, birth weight, parity, infant sex, delivery mode.
- **#4 Subgroup heterogeneity (`R/within_visit_helpers.R` + new
  `transition_subgroup_heterogeneity.Rmd`):** new `transition_subgroups()` and
  `transition_subgroup_scan()` cluster subjects on their per-analyte change (PAM +
  silhouette-chosen k; tertile fallback) for **each antigen × each of the three
  transitions**, reusing `paired_change()`. The new standalone document renders a
  summary table (method, silhouette, subgroup sizes/median Δ, subgroup×arm Fisher p)
  and a per-antigen per-transition distribution figure.
- **#6 IgG vs IgG1 vs IgG+IgG1 (`06_prediction_models.Rmd` + `interaction_and_deming.Rmd`):**
  06 gains **Family 4**, comparing total IgG, IgG1, and IgG+IgG1 as predictors at
  both MatBirth (prior) and InfMon2 (concurrent) on one shared complete-case frame.
  `interaction_and_deming.Rmd` now compares **IgG and IgG1** side-by-side in the
  interaction and Deming slope tables at both visits (`PRED_ANALYTES <- c("IgG","IgG1")`);
  plots stay on IgG1 as the primary visual. The joint IgG+IgG1 model is the domain of
  06 Family 4 (Deming is single-predictor by construction).

## Step 6 — Change #7: ADCD scaffold (inert)

- `parts/05_effector_functions.Rmd` gains an explicitly **inert** scaffold
  (`ADCD_READY <- FALSE`; all chunks `eval=FALSE`) mirroring the three-transition
  pathway for ADCD and the concurrent InfMon2 ADCD↔whole-cell relationship. It is
  design-only until the ADCD transition data are confirmed at the relevant visits.

---

## Latent issue flagged (not fixed here)

`parts/04_antigen_subclass_predictors.Rmd` builds `ANTIGEN_SUBCLASS_FEATS` and
`FCG_FEATS` from **uppercase** tokens (`IGG1`, `FCGR2A`) and matches them
case-sensitively against the post-load **mixed-case** columns (`PT_IgG1`,
`PT_FcgR2a`). On the current data this `intersect()` returns nothing, so Blocks 7–9
silently produce empty/skipped output. Fixing it means switching the additions
constants and the part-04 regexes to mixed case (or making the matches
case-insensitive as `06` already does). It was left untouched because it is outside
the seven requested changes and the recasing cascade can't be knit-verified here.

## Manuscript documents

The two Word summaries (`GaPs_methods_and_results_...docx`,
`statistical_methods_section-PartA...docx`) were **not** edited. They still say
"four defining features" (now three) and use the old Figure numbering; regenerate
them from the updated analysis.

## Data note (casing verdict)

`c_set.RData` carries a mixed-case `analyte` column (`IgG`, `IgG1`, `FcgR2a`) and a
vestigial uppercase pre-combined `feature` column (`PT_IGG1`). `load_serology_data()`
overwrites `feature` with `paste(antigen, analyte)` → mixed case (`PT_IgG1`), which
matches `config/endpoints.R`. TT/DT total IgG and subclasses are present at all
maternal/cord visits, so the analyses are internally consistent — but only because
the loader's casing happens to match the config, which is why the mandatory
`verify_predictor_casing()` guard was added.
