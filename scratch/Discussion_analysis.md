# Discussion section — editorial notes

---

## Part 0 — Reference map of the Discussion as currently written

Paragraph labels (D1–D19) the current order.

| Label | Opening words | Results Section | One-line content |
|---|---|---|---|
| **D1** | "By combining randomised…" | all | High-level summary of every finding |
| **D2** | "This analysis builds directly on GaPs…" | A | Positions against the parent trial; refines "quality" → quantity |
| **D3** | "The comparison between maternal arms…" | A | TT has no pertussis antigens; residual compositional differences (PT ADCD, FHA FcγRIIa, subclass) |
| **D4** | "The predominantly quantitative effect…" | C | Responder phenotype; pre-vaccination serology (CV-AUC 0.88/0.92) |
| **D5** | "Maternal-to-cord changes…" | B | Transfer & decay similar between arms |
| **D6** | "Transfer was nevertheless selective…" | B | IgG1 enrichment; postnatal decay/half-life *(contains placeholder "To add 1–2 sentences")* |
| **D7** | "Clinical factors added little…" | G, I | Clinical block CV-R²; GA/interval/birthweight → PTNA |
| **D8** | "The transferred antibody pool separated…" | D, F | Two antigen-linked pathways; PT-IgG1 vs PRN-IgG1 mechanism |
| **D9** | "Maternal antibody measurements at delivery…" | E, H, I | Maternal-delivery prediction; shared pool; cord adds little; arm modifies prospective not concurrent |
| **D10** | "These findings have practical implications…" | A, C, D, H | Translational: assay panel, maternal-delivery indicator |
| **D11** | "This study has several strengths…" | — | Strengths: trial design, timing |
| **D12** | "Absolute and total-IgG-standardised…" | A, D | Strengths: quantity vs quality; collinearity; FHA not inhibition |
| **D13** | "Generalisability was assessed…" | G, H, I | Strengths: CV, Shapley/LMG, prior-vs-concurrent |
| **D14** | "This study had several limitations…" | — | Limitations: descriptive not causal; sample-size variation; power |
| **D15** | "The analysis addresses the pre-primary…" | B | Limitations: 8-week only; no cord functional assays *(placeholder note on other compartments)* |
| **D16** | "*B. pertussis*-specific total IgG…" | D, F | Limitations: composite readouts; PRN not sole mediator |
| **D17** | "The responder phenotype…" | C | Limitations: responder score not yet a biomarker |
| **D18** | "Finally, GaPs was conducted…" | — | Limitations: single centre; generalisability; not prespecified |
| **D19** | "In conclusion…" | all | Conclusion: parsimonious model |

**Biggest thing to change:** The mechanistic effector-pathway result (Results F, second half — FcγR3b / neutrophil phagocytosis vs complement, and arm-dependent pathway) has not been spelled out in the discussion. 

---

## Edits

### A1. The Discussion paragraph 8 attributes surface killing to complement; the Results demote complement in favour of neutrophil phagocytosis. 
- **Where:** D8 — "PRN IgG1 may contribute more directly to bacterial killing through surface binding and **complement recruitment**."
- **Text issue:** Results F is explicit that complement deposition (ADCD) is not the operative effector. "Complement deposition added no independent variance in any representation"; in the joint effector model "complement deposition did not survive for either readout," and neutrophil phagocytosis (ADNP), coupled to FcγRIIIb, "carries the effector signal most closely tied to both surface functions". The section summary states the operative limb is "Fcγ-receptor 3b engagement and neutrophil phagocytosis — not complement deposition".
- Beyond being inconsistent with Results F, the D8 clause makes a positive mechanistic claim ("contribute… through… complement recruitment") that the study's own mediation analysis specifically failed to support (the serial IgG1→FcγRIIIb→ADCD→outcome path was non-significant).
- **Current:** "PRN IgG1 may contribute more directly to bacterial killing through surface binding and complement recruitment, although association does not establish mediation."
- **Replace with:** "PRN IgG1 may contribute more directly to bacterial killing through surface binding and engagement of neutrophil-mediated effector function (see below), whereas complement deposition did not behave as an independent effector in these data; association nonetheless does not establish mediation."

### A2. Half-life range is cited two different ways, and the study's own estimate sits below both.
- **Where:** Intro (l.47) cites Oguti as "approximately **29–35 days**." D6 (l.318) cites the same meta-analysis as "half-lives of 29–36 day[s]." Results B (l.87) reports the study's *own* estimate as **27.7–29.9 days**, "within the reported 95% CI for four of six antigens."
- **Text issue:** Three different numbers for effectively the same comparison, and the discussion never states the study's own point estimates — it cites only the external range.
- Intro and D6 cite Oguti identically (pick one range and use it in both), and have D6 report the study's own 27.7–29.9-day estimate with the "4 of 6 within CI" qualifier rather than paraphrasing the meta-analysis alone.
- **Current tail:** "…consistent with a previous meta-analysis that reported half-lives of 29–36 day[]s. To add 1-2 sentences"
- **Replace with:** "…with estimated half-lives clustering at 27.7–29.9 days across antigens, within the reported 95% confidence interval for four of six antigens in the individual-participant meta-analysis of Oguti and colleagues (PMID: 34949496). Because half-life varied little by antigen, the cord titre at birth—set by the maternal response and placental transfer—was the principal determinant of how long antibody persisted before the infant primary series, so that maternal response magnitude, rather than differential decay, governed the antibody available at 8 weeks."

### A3. D7's clinical-block CV-R² range does not reconcile with Section I. Also: birth-clinical variables "contributed to" PTNA; Results G and I say clinical variables added nothing.
- **Where:** D7 — "the combined **design and clinical** block explained CV-R² values of **0.14–0.30**."  D7 — "gestational age at birth, vaccine-to-delivery interval and birthweight **contributed to** PT-neutralising activity."
- **Text issue:** Section I (l.292) reports the **design block alone** at CV-R² **0.20 (binding), 0.29 (PTNA), 0.34 (SBA)**, and states adding clinical covariates "did not improve prediction and **slightly reduced** cross-validated fit." A design-plus-clinical range therefore should sit at or just below 0.20–0.34, not 0.14–0.30. The 0.14 lower bound is unexplained by the main text.
- **Text issue:** Results G (l.258): clinical variables "produced no cross-validated gain… Maternal age showed a **borderline** association with PTNA, but the remaining variables contributed little." Results I (l.292): the clinical block "did not improve prediction and slightly reduced fit." Results C (l.122): the vaccination-to-delivery interval "explained little variation."
- **Replace with:** Clinical and pregnancy covariates added little in either role. They did not classify maternal responder status independently of vaccine arm (appendix S.7, Table S7.3), and as a predictive block for infant function they produced no cross-validated gain: the design block, principally maternal vaccine arm, reached CV-R² 0.20–0.34 across outcomes, and adding pre-vaccination and birth clinical covariates did not improve and modestly reduced cross-validated fit (Table S7.4; the combined design-plus-clinical foundation used in the variance analyses was 0.14–0.30). Only maternal age showed a borderline association with PT-neutralising activity; gestational age, vaccination-to-delivery interval and birthweight did not independently predict infant function, and any influence on transfer duration or antibody retention should be regarded as hypothesis-generating (PMIDs: 26797213, 33013843).
- **Replace the clinical sentence in Results Section G:** Clinical variables also added little. Assessed both individually and as pre-vaccination and birth covariate blocks, they produced no cross-validated gain beyond the antigen-specific IgG1 models and modestly reduced fit when added to the design block (Table 6; appendix S.7). Maternal age showed a borderline association with PTNA, but no clinical covariate — including gestational age, vaccination-to-delivery interval or birthweight — independently improved prediction once the matched antibody signal was included.
- **Replace part of the clinical sentence in Results Section I:** Current: "Adding clinical covariates did not improve prediction and slightly reduced cross-validated fit."
- **Edit to:** "Adding pre-vaccination and birth clinical covariates did not improve prediction and slightly reduced cross-validated fit (appendix S.7).

### A4. D5 opening has alternative text.
- **Text issue:** D5 (l.316) — "Maternal-to-cord changes/transfer patterns(?) and early postnatal decay…"
- **Action (edit text):** Resolve the "(?)" and the slashed alternatives to a single clean subject.

### A5. The effector-pathway dissection (FcγRIIIb / neutrophil phagocytosis, not complement; arm-dependent wiring) is essentially absent.
- **What is missing:** Results F establishes (i) surface function is predominantly antibody-quantity-driven; (ii) once quantity is removed by IgG-standardisation, **FcγRIIIb (Fcγ-receptor 3b)** is the dominant per-antibody correlate in Tdap-IPV infants; (iii) among effectors, **neutrophil phagocytosis (ADNP)** is the independent correlate of both SBA and whole-cell binding, while **complement is redundant**; (iv) the FcγRIIIb→readout mapping is **arm-modified** ("shared architecture driven at a different gain"; network M p=0.40, S p=0.02, l.222).
- **Why it belongs:** The Introduction explicitly frames candidate surface mechanisms as "Fcγ-receptor engagement, complement activation, phagocytosis and bacterial killing" (l.45). Results F *answers* which of these operates. Leaving it out means the Discussion poses the intro's mechanistic question and never returns the study's own answer — and, worse, currently gives the wrong answer (A1).

- **What is missing:** Results E shows PRN-driven surface slopes are roughly twice as steep after Tdap-IPV (SBA slope 0.86 vs 0.42, p=0.0014) while the PT→PTNA slope is arm-invariant (p=0.89). D9 notes vaccination "modified some prospective… associations" but does not explain *how*.
- **Why it belongs:** Results F supplies the mechanism (the arm difference is in the FcγRIIIb limb, "same architecture, different gain," not in antibody amount or the IgG1 slope). E + F together are the mechanistic core of the surface story.
- Try out this mechanistic effector paragraph (insert after D8):
- The surface pathway resolved further when antibody quantity was removed. Both SBA and whole-cell binding were predominantly quantity-driven—untransformed pertactin IgG1 alone explained the majority of their variance—but after IgG-standardisation the dominant per-antibody correlate in Tdap-IPV infants was Fcγ-receptor IIIb engagement, and among the three antibody-dependent effector functions it was neutrophil phagocytosis, not complement deposition, that remained independently associated with both readouts. Complement behaved as a redundant correlate of Fcγ-receptor IIIb rather than an independent effector. The maternal-arm difference in surface function was therefore one of pathway gain rather than antibody amount: the antibody–function architecture was shared across arms, but the Fcγ-receptor IIIb→function coupling was engaged after Tdap-IPV and near-silent after tetanus toxoid. This provides the mechanism for the steeper pertactin-associated surface slopes seen after Tdap-IPV and locates the vaccine-sensitive step in a neutrophil-restricted effector limb rather than in the quantity or subclass of transferred antibody.

### A6. Section I / cross-validated pathway paragraphs tighten but add integration
- **What is missing:** Three independent results point the same way. **PT/neutralisation** is arm-independent (E), fully predictable at delivery with no concurrent gain (H/I), and has ~no unique infant variance (I). **PRN/surface** is arm-dependent (E), partly infant-intrinsic (~31% unique at 8 weeks; I), and mechanistically re-wired by arm at the FcγRIIIb limb (F).
- If we state this as one framework, it could convert a list of separate findings into a single interpretive claim and gives the paper a memorable through-line.
- In addition - Divergence arises at maternal production (B); transfer preserves composition and is arm-neutral (B); decay is uniform across antigens and arms so cord titre governs persistence (B); therefore cord and infant samples largely re-measure one shared pool (I), and cord adds no predictive information (I).
- Tighten the discussion section for Results Section I to:
- The layered cross-validated decomposition was not a re-description of the preceding models but the step that placed the entire maternal–infant axis on a single out-of-sample scale, allowing the contribution of each stage to be compared directly and separated from in-sample overfitting. On this common footing, most of the transferable predictive signal was already resolved at maternal delivery: the maternal-production block raised cross-validated R² to 0.59 for SBA, 0.42 for PT-neutralising activity and 0.23 for *B. pertussis*-specific total IgG binding, whereas adding cord serology changed cross-validated fit by no more than 0.024. Cord blood was therefore informative about the mechanism of transfer but largely redundant for prediction once maternal delivery serology was known. Order-independent Shapley/LMG attribution then apportioned the maternal-block variance to a single dominant antigen per outcome—pertactin-specific IgG for the two surface outcomes (approximately 87% for SBA and 75% for binding) and PT-specific IgG for neutralisation (approximately 79%)—reproducing the antigen assignments of the separate multivariable models by a method with a different failure mode, and confirming that the negative conditional FHA coefficients reflected collinearity rather than inhibition.
-
- Partitioning the explained variance into components shared and unique across timepoints recast maternal, cord and infant serology as successive windows on one predominantly shared transferred pool (common component 43–71%; unique maternal share 1–4%) and, in doing so, exposed a dichotomy the per-outcome models could not resolve: PT neutralisation was essentially fully predictable from maternal antibody at delivery, with no gain from concurrent 8-week co-measurement, whereas roughly a third of the variance in *B. pertussis*-specific binding was unique to the infant sample. Read alongside the arm-modification confined to the pertactin surface chain (within-arm SBA R² 0.71 after Tdap-IPV versus 0.20 after tetanus toxoid), this frames PT-directed toxin neutralisation as a largely pure-transfer phenotype—arm-independent and fully predictable from maternal antibody at delivery—and pertactin-driven surface function as transfer-plus-context, arm-dependent and retaining an irreducible infant-intrinsic component. The flat upper rungs of the ladder are thus themselves the result: they quantify how nearly the maternal delivery sample serves as a sufficient statistic for the transferable component of early infant immunity.


### A7. Residual compositional ("quality") differences are listed but not interpreted.
- **What is missing:** Results A shows specific residual per-antibody differences that survive IgG-standardisation — PT-specific ADCD (1.7-fold), FHA-specific FcγRIIa (1.4-fold), and a reproducible subclass pattern (lower PT/FHA/PRN IgG3, higher FHA IgG1; l.63). D3 recites these but stops at "cannot establish whether subclass balance independently caused the residual functional differences."
- **Why it belongs:** The Introduction raises Fc glycosylation and Fc-receptor-binding composition (l.46). The IgG3-lowering / IgG1-favouring pattern is the measurable compositional correlate of that framing and dovetails with the FcγRIIIb finding (C1).
- **Action (reconfigure D3 + link):** Keep the appropriate caution, but connect the surviving FHA-FcγRIIa and the IgG3/IgG1 shift to the FcγRIIIb-neutrophil axis from C1 rather than leaving them as an unexplained list.
) and in D9.

### A8. Placental selection vs shared-pool dominance are reconcilable, not contradictory.
- **Opportunity:** Transfer is demonstrably *selective* by subclass (B/D6, IgG1 enrichment) yet produces *no independent cord-specific predictive signal* (I). The synthesis — selection acts on composition but does not generate information beyond maternal delivery serology — reconciles the two and pre-empts a reviewer seeing D6 and D9 as in tension. One sentence bridging D6 and D9.

### A9. Consider moving the one-sentence "cord adds little / shared pool" claim out of D1's summary so the summary previews rather than pre-empts the Section I argument.
