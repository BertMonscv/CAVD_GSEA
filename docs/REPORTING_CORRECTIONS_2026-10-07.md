# Reporting and provenance corrections

Date: 7 October 2026. Release: v1.0.0.

This dated note corrects reporting and provenance descriptions in the historical analysis archive. The original plans, implementation and archived results remain available under their original commits. The corrections below do not redefine the historical plans or imply that release-stage decisions were made before the analysis.

## Plan, implementation and run chronology

Plan v1 was committed as ed2d35d (tag prereg-frozen) on 23 May 2026. Plan v2 and its amendment were committed as c76f644 (tag prereg-frozen-v2) on 26 May 2026. The analysis implementation pre/run_gsea_v2.R was first committed as 6d6a7cd at 08:00:10 −04:00 on 26 May 2026, after the v2 plan and before the archived run at 18f0f7f (tag gsea-run-v2; 08:28:12 −04:00). The v2 plan tag itself does not contain the script. The numbered plan-section references in historical scripts and prior manuscript wording do not correspond to the unnumbered committed plan documents. Detailed sample exclusions, the log₂ transformation threshold, tie-breaking jitter and the three-category consistency rule are documented in the frozen implementation.

## Gene-set size labels

The historical plan prose transposed the sizes of CUSTOM_MTOR_AUTOPHAGY_LIPOPHAGY_AXIS and CUSTOM_LIPOPHAGY_CORE. The committed membership table, implementation, results and current manuscript consistently contain 20 genes in the axis set and 22 genes in the core set. The analyzed membership is unchanged.

## Historical NES intervals

The historical plotting output included bootstrap-derived NES intervals. The procedure resampled rank values and then assigned the original gene names, breaking the original gene–statistic pairing. Those intervals do not estimate uncertainty in the observed NES and are excluded from the current manuscript figure and supplementary table. Historical files and interval columns are retained for provenance and are not endorsed as inferential results. The independently computed NES estimates, raw P values, within-pool adjusted P values and leading-edge genes are unchanged.

## Multiple-testing scope

The five hypothesis-driven sets were tested within each of the three analysis sets, with Benjamini–Hochberg correction across the five tests in each set. The separate descriptive analysis of all 50 Hallmark sets was performed only in GSE51472, the discovery cohort. CUSTOM_S100A9_RAGE_CALCIFICATION was analyzed separately in each analysis set. The historical v1-to-v2 amendment describes v1 as pooling the descriptive Hallmark analysis with the targeted analysis; the committed v1 plan actually specifies Benjamini–Hochberg correction across the five targeted sets. The executed v2 analysis uses the separate pools reported above.

## Gene-symbol provenance

The archived HGNC notes document verification of 107 unique curated gene symbols without substitutions. Expression matrices use SYMBOL annotations from the recorded platform-specific Bioconductor packages. The implementation does not call a global HGNC symbol-correction function. The plan’s broad statement about conversion of all symbols should therefore be read alongside this implementation record.

## Current figure numbering and source files

Current manuscript Figure 1 contains the cohort design, 15 NES estimates without historical intervals, and the osteogenic leading-edge membership matrix. Historical Figure 2E refers to the earlier GSEA forest plot and is not the current Figure 2, which reports the T1 experiment. Historical Figure 1A corresponds to current Supplementary Figure S10a; 1B to S2; 1C–D to S10b–c; 1E–F to S10d–e; 1G to S10f; and 1H–I to S10g–h. Current assets are provided under figures/current/ with numerical sources under results/current/ and the build entry point scripts/build_current_release.py.

## Current table prose

The four rationale sentences for PRKAA1, SPP1, IBSP and BGLAP were edited for clarity in the current manuscript. All 108 membership records, HGNC identifiers, PubMed identifiers and short citations match the historical committed table. These wording changes do not alter set membership or analysis results.

## Archive scope

This repository and its corresponding Zenodo record document the transcriptomic and GSEA analysis. The microscopy measurements, Fiji macros and single-cell quantification software belong to a separate imaging project. A GSEA DOI does not establish deposition or licensing of that separate project.

## Dynamic KEGG annotations

The historical baseline KEGG panel (current Supplementary Figure S10h) uses an online annotation service. The 7 October 2026 rerun changed pathway labels and order. The archived figure is preserved as the manuscript artifact; the current live-service rerun is separately documented in the reproduction record. Core targeted and Hallmark GSEA results were reproduced exactly.
