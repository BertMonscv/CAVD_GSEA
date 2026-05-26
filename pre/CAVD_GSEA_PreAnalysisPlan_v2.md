# CAVD GSEA Pre-Analysis Plan v2

**Frozen on**: [TO BE FILLED on the day of commit, format YYYY-MM-DD]
**Supersedes**: v1.0 (frozen 2026-05-23, git commit `ed2d35d`, tag `prereg-frozen`). v2 adds four substantive specifications that affect the scientific validity of fgsea output (HGNC symbol normalization, probe collapse rule, `fgseaMultilevel()` with per-call seeding, separated BH pools) and two procedural specifications (commitment of the analysis script at freeze, gene set overlap check before freeze with an explicit Jaccard threshold). It also clarifies the ranking metric. The full list of changes from v1 and the rationale for each is documented in `AMENDMENT_2026-05-26_plan_v1_to_v2.md`.

## Framing

Hypothesis-driven, targeted transcriptomic corroboration of the leptin–mTORC1–autophagy/lipophagy–osteogenesis axis hypothesized from prior literature and validated in vitro in this study. Not an unbiased discovery exercise. All five gene set results will be reported in full regardless of significance.

## Datasets and preprocessing

- **Primary ranking**: GSE51472 (discovery cohort), `logFC` column from limma `topTable()`, sorted in decreasing order
- **Sensitivity ranking**: integrated GSE51472 + GSE12644 after ComBat batch correction
- **Cross-platform validation**: GSE83453, independent NES with 95% CI
- **Probe-to-gene collapse**: retain the probe with the highest interquartile range per gene; same rule applied to all three datasets
- **Gene symbol normalization**: all symbols converted to HGNC current official symbols before gene set matching (e.g., LC3B → MAP1LC3B, ADRP → PLIN2, S6K1 → RPS6KB1). Alias resolution log committed at freeze.

## Pre-specified gene sets (5)

1. Hallmark MTORC1_SIGNALING (MSigDB v7.5 via msigdbr 7.5.1)
2. CUSTOM_LEPTIN_MTOR_SIGNALING (N=22)
3. CUSTOM_MTOR_AUTOPHAGY_LIPOPHAGY_AXIS (N=22)
4. CUSTOM_LIPOPHAGY_CORE (N=20)
5. CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION (N=27)

Final gene set definitions, with a PMID source for each gene in custom sets, are committed as `Supp_Table_S1_v2.tsv` before any fgsea call. Pairwise Jaccard overlap among the four custom gene sets is computed by `verify_gene_set_overlap.R` and logged in `gene_set_overlap_log.tsv` before freeze; any pair with Jaccard > 0.30 is addressed by reassigning shared genes to their most specific set prior to the freeze commit. No modification of gene set composition after the commit.

## Parameters (frozen)

- **Function**: `fgsea::fgseaMultilevel()` (not classic `fgsea()`)
- **Parameters**: `minSize = 8`, `maxSize = 500`, `eps = 0`, `sampleSize = 101`
- **Seed**: `set.seed(42)` called immediately before each `fgseaMultilevel()` invocation
- **Multiple-testing correction**: Benjamini–Hochberg applied **separately** within
  - (a) the 5 pre-specified gene sets — hypothesis-driven pool
  - (b) the 50 Hallmark gene sets — descriptive pool
  These two pools are not merged under any circumstance.
- **Significance threshold**: adjusted P < 0.05
- **Analysis script**: committed as `run_gsea_v2.R` at freeze; the analysis is run from this committed version without subsequent edits, except via dated amendment
- **Software versions**: captured via `sessionInfo()` and committed at run time

## Reporting

- All 5 pre-specified gene sets reported in Supplementary Table S2 (NES, NES 95% CI, raw P, adjusted P, leading-edge genes), regardless of significance
- Cross-cohort NES reported for each pre-specified set across discovery / integrated / GSE83453 (forest plot, Figure 2E)
- Leading-edge genes mapped to their corresponding wet-lab assay in Supplementary Table S5

## Decision rules

- A gene set with fewer than 8 mapped members in a given cohort is reported with a note that statistical power is limited; the raw result is still reported, not dropped.
- Non-significant pre-specified results are reported and interpreted in the Discussion (in particular, biological interpretation of any non-significant lipophagy / leptin–mTOR enrichment), not dropped or reframed as "exploratory".
- No post-hoc modification of gene set composition, ranking metric, fgsea parameters, or BH correction pool. Any change after the freeze commit requires a dated amendment committed before re-running, and is acknowledged in the manuscript.

---

*End of plan.*
