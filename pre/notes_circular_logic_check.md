# Circular logic check — Day 4

**Date:** 2026-05-26
**Input:** `top50_DEGs.txt` (50 genes, prepared Day 0) vs 5 curated gene sets in `gene_sets_final_verified.xlsx`

## Results

| gene_set | set size | top50 overlap | overlap genes | threshold | status |
|---|---:|---:|---|---|---|
| CUSTOM_LEPTIN_MTOR_SIGNALING | 22 | 0 | — | ≤2 | ✓ |
| CUSTOM_MTOR_AUTOPHAGY_LIPOPHAGY_AXIS | 20 | 0 | — | ≤2 | ✓ |
| CUSTOM_LIPOPHAGY_CORE | 22 | 0 | — | ≤2 | ✓ |
| CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION | 27 | 2 | IBSP, SPP1 | 5–8 expected | ⚠ below expected (see note) |
| CUSTOM_S100A9_RAGE_CALCIFICATION | 17 | 2 | MMP9, S100A8 | ≤2 | ✓ |

**Verdict: no circular-logic risk.** All non-osteogenic sets have ≤2 overlap. No replacement needed.

## Note on VIC_OSTEOGENIC lower-than-expected overlap

Expected 5–8 overlap based on assumption that classic osteogenic markers (RUNX2, BGLAP,
ALPL, COL1A1) would appear in top-50 DEGs. They did not — but this reflects the **DEG
signature of this particular dataset**, not a problem with set composition.

**Set composition is intact:** all canonical osteogenic markers are present in the set —
RUNX2, SP7, BGLAP, ALPL, COL1A1/2, IBSP, SPP1, BMP2/4, DLX5, MSX2, PHOSPHO1, ENPP1, ANKH
(17 canonical markers; aliases like OSX=SP7, OCN=BGLAP, BSP=IBSP, OPN=SPP1, TNAP=ALPL
are covered by the official symbols already in the set). The set was not over-diluted.

**Why the top-50 doesn't show classic osteogenic markers:** the top-50 DEG list is
dominated by:
- Immunoglobulin chains (IGH/IGK/IGL family, ~21/50 entries)
- MMPs (MMP1, MMP9, MMP12, MMP13)
- Chemokines (CXCL5, CXCL13)
- B-cell/plasma-cell markers (POU2AF1, SLAMF7, TNFRSF17, CD69, CD52)
- Mast-cell marker (CPA3)

This is an **immune-infiltration + ECM-remodeling** signature, not a classical
osteoblast-differentiation signature. The two osteogenic genes that did make top-50
(IBSP, SPP1) are precisely the genes with dual ECM/mineralization roles — consistent
with the data being driven more by matrix remodeling than by transcriptional osteogenic
reprogramming.

**Implication for GSEA interpretation:** if CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION still
shows significant enrichment despite only 2 leading-edge candidates from top-50, that
enrichment is being driven by the **broader middle-rank signal** (genes ranked 51–N
that are still upregulated but didn't crack the top 50). This is actually a stronger
result biologically — it means the osteogenic program is broadly upregulated, not just
a couple of standout genes.

## Action items

- [x] Document overlap analysis (this file)
- [ ] In Methods Section 2.4, add one sentence: "Pre-specified gene sets were checked
      for overlap with the top-50 DEGs from the discovery dataset; non-osteogenic sets
      showed ≤2 overlapping genes (Supplementary Table SX), confirming the analysis is
      not circular."
- [ ] In Discussion/Limitations, optionally note that the discovery DEG signature is
      dominated by immune infiltration, and that osteogenic program detection therefore
      relies on aggregated mid-rank signal rather than top-of-list genes.
