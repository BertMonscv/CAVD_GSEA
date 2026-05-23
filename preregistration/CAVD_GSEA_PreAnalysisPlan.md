# CAVD GSEA Pre-Analysis Plan (frozen 2026-05-23)

## Datasets and rankings
- Primary ranking: GSE51472 log2FC (limma moderated t)
- Sensitivity ranking: integrated GSE51472+GSE12644 (ComBat)
- Cross-platform validation: GSE83453 (independent NES with 95% CI)

## Pre-specified gene sets (5)
1. Hallmark MTORC1_SIGNALING (MSigDB v7.5 via msigdbr 7.5.1)
2. CUSTOM_LEPTIN_MTOR_SIGNALING (N=22, see Supp Table S1)
3. CUSTOM_MTOR_AUTOPHAGY_LIPOPHAGY_AXIS (N=22, see Supp Table S1)
4. CUSTOM_LIPOPHAGY_CORE (N=20, see Supp Table S1)
5. CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION (N=27, see Supp Table S1)

## Parameters (frozen)
- fgsea: nperm=10000, minSize=8, maxSize=500
- Seed: set.seed(42)
- Multiple-testing correction: Benjamini-Hochberg across all 5 gene sets
- Significance threshold: adjusted P < 0.05
- All 5 results will be reported regardless of significance

## Decision rules
- If a gene set has <8 members in ranked list: report as "insufficient coverage"
- If NES has same sign in all 3 cohorts: report as "directionally consistent"
- No post-hoc gene set modification, no parameter re-tuning
