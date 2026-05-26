# HGNC symbol verification — notes for manuscript

## Summary

All 107 unique gene symbols across the six gene-set sheets (`LEPTIN_MTOR`, `LIPOPHAGY_CORE`, `MTOR_AUTOPHAGY_LIPOPHAGY`, `VIC_OSTEOGENIC`, `S100A9_RAGE_CALCIFICATION`, `ALL_ROWS`) were submitted to the HGNC Multi-symbol checker. HGNC IDs were then mapped back into the workbook via `VLOOKUP` against a new reference sheet `HGNC_REF`.

**Result: all 107 input symbols returned `Approved symbol` matches. No `Previous symbol` matches were returned, and no genuine `Alias symbol`-only substitutions were required. Zero alias→approved replacements were performed on the gene list.**

## Method

1. Submitted the 107 unique symbols to HGNC Multi-symbol checker.
2. Built `HGNC_REF` sheet in the workbook from the checker output, keeping only `Approved symbol` rows (107 rows total).
3. Added `HGNC_ID` column to each gene-set sheet using `=VLOOKUP(A<row>, HGNC_REF!$A$2:$C$108, 3, FALSE)`.
4. Recalculated all formulas — 216 VLOOKUPs evaluated, 0 errors.

## Spurious `Alias symbol` hits in the HGNC output (no action needed)

The HGNC checker returned five duplicate rows where a query symbol coincidentally matches an old alias of an unrelated gene. In every case, the query symbol also has a primary `Approved symbol` match, so the alias row is discarded as a cross-reference artifact, not a real synonym. For full disclosure for the manuscript reviewer:

| Query | Approved-symbol match (used) | Spurious alias-symbol hit (discarded) |
|---|---|---|
| `MTOR` | MTOR (HGNC:3942) — mechanistic target of rapamycin kinase | MTOR also flagged as an alias for HGNC:3942 (self-match duplicate) |
| `RHEB` | RHEB (HGNC:10011) — Ras homolog, mTORC1 binding | `RHEBP1` (HGNC:10010) — RHEB pseudogene 1 |
| `LIPA` | LIPA (HGNC:6617) — lipase A, lysosomal acid type | `SCGB1D1` (HGNC:18395) — secretoglobin family 1D member 1 |
| `SPP1` | SPP1 (HGNC:11255) — secreted phosphoprotein 1 | `CXXC1` (HGNC:24343) — CXXC finger protein 1 |
| `NOS2` | NOS2 (HGNC:7873) — nitric oxide synthase 2 | `NANOS2` (HGNC:23292) — nanos C2HC-type zinc finger 2 |

These five aliases are unrelated to the biology of the queried genes (e.g., the lysosomal acid lipase `LIPA` has nothing to do with the secretoglobin `SCGB1D1` despite sharing a historical alias string). The original symbols `MTOR`, `RHEB`, `LIPA`, `SPP1`, and `NOS2` are the current HGNC-approved symbols and were retained.

## Suggested manuscript wording

> Gene symbols were verified against the HGNC database (Multi-symbol checker, accessed 2026-05-25). All 107 symbols across the analyzed gene sets matched current HGNC-approved symbols and required no renaming; HGNC IDs were attached to each entry for traceability.

## Files

- `gene_sets_final_verified.xlsx` — updated workbook with `HGNC_ID` column on every gene sheet and a new `HGNC_REF` lookup sheet at the end.
- `notes.md` — this file.
