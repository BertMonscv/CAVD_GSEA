# Recovered-source release checks

This record covers the 13 recovered p-S6, VIC identity and co-localization source files listed in `provenance_manifest.json`. Other modules have their own provenance and checks.

- Three recovered Python plotting scripts pass syntax compilation and complete successfully against the retained local data, generating 15 PDF, PNG and TIFF outputs.
- SHA-256 hashes identify the original source and publication copy for each of the 13 files.
- Independent R calculations on the retained local tables reproduced p-S6 group summaries: CTL 279.9 ± 8.4; PM 365.9 ± 49.7; LEP 814.0 ± 33.9; RAPA 194.6 ± 15.0. The three Holm-adjusted paired t-test P values were 0.079385, 0.007832 and 0.006192.
- Independent VIC identity calculations reproduced the retained cell counts and positivity thresholds in Methods S2: vimentin control/stained, 276/185 cells, thresholds 672.12 and 853.90; alpha-SMA control/stained, 132/1613 cells, thresholds 670.33 and 688.66. Positive stained counts at the two thresholds were 185/184 and 1612/1612, respectively.
- LC3B–BODIPY experiment means were independently recalculated, with enrichment calculated per cell before aggregation.
- Fiji image-to-table macros were inspected but not rerun. Python plotting was validated with NumPy 2.5.4, pandas 2.2.3, SciPy 1.18.1, Matplotlib 3.11.2 and Pillow 12.3.0. These are the release-validation versions, not a claim about the historical execution environment. Fontconfig issued cache-permission warnings, but all three scripts exited successfully and wrote their expected outputs.
- Publication adaptations alter path configuration and add explanatory headers. They are not represented as byte-identical copies of the recovered source files.

## Current statistical companions

The independent BODIPY companion was run on the final field table (36 fields, 1,364 cells): 18 experiment values matched exactly; six Holm-adjusted P values differed by less than 1.7e-16. The co-localization companion was run on the final experiment-value table: 12 group summaries and 12 unadjusted paired P values matched to numerical precision (maximum reported comparison difference below 3.6e-15). Inputs were not deposited with the code. Both companions were also checked by an independent reviewer; shuffling input rows preserved the results.
