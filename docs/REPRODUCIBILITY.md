# Reproducing the transcriptomic analyses

## Scope and preserved sources

This repository contains the GEO transcriptomic analysis used in the manuscript: targeted GSEA, descriptive Hallmark GSEA, and the baseline transcriptomic panels. Western blot, cell-culture and microscopy analyses belong to separate source-data packages.

The files under `pre/` and the historical outputs under `results/` are preserved. The analysis implementation was committed at `6d6a7cd`; the saved GSEA run is tagged `gsea-run-v2` (`18f0f7f`). The historical baseline panels were added at `0fe1758`.

`scripts/run_analysis.R` parses the preserved R source and replaces only its input/output path assignments and, for GSEA, `RUN_BOOTSTRAP <- TRUE` with `FALSE` in memory. The frozen script on disk is unchanged. Gene selection, preprocessing, probe collapse, ranking, random seeds, `fgseaMultilevel` parameters, testing pools and consistency rules are unchanged. Source and input SHA-256 hashes are checked before execution. The wrapper writes an override log and session information to the chosen output directory.

The historical bootstrap-derived NES uncertainty columns are excluded from the current manuscript. The wrapper does not calculate them. The archived output schema retains these columns as `NA`; `GSEA_results_no_CI.tsv` omits them. This change does not alter the original NES, P values or leading-edge results.

## 1. Prepare the environment

The verified run used R 4.5.2 and Bioconductor 3.22. `environment/renv.lock` records 150 packages and their versions; `environment/required_packages.tsv` lists the direct dependencies. `environment/validation_sessionInfo.txt` records the environment used for validation.

From the repository root, restore the recorded packages into an isolated library:

```sh
Rscript --vanilla environment/restore.R
```

This creates `environment/library/`. The analysis wrapper uses that library when present. The restore command requires internet access; packages with compiled code also require the build tools appropriate to the operating system. No project startup file is installed or changed.

The validation reported below used the existing installed packages matching the recorded versions. A restore into an empty library has not been tested. The lockfile provides an environment specification, not a claim that every external package archive will remain available.

The GEO retrieval script uses Python 3 and its standard library. It requires no Python packages. Requirements for the current figure assembly are described with its own script.

## 2. Retrieve and verify public GEO inputs

```sh
python3 scripts/fetch_geo_inputs.py
```

The script downloads the gzip-compressed series matrices directly from NCBI GEO, decompresses them, and verifies the uncompressed file size and SHA-256 against `environment/geo_inputs.tsv`. The expected hashes were obtained from the matrices used in the frozen analysis and independently matched to fresh public downloads on 7 October 2026. The local `csv/` directory is excluded from version control.

| Dataset | Expected input | Uncompressed bytes |
| --- | --- | ---: |
| GSE51472 | `GSE51472_series_matrix.txt` | 7,195,403 |
| GSE12644 | `GSE12644_series_matrix.txt` | 12,091,177 |
| GSE83453 | `GSE83453_series_matrix.txt` | 15,933,353 |

Existing valid files are reused. To check without downloading, use `--verify-only`; to repeat the public download, use `--force-download`. `--data-dir DIR` selects another input directory. An optional `--reference-dir DIR` compares downloaded files against an existing collection without modifying the reference files. A mismatch stops acceptance of that file and is recorded in the JSON receipt; it does not replace an existing input.

```sh
python3 scripts/fetch_geo_inputs.py --verify-only
```

`geo_fetch_receipt.json` is written in the input directory unless `--receipt PATH` is supplied. The published validation receipt is under `docs/validation/`.

## 3. Reproduce GSEA

Choose a new output directory, then run:

```sh
Rscript --vanilla scripts/run_analysis.R --output-dir runs/reproduction --mode gsea
```

The wrapper refuses a nonempty selected output subdirectory and never writes into `pre/` or `results/`. To use a different input directory or R library, add `--data-dir DIR` or `--library DIR`.

Outputs appear in `runs/reproduction/gsea/`:

- `Supp_Table_S2_v2.tsv`: the preserved 68-row output schema; deprecated CI fields are `NA`.
- `GSEA_results_no_CI.tsv`: the same results with the three deprecated uncertainty columns removed.
- `coverage_check_v2.tsv`: input and mapped gene-set sizes.
- `cross_cohort_consistency_v2.tsv`: the five pathway consistency classifications.
- `forestplot_data.tsv`: the 15 hypothesis-driven pathway/cohort results in the historical schema.
- `combat_pca_check.pdf`: the integrated-cohort batch-correction assessment.
- `wrapper_overrides.tsv`, `run_provenance.txt` and `sessionInfo.txt`: execution records.

The 68 result rows comprise 15 hypothesis-driven tests, three separate S100A9 supplementary tests and 50 descriptive Hallmark tests in the discovery cohort. The five-test and 50-test pools are adjusted separately. S100A9 is a separate single-test analysis.

On 7 October 2026, the wrapper completed with freshly downloaded inputs. All 68 rows and all 13 retained fields exactly matched the archived TSV as strings, including NES, raw P, adjusted P, coverage and the ordered leading-edge gene lists. The comparison excludes only `NES_SD`, `NES_lower_95CI` and `NES_upper_95CI`. Details and hashes are in `docs/validation/reproduction_2026-10-07.json`.

## 4. Reproduce the baseline panels

```sh
Rscript --vanilla scripts/run_analysis.R --output-dir runs/reproduction --mode baseline
```

This evaluates `pre/draw_figure1.R` with absolute input paths and a separate output directory. It preserves the source's historical output filenames. All required packages are checked before execution, so the source's optional package-install branches are not used.

| Historical artifact | Current manuscript use |
| --- | --- |
| `results/Supp_Table_S2_v2.tsv`, rows with `bh_pool = hallmark_50pool` and `cohort_short = GSE51472` | Supplementary Figure S1, 50 Hallmark results |
| `results/combat_pca_check.pdf` | Additional archived ComBat diagnostic |
| `results/Figure1_panels/Figure1A_PCA_percohort.pdf` | Supplementary Figure S10a |
| `results/Figure1_panels/Figure1B_PCA_combat.pdf` | Supplementary Figure S2 |
| `results/Figure1_panels/Figure1C_volcano_discovery.pdf` | Supplementary Figure S10b |
| `results/Figure1_panels/Figure1D_volcano_validation.pdf` | Supplementary Figure S10c |
| `results/Figure1_panels/Figure1E_markers_discovery.pdf` | Supplementary Figure S10d |
| `results/Figure1_panels/Figure1F_markers_validation.pdf` | Supplementary Figure S10e |
| `results/Figure1_panels/Figure1G_crosscohort_forest.pdf` | Supplementary Figure S10f |
| `results/Figure1_panels/Figure1H_GO.pdf` | Supplementary Figure S10g |
| `results/Figure1_panels/Figure1I_KEGG.pdf` | Supplementary Figure S10h, archived KEGG result |

The GSEA and baseline modes can also run together with `--mode all`. Their outputs use separate `gsea/` and `baseline/` subdirectories. The current main Figure 1 is assembled separately from the frozen results; the historical `Figure1*` filenames above do not refer to that main figure.

### Baseline validation and KEGG version dependence

The 7 October 2026 baseline run reproduced 884 discovery DEGs (542 upregulated, 342 downregulated), 117 validation DEGs, PC1 percentages of 47.3%, 97.4% and 27.7% in the three individual cohorts, and the GO chemotaxis result (adjusted P = 9.430769 × 10⁻³²). Panels A–H and the separate ComBat assessment matched both extracted PDF text and rendered pixels at 72 dpi.

Historical panel I (current Supplementary Figure S10h) calls the live KEGG service through `enrichKEGG`. Its rerun differed from the archived panel: the displayed pathway names and order changed, including the archived `Phagosome` entry and the rerun `Phagocytosis` entry. The archived source does not include a frozen KEGG annotation snapshot. Therefore, the historical KEGG panel is retained as the manuscript artifact; the wrapper's live KEGG result is a dated recomputation, not an exact regeneration of that panel. This difference does not affect the targeted or Hallmark GSEA, which uses the pinned `msigdbr` 7.5.1 collection.

## Validation record

`docs/validation/` contains the public download receipt, core-result comparison, baseline PDF comparison, and sanitized run logs. These records describe what was actually executed. They contain public GEO-derived information and software provenance; no private raw patient records or local user paths are included.
