# Microscopy processing code

This directory contains Fiji macros and Python scripts used for VIC identity characterization, p-S6 immunofluorescence and LC3B–BODIPY co-localization in the CAVD study. It accompanies the main manuscript, Supplementary Methods S2–S4 and S6 and the microscopy summaries in Additional file 3. Code is covered by the repository's MIT license; documentation is covered by CC BY 4.0. The GSEA Zenodo v1.0.0 archive predates this addition and does not contain the microscopy code.

## Scope and provenance

The measurement algorithms and numerical settings are preserved from the recovered analysis files. Publication copies replace computer-specific input and output paths with folder-selection dialogs. The Python scripts already accept portable paths. `provenance_manifest.json` records the original source hash, publication-copy hash, source alias and each adaptation. No source images or participant-identifying information are included.

Original filenames and figure labels inside the code predate the final manuscript. Their current correspondence is:

| Code | Function | Current manuscript location |
| --- | --- | --- |
| `ps6/macros/PS6_scan_v2_1.ijm` | LIF metadata, dynamic range and acquisition checks | Section 2.13; Tables S4–S6 |
| `ps6/macros/PS6_percell.ijm` | DAPI segmentation; perinuclear p-S6 measurements | Section 2.13; Methods S2; Figure 3d |
| `ps6/macros/PS6_suppfig_overlay.ijm` | Segmentation and measurement-region overlays | Figure S6 |
| `ps6/scripts/make_ps6_figure_v2.py` | p-S6 intensity, distributions and control-percentile summaries | Intensity endpoint in Figure 3d; additional analysis views |
| `ps6/scripts/make_supp_figs.py` | Focus, cell-density, ring-area and field-level checks | Figures S7–S9 |
| `colocalization/macros/COLOC_scan_v1.ijm` | Metadata and whole-field intensity histograms | Methods S3.1 and S3.3 |
| `colocalization/macros/COLOC_A_blind_export.ijm` | Coded TIFF export | Methods S3.2 |
| `colocalization/macros/COLOC_A2_outline.ijm` | Interactive blinded cell outlining | Methods S3.2 |
| `colocalization/macros/COLOC_B_measure.ijm` | Per-cell area, intensity, Manders and Pearson measurements | Methods S3.3; Figure 5e–h |
| `colocalization/macros/COLOC_C_representative.ijm` | Selected field and inset display | Figure 5d |
| `colocalization/macros/COLOC_D_supplementary.ijm` | Staining controls and threshold-mask illustration | Figures S4–S5 |

The plotting scripts produce historical stand-alone panels, not the complete final nine-figure layout. Their statistical summaries can be checked against the manuscript. The recovered p-S6 representative-image macro was an earlier cropped-image version and is excluded: the final Figure 3c uses whole fields. This release does not claim that this earlier macro generated the final whole-field layout. Independent BODIPY measurement and statistics are documented in [bodipy/README.md](bodipy/README.md). Co-localization experiment-level statistics are documented in [colocalization/scripts/STATISTICS_README.md](colocalization/scripts/STATISTICS_README.md). The complete final figure assembly is outside this directory's scope.

## Software

The manuscript records Fiji/ImageJ 1.54p, Java 21.0.7 (64-bit) and Bio-Formats 8.5.0 for microscopy processing. Install Fiji with Bio-Formats to use the `.ijm` files. Python plotting requires Python 3, NumPy, pandas, SciPy, Matplotlib and Pillow; dependencies are listed in `requirements.txt`. The requirements file describes package names, not a reconstructed historical lockfile.

Open each macro in Fiji's script editor and run it as ImageJ Macro language. The publication versions prompt for input/output folders. Use a separate empty output folder because measurement macros initialize their output CSVs. Folder dialogs choose a directory; they do not ask for an individual LIF file. Where LIF files have `.lifext` companions, keep both files together.

## VIC identity characterization

`characterization/macros/IF_identity_percell.ijm` is the dedicated v2.4 analysis macro for vimentin and alpha-SMA, with each marker's secondary-only control. It is not a renamed p-S6 macro. The shared nuclear segmentation and 5-micrometre ring geometry follow Methods S2. The macro records all attempted raw fields and technical status/flags. Start with `PROBE=true`; use `PROBE=false` for the complete dataset, with `BATCH=false`. Input subdirectories and filenames identify `VIM`/`SMA` and stained/control samples; names containing `sac` identify the secondary-only control.

The Python script applies the documented acquisition-format restriction before calculating marker-specific thresholds from the appropriate control. It verifies identical FITC exposure within each marker. Thresholds are control mean + 3 SD and control maximum. It reports positivity descriptively for one stained/control coverslip per marker, not three independent biological replicates. Retain both complete CSV exports (including acquisition columns).

```bash
IDENTITY_DATA=/path/to/identity-tables IDENTITY_FIGURES=/path/to/output python microscopy/characterization/scripts/make_FigureS6_MolMed.py
```

The historical output `FigureS6_identity` corresponds to current Figure S3. Its `percell.csv` includes `marker`, `role`, `file`, `series`, `ringMean`, nuclear/ring area and spatial measurements; `perfield.csv` includes status, acquisition width, pixel size and FITC exposure. These files are separate from the identically named p-S6 tables. The retained-format counts are276 vimentin-control cells,185 vimentin-stained cells,132 alpha-SMA-control cells and1613 alpha-SMA-stained cells.

## p-S6 workflow

1. Run `PS6_scan_v2_1.ijm` on the raw LIF directory to generate `PS6_scan.csv` and acquisition logs. The scanner reports both raw and cleared series; this does not make cleared images eligible for quantitative analysis.
2. In `PS6_percell.ijm`, leave `PROBE=true` for an initial small run. Check the log and CSVs; set `PROBE=false` for the complete dataset. Keep `BATCH=false` for the documented Bio-Formats import behavior. Choose a fresh quantification output folder.
3. The macro segments the DAPI channel with Gaussian sigma 2 pixels, Otsu thresholding, hole filling, watershed, nuclear area 30–600 square micrometres and circularity at least 0.35, excluding border-touching nuclei. The measurement ring is the 5-micrometre nuclear dilation minus the nucleus. p-S6 signal is `ringMean - fieldMode` on the raw channel. Two recorded duplicate/misfiled series are excluded by their exact source identifiers; computationally cleared duplicates are excluded separately.
4. Run `PS6_suppfig_overlay.ijm` to show the same segmentation geometry in the fixed, recorded representative fields. The five selected fields and display limits are explicitly recorded in its settings.
5. Put `percell.csv`, `perfield.csv` and `PS6_scan.csv` together in a data folder. Run the plotting scripts with its path and separate output locations:

```bash
PS6_DATA=/path/to/data PS6_FIGURES=/path/to/main-output python microscopy/ps6/scripts/make_ps6_figure_v2.py
PS6_DATA=/path/to/data PS6_FIGURES=/path/to/supp-output python microscopy/ps6/scripts/make_supp_figs.py
```

The main p-S6 inference uses each cell's `net`, the median within each field, then the median of those field medians within each coverslip. Three independent experiments are paired across groups. The three prespecified contrasts are CTL–PM, PM–LEP and LEP–RAPA, with two-sided paired t-tests and Holm adjustment. The plotted p-S6 control-percentile summaries are descriptive distribution analyses.

### p-S6 input columns

- `percell.csv`: `file`, `group`, `rep`, `series`, `seriesName`, `cell`, `ringMean`, `ringMedian`, `ringMax`, `ringArea_um2`, `nucArea_um2`, `fieldMode`, `net`, `pxSize_um`. Groups are `CTL`, `PM`, `LEP`, `RAPA`, and `Sac` (secondary-only control).
- `perfield.csv`: source/series identifiers; `status`; `nCells`; `fieldMode`; `fieldMax`; `dapiNucEdge`; and the acquisition/QC columns generated by the macro. All columns should be retained because supplementary plots use multiple QC measures.
- `PS6_scan.csv`: series metadata and channel statistics generated by the scanner.

`ringMedian_field` and `net_field` in the original field table summarize each ring's median, not the study's ring-mean endpoint. They are retained as legacy measurements but are not substitutes for the primary endpoint. The corrected supplementary plotting script reconstructs the field endpoint from `percell.net`. `focusC1_legacy` is not an exclusion criterion. Primary analysis retains all 149 eligible fields; DAPI focus index at least 0.25 is a sensitivity analysis.

## LC3B–BODIPY workflow

1. Run `COLOC_scan_v1.ijm` for acquisition and whole-field histogram checks. The study's BODIPY and LC3B thresholds were fixed in the pilot dataset; these histograms provide a post hoc check.
2. Run `COLOC_A_blind_export.ijm` on the LIF folder. Select a colocalization working folder. It writes coded TIFFs into `blind/`, initializes `roi/` and saves `KEY_DO_NOT_OPEN.csv`. Seed 20260803 and fixed display ranges are preserved. Keep the key separate from the person outlining cells until outlining is complete.
3. Run `COLOC_A2_outline.ijm`, choosing that same working folder. Outline each eligible cell and press `T` to add it to ROI Manager. The interactive macro saves each field's ROIs. Follow the eligibility criteria in Methods S3.2; the macro does not decide biological eligibility itself.
4. Run `COLOC_B_measure.ijm` in the same working folder. It reads `blind/Fxxx.tif` and `roi/Fxxx.zip` (or a single ROI), then writes `percell_blind.csv`.
5. Decode group, source file and series identifiers with the retained key after outlining. The three experiment identifiers are assigned from the experimental records, not inferred from fluorescent intensities. Group comparisons use the three matched experiment means, as specified in Methods S3.4 and Table S13.
6. `COLOC_C_representative.ijm` and `COLOC_D_supplementary.ijm` read the same blind images/ROIs and use recorded panel selections. They do not rerun the field-ranking selection algorithm. Historical `Figure5F` output labels refer to the co-localization panels now in Figure 5d. Selected fields in the default D4 set are F053, F017 and F004. Recreating the source identifiers requires the original blinded exports and key; relabelled or newly enumerated input directories should not be assumed to reproduce the same codes.

### Co-localization outputs and interpretation

| Column | Definition |
| --- | --- |
| `code`, `roi` | Blinded field identifier and cell ROI index |
| `areaPx`, `areaUm2` | Cell area in pixels and calibrated square micrometres |
| `meanD`, `meanG`, `meanR` | Mean DAPI, BODIPY and LC3B intensities |
| `fracG`, `fracR` | Fractions of cell pixels at or above BODIPY 100 and LC3B 90, respectively, on 8-bit data |
| `M1` | BODIPY intensity within LC3B-positive pixels divided by total cell BODIPY intensity |
| `M2` | LC3B intensity within BODIPY-positive pixels divided by total cell LC3B intensity |
| `pearson` | Pearson coefficient over unthresholded cell pixels |
| `enrich` | Per-cell `M2 / fracG` |
| `nucFrac`, `fracG_inNuc` | Historical exploratory DAPI-dependent measurements; not validated study endpoints and not used in the manuscript conclusions |

The macro uses `-1` as a sentinel for undefined ratios/correlations. Treat these as missing values when summarizing the corresponding endpoint. For enrichment, cells with `fracG=0` are excluded only from that endpoint. Calculate enrichment per cell before averaging within experiment; dividing the experiment-mean M2 by the experiment-mean area fraction gives a different quantity. The main droplet-bearing groups contain 212 cells; nine zero-area cells are excluded from enrichment only. The full study dataset includes 275 cells in 60 fields including staining controls.

M2 is intensity-weighted; its numerator does not additionally require LC3B intensity to exceed 90. An enrichment reference of 1 is descriptive under a simplified spatial-independence model. It is not evidence of autophagic flux.

## Data access and verification limits

Microscopy images, ROIs, coded-image keys and cell/field tables are available from the corresponding authors on reasonable request, as stated in the manuscript. Experiment-level values and summaries for the main co-localization endpoints are in Additional file 3, Table S13; p-S6 summaries are in Tables S4–S6. This code release does not change their access terms.

The three recovered Python plotting scripts were syntax-checked and successfully executed with retained local tables during release preparation. Independent R calculations also reproduced the reported p-S6 and VIC identity summaries. Package versions and execution limits are recorded in VALIDATION.md. The macros were inspected against the reported parameters; a complete new image-to-table run was not performed for this release. Retain the original image identifiers, calibration and acquisition bit depth when reproducing the analyses.

## Version and citation

Use the exact Git commit cited in the manuscript to retrieve this code snapshot. The repository history records this publication copy and later changes. Cite the repository separately from the GSEA-only Zenodo v1.0.0 archive. `FILES.sha256` records every file in this directory except the checksum file itself. Validation scope is documented in `VALIDATION.md` and the assay-specific READMEs.
