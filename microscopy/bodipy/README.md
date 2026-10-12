# Independent BODIPY assay

This directory contains retained Fiji source and a publication reproducibility companion for Supplementary Methods S6, Figure 5a–c and Table S13A. It is separate from the LC3B–BODIPY co-staining assay.

## Contents

- `bodipy_batch.ijm`: retained macro for channel export, manual cell outlines and particle measurement. The published copy is byte-identical to the retained source.
- `make_bodipy_representative_panel.ijm`: retained representative-image assembly macro with portable folder dialogs. Display settings and crop coordinates are unchanged.
- `reproduce_bodipy_statistics.py`: companion written on 12 October 2026 to reproduce the current experiment-level summaries and paired comparisons from recorded field means.
- `input_schema.csv`: column descriptions for the companion input.
- `PROVENANCE.json`: source hashes and adaptations.

The retained measurement macro uses the 0.1 μm² particle-area cutoff reported in Methods S6. The available records do not establish that this exact retained macro version generated every measurement in the final 1,364-cell dataset. The companion reproduces the published aggregation and statistics; it is not presented as the historical image-measurement script. An older field-level Kruskal–Wallis plotting script is excluded because the current analysis uses matched experiment means.

## Fiji measurement workflow

Use Fiji/ImageJ 1.54 with Bio-Formats and the standard ROI Manager and Analyze Particles tools. The assay record specifies Fiji 1.54; it does not identify the exact patch or Bio-Formats version for this independent assay.

1. Install `bodipy_batch.ijm` through **Plugins > Macros > Install**.
2. Run **A Batch prepare single-plane TIFFs**. Select input LIF/TIFF files and an empty output folder. Confirm the actual BODIPY/DAPI channel order. Export preserves raw intensities and calibration. For a stack, the macro chooses the plane with the highest mean BODIPY intensity and records the choice.
3. Run **C Draw cells, script saves**. Outline cells manually, press **T** to add each ROI, and save the matching ROI ZIP. Follow the original experiment's cell-selection records.
4. Run **B Batch quantify all fields** to create `sample_sheet.csv`, fill its recorded experiment and group assignments, then run it again. Supply the frozen raw-intensity threshold from the relevant experiment record. The default `100` is a dialog default, not a documented final-assay threshold.
5. Inspect the exported calibration, bit depth, threshold and particle-area results. The macro appends to `bodipy_rapamycin_cells.csv`; use a fresh results file for a new run to avoid duplicate rows. The statistical companion below consumes field means, not this cell-level export directly.

Keep acquisition and quantification parameters tied to their recorded assay. Do not infer channel order, treatment assignments, thresholds or experiment matching from filenames or macro defaults. The recovered historical cell export has 207 cells in 14 fields and is not included as the final 1,364-cell dataset.

## Representative panel

Run `make_bodipy_representative_panel.ijm` in Fiji and choose BODIPY, matching DAPI and output folders. The retained script expects the three named representative TIFFs and corresponding DAPI files. Their acquisition labels and treatment assignments are part of the associated laboratory records.

The script uses 512 × 512-pixel crops at 0.3251 μm/pixel, one 20 μm scale bar, a BODIPY display range of 490–2900 and a DAPI range of 455–3650. Each channel uses the same display range across all treatment groups. These settings control figure display, not quantitative thresholds. The complete final composite Figure 5 reconstruction source was not found in the retained files; this macro builds its BODIPY representative panel.

## Current statistics

Export the `fields` worksheet of the final assay workbook as UTF-8 CSV, preserving its column names and recorded `batch` labels. See `input_schema.csv`. The expected design comprises 36 fields and 1,364 cells: four fields for each of three groups in each of three matched experiments, B1–B3. The matching was confirmed by the investigators; the script reads those labels and does not derive them from row order.

Run with Python 3.9 or later; no third-party Python packages are required:

```sh
python reproduce_bodipy_statistics.py --fields /path/to/fields.csv --output /path/to/new_results
```

For each endpoint, the script averages four field means within an experiment, then reports mean ± sample SD across the three experiment values. Cell counts document sampling and do not weight experiment means. Three two-sided paired t tests use the recorded B1–B3 matching, with Holm adjustment within each endpoint. Since each test has three pairs, the implementation uses the exact Student t tail for two degrees of freedom.

Outputs are `experiment_values.csv`, `group_summary.csv`, `paired_holm_comparisons.csv` and `run_receipt.json`. The companion was run against the final field-level table: all 18 experiment values and all six adjusted P values matched the current Figure 5 source values to numerical precision.

## Data and license

This release contains code, documentation and an input schema. Raw images, ROIs, source workbooks and measurement data are not included. The manuscript's data-availability statement governs requests for those materials. Table S13 supplies the experiment-level values. Code is provided under the repository's MIT license; documentation follows the repository's stated data/documentation license.
