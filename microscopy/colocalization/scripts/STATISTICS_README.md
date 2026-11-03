# Co-localization statistics companion

`reproduce_coloc_statistics.py` was prepared on 12 October 2026 to reproduce the recorded summaries and comparisons now presented in Figure 5E–H and Table S9B–C. It reads recorded experiment-level values. It does not repeat image segmentation or calculate coefficients from raw images. Current references follow the revised manuscript numbering as of 3 November 2026. The retained source/provenance references to old Figure 5e–h, Table S13B–C and Methods S3 are historical; they do not change the current mapping.

Use Python 3.9 or later; no third-party packages are required:

```sh
python reproduce_coloc_statistics.py --values /path/to/Figure5_plotted_values.csv --output /path/to/new_results
```

The input columns are `panel`, `endpoint`, `group`, `matched_experiment`, `value` and `unit`. Additional columns are allowed. The complete Figure 5 source table may be supplied; the script selects the historical lowercase input labels e–h and ignores other panels. Keep these input values unchanged; they correspond to current panels E–H.

| Panel | Endpoint | Units | Input definition |
|---|---|---|---|
| e | fracG | percent | Mean per-cell BODIPY-positive area fraction within an experiment |
| f | fracR | percent | Mean per-cell LC3B-positive area fraction within an experiment |
| g | M2 | dimensionless fraction | Mean per-cell Manders' M2 within an experiment |
| h | enrich | dimensionless ratio | Mean of valid per-cell M2/BODIPY-area-fraction ratios within an experiment |

Use the exact group labels `PM`, `PM+Lep`, `PM+Lep+Rapa` and the recorded matched experiment IDs `1`, `2`, `3`. Every endpoint must have one value for each group and experiment. Experiment pairing comes from those labels, not input row order.

Each group summary is the mean ± sample SD of three experiment means. Comparisons use two-sided paired t tests (df = 2), with unadjusted P values as specified in Section 2.12 and Methods S5. PM versus PM+Lep+Rapa is the exploratory comparison. The implementation uses the exact Student t tail for two degrees of freedom.

For panel h, the supplied experiment values must already be means of cell-level ratios. Do not calculate them by dividing mean M2 by mean BODIPY area. The original analysis excluded nine undefined ratios where BODIPY-positive cell area was zero; those cells remain in the other endpoints. This script does not infer or reconstruct those cell-level exclusions from experiment means.

The companion was run on the current figure source table. Its 12 summaries and 12 P values matched the supplied current results to numerical precision. Raw images, cell-level measurements and experiment-value inputs are not included in this code release. Table S9 provides the experiment-level values, and the manuscript's data-availability statement governs requests for source materials.
