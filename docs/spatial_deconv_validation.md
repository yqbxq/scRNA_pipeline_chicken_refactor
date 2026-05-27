# Spatial Deconvolution Validation

`07f_deconvolution_validation.R` is opt-in through `--with-deconv-validation`. It compares completed deconvolution method proportions against synthetic truth proportions. If `SPATIAL_VALIDATION_TRUTH_TSV` points to a table with `spot_id`, `cell_type`, and `true_proportion`, that table is used. Otherwise 07f generates a Dirichlet synthetic truth over the completed method spots and cell types as a lightweight local validation path. The manifest records whether `scDesign3` is available so server runs can distinguish full validation-ready environments from local smoke runs.

The script writes:

| File | Purpose |
| --- | --- |
| `validation_manifest.tsv` | Validation status and provenance. |
| `synthetic_truth.tsv` | Synthetic spot truth proportions. |
| `method_summary.tsv` | Per-method validation metrics. |
| `report.md` | Review report for the `spatial_deconv` gate. |

If no completed method outputs are available, the manifest records `skipped_no_method_outputs` without blocking the default ST extension run.
