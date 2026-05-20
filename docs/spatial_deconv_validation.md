# Spatial Deconvolution Validation

`07f_deconvolution_validation.R` is opt-in through `--with-deconv-validation`. It is reserved for scDesign3-based synthetic spot validation after the validation R environment is available.

The script writes:

| File | Purpose |
| --- | --- |
| `validation_manifest.tsv` | Validation status and provenance. |
| `synthetic_truth.tsv` | Synthetic spot truth proportions. |
| `method_summary.tsv` | Per-method validation metrics. |
| `report.md` | Review report for the `spatial_deconv` gate. |

If `scDesign3` or the frozen reference is unavailable, the manifest records `skipped_no_packages` or `skipped_no_reference` without blocking the default ST extension run.
