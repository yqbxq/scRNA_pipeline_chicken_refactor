# Spatial Deconvolution Validation

`07f_deconvolution_validation.R` runs by default in `60_run_spatial_extensions.sh`. It compares completed deconvolution method proportions against truth proportions and writes the question-gate layer consumed by ST-08. `SPATIAL_VALIDATION_MODE` supports `scdesign3`, `external_truth`, and `dirichlet_only`. `dirichlet_only` is a smoke/fallback mode and cannot unlock an I11 PASS; formal PASS requires scDesign3 or external truth plus thresholded metrics.

The script writes:

| File | Purpose |
| --- | --- |
| `validation_manifest.tsv` | Validation status and provenance. |
| `synthetic_generation_summary.tsv` | Truth-generation mode and provenance. |
| `synthetic_sc_metadata.tsv` | Synthetic single-cell metadata when scDesign3 mode runs. |
| `synthetic_spot_counts.rds` / `synthetic_st.rds` | Mixed synthetic spot matrix/object when scDesign3 mode runs. |
| `synthetic_truth.tsv` | Synthetic spot truth proportions. |
| `synthetic_truth_wide.tsv` | Wide truth matrix. |
| `method_prediction_long.tsv` / `method_prediction_wide.tsv` | Method predictions used for validation. |
| `method_summary.tsv` | Per-method validation metrics including RMSE, Pearson, JSD, and dominant accuracy. |
| `celltype_metrics.tsv` / `spot_metrics.tsv` | Fine-grained validation metrics including JSD. |
| `recommended_method_consistency.tsv` | Agreement between 07e recommended method and 07f best validation method. |
| `spatial_question_gate_status.tsv` | I08/I09/I10/I11/I12 deconvolution gates. |
| `report.md` | Review report for the `spatial_deconv` gate. |

If `SPATIAL_VALIDATION_MODE=scdesign3` but `scDesign3` is unavailable, the manifest records `skipped_no_packages`. If the package is available, 07f uses the inherited v09 scDesign3 engine to fit the frozen reference, simulate synthetic cells, mix them into synthetic ST spots, and write the synthetic object contracts. If no completed method outputs are available, it records `skipped_no_method_outputs` without blocking the default ST extension run.
