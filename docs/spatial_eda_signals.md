# Spatial EDA Signals

ST 01 pre-QC writes five triage signals: `high_mito`, `low_feature`, `sparse_tissue`, `sectioning_artifact`, and `mito_qc_unreliable`.

ST 02 post-filter QC adds four fixed signals in `post_filter_triage.tsv`:

| signal | meaning |
| --- | --- |
| `excessive_drop` | Retention is below the configured `excessive_drop_threshold`. |
| `boundary_drop` | Dropped spots form a spatially clustered block. |
| `drift_from_pre_qc` | Pre-QC warned about high mito or low features, but filtering did not remove matching spots. |
| `unbalanced_sections` | Retention differs by more than 25 percentage points across sections. |

These signals are review aids for the `spatial_post_qc` gate and are not inferential tests.
