# Spatial QC Threshold Schema

`config/spatial_qc_thresholds.tsv` controls ST spot-level filtering. A `__DEFAULT__` row is used when a section-specific row is absent.

| column | required | default | meaning |
| --- | --- | --- | --- |
| `section_id` | yes | `__DEFAULT__` | Section id, sample id, or fallback row id. |
| `qc_min_nfeature` | yes | `SPATIAL_QC_MIN_NFEATURE` | Keep spots with at least this many detected features. |
| `qc_max_nfeature` | no | unlimited | Drop likely artifacts with too many detected features. |
| `qc_min_ncount` | yes | `SPATIAL_QC_MIN_NCOUNT` | Keep spots with at least this many counts. |
| `qc_max_ncount` | no | unlimited | Drop likely artifacts with too many counts. |
| `qc_max_mito_pct` | yes | `SPATIAL_QC_MAX_MITO_PCT` | Drop spots above this mitochondrial percentage when mito QC is valid. |
| `mito_set_override` | no | empty | Manual mito gene list or file path used by the raw-object step. |
| `spatial_aware_filter` | no | `false` | When true, detect spatially clustered dropped spots for review. |
| `excessive_drop_threshold` | no | `0.5` | Post-filter triage flags sections below this retention. |

Missing optional columns remain backward-compatible with earlier ST 01 projects.
