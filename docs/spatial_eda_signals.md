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

The `marker_spatial_coherence.tsv` value `morans_i_approx` is a lightweight k-nearest-neighbor Spearman proxy: each spot is compared with the mean expression of its nearest spatial neighbors, then averaged across top HVGs. It is an EDA stability signal, not a formal Moran's I test.

ST 04 sub-region EDA adds six signals in `tables/spatial_04c_subcluster_eda/triage.tsv`:

| signal | meaning |
| --- | --- |
| `small_subcluster` | One or more SNN subclusters fell below `SPATIAL_SUBCLUSTER_SMALL_CLUSTER_FRAC` and were retained with a `_merged_small` suffix. |
| `weak_marker_overlap` | The fraction of subclusters with at least one winning panel marker hit is below `SPATIAL_SUBCLUSTER_TRIAGE_OVERLAP_THRESHOLD`. |
| `low_spatial_coherence` | Mean same-label spatial kNN coherence for `sub_region` is below `SPATIAL_SUBCLUSTER_TRIAGE_COHERENCE_THRESHOLD`. |
| `all_undetermined` | Every assigned sub-region for the layer is `mixed_or_uncertain`. |
| `single_subcluster` | A layer produced only one sub-region label, so it is not informative as a drilldown. |
| `cross_section_imbalance` | With two or more sections, a sub-region's relative section fraction differs by more than `SPATIAL_SUBCLUSTER_TRIAGE_IMBALANCE_THRESHOLD`, default `0.5`. |
