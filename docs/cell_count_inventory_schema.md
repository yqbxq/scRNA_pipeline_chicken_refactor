# Cell Count Inventory Schema

P-R01-A introduces a shared inventory helper for cell-count based evidence gating. The helper is infrastructure only; business stages start consuming it in later PRs.

## Inventory Columns

| Column | Meaning |
|---|---|
| `cluster_id` | Cell type, subtype, region, or subcluster label being evaluated. |
| `sample_id` | Biological sample or section identifier. |
| `group_id` | Comparison group such as `syf` or `f5`. |
| `layer_id` | Optional layer context; blank when the inventory is unlayered. |
| `parent_cluster_id` | Optional parent label used by `merge_to_parent`. |
| `n_cells` | Number of cells/spots in this cluster-sample-group row. |
| `total_umi` | Sum of the selected count column. |
| `median_umi_per_cell` | Median count depth for the row. |
| `fraction_of_sample` | Row cell count divided by all cells from the sample. |
| `fraction_of_group` | Row cell count divided by all cells from the group. |
| `fraction_of_cluster` | Row cell count divided by all cells in the cluster. |

## Helper API

- `build_cell_count_inventory()` builds the inventory from Seurat `meta.data` or a metadata `data.frame`.
- `resolve_thresholds()` resolves env defaults, `cell_count_thresholds.tsv`, and per-comparison overrides.
- `classify_cluster_eligibility()` classifies one `(cluster, comparison)` pair.
- `classify_all_cluster_eligibilities()` produces a `cluster_eligibility.tsv`-ready table.

Threshold priority is:

1. `comparisons.tsv` override columns.
2. `cell_count_thresholds.tsv` exact `celltype`.
3. `cell_count_thresholds.tsv` wildcard `*`.
4. Environment defaults.

## Metadata

`metadata/cell_count_thresholds.tsv` is copied into projects by `init_project.sh`. It contains a wildcard default row plus editable examples.

`comparisons.tsv` gains four optional override columns:

- `min_cells_override`
- `min_samples_override`
- `total_umi_override`
- `single_sample_frac_override`

Use `-` or an empty value to keep the threshold table or env default.
