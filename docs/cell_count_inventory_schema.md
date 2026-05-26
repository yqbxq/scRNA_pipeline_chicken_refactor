# Cell Count Inventory Schema

P-R01-A introduces a shared inventory helper for cell-count based evidence gating. P-R01-B wires the helper into scRNA and spatial 04c subcluster EDA so inventory artifacts are emitted before downstream DE stages consume them.

## Inventory Columns

| Column | Meaning |
|---|---|
| `cluster_id` | Cell type, subtype, region, or subcluster label being evaluated. |
| `sample_id` | Biological sample or section identifier. |
| `group_id` | Comparison group such as `syf` or `f5`. |
| `layer_id` | Optional layer context; blank when the inventory is unlayered. |
| `group_var` | Metadata column used to derive `group_id`; added by 04c callers. |
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

## 04c Manifest Outputs

Both `04c_subcluster_eda.R` and `spatial/04c_subcluster_eda.R` now emit:

| Manifest key | File | Purpose |
|---|---|---|
| `cell_count_inventory_tsv` | `cell_count_inventory.tsv` | Per layer, cluster, sample/section, and comparison group count table. |
| `cluster_eligibility_tsv` | `cluster_eligibility.tsv` | Per `(layer_id, group_var, comparison_id, cluster_id)` evidence tier and recommended action. |
| `inventory_summary_per_comparison_tsv` | `inventory_summary_per_comparison.tsv` | Count of cluster-comparison rows per evidence tier for report rendering. |

The 04c reports include a `Cell-count inventory + evidence tier 建议` section with the summary table and a preview of the eligibility decisions.
