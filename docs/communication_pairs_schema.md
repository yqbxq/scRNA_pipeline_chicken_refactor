# Communication Pairs Schema

`metadata/communication_pairs.tsv` is generated from Tier 1 questions and is the
only input table for 07 communication.

Required columns:

| Column | Meaning |
|---|---|
| `pair_id` | Stable communication task ID. |
| `layer_scope` | Target layer, such as `panorama` or `GC_subcluster`. |
| `sender` / `receiver` | CSV, `*`, or prefix selectors resolved in the chosen cell-type column. |
| `condition_split_var` / `condition_split_values` | Optional condition split, for example `group_id` and `syf,f5`. |
| `tool` | `cellchat`, `nichenet`, `both`, or differential placeholders. |
| `communication_mode` | `baseline`, `condition_split`, or screen/differential mode. |
| `receiver_gene_program_source` | `receiver_marker`, `condition_deg`, or `none`. |
| `baseline_marker_comparison_id` | A01/A02 marker comparison ID for baseline receiver programs. |
| `receiver_deg_comparison_id` | D01/D02/D03 condition DEG comparison ID for split receiver programs. |
| `direction_filter` | `yes` keeps sender->receiver LR rows; `no` keeps full CellChat network. |
| `requires_cell_subtype` | `yes` requires runtime `cell_subtype`; no fallback is allowed. |

07b NicheNet resolves gene programs through
`results/tables/deg/gene_program_registry.tsv`:

- `receiver_gene_program_source=condition_deg` uses `receiver_deg_comparison_id`.
- `receiver_gene_program_source=receiver_marker` uses `baseline_marker_comparison_id`.
- `receiver_gene_program_source=none` skips NicheNet gene-program analysis.

Missing registry rows or missing gene files are recorded as failed/skipped
tasks. 07b no longer falls back to ad hoc `FindMarkers`, and it no longer treats
`pair_id` as a DEG `comparison_id`.

When `requires_cell_subtype=yes`, both CellChat and NicheNet must use the
`cell_subtype` metadata column. Missing `cell_subtype` or unresolved strict
sender/receiver labels is a hard runtime failure.
