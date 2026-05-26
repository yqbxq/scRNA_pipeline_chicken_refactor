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
| `activation_policy` | Runtime policy: `always`, `auto_if_min_cells`, or `derived_from_split`. |
| `min_sender_cells` / `min_receiver_cells` | Positive integer sender/receiver thresholds for auto-gated rows. |
| `min_cells_per_condition` | Minimum `sender_n + receiver_n` after pair subset, condition split, and sender/receiver resolution. |
| `fallback_pair_id` | Concrete Tier2 baseline pair ID to reference when a split is skipped. |
| `derived_from_pair_id` | Concrete split pair IDs used by derived communication rows. |
| `run_baseline_if_split_fails` | Boolean; report baseline fallback for skipped split rows. |
| `requires_all_derived_inputs_pass` | Boolean; derived rows such as `F25` require all split inputs to pass. |
| `receiver_gene_program_source` | `receiver_marker`, `condition_deg`, or `none`. |
| `baseline_marker_comparison_id` | A01/A02 marker comparison ID for baseline receiver programs. |
| `receiver_deg_comparison_id` | D01/D02/D03 condition DEG comparison ID for split receiver programs; must not point to D05 or any `global_context` row. |
| `direction_filter` | `yes` keeps sender->receiver LR rows; `no` keeps full CellChat network. |
| `requires_cell_subtype` | `yes` requires runtime `cell_subtype`; no fallback is allowed. |

07c NicheNet resolves gene programs through
`results/tables/deg/gene_program_registry.tsv`:

- `receiver_gene_program_source=condition_deg` uses `receiver_deg_comparison_id`.
- `receiver_gene_program_source=receiver_marker` uses `baseline_marker_comparison_id`.
- `receiver_gene_program_source=none` is not a valid executable NicheNet gene-program source.

Missing registry rows, invalid roles, missing `top_gene_tsv`, or empty gene
sets are recorded with top-level `status=missing_gene_program` and the specific
cause in `deg_status`. 07c does not run inline `FindMarkers`, does not infer
paths from `pair_id`, and does not fall back across registry path columns.

07 runtime reads only this resolved Tier2 table. Tier1
`metadata/analysis_questions.tsv` policy values are not reinterpreted by
`07a_cellchat.R` or `07c_nichenet.R`.

For `activation_policy=auto_if_min_cells`, 07a/07c compute
`sender_n`, `receiver_n`, and `condition_pair_cell_n = sender_n + receiver_n`
after pair-specific subset, condition split, and sender/receiver role
resolution. The split task runs only when all three thresholds pass.

Skipped split tasks write 0-row schema-valid result TSVs and index rows with
`status=skipped_low_cells`, `success=false`, `result_copied=no`, and separate
fallback fields. Skipped tasks must not write fake method RDS files and must not
copy baseline paths into primary result fields.

Derived communication rows use `activation_policy=derived_from_split`; they are
not executable by 07a/07c. 07e evaluates
`derived_communication_eligibility.tsv`; `F25_GC_internal_diff` is eligible only
when all required expanded `F08` split inputs pass.

When `requires_cell_subtype=yes`, both CellChat and NicheNet must use the
`cell_subtype` metadata column. Missing `cell_subtype` or unresolved strict
sender/receiver labels is a hard runtime failure.

07a/07c additionally require the 04c `cluster_eligibility_tsv` inventory gate.
By default only `primary`, `exploratory`, and `primary_merged` clusters are kept
for communication. The method index tables expose:

| Column | Meaning |
|---|---|
| `inventory_gate_passed` | Whether any cells for the task remained after evidence-tier filtering. |
| `inventory_gate_reason` | Empty on pass, otherwise the filtered evidence tiers or missing eligibility reason. |
| `inventory_gate_eligible_clusters` | Comma-separated cluster labels that passed the communication inventory gate. |

The gate audit tables `cellchat_inventory_gate_log_tsv` and
`nichenet_inventory_gate_log_tsv` record one row per evaluated cluster label.
