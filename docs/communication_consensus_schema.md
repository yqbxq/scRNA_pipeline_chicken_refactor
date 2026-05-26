# Communication Consensus Schema

## 07b LIANA Output

`workflow/04python/07b_liana_consensus.py` writes `results/tables/communication/liana_consensus/liana_consensus_lr.tsv`.

Required columns:

| Column | Meaning |
|---|---|
| `lr_axis_id` | Canonical axis ID: `ligand_complex|receptor_complex|sender->receiver`. |
| `condition_value` | Condition or stage value for the per-condition LIANA run. |
| `source` / `target` | Sender and receiver cell type labels. |
| `ligand` / `receptor` | Original gene or complex names from LIANA. |
| `ligand_human` / `receptor_human` | Human-space names after optional chicken-to-human LUT mapping. |
| `liana_consensus_score` | Rank-aggregate score preference. |
| `n_methods_agreed` | Count of LIANA methods passing their local hit rule. |
| `liana_consensus_hit` | `1` when `n_methods_agreed >= LIANA_MIN_METHODS_AGREED_INSIDE`, otherwise `0`. |
| `status` / `reason` | Execution status and short diagnostic text. |

07d consumes `lr_axis_id` as the join key. 07b does not assign final communication evidence tiers.
