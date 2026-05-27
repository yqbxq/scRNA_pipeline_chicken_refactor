# Spatial Co-Localization Schema

`results/spatial/tables/08_commot/commot_spatial_summary.tsv` is the ST-side
contract consumed by 07d communication consensus.

| Column | Meaning |
|---|---|
| `lr_axis_id` | Canonical ligand-receptor axis ID matching 07d. |
| `section_id` | Spatial section evaluated by COMMOT. |
| `pair_id` | Communication pair ID carried from the 07d consensus table. |
| `condition_value` | Condition/stage label carried from the 07d consensus table. |
| `ligand` / `receptor` | Ligand and receptor symbols passed to COMMOT. |
| `sender` / `receiver` | Sender and receiver labels for the axis. |
| `signal_score` | COMMOT sender-to-receiver score, if available. |
| `spatial_support` | `yes` when the signal score passes `COMMOT_SIGNAL_THRESHOLD`; otherwise `no`. |
| `status` | `ok_commot_run`, `skipped_missing_commot`, `skipped_missing_h5ad`, or another explicit diagnostic status. |
| `reason` | Human-readable status detail. |
