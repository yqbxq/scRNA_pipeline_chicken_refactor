# MultiNicheNet Usage Guide

07c uses `NICHENET_MODE=auto` by default.

| Mode | Behavior |
|---|---|
| `auto` | Use MultiNicheNet only when the runtime packages are available; otherwise fall back to legacy NicheNet. |
| `multinichenet` | Require `multinichenetr` and `SingleCellExperiment`; fail if unavailable. |
| `nichenet_legacy` | Force the existing NicheNet executor. |
| `skip` | Let the 07 stage write a placeholder 07c manifest. |

The 07c index records `multinichenet_mode_used`, `method_evidence_class`, `n_samples_per_group`, and `n_targets_in_receiver_de`. These fields give 07d enough information to separate downstream receiver-program evidence from LIANA ligand-receptor consensus.

Required runtime knobs:

| Variable | Default |
|---|---:|
| `MULTINICHENET_MIN_SAMPLES_PER_GROUP` | `2` |
| `MULTINICHENET_TOP_N_LR` | `250` |
| `MULTINICHENET_TOP_N_TARGETS` | `20` |
| `MULTINICHENET_MIN_CELLS` | `10` |
| `COMMUNICATION_RECEIVER_DE_P_THRESHOLD` | `0.05` |
