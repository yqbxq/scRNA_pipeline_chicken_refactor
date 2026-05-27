# COMMOT Methodology

COMMOT is used here as a spatial support layer, not as a standalone source of
primary communication claims.

The pipeline passes section-level spatial H5AD files and LR candidates from 07d
to COMMOT. The server runtime should provide `commot`, `anndata`, and the normal
`py_spatial` stack. The stage requests cluster-level spatial permutation scoring
with a configurable distance threshold and permutation count.

Outputs are reduced to a pipeline contract:

- `signal_score`: COMMOT-derived sender-to-receiver score when available.
- `spatial_support`: `yes` when `signal_score > COMMOT_SIGNAL_THRESHOLD`.
- `status` and `reason`: explicit runtime or input diagnostics.

COMMOT does not override CellChat's hypothesis-only downgrade. It only supplies
the spatial layer consumed by 07d consensus.
