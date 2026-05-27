# Spatial Communication Pipeline

P-R03-G adds ST stage 08 for COMMOT-backed spatial communication support.

## Stage Layout

| Stage | Script | Output |
|---|---|---|
| 08a | `workflow/05single_script/spatial/08a_spatial_communication_io.R` | `commot_input_manifest.tsv` and `commot_lr_candidates.tsv`. |
| 08f | `workflow/04python/07f_commot_spatial.py` | `commot_lr.tsv`, one row per section and LR candidate. |
| 08b | `workflow/05single_script/spatial/08b_spatial_communication_eda.R` | `commot_spatial_summary.tsv` and a review report. |

`60_run_spatial_extensions.sh` runs this block after the deconvolution gate and
before niche derivation when `COMMOT_ENABLED != no`.

## Runtime Behavior

08f uses Python COMMOT when `commot` and `anndata` are available. If the runtime
is unavailable and `COMMOT_REQUIRE_RUNTIME=no`, the stage writes a schema-valid
skip table so downstream manifests remain explicit. Set
`COMMOT_REQUIRE_RUNTIME=yes` on the server when a missing COMMOT runtime should
block the run.

The COMMOT API used for the real runtime is
`commot.tl.cluster_communication_spatial_permutation`, which accepts an AnnData
object, an LR table, a distance threshold, and a cluster label column.

## 07d Integration

08b writes `results/spatial/tables/08_commot/commot_spatial_summary.tsv`.
07d reads this path through `COMMOT_SPATIAL_SUMMARY_TSV` and fills
`commot_spatial_hit` before assigning evidence tiers. Spatial support can promote
an LR axis only through the 07d evidence matrix; 07e reports the resulting tier.
