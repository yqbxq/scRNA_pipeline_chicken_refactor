# UMAP Decouple Pipeline

P-R02-D moves UMAP calculation out of PCA/integration stages so integration checkpoints stay reusable when only UMAP parameters change.

## Panorama

`03a2_reduce_integrate.R` now computes PCA and integration reductions only. It writes `reduction_candidates.tsv` with the intended `umap_name` for compatibility.

`03a3_compute_umap.R` reads the `03a2` candidate manifest, runs UMAP for each candidate, and writes:

- `results/checkpoints/layers/<layer_id>/umap/*__umap.rds`
- `results/tables/integration/<layer_id>/umap_candidates.tsv`
- `results/manifests/03a3_compute_umap/_manifest.json`

`03b_integration_eda.R` uses the `03a3` UMAP candidate index for UMAP plots while retaining the `03a2` candidate metadata for integration diagnostics.

`03c_cluster.R` clusters on the selected PCA/integration reduction. UMAP is optional metadata for finalized objects and is not used by `FindNeighbors`.

## Spatial

`02a_spatial_integration_eda.R` now computes spatial PCA/integration reductions and diagnostics only. UMAP generation moves to `02a2_compute_umap_spatial.R`, which writes:

- `results/spatial/checkpoints/spatial_panorama_umap.rds`
- `results/spatial/tables/spatial_02a2_compute_umap/spatial_umap_index.tsv`
- `results/spatial/figures/spatial_02a2_compute_umap/*_by_section.png`
- `results/manifests/spatial_02a2_compute_umap/_manifest.json`

`02b_finalize_spatial_clustering.R` continues clustering on PCA/integration reductions instead of UMAP.

## Parameters

Shared UMAP parameters are exported by `workflow/02lib/common.sh` and passed through the environment registry:

- `UMAP_N_NEIGHBORS`
- `UMAP_MIN_DIST`
- `UMAP_SPREAD`
- `UMAP_SEED`
- `UMAP_METRIC`
- `UMAP_LOCAL_CONNECTIVITY`

Changing these parameters invalidates `03a3_compute_umap` under fingerprint checkpointing without forcing PCA/integration recomputation.
