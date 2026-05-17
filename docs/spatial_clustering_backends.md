# Spatial Clustering Backends

Module 03 writes one clustering variant per backend under `results/spatial/checkpoints/clustering_variants/`.

- `b1_seurat_snn`: required R backend. It scans Seurat SNN resolutions and chooses the cluster count closest to `target_clusters`.
- `b2_bayesspace`: optional R backend. Missing or incompatible BayesSpace support is recorded as a failed backend without blocking `b1`.
- `b3_spagcn`: optional Python bridge using `SPAGCN_PY_BIN` and `workflow/04python/cluster_spagcn.py`.
- `b4_stagate`: optional Python bridge using `STAGATE_PY_BIN` and `workflow/04python/cluster_stagate.py`.

`selected_backend.tsv` records `default_choice`, `override_choice`, and `final_choice`. The default is `b1_seurat_snn`; project-specific overrides go in `config/spatial_clustering_override.tsv`.
