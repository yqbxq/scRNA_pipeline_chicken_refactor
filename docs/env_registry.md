# Environment Registry

This file records project-level environment variables that are exported through `workflow/02lib/common.sh` and passed into runtime environments through `workflow/02lib/env_registry.sh`.

| Variable | Default | Purpose |
|---|---:|---|
| `CHECKPOINT_MODE` | `mtime` | Selects checkpoint behavior: `mtime`, `fingerprint`, or `off`. |
| `CHECKPOINT_FINGERPRINT_FALLBACK_ON_ERROR` | `mtime` | Controls fingerprint failure handling: warn and continue with mtime semantics, or `error`. |
| `CHECKPOINT_FINGERPRINT_VERBOSE` | `no` | Prints current and saved fingerprint values during stale checks when set to `yes`. |
| `CHECKPOINT_LARGE_FILE_THRESHOLD_MB` | `100` | Uses large-file mixed digests above this size for heavy binary inputs. |
| `MODULE_CHECKPOINT_FINGERPRINT_VERSION` | `1.0` | Version marker for the checkpoint fingerprint infrastructure. |
| `H5AD_CONTRACT_FILE` | `metadata/h5ad_export_contract.tsv` | Project-level H5AD export contract path. |
| `H5AD_CONTRACT_FAIL_ON` | `error` | Contract violation policy: `error`, `warn`, or `none`. |
| `H5AD_CONTRACT_VERBOSE` | `no` | Prints detailed contract violations when enabled by callers. |
| `MODULE_H5AD_CONTRACT_VERSION` | `1.0` | Version marker for H5AD contract validators. |
| `H5AD_EXPORT_BACKEND` | `zellkonverter` | H5AD writer backend: `zellkonverter`, `sceasy`, `mock`, or `auto`. |
| `H5AD_EXPORT_ASSAY` | `auto` | Assay to mirror into H5AD. |
| `H5AD_EXPORT_LAYERS` | `auto` | Comma-separated Seurat layers/slots to export, or automatic counts/logcounts selection. |
| `H5AD_EXPORT_REDUCTIONS` | `auto` | Comma-separated reductions to export into `obsm`, or all available reductions. |
| `H5AD_EXPORT_COMPRESSION` | `gzip` | Compression mode for real H5AD writers. |
| `H5AD_EXPORT_SPATIAL_PER_SECTION` | `yes` | Writes one spatial H5AD per section when enabled. |
| `H5AD_EXPORT_CONTRACT_FAIL_ON` | inherits `H5AD_CONTRACT_FAIL_ON` | H5AD export-specific contract policy. |
| `MODULE_90_EXPORT_H5AD_VERSION` | `1.0` | Version marker for the 90 H5AD export stage. |
| `MODULE_03A3_VERSION` | `1.0` | Version marker for the independent panorama UMAP stage. |
| `MODULE_SPATIAL_02A2_VERSION` | `1.0` | Version marker for the independent spatial UMAP stage. |
| `UMAP_N_NEIGHBORS` | `30` | Shared UMAP neighbor count for 03a3 and spatial 02a2. |
| `UMAP_MIN_DIST` | `0.3` | Shared UMAP minimum distance. |
| `UMAP_SPREAD` | `1.0` | Shared UMAP spread. |
| `UMAP_SEED` | inherits `RANDOM_SEED` | Shared UMAP random seed. |
| `UMAP_METRIC` | `cosine` | Shared UMAP distance metric. |
| `UMAP_LOCAL_CONNECTIVITY` | `1` | Shared UMAP local connectivity parameter. |
