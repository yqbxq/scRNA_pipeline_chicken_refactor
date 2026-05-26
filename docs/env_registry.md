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
| `H5AD_EXPORT_GATES` | empty | Empty means recommended hard-coded H5AD gate hooks; otherwise comma-separated explicit gate list. |
| `H5AD_EXPORT_ON_FAILURE` | `skip` | H5AD hook failure policy: `skip`, `warn`, or `error`. |
| `H5AD_PYTHON_FALLBACK_TO_RDS` | `yes` | Python H5AD helpers may fall back to legacy RDS/explicit input paths when H5AD is missing. |
| `CELL_COUNT_INVENTORY_THRESHOLD_OVERRIDE_TSV` | `metadata/cell_count_thresholds.tsv` | Per-celltype threshold override table for inventory eligibility. |
| `CELL_COUNT_INVENTORY_MIN_CELLS_PER_SAMPLE` | `20` | Default per-sample minimum cell count. |
| `CELL_COUNT_INVENTORY_MIN_SAMPLES_PER_GROUP` | `2` | Default minimum biological samples per group. |
| `CELL_COUNT_INVENTORY_MIN_TOTAL_UMI` | `10000` | Default minimum total UMI per group. |
| `CELL_COUNT_INVENTORY_MAX_SINGLE_SAMPLE_FRAC` | `0.8` | Default single-sample dominance cutoff. |
| `CELL_COUNT_INVENTORY_MIN_CELLS_SOFT_FLOOR` | `5` | Default soft floor for module-score-only fallback. |
| `CELL_COUNT_INVENTORY_ENABLE_MERGE_TO_PARENT` | `yes` | Allow merge-to-parent evidence tier when parent metadata exists. |
| `CELL_COUNT_INVENTORY_FAIL_ON_MISSING_PARENT` | `no` | Reserved strictness flag for downstream parent merge consumers. |
| `MODULE_INVENTORY_VERSION` | `1.0` | Version marker for inventory helper outputs. |
| `MODULE_04C_VERSION` | `1.2` | Version marker for scRNA 04c subcluster EDA inventory artifacts. |
| `MODULE_SPATIAL_04C_VERSION` | `1.2` | Version marker for spatial 04c subcluster EDA inventory artifacts. |
| `MODULE_05D_VERSION` | `1.1` | Version marker for DEG EDA evidence-tier coverage outputs. |
| `MODULE_07A_CELLCHAT_VERSION` | `1.1` | Version marker for CellChat inventory-gated outputs. |
| `MODULE_07B_LIANA_VERSION` | `0.1` | Version marker for the 07b LIANA+ consensus placeholder. |
| `MODULE_07C_NICHENET_VERSION` | `1.1` | Version marker for NicheNet inventory-gated outputs. |
| `MODULE_07D_CONSENSUS_VERSION` | `0.1` | Version marker for the 07d communication consensus placeholder. |
| `MODULE_07E_EDA_VERSION` | `0.1` | Version marker for the 07e communication EDA layer. |
| `INVENTORY_GATE_REQUIRED` | `yes` | Require 04c `cluster_eligibility_tsv` before communication modules run. |
| `INVENTORY_GATE_ALLOWED_TIERS_COMMUNICATION` | `primary,exploratory,primary_merged` | Evidence tiers allowed to enter CellChat/NicheNet communication analysis. |
| `INVENTORY_GATE_STRICT_MODE` | `no` | When enabled, a cluster must pass all matching eligibility rows to enter communication analysis. |
| `COMMUNICATION_FAIL_ON_NO_PRIMARY` | `no` | Reserved strictness flag for failing communication when no cluster passes inventory gate. |
| `COMMUNICATION_REQUIRE_FULL_PIPELINE` | `no` | When `yes`, 07b/07c/07d/07e failures block the communication stage instead of writing placeholder manifests. |
| `LIANA_CONSENSUS_ENABLED` | `yes` | Enables the 07b LIANA+ placeholder stage; full implementation lands in P-R03-B. |
| `MULTINICHENET_ENABLED` | `auto` | Controls the 07c NicheNet/MultiNicheNet stage. |
| `COMMOT_ENABLED` | `auto` | Reserved for the future COMMOT spatial communication extension. |
| `MODULE_03A3_VERSION` | `1.0` | Version marker for the independent panorama UMAP stage. |
| `MODULE_SPATIAL_02A2_VERSION` | `1.0` | Version marker for the independent spatial UMAP stage. |
| `UMAP_N_NEIGHBORS` | `30` | Shared UMAP neighbor count for 03a3 and spatial 02a2. |
| `UMAP_MIN_DIST` | `0.3` | Shared UMAP minimum distance. |
| `UMAP_SPREAD` | `1.0` | Shared UMAP spread. |
| `UMAP_SEED` | inherits `RANDOM_SEED` | Shared UMAP random seed. |
| `UMAP_METRIC` | `cosine` | Shared UMAP distance metric. |
| `UMAP_LOCAL_CONNECTIVITY` | `1` | Shared UMAP local connectivity parameter. |
