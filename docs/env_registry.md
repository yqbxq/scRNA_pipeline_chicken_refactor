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
| `VELOCITY_INPUT_MODE` | `auto` | `10c_scvelo_dynamical.py` input policy: `auto` tries H5AD then loom, `h5ad` requires H5AD, `loom` keeps legacy loom input. |
| `VELOCITY_H5AD_MODULE` | `GC_subcluster` | scRNA H5AD export module used by H5AD-first velocity when `input_h5ad_path` is not supplied in the velocity reference index. |
| `SPATIAL_NEIGHBORHOOD_H5AD_MODULE` | `spatial_03_region` | Spatial H5AD export module consumed first by ST-06d Squidpy neighborhood. |
| `SPATIAL_NEIGHBORHOOD_GROUP_BY` | `region_label` | H5AD `.obs` label used by ST-06d Squidpy neighborhood. |
| `SPATIAL_REGULATION_H5AD_MODULE` | `spatial_03_region` | Spatial H5AD export module consumed by ST spatial decoupleR/pySCENIC wrappers. |
| `SPATIAL_SVG_H5AD_MODULE` | `spatial_03_region` | Spatial H5AD export module consumed by ST-09a SpatialDE2 sidecar. |
| `SPATIALDE2_ALLOW_PROXY` | `no` | Allows ST-09a to emit exploratory variance-proxy SVG output when SpatialDE2 is unavailable. |
| `SPARKX_ALLOW_PROXY` | `no` | Allows ST-09b to emit exploratory variance-proxy SVG output when SPARK-X is unavailable. |
| `SPARKX_SPATIAL_RDS` | spatial annotated checkpoint | Optional Seurat/RDS input for ST-09b SPARK-X. |
| `SPARKX_SPATIAL_H5AD` | empty | Optional H5AD adapter input for ST-09b proxy fallback. |
| `SVG_REQUIRE_REVIEW` | `no` | When `yes`, stage 60 holds on `spatial_svg` after 09a/09b/09c for manual review. |
| `SVG_LR_SUPPORT_TSV` | `${SPATIAL_TABLE_DIR}/09_svg/09c_consensus/svg_lr_support.tsv` | Optional SVG support table consumed by ST-08e communication consensus. |
| `SVG_ENRICHMENT_HANDOFF_TSV` | `${SPATIAL_TABLE_DIR}/09_svg/09c_consensus/svg_enrichment_handoff.tsv` | Optional SVG gene set handoff consumed by ST-06c enrichment EDA. |
| `SPATIAL_REGULATION_SVG_GENE_SETS` | `${SPATIAL_TABLE_DIR}/09_svg/09c_consensus/svg_gene_sets.tsv` | Optional SVG gene sets recorded by spatial regulation sidecars. |
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
| `MODULE_07A_CELLCHAT_VERSION` | `1.2` | Version marker for CellChat inventory-gated hypothesis-only outputs. |
| `MODULE_07B_LIANA_VERSION` | `1.0` | Version marker for the 07b LIANA+ consensus layer. |
| `MODULE_07C_NICHENET_VERSION` | `1.0` | Version marker for NicheNet/MultiNicheNet inventory-gated outputs. |
| `MODULE_07D_CONSENSUS_VERSION` | `1.0` | Version marker for the 07d communication consensus and evidence tier layer. |
| `MODULE_07E_EDA_VERSION` | `2.0` | Version marker for the 07e communication EDA layer. |
| `INVENTORY_GATE_REQUIRED` | `yes` | Require 04c `cluster_eligibility_tsv` before communication modules run. |
| `INVENTORY_GATE_ALLOWED_TIERS_COMMUNICATION` | `primary,exploratory,primary_merged` | Evidence tiers allowed to enter CellChat/NicheNet communication analysis. |
| `INVENTORY_GATE_STRICT_MODE` | `no` | When enabled, a cluster must pass all matching eligibility rows to enter communication analysis. |
| `COMMUNICATION_FAIL_ON_NO_PRIMARY` | `no` | Reserved strictness flag for failing communication when no cluster passes inventory gate. |
| `COMMUNICATION_REQUIRE_FULL_PIPELINE` | `yes` | When `yes`, 07b/07c/07d/07e failures block the communication stage instead of writing placeholder manifests. |
| `COMMUNICATION_REPORT_FILTER_TIER` | `auto` | Default 07e evidence-tier report filter. `auto` uses `communication_pairs.tsv:evidence_tier_required` when available. |
| `COMMUNICATION_REPORT_TOP_N_PRIMARY` | `30` | Maximum number of filtered axes displayed in the 07e primary-axis preview. |
| `COMMUNICATION_REPORT_HTML` | `yes` | When `yes`, 07e writes a simple HTML companion report next to `report.md`. |
| `CONSENSUS_MIN_METHODS_FOR_PRIMARY` | `2` | Minimum supporting method layers required for a primary 07d communication axis. |
| `CONSENSUS_REQUIRE_NICHENET_FOR_PRIMARY` | `yes` | Require downstream NicheNet/MultiNicheNet support before a 07d axis can be primary. |
| `LIANA_CONSENSUS_ENABLED` | `yes` | Enables the 07b LIANA+ consensus stage. |
| `LIANA_METHODS_LIST` | `cellphonedb,connectome,sca,natmi,logfc,rank_aggregate` | LIANA method set requested by 07b. |
| `LIANA_CONSENSUS_AGGREGATE` | `rank_aggregate` | Consensus score column preference for 07b LIANA outputs. |
| `LIANA_MIN_METHODS_AGREED_INSIDE` | `3` | Minimum method-hit count for `liana_consensus_hit=1`. |
| `LIANA_CELLPHONEDB_PVAL_THRESHOLD` | `0.05` | P-value cutoff for LIANA p-value method hits. |
| `LIANA_RESOURCE_DB` | `consensus` | LIANA ligand-receptor resource name. |
| `LIANA_MIN_CELLS_PER_GROUP` | `100` | Minimum cells per condition for a real LIANA run. |
| `LIANA_CONDITION_COL` | `condition` | H5AD `.obs` column used for per-condition LIANA runs. |
| `ORTHOLOG_CHICKEN_HUMAN_TSV` | `metadata/ortholog_chicken_human.tsv` | Optional chicken-to-human symbol LUT for communication `lr_axis_id` standardization. |
| `MULTINICHENET_ENABLED` | `auto` | Controls the 07c NicheNet/MultiNicheNet stage. |
| `NICHENET_MODE` | `auto` | Selects 07c mode: `auto`, `multinichenet`, `nichenet_legacy`, or `skip`. |
| `MULTINICHENET_MIN_SAMPLES_PER_GROUP` | `2` | Minimum samples per group before 07c can use MultiNicheNet mode. |
| `MULTINICHENET_TOP_N_LR` | `250` | Planned top ligand-receptor cap for MultiNicheNet output. |
| `MULTINICHENET_TOP_N_TARGETS` | `20` | Planned top ligand-target cap for receiver-DE overlap. |
| `MULTINICHENET_MIN_CELLS` | `10` | Planned minimum per-celltype cell count for MultiNicheNet. |
| `COMMUNICATION_RECEIVER_DE_P_THRESHOLD` | `0.05` | Receiver DEG adjusted-p cutoff for downstream target-hit evidence. |
| `COMMOT_ENABLED` | `auto` | Enables the ST 08 COMMOT spatial communication extension. |
| `COMMOT_DISTANCE_THRESHOLD` | `500` | Spatial distance cutoff passed to COMMOT. |
| `COMMOT_SIGNAL_THRESHOLD` | `0` | Minimum COMMOT signal score for `spatial_support=yes` in the 08b summary. |
| `COMMOT_PERMUTATIONS` | `100` | Number of spatial permutations requested for COMMOT cluster-level scoring. |
| `COMMOT_REQUIRE_RUNTIME` | `no` | When `yes`, missing Python COMMOT/anndata runtime is a hard failure instead of a schema-valid skip. |
| `COMMOT_LR_CANDIDATES_TSV` | `${TABLE_DIR}/communication/consensus/method_consensus.tsv` | Source LR axis table used by 08a to prepare COMMOT candidates. |
| `COMMOT_SPATIAL_SUMMARY_TSV` | `${RESULTS_DIR}/spatial/tables/08_commot/commot_spatial_summary.tsv` | 08b spatial support table consumed by 07d consensus. |
| `SPATIAL_VALIDATION_MODE` | `dirichlet_only` | 07f validation mode: `scdesign3`, `external_truth`, or smoke/fallback `dirichlet_only`. |
| `SPATIAL_VALIDATION_RMSE_PASS` / `SPATIAL_VALIDATION_RMSE_WARN` | `0.10` / `0.20` | 07f RMSE thresholds for deconvolution question gates. |
| `SPATIAL_VALIDATION_COR_PASS` / `SPATIAL_VALIDATION_COR_WARN` | `0.80` / `0.60` | 07f correlation thresholds for deconvolution question gates. |
| `SPATIAL_COMMUNICATION_TABLE_DIR` | `${SPATIAL_TABLE_DIR}/08_spatial_communication` | ST-08 integrated evidence-chain tables. |
| `SPATIAL_COMMUNICATION_REPORT_DIR` | `${EDA_REPORT_DIR}/spatial_08_communication` | ST-08 final spatial communication report directory. |
| `MODULE_03A3_VERSION` | `1.0` | Version marker for the independent panorama UMAP stage. |
| `MODULE_SPATIAL_02A2_VERSION` | `1.0` | Version marker for the independent spatial UMAP stage. |
| `UMAP_N_NEIGHBORS` | `30` | Shared UMAP neighbor count for 03a3 and spatial 02a2. |
| `UMAP_MIN_DIST` | `0.3` | Shared UMAP minimum distance. |
| `UMAP_SPREAD` | `1.0` | Shared UMAP spread. |
| `UMAP_SEED` | inherits `RANDOM_SEED` | Shared UMAP random seed. |
| `UMAP_METRIC` | `cosine` | Shared UMAP distance metric. |
| `UMAP_LOCAL_CONNECTIVITY` | `1` | Shared UMAP local connectivity parameter. |
