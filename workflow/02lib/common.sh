#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_FILE="${SCRNA_PIPELINE_CONFIG:-}"

if [[ -n "${CONFIG_FILE}" && ! -f "${CONFIG_FILE}" ]]; then
  echo "缺少配置文件: ${CONFIG_FILE}" >&2
  exit 1
fi

if [[ -n "${CONFIG_FILE}" ]]; then
  # shellcheck disable=SC1090
  source "${CONFIG_FILE}"
fi

choose_pipeline_locale() {
  local requested="${PIPELINE_LOCALE:-}"
  local candidate
  local available
  available="$(locale -a 2>/dev/null || true)"

  for candidate in "${requested}" en_US.utf8 zh_CN.utf8 en_US.UTF-8 zh_CN.UTF-8 C.utf8 C.UTF-8; do
    [[ -n "${candidate}" ]] || continue
    if grep -qi "^${candidate}$" <<< "${available}"; then
      echo "${candidate}"
      return 0
    fi
  done

  echo "[WARN] 未找到已安装的 UTF-8 locale，回退到 C.UTF-8；如 R 报 locale 警告，请安装 en_US.utf8 或 zh_CN.utf8。" >&2
  echo "C.UTF-8"
}

export PIPELINE_LOCALE="${PIPELINE_LOCALE:-$(choose_pipeline_locale)}"
export LANG="${PIPELINE_LOCALE}"
export LC_ALL="${PIPELINE_LOCALE}"
export LC_CTYPE="${PIPELINE_LOCALE}"
export PROJECT_ROOT="${PROJECT_ROOT:-${PIPELINE_ROOT}}"
export DATA_DIR="${DATA_DIR:-${PROJECT_ROOT}/data}"
export RESULTS_DIR="${RESULTS_DIR:-${PROJECT_ROOT}/results}"
export CHECKPOINT_DIR="${CHECKPOINT_DIR:-${RESULTS_DIR}/checkpoints}"
export FIGURE_DIR="${FIGURE_DIR:-${RESULTS_DIR}/figures}"
export TABLE_DIR="${TABLE_DIR:-${RESULTS_DIR}/tables}"
export ST_ENABLED="${ST_ENABLED:-yes}"
export SPATIAL_RESULTS_DIR="${SPATIAL_RESULTS_DIR:-${RESULTS_DIR}/spatial}"
export SPATIAL_CHECKPOINT_DIR="${SPATIAL_CHECKPOINT_DIR:-${SPATIAL_RESULTS_DIR}/checkpoints}"
export SPATIAL_FIGURE_DIR="${SPATIAL_FIGURE_DIR:-${SPATIAL_RESULTS_DIR}/figures}"
export SPATIAL_TABLE_DIR="${SPATIAL_TABLE_DIR:-${SPATIAL_RESULTS_DIR}/tables}"
export SPATIAL_REFERENCE_FREEZE_DIR="${SPATIAL_REFERENCE_FREEZE_DIR:-${RESULTS_DIR}/spatial_reference_frozen}"
export JOINT_RESULTS_DIR="${JOINT_RESULTS_DIR:-${RESULTS_DIR}/joint}"
export JOINT_CHECKPOINT_DIR="${JOINT_CHECKPOINT_DIR:-${JOINT_RESULTS_DIR}/checkpoints}"
export JOINT_FIGURE_DIR="${JOINT_FIGURE_DIR:-${JOINT_RESULTS_DIR}/figures}"
export JOINT_TABLE_DIR="${JOINT_TABLE_DIR:-${JOINT_RESULTS_DIR}/tables}"
export LOG_DIR="${LOG_DIR:-${PROJECT_ROOT}/logs}"
export ENV_DIR="${ENV_DIR:-${PROJECT_ROOT}/envs}"
export RESOURCE_DIR="${RESOURCE_DIR:-${PROJECT_ROOT}/resources}"
export MANIFEST_DIR="${MANIFEST_DIR:-${RESULTS_DIR}/manifests}"
export REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_ROOT}/reference}"
export REFERENCE_ROOT="${REFERENCE_ROOT:-$(dirname "${REFERENCE_DIR}")}"
export REFERENCE_VERSION="${REFERENCE_VERSION:-ensembl_release112}"
export FASTQ_DIR="${FASTQ_DIR:-${PROJECT_ROOT}/fastq}"
export GENOME_FASTA_GZ="${GENOME_FASTA_GZ:-${REFERENCE_DIR}/genome.fa}"
export REFERENCE_GTF="${REFERENCE_GTF:-${REFERENCE_DIR}/genes.gtf}"
export CLEAN_GTF="${CLEAN_GTF:-${REFERENCE_DIR}/genes.clean.gtf}"
export PAPER_MATRIX_SOURCE="${PAPER_MATRIX_SOURCE:-}"
export CELLRANGER_OUT_DIR="${CELLRANGER_OUT_DIR:-}"
export SRA_TOOLS_DIR="${SRA_TOOLS_DIR:-}"
export MAIN_THREADS="${MAIN_THREADS:-8}"
export SAMPLE_NAMES="${SAMPLE_NAMES:-}"
export RAW_SAMPLES="${RAW_SAMPLES:-}"
export RAW_GROUP_1_NAME="${RAW_GROUP_1_NAME:-Group1}"
export RAW_GROUP_1_SAMPLES="${RAW_GROUP_1_SAMPLES:-}"
export RAW_GROUP_2_NAME="${RAW_GROUP_2_NAME:-Group2}"
export RAW_GROUP_2_SAMPLES="${RAW_GROUP_2_SAMPLES:-}"
export ANALYSIS_GROUP_1_NAME="${ANALYSIS_GROUP_1_NAME:-}"
export ANALYSIS_GROUP_1_SAMPLES="${ANALYSIS_GROUP_1_SAMPLES:-}"
export ANALYSIS_GROUP_2_NAME="${ANALYSIS_GROUP_2_NAME:-}"
export ANALYSIS_GROUP_2_SAMPLES="${ANALYSIS_GROUP_2_SAMPLES:-}"
export DEG_IDENT_1="${DEG_IDENT_1:-}"
export DEG_IDENT_2="${DEG_IDENT_2:-}"
export RANDOM_SEED="${RANDOM_SEED:-42}"
export QC_MIN_NFEATURE="${QC_MIN_NFEATURE:-200}"
export QC_MIN_NCOUNT="${QC_MIN_NCOUNT:-1000}"
export QC_MIN_LOG10UMI="${QC_MIN_LOG10UMI:-0.7}"
export QC_MAX_MITO_PCT="${QC_MAX_MITO_PCT:-20}"
export MODULE_02_VERSION="${MODULE_02_VERSION:-1.0}"
export NORMALIZATION_METHODS="${NORMALIZATION_METHODS:-lognorm}"
export INTEGRATION_MODES="${INTEGRATION_MODES:-${INTEGRATION_MODE:-harmony}}"
export VARS_TO_REGRESS_DEFAULT="${VARS_TO_REGRESS_DEFAULT:-}"
export RES_FINE_STEP="${RES_FINE_STEP:-0.005}"
export PCA_DIMS_PANORAMA="${PCA_DIMS_PANORAMA:-${PCA_DIMS:-1:30}}"
export PCA_DIMS_SUBCLUSTER="${PCA_DIMS_SUBCLUSTER:-1:20}"
export SCDESIGN3_N_SIM="${SCDESIGN3_N_SIM:-5}"
export SCDESIGN3_FAMILY="${SCDESIGN3_FAMILY:-nb}"
export SCDESIGN3_N_CORES="${SCDESIGN3_N_CORES:-${MAIN_THREADS:-4}}"
export SCDESIGN3_MAX_CELLS_PER_LABEL="${SCDESIGN3_MAX_CELLS_PER_LABEL:-2000}"
export SCDESIGN3_N_HVG="${SCDESIGN3_N_HVG:-2000}"
export SCDESIGN3_N_PCS="${SCDESIGN3_N_PCS:-30}"
export PANORAMA_LAYER_ID="${PANORAMA_LAYER_ID:-panorama}"
export ANNOTATION_HUB_PATH_CLUSTERED="${ANNOTATION_HUB_PATH_CLUSTERED:-${CHECKPOINT_DIR}/02_after_clustering.rds}"
export MODULE_03_VERSION="${MODULE_03_VERSION:-1.0}"
export MODULE_04_VERSION="${MODULE_04_VERSION:-1.0}"
export MODULE_05_VERSION="${MODULE_05_VERSION:-1.0}"
export MAX_INTEGRATION_CANDIDATES_PER_LAYER="${MAX_INTEGRATION_CANDIDATES_PER_LAYER:-4}"
export USE_EXISTING_SIF="${USE_EXISTING_SIF:-no}"
export EXISTING_R_SIF="${EXISTING_R_SIF:-}"
export R_LIBS_MAIN="${R_LIBS_MAIN:-}"
export R_LIBS_SCENIC="${R_LIBS_SCENIC:-}"
export R_LIBS_INTERACTION="${R_LIBS_INTERACTION:-}"
export R_LIBS_DECOUPLER="${R_LIBS_DECOUPLER:-}"
export R_MAIN_ENV_PREFIX="${R_MAIN_ENV_PREFIX:-${ENV_DIR}/conda/r_main}"
export R_SCENIC_ENV_PREFIX="${R_SCENIC_ENV_PREFIX:-${ENV_DIR}/conda/r_scenic}"
export R_INTERACTION_ENV_PREFIX="${R_INTERACTION_ENV_PREFIX:-${ENV_DIR}/conda/r_interaction}"
export R_DECOUPLER_ENV_PREFIX="${R_DECOUPLER_ENV_PREFIX:-${ENV_DIR}/conda/r_decoupler}"
export R_LEGACY_ENV_PREFIX="${R_LEGACY_ENV_PREFIX:-${ENV_DIR}/conda/r_legacy}"
export R_SPATIAL_ENV_PREFIX="${R_SPATIAL_ENV_PREFIX:-${ENV_DIR}/conda/r_spatial}"
export PY_SPATIAL_ENV_PREFIX="${PY_SPATIAL_ENV_PREFIX:-${ENV_DIR}/conda/py_spatial}"
export PY_SPATIAL_LEGACY_ENV_PREFIX="${PY_SPATIAL_LEGACY_ENV_PREFIX:-${ENV_DIR}/conda/py_spatial_legacy}"
export PY_CELL2LOCATION_ENV_PREFIX="${PY_CELL2LOCATION_ENV_PREFIX:-${ENV_DIR}/conda/py_cell2location}"
export PY_SAW_ENV_PREFIX="${PY_SAW_ENV_PREFIX:-${ENV_DIR}/conda/py_saw}"
export R_VALIDATION_ENV_PREFIX="${R_VALIDATION_ENV_PREFIX:-${ENV_DIR}/conda/r_validation}"
export PYSCENIC_ENV_PREFIX="${PYSCENIC_ENV_PREFIX:-${ENV_DIR}/conda/pyscenic}"
export VELOCITY_ENV_PREFIX="${VELOCITY_ENV_PREFIX:-${ENV_DIR}/conda/velocity}"
export SCVELO_ENV_PREFIX="${SCVELO_ENV_PREFIX:-${ENV_DIR}/conda/scvelo}"
export VELOCITY_DIR="${VELOCITY_DIR:-${RESULTS_DIR}/velocity}"
export VELOCITY_INPUT_DIR="${VELOCITY_INPUT_DIR:-${VELOCITY_DIR}/input}"
export VELOCITY_LOOM_DIR="${VELOCITY_LOOM_DIR:-${VELOCITY_DIR}/loom}"
export VELOCITY_OUTPUT_DIR="${VELOCITY_OUTPUT_DIR:-${VELOCITY_DIR}/output}"
export SCVELO_THREADS="${SCVELO_THREADS:-${MAIN_THREADS:-8}}"
export MODULE_10_VERSION="${MODULE_10_VERSION:-1.0}"
export VELOCITY_SAMPLE_IDS="${VELOCITY_SAMPLE_IDS:-}"
export VELOCITY_BAM_PATTERN="${VELOCITY_BAM_PATTERN:-anno_decon_sorted.bam}"
export VELOCITY_H5_PATTERN="${VELOCITY_H5_PATTERN:-filtered_feature_bc_matrix.h5}"
export VELOCITY_GTF="${VELOCITY_GTF:-${CLEAN_GTF}}"
export VELOCYTO_REPEAT_MASK_GTF="${VELOCYTO_REPEAT_MASK_GTF:-}"
export VELOCYTO_THREADS="${VELOCYTO_THREADS:-${MAIN_THREADS:-8}}"
export VELOCITY_MIN_REFERENCE_CELLS="${VELOCITY_MIN_REFERENCE_CELLS:-10}"
export SCVELO_MIN_SHARED_COUNTS="${SCVELO_MIN_SHARED_COUNTS:-20}"
export SCVELO_TOP_GENES="${SCVELO_TOP_GENES:-2000}"
export SCVELO_N_PCS="${SCVELO_N_PCS:-30}"
export SCVELO_N_NEIGHBORS="${SCVELO_N_NEIGHBORS:-30}"
export VELOCITY_DRIVER_TOP_N="${VELOCITY_DRIVER_TOP_N:-200}"
export CELLRANK_MIN_CELLS="${CELLRANK_MIN_CELLS:-200}"
export CELLRANK_MIN_VELOCITY_CONFIDENCE="${CELLRANK_MIN_VELOCITY_CONFIDENCE:-0.05}"
export CELLRANK_N_STATES="${CELLRANK_N_STATES:-6}"
export SCENIC_DB_DIR="${SCENIC_DB_DIR:-${RESOURCE_DIR}/scenic_db}"
export SCENIC_INPUT_DIR="${SCENIC_INPUT_DIR:-${RESULTS_DIR}/scenic_input}"
export SCENIC_OUTPUT_DIR="${SCENIC_OUTPUT_DIR:-${RESULTS_DIR}/scenic_output}"
export SCENIC_TF_LIST="${SCENIC_TF_LIST:-${SCENIC_DB_DIR}/hs_hgnc_tfs.txt}"
export SCENIC_MOTIF_ANN="${SCENIC_MOTIF_ANN:-${SCENIC_DB_DIR}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl}"
export SCENIC_DB_500BP="${SCENIC_DB_500BP:-${SCENIC_DB_DIR}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather}"
export SCENIC_DB_10KB="${SCENIC_DB_10KB:-${SCENIC_DB_DIR}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather}"
export SCENIC_RESOURCE_MANIFEST="${SCENIC_RESOURCE_MANIFEST:-${SCENIC_DB_DIR}/scenic_resources_manifest.json}"
export SCENIC_ORTHOLOG_MAP_FILE="${SCENIC_ORTHOLOG_MAP_FILE:-}"
export DECOUPLER_RESOURCE_DIR="${DECOUPLER_RESOURCE_DIR:-${RESOURCE_DIR}/decoupler}"
export ENSEMBL_MIRROR="${ENSEMBL_MIRROR:-asia}"
export ORTHOLOG_CACHE_DIR="${ORTHOLOG_CACHE_DIR:-${RESULTS_DIR}/ortholog_cache}"
export ORTHOLOG_MANIFEST="${ORTHOLOG_MANIFEST:-${ORTHOLOG_CACHE_DIR}/_manifest.json}"
export CELL_CYCLE_GENES_RDS="${CELL_CYCLE_GENES_RDS:-${ORTHOLOG_CACHE_DIR}/chicken_cc_genes.rds}"
export PROJECT_CONFIG_DIR="${PROJECT_CONFIG_DIR:-${PROJECT_ROOT}/config}"
export RAW_DIR="${RAW_DIR:-${PROJECT_ROOT}/raw}"
export RAW_RECEIVED_DIR="${RAW_RECEIVED_DIR:-${RAW_DIR}/received}"
export METADATA_DIR="${METADATA_DIR:-${PROJECT_ROOT}/metadata}"
export REPORT_DIR="${REPORT_DIR:-${PROJECT_ROOT}/reports}"
export INTAKE_REPORT_DIR="${INTAKE_REPORT_DIR:-${REPORT_DIR}/intake}"
export EDA_REPORT_DIR="${EDA_REPORT_DIR:-${REPORT_DIR}/eda}"
export STATUS_DIR="${STATUS_DIR:-${PROJECT_ROOT}/status}"
export WORKFLOW_STATUS_FILE="${WORKFLOW_STATUS_FILE:-${STATUS_DIR}/workflow_status.json}"
export SAMPLE_SHEET="${SAMPLE_SHEET:-${METADATA_DIR}/samples.tsv}"
export CANONICAL_SAMPLE_SHEET="${CANONICAL_SAMPLE_SHEET:-${METADATA_DIR}/samples.canonical.tsv}"
export COMPARISON_SHEET="${COMPARISON_SHEET:-${METADATA_DIR}/comparisons.tsv}"
export SECTION_SHEET="${SECTION_SHEET:-${METADATA_DIR}/sections.tsv}"
export SPATIAL_REFERENCE_INVENTORY_FILE="${SPATIAL_REFERENCE_INVENTORY_FILE:-${METADATA_DIR}/spatial_reference_inventory.tsv}"
export COMMUNICATION_PAIRS_SHEET="${COMMUNICATION_PAIRS_SHEET:-${METADATA_DIR}/communication_pairs.tsv}"
export TRAJECTORY_PAIRS_SHEET="${TRAJECTORY_PAIRS_SHEET:-${METADATA_DIR}/trajectory_pairs.tsv}"
export SCENIC_TARGETS_SHEET="${SCENIC_TARGETS_SHEET:-${METADATA_DIR}/scenic_targets.tsv}"
export ENRICHMENT_TARGETS_SHEET="${ENRICHMENT_TARGETS_SHEET:-${METADATA_DIR}/enrichment_targets.tsv}"
export GENE_PROGRAM_TARGETS_SHEET="${GENE_PROGRAM_TARGETS_SHEET:-${METADATA_DIR}/gene_program_targets.tsv}"
export DECONV_PAIRS_SHEET="${DECONV_PAIRS_SHEET:-${METADATA_DIR}/deconv_pairs.tsv}"
export SPATIAL_PAIRS_SHEET="${SPATIAL_PAIRS_SHEET:-${METADATA_DIR}/spatial_pairs.tsv}"
export SCDESIGN3_QUESTION_MAP="${SCDESIGN3_QUESTION_MAP:-${METADATA_DIR}/scdesign3_question_map.tsv}"
export SCDESIGN3_TARGETS_SHEET="${SCDESIGN3_TARGETS_SHEET:-${METADATA_DIR}/scdesign3_targets.tsv}"
export SCDESIGN3_SIMULATION_DESIGNS_SHEET="${SCDESIGN3_SIMULATION_DESIGNS_SHEET:-${METADATA_DIR}/scdesign3_simulation_designs.tsv}"
export SCDESIGN3_THRESHOLDS_SHEET="${SCDESIGN3_THRESHOLDS_SHEET:-${METADATA_DIR}/scdesign3_thresholds.tsv}"
export DELIVERY_MANIFEST="${DELIVERY_MANIFEST:-${METADATA_DIR}/delivery_manifest.tsv}"
export RECEIVED_FILES_MANIFEST="${RECEIVED_FILES_MANIFEST:-${METADATA_DIR}/received_files_manifest.tsv}"
export INPUT_INVENTORY_FILE="${INPUT_INVENTORY_FILE:-${INTAKE_REPORT_DIR}/input_inventory.tsv}"
export SPATIAL_INPUT_INVENTORY_FILE="${SPATIAL_INPUT_INVENTORY_FILE:-${INTAKE_REPORT_DIR}/spatial_input_inventory.tsv}"
export BRANCH_READINESS_FILE="${BRANCH_READINESS_FILE:-${INTAKE_REPORT_DIR}/branch_readiness.tsv}"
export INTAKE_SUMMARY_FILE="${INTAKE_SUMMARY_FILE:-${INTAKE_REPORT_DIR}/intake_summary.md}"
export WAIVER_FILE="${WAIVER_FILE:-${PROJECT_CONFIG_DIR}/waivers.tsv}"
export QC_THRESHOLD_FILE="${QC_THRESHOLD_FILE:-${PROJECT_CONFIG_DIR}/qc_thresholds.tsv}"
export SPATIAL_QC_THRESHOLD_FILE="${SPATIAL_QC_THRESHOLD_FILE:-${PROJECT_CONFIG_DIR}/spatial_qc_thresholds.tsv}"
export EDA_GATE_FILE="${EDA_GATE_FILE:-${PROJECT_CONFIG_DIR}/eda_gates.tsv}"
export OBJECT_LAYER_CONFIG_FILE="${OBJECT_LAYER_CONFIG_FILE:-${PROJECT_CONFIG_DIR}/object_layers.tsv}"
export SPATIAL_OBJECT_LAYER_FILE="${SPATIAL_OBJECT_LAYER_FILE:-${PROJECT_CONFIG_DIR}/spatial_object_layers.tsv}"
export SPATIAL_INTAKE_CONTRACT_FILE="${SPATIAL_INTAKE_CONTRACT_FILE:-${PROJECT_CONFIG_DIR}/spatial_intake_contract.tsv}"
export MARKER_PANEL_DIR="${MARKER_PANEL_DIR:-${PROJECT_CONFIG_DIR}/marker_panels}"
export MITO_GENE_LIST_FILE="${MITO_GENE_LIST_FILE:-${PROJECT_CONFIG_DIR}/mito_gene_list.txt}"
export AMBIENT_REPORT_DIR="${AMBIENT_REPORT_DIR:-${EDA_REPORT_DIR}/ambient}"
export PRE_QC_REPORT_DIR="${PRE_QC_REPORT_DIR:-${EDA_REPORT_DIR}/pre_qc}"
export POST_QC_REPORT_DIR="${POST_QC_REPORT_DIR:-${EDA_REPORT_DIR}/post_qc}"
export INTEGRATION_REPORT_DIR="${INTEGRATION_REPORT_DIR:-${EDA_REPORT_DIR}/integration}"
export ANNOTATION_REPORT_DIR="${ANNOTATION_REPORT_DIR:-${EDA_REPORT_DIR}/annotation}"
export SUBCLUSTER_REPORT_DIR="${SUBCLUSTER_REPORT_DIR:-${EDA_REPORT_DIR}/subcluster}"
export DEG_REPORT_DIR="${DEG_REPORT_DIR:-${EDA_REPORT_DIR}/deg}"
export ENRICHMENT_REPORT_DIR="${ENRICHMENT_REPORT_DIR:-${EDA_REPORT_DIR}/enrichment}"
export COMMUNICATION_REPORT_DIR="${COMMUNICATION_REPORT_DIR:-${EDA_REPORT_DIR}/communication}"
export REGULATION_REPORT_DIR="${REGULATION_REPORT_DIR:-${EDA_REPORT_DIR}/regulation}"
export SELECTED_INTEGRATION_FILE="${SELECTED_INTEGRATION_FILE:-${INTEGRATION_REPORT_DIR}/panorama/selected_integration.txt}"
export LAYER_STATUS_FILE="${LAYER_STATUS_FILE:-${TABLE_DIR}/layer_status.tsv}"
export SUBCLUSTER_REVIEW_SUMMARY_FILE="${SUBCLUSTER_REVIEW_SUMMARY_FILE:-${TABLE_DIR}/subcluster/subcluster_review_summary.tsv}"
export SUBCLUSTER_CANDIDATE_LAYERS_COUNT_FILE="${SUBCLUSTER_CANDIDATE_LAYERS_COUNT_FILE:-${TABLE_DIR}/subcluster/candidate_layers_count.txt}"
export DNBC4TOOLS_OUT_DIR="${DNBC4TOOLS_OUT_DIR:-${DATA_DIR}/dnbc4tools_out}"
export STAR_INDEX_DIR="${STAR_INDEX_DIR:-${REFERENCE_DIR}/star_index}"
export STAR_THREADS="${STAR_THREADS:-${MAIN_THREADS:-8}}"
export DNBC4TOOLS_THREADS="${DNBC4TOOLS_THREADS:-${MAIN_THREADS:-8}}"
export STAR_SA_INDEX_NBASES="${STAR_SA_INDEX_NBASES:-13}"
export INPUT_STANDARDIZE_MODE="${INPUT_STANDARDIZE_MODE:-symlink}"
export ANNOTATION_HUB_PATH="${ANNOTATION_HUB_PATH:-${CHECKPOINT_DIR}/03_after_annotation.rds}"
export MIN_BIOLOGICAL_REPLICATES="${MIN_BIOLOGICAL_REPLICATES:-2}"
export DEG_MIN_CELLS_PER_GROUP="${DEG_MIN_CELLS_PER_GROUP:-3}"
export DEG_LOGFC_THRESHOLD="${DEG_LOGFC_THRESHOLD:-0}"
export DEG_ALPHA="${DEG_ALPHA:-0.05}"
export MODULE_06_VERSION="${MODULE_06_VERSION:-1.0}"
export ENRICHMENT_SPECIES_STRATEGY="${ENRICHMENT_SPECIES_STRATEGY:-chicken_primary}"
export ENRICHMENT_PVALUE_CUTOFF="${ENRICHMENT_PVALUE_CUTOFF:-0.05}"
export ENRICHMENT_QVALUE_CUTOFF="${ENRICHMENT_QVALUE_CUTOFF:-0.2}"
export ENRICHMENT_TOP_N_GENES="${ENRICHMENT_TOP_N_GENES:-100}"
export ENRICHMENT_MIN_INPUT_GENES="${ENRICHMENT_MIN_INPUT_GENES:-5}"
export ENRICHMENT_MIN_GS_SIZE="${ENRICHMENT_MIN_GS_SIZE:-10}"
export ENRICHMENT_MAX_GS_SIZE="${ENRICHMENT_MAX_GS_SIZE:-500}"
export ENRICHMENT_KEGG_TIMEOUT_SEC="${ENRICHMENT_KEGG_TIMEOUT_SEC:-60}"
export MODULE_07_VERSION="${MODULE_07_VERSION:-1.0}"
export COMMUNICATION_CELL_TYPE_COL="${COMMUNICATION_CELL_TYPE_COL:-cell_subtype}"
export COMMUNICATION_MIN_CELLS_PER_CELLTYPE="${COMMUNICATION_MIN_CELLS_PER_CELLTYPE:-20}"
export COMMUNICATION_ORTHOLOG_MIN_COVERAGE="${COMMUNICATION_ORTHOLOG_MIN_COVERAGE:-0.30}"
export CELLCHAT_MIN_CELLS_PER_GROUP="${CELLCHAT_MIN_CELLS_PER_GROUP:-20}"
export CELLCHAT_PROB_CUTOFF="${CELLCHAT_PROB_CUTOFF:-0.05}"
export NICHENET_RESOURCE_DIR="${NICHENET_RESOURCE_DIR:-${RESOURCE_DIR}/nichenet}"
export NICHENET_EXPRESSION_PCT="${NICHENET_EXPRESSION_PCT:-0.10}"
export NICHENET_TOP_LIGAND_N="${NICHENET_TOP_LIGAND_N:-20}"
export NICHENET_TOP_TARGET_N="${NICHENET_TOP_TARGET_N:-200}"
export MODULE_08_VERSION="${MODULE_08_VERSION:-1.0}"
export REGULATION_LAYERS="${REGULATION_LAYERS:-${PANORAMA_LAYER_ID:-panorama}}"
export SCENIC_MODULE_TOP_N="${SCENIC_MODULE_TOP_N:-50}"
export SCENIC_MODULE_MIN_GENES="${SCENIC_MODULE_MIN_GENES:-10}"
export SCENIC_REGULON_MIN_TARGETS="${SCENIC_REGULON_MIN_TARGETS:-5}"
export SCENIC_NES_THRESHOLD="${SCENIC_NES_THRESHOLD:-3}"
export SCENIC_AUC_MAX_RANK_FRACTION="${SCENIC_AUC_MAX_RANK_FRACTION:-0.05}"
export SCENIC_THREADS="${SCENIC_THREADS:-8}"
export SCENIC_ORTHOLOG_MIN_COVERAGE="${SCENIC_ORTHOLOG_MIN_COVERAGE:-0.30}"
export DECOUPLER_TF_METHOD="${DECOUPLER_TF_METHOD:-wmean}"
export DECOUPLER_PATHWAY_METHOD="${DECOUPLER_PATHWAY_METHOD:-ulm}"
export DECOUPLER_TOP_TF_N="${DECOUPLER_TOP_TF_N:-25}"
export DECOUPLER_MIN_TARGETS="${DECOUPLER_MIN_TARGETS:-5}"
export DECOUPLER_CONFIDENCE_LEVELS="${DECOUPLER_CONFIDENCE_LEVELS:-A,B,C}"
export DECOUPLER_USE_CACHE="${DECOUPLER_USE_CACHE:-yes}"
export DECOUPLER_ACTIVITY_LEVEL="${DECOUPLER_ACTIVITY_LEVEL:-group_average}"
export DECOUPLER_GROUP_COL="${DECOUPLER_GROUP_COL:-auto}"
export ALLOW_REGULATION_EMPTY_REPORT="${ALLOW_REGULATION_EMPTY_REPORT:-no}"
export TRIAGE_FRAC_BELOW_CUTOFF="${TRIAGE_FRAC_BELOW_CUTOFF:-0.35}"
export TRIAGE_FRAC_ABOVE_MITO="${TRIAGE_FRAC_ABOVE_MITO:-0.25}"
export TRIAGE_DENSITY_PEAKS="${TRIAGE_DENSITY_PEAKS:-2}"
export AMBIENT_PRIMARY_METHOD="${AMBIENT_PRIMARY_METHOD:-soupx}"
export AMBIENT_FALLBACK_METHOD="${AMBIENT_FALLBACK_METHOD:-decontx}"
export AMBIENT_APPLY_POLICY="${AMBIENT_APPLY_POLICY:-manual}"
export AMBIENT_MIN_CELLS="${AMBIENT_MIN_CELLS:-50}"
export AMBIENT_CLUSTER_DIMS="${AMBIENT_CLUSTER_DIMS:-1:20}"
export AMBIENT_CLUSTER_RESOLUTION="${AMBIENT_CLUSTER_RESOLUTION:-0.4}"
export AMBIENT_MARKER_TOP_N="${AMBIENT_MARKER_TOP_N:-3}"
export AMBIENT_RECOMMEND_MIN_CONTAMINATION="${AMBIENT_RECOMMEND_MIN_CONTAMINATION:-0.05}"
export CELLBENDER_MODE="${CELLBENDER_MODE:-stub}"
export CELLBENDER_FPR="${CELLBENDER_FPR:-0.01}"
export CELLBENDER_CUDA="${CELLBENDER_CUDA:-yes}"
export CELLBENDER_EXTRA_ARGS="${CELLBENDER_EXTRA_ARGS:-}"
export DOUBLET_RATE_PER_1K="${DOUBLET_RATE_PER_1K:-${DOUBLET_RATE:-0.008}}"
export DOUBLET_RATE="${DOUBLET_RATE:-0.008}"
export DOUBLET_PRIMARY_CALLER="${DOUBLET_PRIMARY_CALLER:-scDblFinder}"
export DOUBLET_SECONDARY_CALLER="${DOUBLET_SECONDARY_CALLER:-DoubletFinder}"
export DOUBLET_SECONDARY_ENABLED="${DOUBLET_SECONDARY_ENABLED:-yes}"
export DOUBLET_MIN_CELLS="${DOUBLET_MIN_CELLS:-50}"
export DOUBLET_DIMS="${DOUBLET_DIMS:-1:20}"
export HVG_NFEATURES="${HVG_NFEATURES:-2000}"
export PCA_DIMS="${PCA_DIMS:-1:30}"
export TARGET_CLUSTERS="${TARGET_CLUSTERS:-15}"
export RES_RANGE="${RES_RANGE:-0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60}"
export TRAJECTORY_START="${TRAJECTORY_START:-}"
export TRADESEQ_KNOTS="${TRADESEQ_KNOTS:-6}"
export MODULE_09_VERSION="${MODULE_09_VERSION:-1.0}"
export TRAJECTORY_DIR="${TRAJECTORY_DIR:-${RESULTS_DIR}/trajectory}"
export TRAJECTORY_REPORT_DIR="${TRAJECTORY_REPORT_DIR:-${EDA_REPORT_DIR}/trajectory}"
export VELOCITY_REPORT_DIR="${VELOCITY_REPORT_DIR:-${EDA_REPORT_DIR}/velocity}"
export SPATIAL_PRE_QC_REPORT_DIR="${SPATIAL_PRE_QC_REPORT_DIR:-${EDA_REPORT_DIR}/spatial_pre_qc}"
export SPATIAL_POST_QC_REPORT_DIR="${SPATIAL_POST_QC_REPORT_DIR:-${EDA_REPORT_DIR}/spatial_post_qc}"
export SPATIAL_INTEGRATION_REPORT_DIR="${SPATIAL_INTEGRATION_REPORT_DIR:-${EDA_REPORT_DIR}/spatial_integration}"
export SPATIAL_REGION_ANNOTATION_REPORT_DIR="${SPATIAL_REGION_ANNOTATION_REPORT_DIR:-${EDA_REPORT_DIR}/spatial_region_annotation}"
export SPATIAL_DECONV_REPORT_DIR="${SPATIAL_DECONV_REPORT_DIR:-${EDA_REPORT_DIR}/spatial_deconv}"
export JOINT_TRAJECTORY_REPORT_DIR="${JOINT_TRAJECTORY_REPORT_DIR:-${EDA_REPORT_DIR}/joint_trajectory}"
export JOINT_COMMUNICATION_REPORT_DIR="${JOINT_COMMUNICATION_REPORT_DIR:-${EDA_REPORT_DIR}/joint_communication}"
export SPATIAL_QC_MIN_NFEATURE="${SPATIAL_QC_MIN_NFEATURE:-200}"
export SPATIAL_QC_MIN_NCOUNT="${SPATIAL_QC_MIN_NCOUNT:-500}"
export SPATIAL_QC_MAX_MITO_PCT="${SPATIAL_QC_MAX_MITO_PCT:-20}"
export SPATIAL_SCT_VST_FLAVOR="${SPATIAL_SCT_VST_FLAVOR:-v2}"
export SPATIAL_DEFAULT_CLUSTER_RESOLUTION="${SPATIAL_DEFAULT_CLUSTER_RESOLUTION:-0.6}"
export SPATIAL_CLUSTER_RESOLUTIONS="${SPATIAL_CLUSTER_RESOLUTIONS:-0.4,0.6,0.8}"
export SPATIAL_INTEGRATION_MODE="${SPATIAL_INTEGRATION_MODE:-none}"
export SPATIAL_DECONV_PRIMARY="${SPATIAL_DECONV_PRIMARY:-rctd}"
export SPATIAL_DECONV_FALLBACK="${SPATIAL_DECONV_FALLBACK:-seurat_transfer,card,cell2location}"
export SPATIAL_SVG_METHOD="${SPATIAL_SVG_METHOD:-spatialde}"
export TRAJECTORY_HVG_NFEATURES="${TRAJECTORY_HVG_NFEATURES:-${HVG_NFEATURES:-2000}}"
export TRAJECTORY_PCA_DIMS="${TRAJECTORY_PCA_DIMS:-1:30}"
export TRAJECTORY_UMAP_N_NEIGHBORS="${TRAJECTORY_UMAP_N_NEIGHBORS:-30}"
export TRAJECTORY_SPLIT_MIN_CELLS="${TRAJECTORY_SPLIT_MIN_CELLS:-50}"
export TRAJECTORY_BALANCE_WARN_FRACTION="${TRAJECTORY_BALANCE_WARN_FRACTION:-0.30}"
export TRAJECTORY_CONSENSUS_MIN_METHODS="${TRAJECTORY_CONSENSUS_MIN_METHODS:-2}"
export TRAJECTORY_CONSENSUS_CONFLICT_RHO="${TRAJECTORY_CONSENSUS_CONFLICT_RHO:-0.30}"

export PATH="${SRA_TOOLS_DIR:+${SRA_TOOLS_DIR}:}${PATH}"

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

warn() {
  echo "[WARN] $*" >&2
}

ensure_dir() {
  local path
  for path in "$@"; do
    [[ -n "${path}" ]] || continue
    mkdir -p "${path}"
  done
}

ensure_eda_control_files() {
  ensure_dir \
    "${PROJECT_CONFIG_DIR}" \
    "${MARKER_PANEL_DIR}" \
    "${AMBIENT_REPORT_DIR}" \
    "${EDA_REPORT_DIR}" \
    "${PRE_QC_REPORT_DIR}" \
    "${POST_QC_REPORT_DIR}" \
    "${INTEGRATION_REPORT_DIR}" \
    "${ANNOTATION_REPORT_DIR}" \
    "${SUBCLUSTER_REPORT_DIR}" \
    "${DEG_REPORT_DIR}" \
    "${ENRICHMENT_REPORT_DIR}" \
    "${COMMUNICATION_REPORT_DIR}" \
    "${REGULATION_REPORT_DIR}" \
    "${TRAJECTORY_REPORT_DIR}" \
    "${VELOCITY_REPORT_DIR}" \
    "${SPATIAL_PRE_QC_REPORT_DIR}" \
    "${SPATIAL_POST_QC_REPORT_DIR}" \
    "${SPATIAL_INTEGRATION_REPORT_DIR}" \
    "${SPATIAL_REGION_ANNOTATION_REPORT_DIR}" \
    "${SPATIAL_DECONV_REPORT_DIR}" \
    "${JOINT_TRAJECTORY_REPORT_DIR}" \
    "${JOINT_COMMUNICATION_REPORT_DIR}"

  if [[ ! -s "${EDA_GATE_FILE}" ]]; then
    {
      printf 'gate_id\tstatus\tapproved_by\tapproved_at\tgate_profile\tnotes\n'
      printf 'ambient\tpending\t\t\tdefault\t\n'
      printf 'pre_qc\tpending\t\t\tdefault\t\n'
      printf 'post_qc\tpending\t\t\tdefault\t\n'
      printf 'integration\tpending\t\t\tdefault\t\n'
      printf 'annotation\tpending\t\t\tdefault\t\n'
      printf 'subcluster\tpending\t\t\tdefault\t\n'
      printf 'deg\tpending\t\t\tdefault\t\n'
      printf 'enrichment\tpending\t\t\tdefault\t\n'
      printf 'communication\tpending\t\t\tdefault\t\n'
      printf 'regulation\tpending\t\t\tdefault\t\n'
      printf 'trajectory_inputs\tpending\t\t\tdefault\t\n'
      printf 'trajectory_methods\tpending\t\t\tdefault\t\n'
      printf 'trajectory_finalize\tpending\t\t\tdefault\t\n'
      printf 'velocity_inputs\tpending\t\t\tdefault\t\n'
      printf 'velocity_finalize\tpending\t\t\tdefault\t\n'
      printf 'scdesign3_targets\tpending\t\t\tdefault\t\n'
      printf 'scdesign3_validated\tpending\t\t\tdefault\t\n'
      printf 'spatial_pre_qc\tpending\t\t\tdefault\t\n'
      printf 'spatial_post_qc\tpending\t\t\tdefault\t\n'
      printf 'spatial_integration\tpending\t\t\tdefault\t\n'
      printf 'spatial_region_annotation\tpending\t\t\tdefault\t\n'
      printf 'spatial_deconv\tpending\t\t\tdefault\t\n'
      printf 'joint_trajectory\tpending\t\t\tdefault\t\n'
      printf 'joint_communication\tpending\t\t\tdefault\t\n'
    } > "${EDA_GATE_FILE}"
  fi

  local gate_id
  for gate_id in ambient pre_qc post_qc integration annotation subcluster deg enrichment communication regulation trajectory_inputs trajectory_methods trajectory_finalize velocity_inputs velocity_finalize scdesign3_targets scdesign3_validated spatial_pre_qc spatial_post_qc spatial_integration spatial_region_annotation spatial_deconv joint_trajectory joint_communication; do
    if ! awk -F '\t' -v gate="${gate_id}" 'NR > 1 && $1 == gate { found = 1 } END { exit(found ? 0 : 1) }' "${EDA_GATE_FILE}" >/dev/null 2>&1; then
      printf '%s\tpending\t\t\tdefault\t\n' "${gate_id}" >> "${EDA_GATE_FILE}"
    fi
  done

  if [[ ! -s "${QC_THRESHOLD_FILE}" ]]; then
    cat > "${QC_THRESHOLD_FILE}" <<'EOF'
sample_id	qc_min_nfeature	qc_min_ncount	qc_min_log10umi	qc_max_mito_pct
EOF
  fi
}

ensure_object_layer_config_file() {
  ensure_dir "${PROJECT_CONFIG_DIR}"
  [[ -s "${OBJECT_LAYER_CONFIG_FILE}" ]] && return 0

  {
    printf '# object_layers.tsv\n'
    printf '# - panorama 行是根层，默认从全部 post-QC 细胞建模\n'
    printf '# - subcluster 行从 parent_layer 的 metadata 中按 sample_include 或 selection_column/selection_values 取细胞\n'
    printf '# - 默认示例保留 subcluster_1 / subcluster_2；请按 sample_include 或 panorama 注释结果填好后再启用\n'
    printf 'layer_id\tlayer_role\tenabled\tparent_layer\tsample_include\tsample_exclude\tselection_column\tselection_values\trebuild_normalization\thvg_nfeatures\tpca_dims\ttarget_clusters\tres_range\tres_fine_step\tnormalization_methods\tintegration_mode\tvars_to_regress\tdescription\n'
    printf 'panorama\tpanorama\tyes\t\t\t\t\t\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${HVG_NFEATURES:-2000}" \
      "${PCA_DIMS_PANORAMA:-${PCA_DIMS:-1:30}}" \
      "${TARGET_CLUSTERS:-15}" \
      "${RES_RANGE:-0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60}" \
      "${RES_FINE_STEP:-0.005}" \
      "${NORMALIZATION_METHODS:-lognorm}" \
      "${INTEGRATION_MODES:-${INTEGRATION_MODE:-harmony}}" \
      "${VARS_TO_REGRESS_DEFAULT:-}" \
      'Root panorama object built from all post-QC cells.'
    printf 'subcluster_1\tsubcluster\tno\tpanorama\t\t\t\t\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${HVG_NFEATURES:-2000}" \
      "${PCA_DIMS_SUBCLUSTER:-1:20}" \
      "8" \
      "0.10,0.15,0.20,0.25,0.30,0.35,0.40" \
      "${RES_FINE_STEP:-0.005}" \
      "" \
      "" \
      "" \
      'Fill sample_include or selection_column/selection_values from panorama results before enabling this subcluster.'
    printf 'subcluster_2\tsubcluster\tno\tpanorama\t\t\t\t\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${HVG_NFEATURES:-2000}" \
      "${PCA_DIMS_SUBCLUSTER:-1:20}" \
      "6" \
      "0.10,0.15,0.20,0.25,0.30,0.35,0.40" \
      "${RES_FINE_STEP:-0.005}" \
      "" \
      "" \
      "" \
      'Fill sample_include or selection_column/selection_values from panorama results before enabling this subcluster.'
  } > "${OBJECT_LAYER_CONFIG_FILE}"
}

ensure_spatial_config_files() {
  ensure_dir "${PROJECT_CONFIG_DIR}" "${MARKER_PANEL_DIR}" "${METADATA_DIR}"

  if [[ ! -s "${SPATIAL_QC_THRESHOLD_FILE}" && -f "${PIPELINE_ROOT}/config/spatial_qc_thresholds.tsv.template" ]]; then
    cp -f "${PIPELINE_ROOT}/config/spatial_qc_thresholds.tsv.template" "${SPATIAL_QC_THRESHOLD_FILE}"
  fi
  if [[ ! -s "${SPATIAL_OBJECT_LAYER_FILE}" && -f "${PIPELINE_ROOT}/config/spatial_object_layers.tsv.template" ]]; then
    cp -f "${PIPELINE_ROOT}/config/spatial_object_layers.tsv.template" "${SPATIAL_OBJECT_LAYER_FILE}"
  fi
  if [[ ! -s "${SPATIAL_INTAKE_CONTRACT_FILE}" && -f "${PIPELINE_ROOT}/config/spatial_intake_contract.tsv.template" ]]; then
    cp -f "${PIPELINE_ROOT}/config/spatial_intake_contract.tsv.template" "${SPATIAL_INTAKE_CONTRACT_FILE}"
  fi
  if [[ ! -s "${MARKER_PANEL_DIR}/spatial_region_panel.tsv.template" && -f "${PIPELINE_ROOT}/config/marker_panels/spatial_region_panel.tsv.template" ]]; then
    cp -f "${PIPELINE_ROOT}/config/marker_panels/spatial_region_panel.tsv.template" "${MARKER_PANEL_DIR}/spatial_region_panel.tsv.template"
  fi
  if [[ ! -s "${SECTION_SHEET}" && -f "${PIPELINE_ROOT}/metadata/templates/sections.tsv.template" ]]; then
    cp -f "${PIPELINE_ROOT}/metadata/templates/sections.tsv.template" "${SECTION_SHEET}"
  fi
  if [[ ! -s "${SPATIAL_REFERENCE_INVENTORY_FILE}" && -f "${PIPELINE_ROOT}/metadata/spatial_reference_inventory.tsv" ]]; then
    cp -f "${PIPELINE_ROOT}/metadata/spatial_reference_inventory.tsv" "${SPATIAL_REFERENCE_INVENTORY_FILE}"
  fi
}

ensure_marker_panel_dir() {
  ensure_dir "${MARKER_PANEL_DIR}"

  if [[ -f "${PIPELINE_ROOT}/config/marker_panels/README.md" && ! -f "${MARKER_PANEL_DIR}/README.md" ]]; then
    cp -f "${PIPELINE_ROOT}/config/marker_panels/README.md" "${MARKER_PANEL_DIR}/README.md"
  fi
  if [[ -f "${PIPELINE_ROOT}/config/marker_panels/example_panel.tsv.template" && ! -f "${MARKER_PANEL_DIR}/example_panel.tsv.template" ]]; then
    cp -f "${PIPELINE_ROOT}/config/marker_panels/example_panel.tsv.template" "${MARKER_PANEL_DIR}/example_panel.tsv.template"
  fi
}

ensure_mito_gene_list_file() {
  ensure_dir "${PROJECT_CONFIG_DIR}"

  if [[ -s "${MITO_GENE_LIST_FILE}" ]]; then
    return 0
  fi

  if [[ -f "${PIPELINE_ROOT}/config/mito_gene_list.txt.template" ]]; then
    cp -f "${PIPELINE_ROOT}/config/mito_gene_list.txt.template" "${MITO_GENE_LIST_FILE}"
    return 0
  fi

  cat > "${MITO_GENE_LIST_FILE}" <<'EOF'
# One mito gene identifier per line.
# Use gene symbols or feature IDs that exactly match your matrix rownames.
# Leave this file unchanged when GTF-based mito detection is sufficient.
EOF
}

now_iso() {
  date --iso-8601=seconds
}

detect_python() {
  if command -v python3 >/dev/null 2>&1; then
    echo "python3"
    return
  fi

  local candidate=""
  for candidate in \
    "${R_MAIN_ENV_PREFIX:-}/bin/python" \
    "${R_SCENIC_ENV_PREFIX:-}/bin/python" \
    "${PY_SPATIAL_ENV_PREFIX:-}/bin/python" \
    "${VELOCITY_ENV_PREFIX:-}/bin/python"; do
    if [[ -n "${candidate}" && -x "${candidate}" ]]; then
      echo "${candidate}"
      return
    fi
  done

  if command -v python >/dev/null 2>&1 && python - <<'PY' >/dev/null 2>&1
import sys
sys.exit(0 if sys.version_info[0] >= 3 else 1)
PY
  then
    echo "python"
    return
  fi

  die "系统中缺少可用的 Python 3 解释器，无法写入 workflow 状态文件。"
}

bool_string() {
  case "${1:-false}" in
    true|TRUE|1|yes|YES|on|ON)
      echo "true"
      ;;
    *)
      echo "false"
      ;;
  esac
}

tsv_get_col_index() {
  local file_path="$1"
  local column_name="$2"
  awk -F '\t' -v target="${column_name}" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == target) {
          print i
          found = 1
          exit 0
        }
      }
      exit(found ? 0 : 1)
    }
  ' "${file_path}"
}

tsv_require_columns() {
  local file_path="$1"
  shift
  local column_name
  [[ -f "${file_path}" ]] || die "缺少 TSV 文件: ${file_path}"
  for column_name in "$@"; do
    if ! tsv_get_col_index "${file_path}" "${column_name}" >/dev/null; then
      die "TSV 文件 ${file_path} 缺少必需列: ${column_name}"
    fi
  done
}

ensure_metadata_fresh() {
  local questions_file="${ANALYSIS_QUESTIONS_FILE:-${METADATA_DIR}/analysis_questions.tsv}"
  local generated_tables=(
    "${COMPARISON_SHEET}"
    "${COMMUNICATION_PAIRS_SHEET}"
    "${TRAJECTORY_PAIRS_SHEET}"
    "${SCENIC_TARGETS_SHEET}"
    "${ENRICHMENT_TARGETS_SHEET}"
    "${GENE_PROGRAM_TARGETS_SHEET}"
    "${DECONV_PAIRS_SHEET}"
    "${SPATIAL_PAIRS_SHEET}"
    "${SCDESIGN3_QUESTION_MAP}"
    "${SCDESIGN3_TARGETS_SHEET}"
    "${SCDESIGN3_SIMULATION_DESIGNS_SHEET}"
    "${SCDESIGN3_THRESHOLDS_SHEET}"
  )
  local generator="${PIPELINE_ROOT}/workflow/03stages/95_run_metadata_generator.sh"
  local validator="${PIPELINE_ROOT}/workflow/03stages/96_validate_metadata.sh"
  local table_path
  local needs_generate=0

  [[ -s "${questions_file}" ]] || die "缺少 analysis questions 表: ${questions_file}"
  [[ -x "${generator}" ]] || die "缺少可执行 metadata generator: ${generator}"
  [[ -x "${validator}" ]] || die "缺少可执行 metadata validator: ${validator}"

  for table_path in "${generated_tables[@]}"; do
    if [[ ! -s "${table_path}" || "${questions_file}" -nt "${table_path}" ]]; then
      needs_generate=1
      break
    fi
  done

  if [[ "${needs_generate}" == "1" ]]; then
    echo "[metadata] analysis_questions.tsv 更新或 Tier 2 表缺失，重新生成 metadata。"
    RUN_METADATA_VALIDATION_AFTER_GENERATE=no bash "${generator}"
  else
    echo "[metadata] Tier 2 metadata 已是最新。"
  fi

  bash "${validator}"
}

prepare_project_state_dirs() {
  ensure_dir \
    "${PROJECT_CONFIG_DIR}" \
    "${RAW_RECEIVED_DIR}" \
    "${METADATA_DIR}" \
    "${REPORT_DIR}" \
    "${INTAKE_REPORT_DIR}" \
    "${EDA_REPORT_DIR}" \
    "${PRE_QC_REPORT_DIR}" \
    "${POST_QC_REPORT_DIR}" \
    "${INTEGRATION_REPORT_DIR}" \
    "${ANNOTATION_REPORT_DIR}" \
    "${SUBCLUSTER_REPORT_DIR}" \
    "${DEG_REPORT_DIR}" \
    "${ENRICHMENT_REPORT_DIR}" \
    "${COMMUNICATION_REPORT_DIR}" \
    "${REGULATION_REPORT_DIR}" \
    "${TRAJECTORY_REPORT_DIR}" \
    "${VELOCITY_REPORT_DIR}" \
    "${SPATIAL_RESULTS_DIR}" \
    "${SPATIAL_CHECKPOINT_DIR}" \
    "${SPATIAL_FIGURE_DIR}" \
    "${SPATIAL_TABLE_DIR}" \
    "${SPATIAL_REFERENCE_FREEZE_DIR}" \
    "${JOINT_RESULTS_DIR}" \
    "${JOINT_CHECKPOINT_DIR}" \
    "${JOINT_FIGURE_DIR}" \
    "${JOINT_TABLE_DIR}" \
    "${SPATIAL_PRE_QC_REPORT_DIR}" \
    "${SPATIAL_POST_QC_REPORT_DIR}" \
    "${SPATIAL_INTEGRATION_REPORT_DIR}" \
    "${SPATIAL_REGION_ANNOTATION_REPORT_DIR}" \
    "${SPATIAL_DECONV_REPORT_DIR}" \
    "${JOINT_TRAJECTORY_REPORT_DIR}" \
    "${JOINT_COMMUNICATION_REPORT_DIR}" \
    "${STATUS_DIR}" \
    "${MANIFEST_DIR}" \
    "${LOG_DIR}"
  ensure_eda_control_files
  ensure_object_layer_config_file
  ensure_spatial_config_files
  ensure_marker_panel_dir
  ensure_mito_gene_list_file
}

LIB_DIR="${SCRIPT_DIR}"
for lib_file in gate.sh stage.sh env_registry.sh executor.sh; do
  # shellcheck disable=SC1090
  source "${LIB_DIR}/${lib_file}"
done
