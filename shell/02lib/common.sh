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
export SCDESIGN3_N_SIM="${SCDESIGN3_N_SIM:-1}"
export SCDESIGN3_FAMILY="${SCDESIGN3_FAMILY:-nb}"
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
export R_MAIN_ENV_PREFIX="${R_MAIN_ENV_PREFIX:-${ENV_DIR}/conda/r_main}"
export R_SCENIC_ENV_PREFIX="${R_SCENIC_ENV_PREFIX:-${ENV_DIR}/conda/r_scenic}"
export PYSCENIC_ENV_PREFIX="${PYSCENIC_ENV_PREFIX:-${ENV_DIR}/conda/pyscenic}"
export SCVELO_ENV_PREFIX="${SCVELO_ENV_PREFIX:-${ENV_DIR}/conda/scvelo}"
export SCENIC_TF_LIST="${SCENIC_TF_LIST:-}"
export SCENIC_MOTIF_ANN="${SCENIC_MOTIF_ANN:-}"
export SCENIC_DB_500BP="${SCENIC_DB_500BP:-}"
export SCENIC_DB_10KB="${SCENIC_DB_10KB:-}"
export SCENIC_ORTHOLOG_MAP_FILE="${SCENIC_ORTHOLOG_MAP_FILE:-}"
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
export DELIVERY_MANIFEST="${DELIVERY_MANIFEST:-${METADATA_DIR}/delivery_manifest.tsv}"
export RECEIVED_FILES_MANIFEST="${RECEIVED_FILES_MANIFEST:-${METADATA_DIR}/received_files_manifest.tsv}"
export INPUT_INVENTORY_FILE="${INPUT_INVENTORY_FILE:-${INTAKE_REPORT_DIR}/input_inventory.tsv}"
export BRANCH_READINESS_FILE="${BRANCH_READINESS_FILE:-${INTAKE_REPORT_DIR}/branch_readiness.tsv}"
export INTAKE_SUMMARY_FILE="${INTAKE_SUMMARY_FILE:-${INTAKE_REPORT_DIR}/intake_summary.md}"
export WAIVER_FILE="${WAIVER_FILE:-${PROJECT_CONFIG_DIR}/waivers.tsv}"
export QC_THRESHOLD_FILE="${QC_THRESHOLD_FILE:-${PROJECT_CONFIG_DIR}/qc_thresholds.tsv}"
export EDA_GATE_FILE="${EDA_GATE_FILE:-${PROJECT_CONFIG_DIR}/eda_gates.tsv}"
export OBJECT_LAYER_CONFIG_FILE="${OBJECT_LAYER_CONFIG_FILE:-${PROJECT_CONFIG_DIR}/object_layers.tsv}"
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

export PATH="${SRA_TOOLS_DIR:+${SRA_TOOLS_DIR}:}${PATH}"

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

warn() {
  echo "[WARN] $*" >&2
}

ensure_dir() {
  mkdir -p "$@"
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
    "${ENRICHMENT_REPORT_DIR}"

  if [[ ! -s "${EDA_GATE_FILE}" ]]; then
    {
      printf 'gate_id\tstatus\tapproved_by\tnotes\n'
      printf 'pre_qc\tpending\t\t\n'
      printf 'post_qc\tpending\t\t\n'
      printf 'integration\tpending\t\t\n'
      printf 'annotation\tpending\t\t\n'
      printf 'subcluster\tpending\t\t\n'
      printf 'deg\tpending\t\t\n'
    } > "${EDA_GATE_FILE}"
  fi

  local gate_id
  for gate_id in pre_qc post_qc integration annotation subcluster deg; do
    if ! awk -F '\t' -v gate="${gate_id}" 'NR > 1 && $1 == gate { found = 1 } END { exit(found ? 0 : 1) }' "${EDA_GATE_FILE}" >/dev/null 2>&1; then
      printf '%s\tpending\t\t\n' "${gate_id}" >> "${EDA_GATE_FILE}"
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
    "${STATUS_DIR}" \
    "${MANIFEST_DIR}" \
    "${LOG_DIR}"
  ensure_eda_control_files
  ensure_object_layer_config_file
  ensure_marker_panel_dir
  ensure_mito_gene_list_file
}

LIB_DIR="${SCRIPT_DIR}"
for lib_file in gate.sh stage.sh env_registry.sh executor.sh; do
  # shellcheck disable=SC1090
  source "${LIB_DIR}/${lib_file}"
done
