detect_container_runtime() {
  if command -v apptainer >/dev/null 2>&1; then
    echo "apptainer"
    return
  fi
  if command -v singularity >/dev/null 2>&1; then
    echo "singularity"
    return
  fi
  echo ""
}

detect_conda_frontend() {
  if command -v micromamba >/dev/null 2>&1; then
    echo "micromamba"
    return
  fi
  if command -v mamba >/dev/null 2>&1; then
    echo "mamba"
    return
  fi
  if command -v conda >/dev/null 2>&1; then
    echo "conda"
    return
  fi
  echo ""
}

CONTAINER_RUNTIME="${CONTAINER_RUNTIME:-$(detect_container_runtime)}"
CONDA_FRONTEND="${CONDA_FRONTEND:-$(detect_conda_frontend)}"

declare -ag CONTAINER_BIND_ARGS=()
declare -ag BIND_PATHS=(
  "${PROJECT_ROOT}"
  "${PAPER_MATRIX_SOURCE:-}"
  "${CELLRANGER_OUT_DIR:-}"
  "${DNBC4TOOLS_OUT_DIR:-}"
  "${ORTHOLOG_CACHE_DIR:-}"
  "${REFERENCE_DIR:-}"
  "${ENV_DIR:-}"
)

if [[ -n "${EXISTING_R_SIF:-}" ]]; then
  BIND_PATHS+=("$(dirname "${EXISTING_R_SIF}")")
fi

for bind_path in "${BIND_PATHS[@]}"; do
  [[ -n "${bind_path}" ]] || continue
  if [[ -e "${bind_path}" ]]; then
    CONTAINER_BIND_ARGS+=(--bind "${bind_path}:${bind_path}")
  fi
done

run_r_main() {
  local script_path="$1"
  shift || true

  if [[ "${USE_EXISTING_SIF}" == "yes" ]]; then
    [[ -n "${CONTAINER_RUNTIME}" ]] || die "USE_EXISTING_SIF=yes，但系统中没有 apptainer/singularity。"
    [[ -f "${EXISTING_R_SIF}" ]] || die "找不到旧 R 容器: ${EXISTING_R_SIF}"

    "${CONTAINER_RUNTIME}" exec "${CONTAINER_BIND_ARGS[@]}" "${EXISTING_R_SIF}" \
      env \
      PIPELINE_ROOT="${PIPELINE_ROOT}" \
      R_LIBS_USER="${R_LIBS_MAIN}" \
      PROJECT_ROOT="${PROJECT_ROOT}" \
      DATA_DIR="${DATA_DIR}" \
      RESULTS_DIR="${RESULTS_DIR}" \
      CHECKPOINT_DIR="${CHECKPOINT_DIR}" \
      FIGURE_DIR="${FIGURE_DIR}" \
      TABLE_DIR="${TABLE_DIR}" \
      LOG_DIR="${LOG_DIR}" \
      ORTHOLOG_CACHE_DIR="${ORTHOLOG_CACHE_DIR}" \
      CELL_CYCLE_GENES_RDS="${CELL_CYCLE_GENES_RDS}" \
      REFERENCE_DIR="${REFERENCE_DIR}" \
      REFERENCE_GTF="${REFERENCE_GTF}" \
      CLEAN_GTF="${CLEAN_GTF}" \
      GENOME_FASTA_GZ="${GENOME_FASTA_GZ}" \
      CELLRANGER_OUT_DIR="${CELLRANGER_OUT_DIR}" \
      DNBC4TOOLS_OUT_DIR="${DNBC4TOOLS_OUT_DIR}" \
      STAR_INDEX_DIR="${STAR_INDEX_DIR}" \
      STAR_THREADS="${STAR_THREADS}" \
      DNBC4TOOLS_THREADS="${DNBC4TOOLS_THREADS}" \
      STAR_SA_INDEX_NBASES="${STAR_SA_INDEX_NBASES}" \
      VELOCITY_INPUT_DIR="${VELOCITY_INPUT_DIR}" \
      VELOCITY_LOOM_DIR="${VELOCITY_LOOM_DIR}" \
      VELOCITY_OUTPUT_DIR="${VELOCITY_OUTPUT_DIR}" \
      SCENIC_INPUT_DIR="${SCENIC_INPUT_DIR}" \
      SCENIC_OUTPUT_DIR="${SCENIC_OUTPUT_DIR}" \
      SCENIC_ORTHOLOG_MAP_FILE="${SCENIC_ORTHOLOG_MAP_FILE:-}" \
      RESOURCE_DIR="${RESOURCE_DIR}" \
      RAW_GROUP_1_NAME="${RAW_GROUP_1_NAME}" \
      RAW_GROUP_1_SAMPLES="${RAW_GROUP_1_SAMPLES}" \
      RAW_GROUP_2_NAME="${RAW_GROUP_2_NAME}" \
      RAW_GROUP_2_SAMPLES="${RAW_GROUP_2_SAMPLES}" \
      RAW_SAMPLES="${RAW_SAMPLES}" \
      DEG_IDENT_1="${DEG_IDENT_1}" \
      DEG_IDENT_2="${DEG_IDENT_2}" \
      RANDOM_SEED="${RANDOM_SEED}" \
      SAMPLE_NAMES="${SAMPLE_NAMES}" \
      OBJECT_LAYER_CONFIG_FILE="${OBJECT_LAYER_CONFIG_FILE}" \
      QC_MIN_NFEATURE="${QC_MIN_NFEATURE}" \
      QC_MIN_NCOUNT="${QC_MIN_NCOUNT}" \
      QC_MIN_LOG10UMI="${QC_MIN_LOG10UMI}" \
      QC_MAX_MITO_PCT="${QC_MAX_MITO_PCT}" \
      AMBIENT_PRIMARY_METHOD="${AMBIENT_PRIMARY_METHOD}" \
      AMBIENT_FALLBACK_METHOD="${AMBIENT_FALLBACK_METHOD}" \
      AMBIENT_APPLY_POLICY="${AMBIENT_APPLY_POLICY}" \
      AMBIENT_MIN_CELLS="${AMBIENT_MIN_CELLS}" \
      AMBIENT_CLUSTER_DIMS="${AMBIENT_CLUSTER_DIMS}" \
      AMBIENT_CLUSTER_RESOLUTION="${AMBIENT_CLUSTER_RESOLUTION}" \
      AMBIENT_MARKER_TOP_N="${AMBIENT_MARKER_TOP_N}" \
      AMBIENT_RECOMMEND_MIN_CONTAMINATION="${AMBIENT_RECOMMEND_MIN_CONTAMINATION}" \
      CELLBENDER_MODE="${CELLBENDER_MODE}" \
      CELLBENDER_FPR="${CELLBENDER_FPR}" \
      CELLBENDER_CUDA="${CELLBENDER_CUDA}" \
      CELLBENDER_EXTRA_ARGS="${CELLBENDER_EXTRA_ARGS}" \
      DOUBLET_RATE="${DOUBLET_RATE}" \
      DOUBLET_RATE_PER_1K="${DOUBLET_RATE_PER_1K}" \
      DOUBLET_PRIMARY_CALLER="${DOUBLET_PRIMARY_CALLER}" \
      DOUBLET_SECONDARY_CALLER="${DOUBLET_SECONDARY_CALLER}" \
      DOUBLET_SECONDARY_ENABLED="${DOUBLET_SECONDARY_ENABLED}" \
      DOUBLET_MIN_CELLS="${DOUBLET_MIN_CELLS}" \
      DOUBLET_DIMS="${DOUBLET_DIMS}" \
      HVG_NFEATURES="${HVG_NFEATURES}" \
      PCA_DIMS="${PCA_DIMS}" \
      TARGET_CLUSTERS="${TARGET_CLUSTERS}" \
      RES_RANGE="${RES_RANGE}" \
      TRAJECTORY_START="${TRAJECTORY_START}" \
      TRADESEQ_KNOTS="${TRADESEQ_KNOTS}" \
      ENSEMBL_MIRROR="${ENSEMBL_MIRROR}" \
      MAIN_THREADS="${MAIN_THREADS}" \
      Rscript "${script_path}" "$@"
    return
  fi

  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda，无法运行主流程 R 环境。"
  [[ -d "${R_MAIN_ENV_PREFIX}" ]] || die "主流程 R 环境不存在: ${R_MAIN_ENV_PREFIX}"

  env \
    PIPELINE_ROOT="${PIPELINE_ROOT}" \
    R_LIBS_USER="${R_LIBS_MAIN}" \
    PROJECT_ROOT="${PROJECT_ROOT}" \
    DATA_DIR="${DATA_DIR}" \
    RESULTS_DIR="${RESULTS_DIR}" \
    CHECKPOINT_DIR="${CHECKPOINT_DIR}" \
    FIGURE_DIR="${FIGURE_DIR}" \
    TABLE_DIR="${TABLE_DIR}" \
    LOG_DIR="${LOG_DIR}" \
    ORTHOLOG_CACHE_DIR="${ORTHOLOG_CACHE_DIR}" \
    CELL_CYCLE_GENES_RDS="${CELL_CYCLE_GENES_RDS}" \
    REFERENCE_DIR="${REFERENCE_DIR}" \
    REFERENCE_GTF="${REFERENCE_GTF}" \
    CLEAN_GTF="${CLEAN_GTF}" \
    GENOME_FASTA_GZ="${GENOME_FASTA_GZ}" \
    CELLRANGER_OUT_DIR="${CELLRANGER_OUT_DIR}" \
    DNBC4TOOLS_OUT_DIR="${DNBC4TOOLS_OUT_DIR}" \
    STAR_INDEX_DIR="${STAR_INDEX_DIR}" \
    STAR_THREADS="${STAR_THREADS}" \
    DNBC4TOOLS_THREADS="${DNBC4TOOLS_THREADS}" \
    STAR_SA_INDEX_NBASES="${STAR_SA_INDEX_NBASES}" \
    VELOCITY_INPUT_DIR="${VELOCITY_INPUT_DIR}" \
    VELOCITY_LOOM_DIR="${VELOCITY_LOOM_DIR}" \
    VELOCITY_OUTPUT_DIR="${VELOCITY_OUTPUT_DIR}" \
    SCENIC_INPUT_DIR="${SCENIC_INPUT_DIR}" \
    SCENIC_OUTPUT_DIR="${SCENIC_OUTPUT_DIR}" \
    SCENIC_ORTHOLOG_MAP_FILE="${SCENIC_ORTHOLOG_MAP_FILE:-}" \
    RESOURCE_DIR="${RESOURCE_DIR}" \
    RAW_GROUP_1_NAME="${RAW_GROUP_1_NAME}" \
    RAW_GROUP_1_SAMPLES="${RAW_GROUP_1_SAMPLES}" \
    RAW_GROUP_2_NAME="${RAW_GROUP_2_NAME}" \
    RAW_GROUP_2_SAMPLES="${RAW_GROUP_2_SAMPLES}" \
    RAW_SAMPLES="${RAW_SAMPLES}" \
    DEG_IDENT_1="${DEG_IDENT_1}" \
    DEG_IDENT_2="${DEG_IDENT_2}" \
    RANDOM_SEED="${RANDOM_SEED}" \
    SAMPLE_NAMES="${SAMPLE_NAMES}" \
    OBJECT_LAYER_CONFIG_FILE="${OBJECT_LAYER_CONFIG_FILE}" \
    QC_MIN_NFEATURE="${QC_MIN_NFEATURE}" \
    QC_MIN_NCOUNT="${QC_MIN_NCOUNT}" \
    QC_MIN_LOG10UMI="${QC_MIN_LOG10UMI}" \
    QC_MAX_MITO_PCT="${QC_MAX_MITO_PCT}" \
    AMBIENT_PRIMARY_METHOD="${AMBIENT_PRIMARY_METHOD}" \
    AMBIENT_FALLBACK_METHOD="${AMBIENT_FALLBACK_METHOD}" \
    AMBIENT_APPLY_POLICY="${AMBIENT_APPLY_POLICY}" \
    AMBIENT_MIN_CELLS="${AMBIENT_MIN_CELLS}" \
    AMBIENT_CLUSTER_DIMS="${AMBIENT_CLUSTER_DIMS}" \
    AMBIENT_CLUSTER_RESOLUTION="${AMBIENT_CLUSTER_RESOLUTION}" \
    AMBIENT_MARKER_TOP_N="${AMBIENT_MARKER_TOP_N}" \
    AMBIENT_RECOMMEND_MIN_CONTAMINATION="${AMBIENT_RECOMMEND_MIN_CONTAMINATION}" \
    CELLBENDER_MODE="${CELLBENDER_MODE}" \
    CELLBENDER_FPR="${CELLBENDER_FPR}" \
    CELLBENDER_CUDA="${CELLBENDER_CUDA}" \
    CELLBENDER_EXTRA_ARGS="${CELLBENDER_EXTRA_ARGS}" \
    DOUBLET_RATE="${DOUBLET_RATE}" \
    DOUBLET_RATE_PER_1K="${DOUBLET_RATE_PER_1K}" \
    DOUBLET_PRIMARY_CALLER="${DOUBLET_PRIMARY_CALLER}" \
    DOUBLET_SECONDARY_CALLER="${DOUBLET_SECONDARY_CALLER}" \
    DOUBLET_SECONDARY_ENABLED="${DOUBLET_SECONDARY_ENABLED}" \
    DOUBLET_MIN_CELLS="${DOUBLET_MIN_CELLS}" \
    DOUBLET_DIMS="${DOUBLET_DIMS}" \
    HVG_NFEATURES="${HVG_NFEATURES}" \
    PCA_DIMS="${PCA_DIMS}" \
    TARGET_CLUSTERS="${TARGET_CLUSTERS}" \
    RES_RANGE="${RES_RANGE}" \
    TRAJECTORY_START="${TRAJECTORY_START}" \
    TRADESEQ_KNOTS="${TRADESEQ_KNOTS}" \
    ENSEMBL_MIRROR="${ENSEMBL_MIRROR}" \
    MAIN_THREADS="${MAIN_THREADS}" \
    "${CONDA_FRONTEND}" run -p "${R_MAIN_ENV_PREFIX}" Rscript "${script_path}" "$@"
}

run_r_scenic() {
  local script_path="$1"
  shift || true

  if [[ "${USE_EXISTING_SIF}" == "yes" ]]; then
    [[ -n "${CONTAINER_RUNTIME}" ]] || die "USE_EXISTING_SIF=yes，但系统中没有 apptainer/singularity。"
    [[ -f "${EXISTING_R_SIF}" ]] || die "找不到旧 R 容器: ${EXISTING_R_SIF}"

    "${CONTAINER_RUNTIME}" exec "${CONTAINER_BIND_ARGS[@]}" "${EXISTING_R_SIF}" \
      env \
      PIPELINE_ROOT="${PIPELINE_ROOT}" \
      R_LIBS_USER="${R_LIBS_SCENIC}" \
      PROJECT_ROOT="${PROJECT_ROOT}" \
      DATA_DIR="${DATA_DIR}" \
      RESULTS_DIR="${RESULTS_DIR}" \
      CHECKPOINT_DIR="${CHECKPOINT_DIR}" \
      FIGURE_DIR="${FIGURE_DIR}" \
      TABLE_DIR="${TABLE_DIR}" \
      LOG_DIR="${LOG_DIR}" \
      SCENIC_INPUT_DIR="${SCENIC_INPUT_DIR}" \
      SCENIC_OUTPUT_DIR="${SCENIC_OUTPUT_DIR}" \
      SCENIC_TF_LIST="${SCENIC_TF_LIST}" \
      SCENIC_MOTIF_ANN="${SCENIC_MOTIF_ANN}" \
      SCENIC_DB_500BP="${SCENIC_DB_500BP}" \
      SCENIC_DB_10KB="${SCENIC_DB_10KB}" \
      RANDOM_SEED="${RANDOM_SEED}" \
      Rscript "${script_path}" "$@"
    return
  fi

  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda，无法运行 SCENIC R 环境。"
  [[ -d "${R_SCENIC_ENV_PREFIX}" ]] || die "SCENIC R 环境不存在: ${R_SCENIC_ENV_PREFIX}"

  env \
    PIPELINE_ROOT="${PIPELINE_ROOT}" \
    R_LIBS_USER="${R_LIBS_SCENIC}" \
    PROJECT_ROOT="${PROJECT_ROOT}" \
    DATA_DIR="${DATA_DIR}" \
    RESULTS_DIR="${RESULTS_DIR}" \
    CHECKPOINT_DIR="${CHECKPOINT_DIR}" \
    FIGURE_DIR="${FIGURE_DIR}" \
    TABLE_DIR="${TABLE_DIR}" \
    LOG_DIR="${LOG_DIR}" \
    SCENIC_INPUT_DIR="${SCENIC_INPUT_DIR}" \
    SCENIC_OUTPUT_DIR="${SCENIC_OUTPUT_DIR}" \
    SCENIC_TF_LIST="${SCENIC_TF_LIST}" \
    SCENIC_MOTIF_ANN="${SCENIC_MOTIF_ANN}" \
    SCENIC_DB_500BP="${SCENIC_DB_500BP}" \
    SCENIC_DB_10KB="${SCENIC_DB_10KB}" \
    RANDOM_SEED="${RANDOM_SEED}" \
    "${CONDA_FRONTEND}" run -p "${R_SCENIC_ENV_PREFIX}" Rscript "${script_path}" "$@"
}

run_in_conda_prefix() {
  local prefix="$1"
  shift
  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda。"
  "${CONDA_FRONTEND}" run -p "${prefix}" "$@"
}

create_or_update_conda_env() {
  local prefix="$1"
  local yaml_path="$2"
  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda。"

  if [[ -d "${prefix}" ]]; then
    "${CONDA_FRONTEND}" env update -y -p "${prefix}" -f "${yaml_path}" --prune
  else
    "${CONDA_FRONTEND}" env create -y -p "${prefix}" -f "${yaml_path}"
  fi
}

run_pyscenic() {
  [[ -d "${PYSCENIC_ENV_PREFIX}" ]] || die "pySCENIC 环境不存在: ${PYSCENIC_ENV_PREFIX}"
  run_in_conda_prefix "${PYSCENIC_ENV_PREFIX}" "$@"
}

run_r_legacy() {
  [[ -d "${R_LEGACY_ENV_PREFIX}" ]] || die "legacy R 环境不存在: ${R_LEGACY_ENV_PREFIX}"
  run_in_conda_prefix "${R_LEGACY_ENV_PREFIX}" Rscript "$@"
}

run_r_interaction() {
  [[ -d "${R_INTERACTION_ENV_PREFIX}" ]] || die "通讯分析 R 环境不存在: ${R_INTERACTION_ENV_PREFIX}"
  run_in_conda_prefix "${R_INTERACTION_ENV_PREFIX}" Rscript "$@"
}

run_r_spatial() {
  [[ -d "${R_SPATIAL_ENV_PREFIX}" ]] || die "空间分析 R 环境不存在: ${R_SPATIAL_ENV_PREFIX}"
  run_in_conda_prefix "${R_SPATIAL_ENV_PREFIX}" Rscript "$@"
}

run_py_spatial() {
  [[ -d "${PY_SPATIAL_ENV_PREFIX}" ]] || die "空间分析 Python 环境不存在: ${PY_SPATIAL_ENV_PREFIX}"
  run_in_conda_prefix "${PY_SPATIAL_ENV_PREFIX}" "$@"
}

run_py_spatial_legacy() {
  [[ -d "${PY_SPATIAL_LEGACY_ENV_PREFIX}" ]] || die "空间 legacy Python 环境不存在: ${PY_SPATIAL_LEGACY_ENV_PREFIX}"
  run_in_conda_prefix "${PY_SPATIAL_LEGACY_ENV_PREFIX}" "$@"
}

run_py_cell2location() {
  [[ -d "${PY_CELL2LOCATION_ENV_PREFIX}" ]] || die "cell2location Python 环境不存在: ${PY_CELL2LOCATION_ENV_PREFIX}"
  run_in_conda_prefix "${PY_CELL2LOCATION_ENV_PREFIX}" "$@"
}

run_r_validation() {
  [[ -d "${R_VALIDATION_ENV_PREFIX}" ]] || die "验证分析 R 环境不存在: ${R_VALIDATION_ENV_PREFIX}"
  run_in_conda_prefix "${R_VALIDATION_ENV_PREFIX}" Rscript "$@"
}
