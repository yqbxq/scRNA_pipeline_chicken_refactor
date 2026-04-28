#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PIPELINE_CONFIG_DIR="${PIPELINE_ROOT}/config"
GLOBAL_CONFIG="${PIPELINE_ROOT}/config/server_config.sh"
ROOT_CONFIG_SHIM="${PIPELINE_ROOT}/server_config.sh"

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

usage() {
  cat <<'EOF'
用法：
  bash workflow/00_init_project.sh \
    --project-root /home/user_test/projects/chicken_case01 \
    --fastq-source-dir /path/to/fastq \
    --genome-fasta /path/to/genome.fa.gz \
    --reference-gtf /path/to/genes.gtf \
    --sample-names Ctrl_1,Ctrl_2,Treat_1,Treat_2 \
    --group1-name Ctrl \
    --group1-samples Ctrl_1,Ctrl_2 \
    --group2-name Treat \
    --group2-samples Treat_1,Treat_2

常用可选项：
  --cellranger-out-source-dir /path/to/existing/cellranger_out
  --matrix-source-dir /path/to/existing/matrix_root
  --paper-matrix-source /home/user_test/scRNA_sichuan/data
  --cellranger-bin /Bio/bin/cel1ranger
  --sra-tools-dir /conda_envs/softwares/sratoolkit.3.0.7-ubuntu64/bin

说明：
  1. 该脚本会在新项目目录下自动创建标准目录结构。
  2. FASTQ、参考文件、已有 cellranger_out、已有 matrix 会用软链接接入新项目。
  3. 会生成项目专属 config/project_config.sh。
  4. 会把全局 config/server_config.sh 自动切换到这个新项目。
EOF
}

link_file_into_dir() {
  local src="$1"
  local dst_dir="$2"
  [[ -f "${src}" ]] || die "文件不存在: ${src}"
  mkdir -p "${dst_dir}"
  ln -sfn "${src}" "${dst_dir}/$(basename "${src}")"
}

link_fastq_dir() {
  local src_dir="$1"
  local dst_dir="$2"
  [[ -d "${src_dir}" ]] || die "FASTQ 目录不存在: ${src_dir}"
  mkdir -p "${dst_dir}"

  local found=0
  while IFS= read -r -d '' file_path; do
    ln -sfn "${file_path}" "${dst_dir}/$(basename "${file_path}")"
    found=1
  done < <(find "${src_dir}" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) -print0)

  [[ "${found}" == "1" ]] || die "在 ${src_dir} 中没有找到 FASTQ 文件"
}

link_sample_dirs() {
  local src_root="$1"
  local dst_root="$2"
  local samples_csv="$3"
  local sample_array=()
  IFS=',' read -r -a sample_array <<< "${samples_csv}"

  mkdir -p "${dst_root}"
  for sample in "${sample_array[@]}"; do
    [[ -n "${sample}" ]] || continue
    [[ -d "${src_root}/${sample}" ]] || die "缺少样本目录: ${src_root}/${sample}"
    ln -sfn "${src_root}/${sample}" "${dst_root}/${sample}"
  done
}

write_export() {
  local name="$1"
  local value="$2"
  printf 'export %s=%q\n' "${name}" "${value}"
}

PROJECT_ROOT=""
FASTQ_SOURCE_DIR=""
GENOME_FASTA_GZ_SRC=""
REFERENCE_GTF_SRC=""
CELLRANGER_OUT_SOURCE_DIR=""
MATRIX_SOURCE_DIR=""
SAMPLE_NAMES=""
GROUP1_NAME=""
GROUP1_SAMPLES=""
GROUP2_NAME=""
GROUP2_SAMPLES=""
RAW_SAMPLES=""
RAW_GROUP1_NAME=""
RAW_GROUP1_SAMPLES=""
RAW_GROUP2_NAME=""
RAW_GROUP2_SAMPLES=""
PAPER_MATRIX_SOURCE="/home/user_test/scRNA_sichuan/data"
CELLRANGER_BIN="/Bio/bin/cel1ranger"
SRA_TOOLS_DIR="/conda_envs/softwares/sratoolkit.3.0.7-ubuntu64/bin"
MAIN_THREADS="16"
VELOCYTO_THREADS="8"
SCVELO_THREADS="16"
SCENIC_THREADS="16"
CELLRANGER_THREADS="96"
CELLRANGER_MEM_GB="128"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root)
      PROJECT_ROOT="$2"
      shift 2
      ;;
    --fastq-source-dir)
      FASTQ_SOURCE_DIR="$2"
      shift 2
      ;;
      --genome-fasta-gz|--genome-fasta)
        GENOME_FASTA_GZ_SRC="$2"
        shift 2
        ;;
    --reference-gtf)
      REFERENCE_GTF_SRC="$2"
      shift 2
      ;;
    --cellranger-out-source-dir)
      CELLRANGER_OUT_SOURCE_DIR="$2"
      shift 2
      ;;
    --matrix-source-dir)
      MATRIX_SOURCE_DIR="$2"
      shift 2
      ;;
    --sample-names)
      SAMPLE_NAMES="$2"
      shift 2
      ;;
    --group1-name)
      GROUP1_NAME="$2"
      shift 2
      ;;
    --group1-samples)
      GROUP1_SAMPLES="$2"
      shift 2
      ;;
    --group2-name)
      GROUP2_NAME="$2"
      shift 2
      ;;
    --group2-samples)
      GROUP2_SAMPLES="$2"
      shift 2
      ;;
    --raw-samples)
      RAW_SAMPLES="$2"
      shift 2
      ;;
    --raw-group1-name)
      RAW_GROUP1_NAME="$2"
      shift 2
      ;;
    --raw-group1-samples)
      RAW_GROUP1_SAMPLES="$2"
      shift 2
      ;;
    --raw-group2-name)
      RAW_GROUP2_NAME="$2"
      shift 2
      ;;
    --raw-group2-samples)
      RAW_GROUP2_SAMPLES="$2"
      shift 2
      ;;
    --paper-matrix-source)
      PAPER_MATRIX_SOURCE="$2"
      shift 2
      ;;
    --cellranger-bin)
      CELLRANGER_BIN="$2"
      shift 2
      ;;
    --sra-tools-dir)
      SRA_TOOLS_DIR="$2"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      die "未知参数: $1"
      ;;
  esac
done

[[ -n "${PROJECT_ROOT}" ]] || die "必须提供 --project-root"
[[ -n "${SAMPLE_NAMES}" ]] || die "必须提供 --sample-names"
[[ -n "${GROUP1_NAME}" ]] || die "必须提供 --group1-name"
[[ -n "${GROUP1_SAMPLES}" ]] || die "必须提供 --group1-samples"
[[ -n "${GROUP2_NAME}" ]] || die "必须提供 --group2-name"
[[ -n "${GROUP2_SAMPLES}" ]] || die "必须提供 --group2-samples"

if [[ -n "${FASTQ_SOURCE_DIR}" || -z "${CELLRANGER_OUT_SOURCE_DIR}" ]]; then
  [[ -n "${GENOME_FASTA_GZ_SRC}" ]] || die "从 FASTQ 开始时必须提供 --genome-fasta"
  [[ -n "${REFERENCE_GTF_SRC}" ]] || die "从 FASTQ 开始时必须提供 --reference-gtf"
fi

if [[ -z "${RAW_SAMPLES}" ]]; then
  RAW_SAMPLES="${SAMPLE_NAMES}"
fi
if [[ -z "${RAW_GROUP1_NAME}" ]]; then
  RAW_GROUP1_NAME="${GROUP1_NAME}"
fi
if [[ -z "${RAW_GROUP1_SAMPLES}" ]]; then
  RAW_GROUP1_SAMPLES="${GROUP1_SAMPLES}"
fi
if [[ -z "${RAW_GROUP2_NAME}" ]]; then
  RAW_GROUP2_NAME="${GROUP2_NAME}"
fi
if [[ -z "${RAW_GROUP2_SAMPLES}" ]]; then
  RAW_GROUP2_SAMPLES="${GROUP2_SAMPLES}"
fi

mkdir -p "$(dirname "${PROJECT_ROOT}")"
PROJECT_ROOT="$(cd "$(dirname "${PROJECT_ROOT}")" && pwd)/$(basename "${PROJECT_ROOT}")"
PROJECT_CONFIG_DIR="${PROJECT_ROOT}/config"
PROJECT_MITO_GENE_LIST_FILE="${PROJECT_CONFIG_DIR}/mito_gene_list.txt"
PROJECT_MARKER_PANEL_DIR="${PROJECT_CONFIG_DIR}/marker_panels"
PROJECT_CONFIG_FILE="${PROJECT_CONFIG_DIR}/project_config.sh"
PROJECT_RAW_DIR="${PROJECT_ROOT}/raw"
PROJECT_RAW_RECEIVED_DIR="${PROJECT_RAW_DIR}/received"
PROJECT_METADATA_DIR="${PROJECT_ROOT}/metadata"
PROJECT_REPORT_DIR="${PROJECT_ROOT}/reports"
PROJECT_INTAKE_REPORT_DIR="${PROJECT_REPORT_DIR}/intake"
PROJECT_STATUS_DIR="${PROJECT_ROOT}/status"
PROJECT_STATUS_FILE="${PROJECT_STATUS_DIR}/workflow_status.json"

DEFAULT_INPUT_MODE="matrix"
if [[ -n "${MATRIX_SOURCE_DIR}" ]]; then
  DEFAULT_INPUT_MODE="matrix"
elif [[ -n "${CELLRANGER_OUT_SOURCE_DIR}" ]]; then
  DEFAULT_INPUT_MODE="cellranger_out"
elif [[ -n "${FASTQ_SOURCE_DIR}" ]]; then
  DEFAULT_INPUT_MODE="fastq"
fi

DEFAULT_RUN_VELOCITY="no"
if [[ -n "${FASTQ_SOURCE_DIR}" || -n "${CELLRANGER_OUT_SOURCE_DIR}" ]]; then
  DEFAULT_RUN_VELOCITY="yes"
fi

mkdir -p \
  "${PROJECT_ROOT}" \
  "${PROJECT_CONFIG_DIR}" \
  "${PROJECT_MARKER_PANEL_DIR}" \
  "${PROJECT_RAW_RECEIVED_DIR}" \
  "${PROJECT_METADATA_DIR}" \
  "${PROJECT_REPORT_DIR}" \
  "${PROJECT_INTAKE_REPORT_DIR}" \
  "${PROJECT_REPORT_DIR}/eda/pre_qc" \
  "${PROJECT_REPORT_DIR}/eda/post_qc" \
  "${PROJECT_REPORT_DIR}/eda/integration" \
  "${PROJECT_REPORT_DIR}/eda/annotation" \
  "${PROJECT_REPORT_DIR}/eda/deg" \
  "${PROJECT_STATUS_DIR}" \
  "${PROJECT_ROOT}/data" \
  "${PROJECT_ROOT}/results/checkpoints" \
  "${PROJECT_ROOT}/results/figures" \
  "${PROJECT_ROOT}/results/tables" \
  "${PROJECT_ROOT}/results/velocity/input" \
  "${PROJECT_ROOT}/results/velocity/loom" \
  "${PROJECT_ROOT}/results/velocity/output" \
  "${PROJECT_ROOT}/results/scenic_input" \
  "${PROJECT_ROOT}/results/scenic_output" \
  "${PROJECT_ROOT}/logs" \
  "${PROJECT_ROOT}/envs" \
  "${PROJECT_ROOT}/resources/scenic_db" \
  "${PROJECT_ROOT}/fastq" \
  "${PROJECT_ROOT}/reference" \
  "${PROJECT_ROOT}/cellranger_out"

if [[ -n "${FASTQ_SOURCE_DIR}" ]]; then
  link_fastq_dir "${FASTQ_SOURCE_DIR}" "${PROJECT_ROOT}/fastq"
fi

if [[ -n "${GENOME_FASTA_GZ_SRC}" ]]; then
  link_file_into_dir "${GENOME_FASTA_GZ_SRC}" "${PROJECT_ROOT}/reference"
fi
if [[ -n "${REFERENCE_GTF_SRC}" ]]; then
  link_file_into_dir "${REFERENCE_GTF_SRC}" "${PROJECT_ROOT}/reference"
fi

if [[ -n "${CELLRANGER_OUT_SOURCE_DIR}" ]]; then
  link_sample_dirs "${CELLRANGER_OUT_SOURCE_DIR}" "${PROJECT_ROOT}/cellranger_out" "${RAW_SAMPLES}"
fi

if [[ -n "${MATRIX_SOURCE_DIR}" ]]; then
  link_sample_dirs "${MATRIX_SOURCE_DIR}" "${PROJECT_ROOT}/data" "${SAMPLE_NAMES}"
fi

GENOME_FASTA_GZ_LOCAL="${PROJECT_ROOT}/reference/$(basename "${GENOME_FASTA_GZ_SRC}")"
REFERENCE_GTF_LOCAL="${PROJECT_ROOT}/reference/$(basename "${REFERENCE_GTF_SRC}")"
gtf_name="$(basename "${REFERENCE_GTF_LOCAL}")"
gtf_stem="${gtf_name%.gtf}"
if [[ "${gtf_stem}" == "${gtf_name}" ]]; then
  gtf_stem="${gtf_name}"
fi
CLEAN_GTF_LOCAL="${PROJECT_ROOT}/reference/${gtf_stem}_clean.gtf"
CELLRANGER_REF_LOCAL="${PROJECT_ROOT}/reference/cellranger_ref"

{
  echo '#!/usr/bin/env bash'
  echo
  echo '# =============================== [TUNABLE] ==============================='
  echo '# 当前激活的正式项目根目录。'
  write_export PROJECT_ROOT "${PROJECT_ROOT}"
  echo
  echo '# 论文验证数据目录。正式项目通常不会用到，仅保留作为验证入口。'
  write_export PAPER_MATRIX_SOURCE "${PAPER_MATRIX_SOURCE}"
  echo
  echo '# 旧复现容器。正式项目默认不再使用。'
  write_export EXISTING_R_SIF ""
  echo
  echo '# 服务器已有软件路径。'
  write_export CELLRANGER_BIN "${CELLRANGER_BIN}"
  write_export SRA_TOOLS_DIR "${SRA_TOOLS_DIR}"
  echo
  echo '# 正式项目统一采用项目目录内的标准输入输出位置。'
  write_export FASTQ_DIR "${PROJECT_ROOT}/fastq"
  write_export CELLRANGER_OUT_DIR "${PROJECT_ROOT}/cellranger_out"
  write_export REFERENCE_DIR "${PROJECT_ROOT}/reference"
  write_export GENOME_FASTA_GZ "${GENOME_FASTA_GZ_LOCAL}"
  write_export REFERENCE_GTF "${REFERENCE_GTF_LOCAL}"
  write_export CLEAN_GTF "${CLEAN_GTF_LOCAL}"
  write_export CELLRANGER_REF_DIR "${CELLRANGER_REF_LOCAL}"
  echo
  echo '# 运行模式开关。正式项目默认使用独立 conda 环境。'
  write_export USE_EXISTING_SIF "no"
  write_export AUTHOR_INPUT_MODE "symlink"
  write_export RANDOM_SEED "42"
  write_export INTEGRATION_MODE "harmony"
  echo
  echo '# 计算资源参数。'
  write_export MAIN_THREADS "${MAIN_THREADS}"
  write_export VELOCYTO_THREADS "${VELOCYTO_THREADS}"
  write_export SCVELO_THREADS "${SCVELO_THREADS}"
  write_export SCENIC_THREADS "${SCENIC_THREADS}"
  write_export CELLRANGER_THREADS "${CELLRANGER_THREADS}"
  write_export CELLRANGER_MEM_GB "${CELLRANGER_MEM_GB}"
  echo
  echo '# 主流程输入样本名。'
  echo '# 这里写的是 data/ 下每个样本目录名。'
  write_export SAMPLE_NAMES "${SAMPLE_NAMES}"
  echo
  echo '# 主流程分组定义。'
  echo '# 这里决定 DEG 比较分组。组成员应写样本目录名，而不是原始 FASTQ 文件名。'
  write_export ANALYSIS_GROUP_1_NAME "${GROUP1_NAME}"
  write_export ANALYSIS_GROUP_1_SAMPLES "${GROUP1_SAMPLES}"
  write_export ANALYSIS_GROUP_2_NAME "${GROUP2_NAME}"
  write_export ANALYSIS_GROUP_2_SAMPLES "${GROUP2_SAMPLES}"
  echo
  echo '# 差异分析对比组。应填写上面 ANALYSIS_GROUP_*_NAME 的名字。'
  write_export DEG_IDENT_1 "${GROUP1_NAME}"
  write_export DEG_IDENT_2 "${GROUP2_NAME}"
  write_export DEG_MIN_CELLS_PER_GROUP "3"
  write_export DEG_LOGFC_THRESHOLD "0"
  write_export DEG_ALPHA "0.05"
  echo
  echo '# RNA velocity 使用的 raw 样本。'
  echo '# 正式项目里通常与 SAMPLE_NAMES 保持一致；论文验证模式才会出现 SAMPLE_NAMES 与 RAW_SAMPLES 不同。'
  write_export RAW_GROUP_1_NAME "${RAW_GROUP1_NAME}"
  write_export RAW_GROUP_1_SAMPLES "${RAW_GROUP1_SAMPLES}"
  write_export RAW_GROUP_2_NAME "${RAW_GROUP2_NAME}"
  write_export RAW_GROUP_2_SAMPLES "${RAW_GROUP2_SAMPLES}"
  write_export RAW_SAMPLES "${RAW_SAMPLES}"
  echo
  echo '# 主流程参数。'
  write_export QC_MIN_NFEATURE "200"
  write_export QC_MIN_NCOUNT "1000"
  write_export QC_MIN_LOG10UMI "0.7"
  write_export QC_MAX_MITO_PCT "20"
  write_export AMBIENT_PRIMARY_METHOD "soupx"
  write_export AMBIENT_FALLBACK_METHOD "decontx"
  write_export AMBIENT_APPLY_POLICY "manual"
  write_export AMBIENT_MIN_CELLS "50"
  write_export AMBIENT_CLUSTER_DIMS "1:20"
  write_export AMBIENT_CLUSTER_RESOLUTION "0.4"
  write_export AMBIENT_MARKER_TOP_N "3"
  write_export AMBIENT_RECOMMEND_MIN_CONTAMINATION "0.05"
  write_export CELLBENDER_MODE "stub"
  write_export CELLBENDER_FPR "0.01"
  write_export CELLBENDER_CUDA "yes"
  write_export CELLBENDER_EXTRA_ARGS ""
  write_export DOUBLET_RATE "0.008"
  write_export DOUBLET_RATE_PER_1K "0.008"
  write_export DOUBLET_PRIMARY_CALLER "scDblFinder"
  write_export DOUBLET_SECONDARY_CALLER "DoubletFinder"
  write_export DOUBLET_SECONDARY_ENABLED "yes"
  write_export DOUBLET_MIN_CELLS "50"
  write_export DOUBLET_DIMS "1:20"
  write_export HVG_NFEATURES "2000"
  write_export PCA_DIMS "1:30"
  write_export TARGET_CLUSTERS "15"
  write_export RES_RANGE "0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60"
  write_export RES_FINE_STEP "0.005"
  write_export TRAJECTORY_START ""
  write_export TRADESEQ_KNOTS "6"
  echo
  echo '# 外部镜像与资源。'
  write_export ENSEMBL_MIRROR "asia"
  echo
  echo '# ================================ [FIXED] ================================'
  echo 'export DATA_DIR="${PROJECT_ROOT}/data"'
  echo 'export RESULTS_DIR="${PROJECT_ROOT}/results"'
  echo 'export ORTHOLOG_CACHE_DIR="${RESULTS_DIR}/ortholog_cache"'
  echo 'export CHECKPOINT_DIR="${RESULTS_DIR}/checkpoints"'
  echo 'export FIGURE_DIR="${RESULTS_DIR}/figures"'
  echo 'export TABLE_DIR="${RESULTS_DIR}/tables"'
  echo 'export LOG_DIR="${PROJECT_ROOT}/logs"'
  echo 'export ENV_DIR="${PROJECT_ROOT}/envs"'
  echo 'export RESOURCE_DIR="${PROJECT_ROOT}/resources"'
  echo 'export PROJECT_CONFIG_DIR="${PROJECT_ROOT}/config"'
  echo 'export RAW_DIR="${PROJECT_ROOT}/raw"'
  echo 'export RAW_RECEIVED_DIR="${RAW_DIR}/received"'
  echo 'export METADATA_DIR="${PROJECT_ROOT}/metadata"'
  echo 'export REPORT_DIR="${PROJECT_ROOT}/reports"'
  echo 'export INTAKE_REPORT_DIR="${REPORT_DIR}/intake"'
  echo 'export EDA_REPORT_DIR="${REPORT_DIR}/eda"'
  echo 'export STATUS_DIR="${PROJECT_ROOT}/status"'
  echo 'export WORKFLOW_STATUS_FILE="${STATUS_DIR}/workflow_status.json"'
  echo 'export SHARED_ENV_DIR="${SHARED_ENV_DIR:-/home/user_test/syf_f5/01shared_resources/envs/scRNA_pipeline_chicken}"'
  echo 'export PIPELINE_ENV_DIR="${PIPELINE_ENV_DIR:-${SHARED_ENV_DIR}}"'
  echo 'export PIPELINE_RESOURCE_DIR="${PIPELINE_ROOT}/resources"'
  echo 'export PROJECT_SCENIC_DB_DIR="${RESOURCE_DIR}/scenic_db"'
  echo 'export PIPELINE_SCENIC_DB_DIR="${PIPELINE_RESOURCE_DIR}/scenic_db"'
  echo 'export VELOCITY_DIR="${RESULTS_DIR}/velocity"'
  echo 'export VELOCITY_INPUT_DIR="${VELOCITY_DIR}/input"'
  echo 'export VELOCITY_LOOM_DIR="${VELOCITY_DIR}/loom"'
  echo 'export VELOCITY_OUTPUT_DIR="${VELOCITY_DIR}/output"'
  echo 'export SCENIC_INPUT_DIR="${RESULTS_DIR}/scenic_input"'
  echo 'export SCENIC_OUTPUT_DIR="${RESULTS_DIR}/scenic_output"'
  echo 'export SCENIC_DB_DIR="${SCENIC_DB_DIR_OVERRIDE:-${PROJECT_SCENIC_DB_DIR}}"'
  echo 'export SAMPLE_SHEET="${METADATA_DIR}/samples.tsv"'
  echo 'export CANONICAL_SAMPLE_SHEET="${METADATA_DIR}/samples.canonical.tsv"'
  echo 'export COMPARISON_SHEET="${METADATA_DIR}/comparisons.tsv"'
  echo 'export COMMUNICATION_PAIRS_SHEET="${METADATA_DIR}/communication_pairs.tsv"'
  echo 'export DELIVERY_MANIFEST="${METADATA_DIR}/delivery_manifest.tsv"'
  echo 'export RECEIVED_FILES_MANIFEST="${METADATA_DIR}/received_files_manifest.tsv"'
  echo 'export INPUT_INVENTORY_FILE="${INTAKE_REPORT_DIR}/input_inventory.tsv"'
  echo 'export BRANCH_READINESS_FILE="${INTAKE_REPORT_DIR}/branch_readiness.tsv"'
  echo 'export INTAKE_SUMMARY_FILE="${INTAKE_REPORT_DIR}/intake_summary.md"'
  echo 'export WAIVER_FILE="${PROJECT_CONFIG_DIR}/waivers.tsv"'
  echo 'export QC_THRESHOLD_FILE="${PROJECT_CONFIG_DIR}/qc_thresholds.tsv"'
  echo 'export EDA_GATE_FILE="${PROJECT_CONFIG_DIR}/eda_gates.tsv"'
  echo 'export OBJECT_LAYER_CONFIG_FILE="${PROJECT_CONFIG_DIR}/object_layers.tsv"'
  echo 'export MARKER_PANEL_DIR="${PROJECT_CONFIG_DIR}/marker_panels"'
  echo 'export MITO_GENE_LIST_FILE="${PROJECT_CONFIG_DIR}/mito_gene_list.txt"'
  echo 'export MIN_BIOLOGICAL_REPLICATES="2"'
  echo 'export AMBIENT_REPORT_DIR="${EDA_REPORT_DIR}/ambient"'
  echo 'export PRE_QC_REPORT_DIR="${EDA_REPORT_DIR}/pre_qc"'
  echo 'export POST_QC_REPORT_DIR="${EDA_REPORT_DIR}/post_qc"'
  echo 'export INTEGRATION_REPORT_DIR="${EDA_REPORT_DIR}/integration"'
  echo 'export ANNOTATION_REPORT_DIR="${EDA_REPORT_DIR}/annotation"'
  echo 'export DEG_REPORT_DIR="${EDA_REPORT_DIR}/deg"'
  echo 'export ENRICHMENT_REPORT_DIR="${EDA_REPORT_DIR}/enrichment"'
  echo 'export COMMUNICATION_REPORT_DIR="${EDA_REPORT_DIR}/communication"'
  echo 'export NICHENET_RESOURCE_DIR="${RESOURCE_DIR}/nichenet"'
  echo
  echo 'export R_MAIN_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/r_main"'
  echo 'export R_SCENIC_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/r_scenic"'
  echo 'export R_LIBS_MAIN="${PIPELINE_ENV_DIR}/R_libs_main"'
  echo 'export R_LIBS_SCENIC="${PIPELINE_ENV_DIR}/R_libs_scenic"'
  echo 'export R_LIBS_INTERACTION=""'
  echo 'export VELOCITY_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/velocity"'
  echo 'export PYSCENIC_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/pyscenic"'
  echo 'export R_INTERACTION_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/r_interaction"'
  echo 'export R_SPATIAL_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/r_spatial"'
  echo 'export PY_SPATIAL_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/py_spatial"'
  echo 'export PY_SPATIAL_LEGACY_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/py_spatial_legacy"'
  echo 'export PY_CELL2LOCATION_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/py_cell2location"'
  echo 'export R_VALIDATION_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/r_validation"'
  echo
  echo 'export SCENIC_TF_LIST="${SCENIC_DB_DIR}/hs_hgnc_tfs.txt"'
  echo 'export SCENIC_MOTIF_ANN="${SCENIC_DB_DIR}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"'
  echo 'export SCENIC_DB_500BP="${SCENIC_DB_DIR}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"'
  echo 'export SCENIC_DB_10KB="${SCENIC_DB_DIR}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"'
} > "${PROJECT_CONFIG_FILE}"

if [[ -f "${PIPELINE_ROOT}/config/marker_panels/README.md" ]]; then
  cp -f "${PIPELINE_ROOT}/config/marker_panels/README.md" "${PROJECT_MARKER_PANEL_DIR}/README.md"
fi
if [[ -f "${PIPELINE_ROOT}/config/marker_panels/example_panel.tsv.template" ]]; then
  cp -f "${PIPELINE_ROOT}/config/marker_panels/example_panel.tsv.template" "${PROJECT_MARKER_PANEL_DIR}/example_panel.tsv.template"
fi
if [[ -f "${PIPELINE_ROOT}/config/mito_gene_list.txt.template" ]]; then
  cp -f "${PIPELINE_ROOT}/config/mito_gene_list.txt.template" "${PROJECT_MITO_GENE_LIST_FILE}"
fi

chmod +x "${PROJECT_CONFIG_FILE}"

IFS=',' read -r -a sample_names_array <<< "${SAMPLE_NAMES}"
IFS=',' read -r -a group1_samples_array <<< "${GROUP1_SAMPLES}"
IFS=',' read -r -a group2_samples_array <<< "${GROUP2_SAMPLES}"
IFS=',' read -r -a raw_samples_array <<< "${RAW_SAMPLES}"
IFS=',' read -r -a raw_group1_samples_array <<< "${RAW_GROUP1_SAMPLES}"
IFS=',' read -r -a raw_group2_samples_array <<< "${RAW_GROUP2_SAMPLES}"

main_sample_run_velocity="${DEFAULT_RUN_VELOCITY}"
if [[ "${#raw_samples_array[@]}" -ne "${#sample_names_array[@]}" ]]; then
  main_sample_run_velocity="no"
else
  for idx in "${!sample_names_array[@]}"; do
    if [[ "${sample_names_array[$idx]}" != "${raw_samples_array[$idx]}" ]]; then
      main_sample_run_velocity="no"
      break
    fi
  done
fi

{
  printf 'sample_id\tcondition\tbiological_replicate\ttechnical_replicate\tbatch\tinput_mode\tinput_source\tsource_path\tplatform\tgene_id_type\treference_version\ttimepoint\ttissue\tchemistry\trun_main\trun_velocity\trun_scenic\n'
  for sample in "${sample_names_array[@]}"; do
    [[ -n "${sample}" ]] || continue
    condition="unassigned"
    for group_sample in "${group1_samples_array[@]}"; do
      if [[ "${sample}" == "${group_sample}" ]]; then
        condition="${GROUP1_NAME}"
        break
      fi
    done
    if [[ "${condition}" == "unassigned" ]]; then
      for group_sample in "${group2_samples_array[@]}"; do
        if [[ "${sample}" == "${group_sample}" ]]; then
          condition="${GROUP2_NAME}"
          break
        fi
      done
    fi

    source_path=""
    platform="generic_mex"
    case "${DEFAULT_INPUT_MODE}" in
      fastq)
        source_path="${FASTQ_SOURCE_DIR}"
        platform="10x_cellranger"
        ;;
      cellranger_out)
        source_path="${CELLRANGER_OUT_SOURCE_DIR}/${sample}"
        platform="10x_cellranger"
        ;;
      matrix)
        if [[ -n "${MATRIX_SOURCE_DIR}" ]]; then
          source_path="${MATRIX_SOURCE_DIR}/${sample}"
        elif [[ -n "${PAPER_MATRIX_SOURCE}" ]]; then
          source_path="${PAPER_MATRIX_SOURCE}/${sample}"
        fi
        ;;
      *)
        ;;
    esac

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${sample}" \
      "${condition}" \
      "${sample}" \
      "1" \
      "default" \
      "${DEFAULT_INPUT_MODE}" \
      "init_project" \
      "${source_path}" \
      "${platform}" \
      "auto" \
      "" \
      "" \
      "" \
      "" \
      "yes" \
      "${main_sample_run_velocity}" \
      "yes"
  done

  for raw_sample in "${raw_samples_array[@]}"; do
    [[ -n "${raw_sample}" ]] || continue

    skip_raw_row="0"
    for sample in "${sample_names_array[@]}"; do
      if [[ "${raw_sample}" == "${sample}" ]]; then
        skip_raw_row="1"
        break
      fi
    done
    [[ "${skip_raw_row}" == "1" ]] && continue

    condition="unassigned"
    for group_sample in "${raw_group1_samples_array[@]}"; do
      if [[ "${raw_sample}" == "${group_sample}" ]]; then
        condition="${RAW_GROUP1_NAME}"
        break
      fi
    done
    if [[ "${condition}" == "unassigned" ]]; then
      for group_sample in "${raw_group2_samples_array[@]}"; do
        if [[ "${raw_sample}" == "${group_sample}" ]]; then
          condition="${RAW_GROUP2_NAME}"
          break
        fi
      done
    fi

    raw_source_path=""
    if [[ -n "${CELLRANGER_OUT_SOURCE_DIR}" ]]; then
      raw_source_path="${CELLRANGER_OUT_SOURCE_DIR}/${raw_sample}"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${raw_sample}" \
      "${condition}" \
      "${raw_sample}" \
      "1" \
      "default" \
      "cellranger_out" \
      "init_project_raw" \
      "${raw_source_path}" \
      "10x_cellranger" \
      "auto" \
      "" \
      "" \
      "" \
      "" \
      "no" \
      "yes" \
      "no"
  done
} > "${PROJECT_METADATA_DIR}/samples.tsv"

{
  printf 'comparison_id\tident_1\tident_2\tenabled\tgroup_var\tbatch_var\tlayer_scope\tmin_biological_replicates\tsubset_column\tsubset_value\tforce_exploratory\tmin_cells_per_group\tlogfc_threshold\n'
  printf '%s_vs_%s\t%s\t%s\tyes\tgroup_id\tbatch\t*\t2\t\tno\t3\t0\n' "${GROUP1_NAME}" "${GROUP2_NAME}" "${GROUP1_NAME}" "${GROUP2_NAME}"
} > "${PROJECT_METADATA_DIR}/comparisons.tsv"

{
  printf 'check_id\tscope\treason\tapproved_by\n'
} > "${PROJECT_CONFIG_DIR}/waivers.tsv"

{
  printf 'sample_id\tqc_min_nfeature\tqc_min_ncount\tqc_min_log10umi\tqc_max_mito_pct\n'
} > "${PROJECT_CONFIG_DIR}/qc_thresholds.tsv"

{
  printf 'gate_id\tstatus\tapproved_by\tnotes\n'
  printf 'pre_qc\tpending\t\t\n'
  printf 'post_qc\tpending\t\t\n'
  printf 'integration\tpending\t\t\n'
  printf 'annotation\tpending\t\t\n'
  printf 'deg\tpending\t\t\n'
} > "${PROJECT_CONFIG_DIR}/eda_gates.tsv"

{
  printf '# object_layers.tsv\n'
  printf '# - panorama 行是根层，默认从全部 post-QC 细胞建模\n'
  printf '# - subcluster 行从 parent_layer 的 metadata 中按 selection_column/selection_values 取细胞\n'
  printf '# - 默认示例保留 subcluster_1 / subcluster_2，但 selection_values 需要在审阅 panorama 后填写\n'
  printf 'layer_id\tlayer_role\tenabled\tparent_layer\tsample_include\tsample_exclude\tselection_column\tselection_values\trebuild_normalization\thvg_nfeatures\tpca_dims\ttarget_clusters\tres_range\tres_fine_step\tintegration_mode\tdescription\n'
  printf 'panorama\tpanorama\tyes\t\t\t\t\t\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "2000" \
    "1:30" \
    "15" \
    "0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60" \
    "0.005" \
    "harmony" \
    'Root panorama object built from all post-QC cells.'
  printf 'subcluster_1\tsubcluster\tyes\tpanorama\t\t\tpanorama_cluster\t\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "2000" \
    "1:20" \
    "8" \
    "0.10,0.15,0.20,0.25,0.30,0.35,0.40" \
    "0.005" \
    "harmony" \
    'Edit selection_column and selection_values to define this subcluster from panorama results.'
  printf 'subcluster_2\tsubcluster\tyes\tpanorama\t\t\tpanorama_cluster\t\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "2000" \
    "1:20" \
    "6" \
    "0.10,0.15,0.20,0.25,0.30,0.35,0.40" \
    "0.005" \
    "harmony" \
    'Edit selection_column and selection_values to define this subcluster from panorama results.'
} > "${PROJECT_CONFIG_DIR}/object_layers.tsv"

{
  printf 'delivery_id\tsource_path\tlink_target\tregistered_at\tmode\n'
} > "${PROJECT_METADATA_DIR}/delivery_manifest.tsv"

{
  printf 'delivery_id\trelative_path\tfile_type\tregistered_at\n'
} > "${PROJECT_METADATA_DIR}/received_files_manifest.tsv"

cat > "${PROJECT_STATUS_FILE}" <<EOF
{
  "updated_at": "$(date --iso-8601=seconds)",
  "project_root": "${PROJECT_ROOT}",
  "config_file": "${PROJECT_CONFIG_FILE}",
  "current_stage": "project_initialized",
  "next_step": "bash ${PIPELINE_ROOT}/workflow/01_register_delivery.sh --help",
  "project_input_mode": "${DEFAULT_INPUT_MODE}",
  "paths": {
    "sample_sheet": "${PROJECT_METADATA_DIR}/samples.tsv",
    "canonical_sample_sheet": "${PROJECT_METADATA_DIR}/samples.canonical.tsv",
    "comparison_sheet": "${PROJECT_METADATA_DIR}/comparisons.tsv",
    "object_layer_config": "${PROJECT_CONFIG_DIR}/object_layers.tsv",
    "input_inventory": "${PROJECT_INTAKE_REPORT_DIR}/input_inventory.tsv",
    "branch_readiness": "${PROJECT_INTAKE_REPORT_DIR}/branch_readiness.tsv",
    "intake_summary": "${PROJECT_INTAKE_REPORT_DIR}/intake_summary.md",
    "status_dir": "${PROJECT_STATUS_DIR}"
  },
  "annotation_hub": {
    "path": "${PROJECT_ROOT}/results/checkpoints/03_after_annotation.rds",
    "exists": false,
    "produced_by": "20_run_main_pipeline.sh / r/03_annotation.R"
  },
  "status": {
    "metadata_valid": false,
    "standardized_inputs": false,
    "input_eda_complete": false,
    "ambient_branch_complete": false,
    "pre_qc_eda_complete": false,
    "pre_qc_gate_passed": false,
    "post_qc_eda_complete": false,
    "post_qc_gate_passed": false,
    "integration_eda_complete": false,
    "integration_gate_passed": false,
    "annotation_eda_complete": false,
    "annotation_gate_passed": false,
    "deg_gate_passed": false,
    "05_deg_completed": false,
    "main_upstream_ready": false,
    "velocity_upstream_ready": false,
    "ambient_upstream_ready": false,
    "scenic_upstream_ready": false,
    "de_replicate_ready": false,
    "main_ready": false,
    "deg_ready": false,
    "enrichment_ready": false,
    "trajectory_ready": false,
    "velocity_reference_ready": false,
    "scenic_export_ready": false
  },
  "waivers": []
}
EOF

mkdir -p "${PIPELINE_CONFIG_DIR}"
{
  echo '#!/usr/bin/env bash'
  echo
  echo '# 由 workflow/00_init_project.sh 自动更新。'
  echo '# 当前全局激活项目配置如下：'
  printf 'source %q\n' "${PROJECT_CONFIG_FILE}"
} > "${GLOBAL_CONFIG}"

chmod +x "${GLOBAL_CONFIG}"

{
  echo '#!/usr/bin/env bash'
  echo
  echo '# 兼容入口：真实配置已经迁移到 config/server_config.sh。'
  printf 'source %q\n' "${GLOBAL_CONFIG}"
} > "${ROOT_CONFIG_SHIM}"

chmod +x "${ROOT_CONFIG_SHIM}"

echo "项目框架已创建: ${PROJECT_ROOT}"
echo "项目配置已写入: ${PROJECT_CONFIG_FILE}"
echo "全局配置已切换到: ${PROJECT_CONFIG_FILE}"
echo
echo "下一步建议执行："
echo "1. 检查配置与 metadata: ${PROJECT_CONFIG_FILE} / ${PROJECT_METADATA_DIR}/samples.tsv"
echo "2. 登记交付物: bash ${PIPELINE_ROOT}/workflow/01_register_delivery.sh --help"
echo "3. 校验 metadata: bash ${PIPELINE_ROOT}/workflow/02_validate_metadata.sh"
echo "4. 审计输入: bash ${PIPELINE_ROOT}/workflow/03_audit_inputs.sh"
echo "5. 从 FASTQ 开始时先跑: bash ${PIPELINE_ROOT}/workflow/10_run_cellranger_from_fastq.sh"
echo "6. 标准化主流程输入: bash ${PIPELINE_ROOT}/workflow/04_standardize_inputs.sh"
echo "7. 生成 intake summary: bash ${PIPELINE_ROOT}/workflow/05_input_summary.sh"
echo "8. 安装环境: bash ${PIPELINE_ROOT}/workflow/01_install_envs.sh"
echo "9. 跑主流程: bash ${PIPELINE_ROOT}/workflow/20_run_main_pipeline.sh"
