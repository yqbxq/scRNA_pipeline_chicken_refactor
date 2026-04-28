#!/usr/bin/env bash

# 注意：
# 这个文件当前是“激活中的项目配置”。
# 在论文验证示例中，SAMPLE_NAMES=HF_GCs,PHF_GCs，而 RAW_SAMPLES=5 个 SRR；
# 这是因为作者直接提供了聚合后的 matrix，但 RNA velocity 仍然依赖原始 Cell Ranger 样本。
# 正式项目模式下，通常应先运行 workflow/00_init_project.sh 生成新项目配置，
# 此时 SAMPLE_NAMES 与 RAW_SAMPLES 一般保持一致，都是你自己的样本名。

# =============================== [TUNABLE] ===============================
# Linux 服务器上的正式项目根目录。
export PROJECT_ROOT="/home/user_test/scRNA_pipeline_chicken"

# 论文验证数据来源。
# 只有在“论文验证模式”下才需要用到。
export PAPER_MATRIX_SOURCE="/home/user_test/shared_resources/paper_inputs/chicken_sichuan/data"

# 旧复现工程中的资源。
# 当前仅作为可复用参考，不再作为默认运行时。
export EXISTING_R_SIF=""
export CELLRANGER_OUT_DIR="/home/user_test/shared_resources/paper_inputs/chicken_sichuan/cellranger_out"
export REFERENCE_DIR="/home/user_test/shared_resources/reference/chicken"

# 服务器上已有的软件路径。
export CELLRANGER_BIN="/Bio/bin/cel1ranger"
export SRA_TOOLS_DIR="/conda_envs/softwares/sratoolkit.3.0.7-ubuntu64/bin"

# 从原始 FASTQ 开始时需要的输入。
export FASTQ_DIR="/home/user_test/shared_resources/paper_inputs/chicken_sichuan/fastq"
export GENOME_FASTA_GZ="${REFERENCE_DIR}/genome.fa"
export REFERENCE_GTF="${REFERENCE_DIR}/GRCg7b_genomic.gtf"
export CLEAN_GTF="${REFERENCE_DIR}/GRCg7b_genomic_clean.gtf"
export CELLRANGER_REF_DIR="${REFERENCE_DIR}/GRCg7b"

# 论文验证数据的 raw 样本分组。
# 换成你自己的项目时，把名字和样本列表改掉即可。
export RAW_GROUP_1_NAME="HF_GCs"
export RAW_GROUP_1_SAMPLES="SRR32621893,SRR32621894,SRR32621895"
export RAW_GROUP_2_NAME="PHF_GCs"
export RAW_GROUP_2_SAMPLES="SRR32621891,SRR32621892"
export RAW_SAMPLES="SRR32621891,SRR32621892,SRR32621893,SRR32621894,SRR32621895"

# 运行模式开关。
# 正式服务器环境默认不再复用旧 sif。
export USE_EXISTING_SIF="no"
export AUTHOR_INPUT_MODE="symlink"
export INPUT_STANDARDIZE_MODE="symlink"
export RANDOM_SEED="42"
export INTEGRATION_MODE="harmony"

# 计算资源参数。
export MAIN_THREADS="16"
export VELOCYTO_THREADS="8"
export SCVELO_THREADS="16"
export SCENIC_THREADS="16"
export CELLRANGER_THREADS="96"
export CELLRANGER_MEM_GB="128"

# 主流程样本名。
# 论文验证模式下当前默认是 HF_GCs,PHF_GCs。
# 以后换项目时，这里改成你自己的 matrix 目录名。
export SAMPLE_NAMES="HF_GCs,PHF_GCs"

# 主流程分组定义。
# 论文验证模式下成员就是 HF_GCs / PHF_GCs 两个矩阵目录名；
# 正式项目模式下成员应写成你的样本目录名，例如 Ctrl_1,Ctrl_2 / Treat_1,Treat_2。
export ANALYSIS_GROUP_1_NAME="HF_GCs"
export ANALYSIS_GROUP_1_SAMPLES="HF_GCs"
export ANALYSIS_GROUP_2_NAME="PHF_GCs"
export ANALYSIS_GROUP_2_SAMPLES="PHF_GCs"

# 差异分析对比组。
# 这里应填写上面 ANALYSIS_GROUP_*_NAME 的组名。
export DEG_IDENT_1="HF_GCs"
export DEG_IDENT_2="PHF_GCs"
export MIN_BIOLOGICAL_REPLICATES="2"

# 主流程参数。
export QC_MIN_NFEATURE="200"
export QC_MIN_NCOUNT="1000"
export QC_MIN_LOG10UMI="0.7"
export QC_MAX_MITO_PCT="20"
export AMBIENT_PRIMARY_METHOD="soupx"
export AMBIENT_FALLBACK_METHOD="decontx"
export AMBIENT_APPLY_POLICY="manual"
export AMBIENT_MIN_CELLS="50"
export AMBIENT_CLUSTER_DIMS="1:20"
export AMBIENT_CLUSTER_RESOLUTION="0.4"
export AMBIENT_MARKER_TOP_N="3"
export AMBIENT_RECOMMEND_MIN_CONTAMINATION="0.05"
export CELLBENDER_MODE="stub"
export CELLBENDER_FPR="0.01"
export CELLBENDER_CUDA="yes"
export CELLBENDER_EXTRA_ARGS=""
export DOUBLET_RATE="0.008"
export DOUBLET_RATE_PER_1K="0.008"
export DOUBLET_PRIMARY_CALLER="scDblFinder"
export DOUBLET_SECONDARY_CALLER="DoubletFinder"
export DOUBLET_SECONDARY_ENABLED="yes"
export DOUBLET_MIN_CELLS="50"
export DOUBLET_DIMS="1:20"
export HVG_NFEATURES="2000"
export PCA_DIMS="1:30"
export TARGET_CLUSTERS="15"
export RES_RANGE="0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60"
export RES_FINE_STEP="0.005"
export TRAJECTORY_START=""
export TRAJECTORY_COARSE_LABEL="cell_type"
export TRAJECTORY_FINE_LABEL="seurat_clusters"
export TRAJECTORY_FINE_START_CLUSTER=""
export TRAJECTORY_FINE_TOP_N="4"
export TRAJECTORY_TRADESEQ_LABEL="cell_type"
export TRADESEQ_KNOTS="6"

# 外部镜像与资源。
export ENSEMBL_MIRROR="asia"
export SCENIC_ORTHOLOG_MODE="strict"
export SCENIC_ORTHOLOG_MAP_FILE=""

# ================================ [FIXED] ================================
export PROJECT_CONFIG_DIR="${PROJECT_ROOT}/config"
export DATA_DIR="${PROJECT_ROOT}/data"
export RESULTS_DIR="${PROJECT_ROOT}/results"
export ORTHOLOG_CACHE_DIR="${RESULTS_DIR}/ortholog_cache"
export CHECKPOINT_DIR="${RESULTS_DIR}/checkpoints"
export FIGURE_DIR="${RESULTS_DIR}/figures"
export TABLE_DIR="${RESULTS_DIR}/tables"
export LOG_DIR="${PROJECT_ROOT}/logs"
export ENV_DIR="${PROJECT_ROOT}/envs"
export RESOURCE_DIR="${PROJECT_ROOT}/resources"
export RAW_DIR="${PROJECT_ROOT}/raw"
export RAW_RECEIVED_DIR="${RAW_DIR}/received"
export METADATA_DIR="${PROJECT_ROOT}/metadata"
export REPORT_DIR="${PROJECT_ROOT}/reports"
export INTAKE_REPORT_DIR="${REPORT_DIR}/intake"
export EDA_REPORT_DIR="${REPORT_DIR}/eda"
export STATUS_DIR="${PROJECT_ROOT}/status"
export WORKFLOW_STATUS_FILE="${STATUS_DIR}/workflow_status.json"
export SHARED_ENV_DIR="${SHARED_ENV_DIR:-/home/user_test/syf_f5/01shared_resources/envs/scRNA_pipeline_chicken}"
export PIPELINE_ENV_DIR="${PIPELINE_ENV_DIR:-${SHARED_ENV_DIR}}"
export PIPELINE_RESOURCE_DIR="${PIPELINE_ROOT}/resources"
export PROJECT_SCENIC_DB_DIR="${RESOURCE_DIR}/scenic_db"
export PIPELINE_SCENIC_DB_DIR="${PIPELINE_RESOURCE_DIR}/scenic_db"
export VELOCITY_DIR="${RESULTS_DIR}/velocity"
export VELOCITY_INPUT_DIR="${VELOCITY_DIR}/input"
export VELOCITY_LOOM_DIR="${VELOCITY_DIR}/loom"
export VELOCITY_OUTPUT_DIR="${VELOCITY_DIR}/output"
export SCENIC_INPUT_DIR="${RESULTS_DIR}/scenic_input"
export SCENIC_OUTPUT_DIR="${RESULTS_DIR}/scenic_output"
export SCENIC_DB_DIR="${SCENIC_DB_DIR_OVERRIDE:-${PROJECT_SCENIC_DB_DIR}}"
export SAMPLE_SHEET="${METADATA_DIR}/samples.tsv"
export CANONICAL_SAMPLE_SHEET="${METADATA_DIR}/samples.canonical.tsv"
export COMPARISON_SHEET="${METADATA_DIR}/comparisons.tsv"
export COMMUNICATION_PAIRS_SHEET="${METADATA_DIR}/communication_pairs.tsv"
export DELIVERY_MANIFEST="${METADATA_DIR}/delivery_manifest.tsv"
export RECEIVED_FILES_MANIFEST="${METADATA_DIR}/received_files_manifest.tsv"
export INPUT_INVENTORY_FILE="${INTAKE_REPORT_DIR}/input_inventory.tsv"
export BRANCH_READINESS_FILE="${INTAKE_REPORT_DIR}/branch_readiness.tsv"
export INTAKE_SUMMARY_FILE="${INTAKE_REPORT_DIR}/intake_summary.md"
export WAIVER_FILE="${PROJECT_CONFIG_DIR}/waivers.tsv"
export QC_THRESHOLD_FILE="${PROJECT_CONFIG_DIR}/qc_thresholds.tsv"
export EDA_GATE_FILE="${PROJECT_CONFIG_DIR}/eda_gates.tsv"
export OBJECT_LAYER_CONFIG_FILE="${PROJECT_CONFIG_DIR}/object_layers.tsv"
export MARKER_PANEL_DIR="${PROJECT_CONFIG_DIR}/marker_panels"
export AMBIENT_REPORT_DIR="${EDA_REPORT_DIR}/ambient"
export PRE_QC_REPORT_DIR="${EDA_REPORT_DIR}/pre_qc"
export POST_QC_REPORT_DIR="${EDA_REPORT_DIR}/post_qc"
export INTEGRATION_REPORT_DIR="${EDA_REPORT_DIR}/integration"
export ANNOTATION_REPORT_DIR="${EDA_REPORT_DIR}/annotation"
export SUBCLUSTER_REPORT_DIR="${EDA_REPORT_DIR}/subcluster"
export DEG_REPORT_DIR="${EDA_REPORT_DIR}/deg"
export ENRICHMENT_REPORT_DIR="${EDA_REPORT_DIR}/enrichment"
export COMMUNICATION_REPORT_DIR="${EDA_REPORT_DIR}/communication"
export NICHENET_RESOURCE_DIR="${RESOURCE_DIR}/nichenet"

# 共享的正式服务器环境路径。
export R_MAIN_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/r_main"
export R_SCENIC_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/r_scenic"
export R_LIBS_MAIN="${PIPELINE_ENV_DIR}/R_libs_main"
export R_LIBS_SCENIC="${PIPELINE_ENV_DIR}/R_libs_scenic"
export R_LIBS_INTERACTION="${R_LIBS_INTERACTION:-}"
export VELOCITY_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/velocity"
export PYSCENIC_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/pyscenic"
export R_LEGACY_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/r_heart_legacy"
export R_INTERACTION_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/r_interaction"
export R_SPATIAL_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/r_spatial"
export PY_SPATIAL_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/py_spatial"
export PY_SPATIAL_LEGACY_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/py_spatial_legacy"
export PY_CELL2LOCATION_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/py_cell2location"
export R_VALIDATION_ENV_PREFIX="${PIPELINE_ENV_DIR}/conda/r_validation"

export SCENIC_TF_LIST="${SCENIC_DB_DIR}/hs_hgnc_tfs.txt"
export SCENIC_MOTIF_ANN="${SCENIC_DB_DIR}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
export SCENIC_DB_500BP="${SCENIC_DB_DIR}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
export SCENIC_DB_10KB="${SCENIC_DB_DIR}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
