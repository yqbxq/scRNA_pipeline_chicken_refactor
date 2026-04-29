#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/scrna_intake_smoke.XXXXXX")"
trap 'rm -rf "${TMP_ROOT}"' EXIT

PROJECT_ROOT="${TMP_ROOT}/project"
FASTQ_DIR="${TMP_ROOT}/fastq"
CELLRANGER_OUT_DIR="${TMP_ROOT}/cellranger_out"
DNB_SOURCE_DIR="${TMP_ROOT}/dnbelab"
NO_RAW_SOURCE_DIR="${TMP_ROOT}/generic_no_raw"
REFERENCE_DIR="${TMP_ROOT}/reference"
METADATA_DIR="${PROJECT_ROOT}/metadata"
REPORT_DIR="${PROJECT_ROOT}/reports"
INTAKE_REPORT_DIR="${REPORT_DIR}/intake"
STATUS_DIR="${PROJECT_ROOT}/status"
DATA_DIR="${PROJECT_ROOT}/data"
RESULTS_DIR="${PROJECT_ROOT}/results"
CHECKPOINT_DIR="${RESULTS_DIR}/checkpoints"
FIGURE_DIR="${RESULTS_DIR}/figures"
TABLE_DIR="${RESULTS_DIR}/tables"
LOG_DIR="${PROJECT_ROOT}/logs"
CONFIG_FILE="${TMP_ROOT}/config.sh"

mkdir -p \
  "${FASTQ_DIR}" \
  "${CELLRANGER_OUT_DIR}/CR_S1/outs" \
  "${DNB_SOURCE_DIR}/DNB_S1" \
  "${NO_RAW_SOURCE_DIR}/NO_RAW_S1" \
  "${REFERENCE_DIR}" \
  "${METADATA_DIR}" \
  "${INTAKE_REPORT_DIR}" \
  "${STATUS_DIR}" \
  "${DATA_DIR}" \
  "${CHECKPOINT_DIR}" \
  "${FIGURE_DIR}" \
  "${TABLE_DIR}" \
  "${LOG_DIR}" \
  "${PROJECT_ROOT}/config"

create_mex_dir() {
  local target_dir="$1"
  local feature_1="$2"
  local feature_2="$3"
  mkdir -p "${target_dir}"
  cat > "${target_dir}/matrix.mtx" <<'EOF'
%%MatrixMarket matrix coordinate integer general
%
2 2 4
1 1 3
1 2 1
2 1 2
2 2 4
EOF
  gzip -f "${target_dir}/matrix.mtx"
  cat > "${target_dir}/barcodes.tsv" <<'EOF'
cellA
cellB
EOF
  gzip -f "${target_dir}/barcodes.tsv"
  cat > "${target_dir}/features.tsv" <<EOF
${feature_1}	${feature_1}	Gene Expression
${feature_2}	${feature_2}	Gene Expression
EOF
  gzip -f "${target_dir}/features.tsv"
}

create_mex_dir "${CELLRANGER_OUT_DIR}/CR_S1/outs/filtered_feature_bc_matrix" "MT-CO1" "RPLP0"
create_mex_dir "${CELLRANGER_OUT_DIR}/CR_S1/outs/raw_feature_bc_matrix" "MT-CO1" "RPLP0"
create_mex_dir "${DNB_SOURCE_DIR}/DNB_S1/filter_matrix" "ENSGALG00000000001" "ENSGALG00000000002"
create_mex_dir "${DNB_SOURCE_DIR}/DNB_S1/raw_matrix" "ENSGALG00000000001" "ENSGALG00000000002"
create_mex_dir "${NO_RAW_SOURCE_DIR}/NO_RAW_S1/filtered_feature_bc_matrix" "GENE1" "GENE2"

cat > "${CELLRANGER_OUT_DIR}/CR_S1/outs/metrics_summary.csv" <<'EOF'
Estimated Number of Cells,Mean Reads per Cell,Median Genes per Cell,Sequencing Saturation,Fraction Reads in Cells
1000,25000,900,68%,72%
EOF
touch "${CELLRANGER_OUT_DIR}/CR_S1/outs/possorted_genome_bam.bam"

cat > "${DNB_SOURCE_DIR}/DNB_S1/metrics_summary.xls" <<'EOF'
Estimated Number of Cells	Mean Reads per Cell	Median Genes per Cell	Sequencing Saturation	Fraction Reads in Cells
800	18000	700	55%	61%
EOF
touch "${DNB_SOURCE_DIR}/DNB_S1/anno_decon_sorted.bam"

cat > "${METADATA_DIR}/samples.tsv" <<EOF
sample_id	condition	biological_replicate	technical_replicate	batch	input_mode	input_source	source_path	platform	gene_id_type	reference_version	timepoint	tissue	chemistry	run_main	run_velocity	run_scenic
CR_S1	Ctrl	Ctrl_1	1	batchA	cellranger_out	smoke	${CELLRANGER_OUT_DIR}/CR_S1	10x_cellranger	auto	GRCg7b	E12	ovary	v3	yes	yes	yes
DNB_S1	Treat	Treat_1	1	batchB	matrix	smoke	${DNB_SOURCE_DIR}/DNB_S1	dnbelab_c	auto	GRCg7b	E12	ovary	DNB	yes	no	yes
NO_RAW_S1	Treat	Treat_2	1	batchC	matrix	smoke	${NO_RAW_SOURCE_DIR}/NO_RAW_S1	generic_mex	symbol	GRCg7b	E12	ovary	v3	yes	no	yes
EOF

cat > "${METADATA_DIR}/comparisons.tsv" <<'EOF'
comparison_id	ident_1	ident_2	enabled
Ctrl_vs_Treat	Ctrl	Treat	yes
EOF

cat > "${PROJECT_ROOT}/config/waivers.tsv" <<'EOF'
check_id	scope	reason	approved_by
EOF

cat > "${CONFIG_FILE}" <<EOF
#!/usr/bin/env bash
export PROJECT_ROOT="${PROJECT_ROOT}"
export FASTQ_DIR="${FASTQ_DIR}"
export CELLRANGER_OUT_DIR="${CELLRANGER_OUT_DIR}"
export REFERENCE_DIR="${REFERENCE_DIR}"
export SRA_TOOLS_DIR=""
export CELLRANGER_BIN=""
export PAPER_MATRIX_SOURCE=""
export REFERENCE_GTF="${REFERENCE_DIR}/dummy.gtf"
export CLEAN_GTF=""
export GENOME_FASTA_GZ="${REFERENCE_DIR}/dummy.fa"
export CELLRANGER_REF_DIR="${REFERENCE_DIR}/dummy_ref"
export SAMPLE_NAMES="CR_S1,DNB_S1,NO_RAW_S1"
export ANALYSIS_GROUP_1_NAME="Ctrl"
export ANALYSIS_GROUP_1_SAMPLES="CR_S1"
export ANALYSIS_GROUP_2_NAME="Treat"
export ANALYSIS_GROUP_2_SAMPLES="DNB_S1,NO_RAW_S1"
export RAW_GROUP_1_NAME="Ctrl"
export RAW_GROUP_1_SAMPLES="CR_S1"
export RAW_GROUP_2_NAME="Treat"
export RAW_GROUP_2_SAMPLES="DNB_S1"
export RAW_SAMPLES="CR_S1,DNB_S1"
export DEG_IDENT_1="Ctrl"
export DEG_IDENT_2="Treat"
export RANDOM_SEED="42"
export INTEGRATION_MODE="harmony"
export MAIN_THREADS="1"
export PROJECT_CONFIG_DIR="${PROJECT_ROOT}/config"
export DATA_DIR="${DATA_DIR}"
export RESULTS_DIR="${RESULTS_DIR}"
export CHECKPOINT_DIR="${CHECKPOINT_DIR}"
export FIGURE_DIR="${FIGURE_DIR}"
export TABLE_DIR="${TABLE_DIR}"
export LOG_DIR="${LOG_DIR}"
export RAW_DIR="${PROJECT_ROOT}/raw"
export RAW_RECEIVED_DIR="${PROJECT_ROOT}/raw/received"
export METADATA_DIR="${METADATA_DIR}"
export REPORT_DIR="${REPORT_DIR}"
export INTAKE_REPORT_DIR="${INTAKE_REPORT_DIR}"
export EDA_REPORT_DIR="${REPORT_DIR}/eda"
export STATUS_DIR="${STATUS_DIR}"
export WORKFLOW_STATUS_FILE="${STATUS_DIR}/workflow_status.json"
export SAMPLE_SHEET="${METADATA_DIR}/samples.tsv"
export CANONICAL_SAMPLE_SHEET="${METADATA_DIR}/samples.canonical.tsv"
export COMPARISON_SHEET="${METADATA_DIR}/comparisons.tsv"
export DELIVERY_MANIFEST="${METADATA_DIR}/delivery_manifest.tsv"
export RECEIVED_FILES_MANIFEST="${METADATA_DIR}/received_files_manifest.tsv"
export INPUT_INVENTORY_FILE="${INTAKE_REPORT_DIR}/input_inventory.tsv"
export BRANCH_READINESS_FILE="${INTAKE_REPORT_DIR}/branch_readiness.tsv"
export INTAKE_SUMMARY_FILE="${INTAKE_REPORT_DIR}/intake_summary.md"
export WAIVER_FILE="${PROJECT_ROOT}/config/waivers.tsv"
export QC_THRESHOLD_FILE="${PROJECT_ROOT}/config/qc_thresholds.tsv"
export EDA_GATE_FILE="${PROJECT_ROOT}/config/eda_gates.tsv"
export PRE_QC_REPORT_DIR="${REPORT_DIR}/eda/pre_qc"
export POST_QC_REPORT_DIR="${REPORT_DIR}/eda/post_qc"
export INTEGRATION_REPORT_DIR="${REPORT_DIR}/eda/integration"
export ANNOTATION_REPORT_DIR="${REPORT_DIR}/eda/annotation"
export INPUT_STANDARDIZE_MODE="symlink"
export USE_EXISTING_SIF="no"
EOF

SCRNA_PIPELINE_CONFIG="${CONFIG_FILE}" bash "${PIPELINE_ROOT}/workflow/03stages/validate_metadata.sh"
SCRNA_PIPELINE_CONFIG="${CONFIG_FILE}" bash "${PIPELINE_ROOT}/workflow/03stages/audit_inputs.sh"
SCRNA_PIPELINE_CONFIG="${CONFIG_FILE}" bash "${PIPELINE_ROOT}/workflow/03stages/standardize_inputs.sh"
SCRNA_PIPELINE_CONFIG="${CONFIG_FILE}" bash "${PIPELINE_ROOT}/workflow/03stages/input_summary.sh"

python3 - <<PY
import csv
from pathlib import Path

inventory = Path("${INTAKE_REPORT_DIR}/input_inventory.tsv")
readiness = Path("${INTAKE_REPORT_DIR}/branch_readiness.tsv")
summary = Path("${INTAKE_REPORT_DIR}/intake_summary.md")

with inventory.open("r", encoding="utf-8", newline="") as handle:
    rows = {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

cr = rows["CR_S1"]
dnb = rows["DNB_S1"]
no_raw = rows["NO_RAW_S1"]

assert cr["platform_resolved"] == "10x_cellranger", cr
assert cr["filtered_matrix_dir"].endswith("outs/filtered_feature_bc_matrix"), cr
assert "raw_feature_bc_matrix" in cr["raw_matrix_dir"], cr
assert cr["metrics_path"].endswith("metrics_summary.csv"), cr
assert cr["bam_path"].endswith("possorted_genome_bam.bam"), cr

assert dnb["platform_resolved"] == "dnbelab_c", dnb
assert dnb["filtered_matrix_dir"].endswith("filter_matrix"), dnb
assert dnb["raw_matrix_dir"].endswith("raw_matrix"), dnb
assert dnb["metrics_path"].endswith("metrics_summary.xls"), dnb
assert dnb["bam_path"].endswith("anno_decon_sorted.bam"), dnb

assert no_raw["platform_resolved"] == "generic_mex", no_raw
assert no_raw["filtered_matrix_dir"].endswith("filtered_feature_bc_matrix"), no_raw
assert no_raw["raw_matrix_dir"] == "", no_raw

with readiness.open("r", encoding="utf-8", newline="") as handle:
    readiness_rows = {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

assert readiness_rows["CR_S1"]["ambient_upstream_ready"] == "true", readiness_rows["CR_S1"]
assert readiness_rows["DNB_S1"]["ambient_upstream_ready"] == "true", readiness_rows["DNB_S1"]
assert readiness_rows["CR_S1"]["ambient_soupx_ready"] == "true", readiness_rows["CR_S1"]
assert readiness_rows["CR_S1"]["ambient_preferred_method"] == "soupx", readiness_rows["CR_S1"]
assert readiness_rows["DNB_S1"]["ambient_soupx_ready"] == "true", readiness_rows["DNB_S1"]
assert readiness_rows["NO_RAW_S1"]["ambient_soupx_ready"] == "false", readiness_rows["NO_RAW_S1"]
assert readiness_rows["NO_RAW_S1"]["ambient_decontx_ready"] == "true", readiness_rows["NO_RAW_S1"]
assert readiness_rows["NO_RAW_S1"]["ambient_preferred_method"] == "decontx", readiness_rows["NO_RAW_S1"]
assert "soupx_requires_raw_matrix" in readiness_rows["NO_RAW_S1"]["ambient_notes"], readiness_rows["NO_RAW_S1"]

text = summary.read_text(encoding="utf-8")
assert "filtered_matrix_dir" in text, text
assert "raw_matrix_dir" in text, text
assert "ambient: preferred=" in text, text
assert "CR_S1" in text and "DNB_S1" in text and "NO_RAW_S1" in text, text
PY

test -e "${DATA_DIR}/CR_S1"
test -e "${DATA_DIR}/DNB_S1"
test -e "${DATA_DIR}/NO_RAW_S1"

echo "smoke_intake_contract: PASS"
