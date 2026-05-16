#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/spatial_intake_smoke.XXXXXX")"
trap 'rm -rf "${TMP_ROOT}"' EXIT

if ! command -v Rscript >/dev/null 2>&1; then
  echo "skip: Rscript is not available"
  exit 0
fi

if ! Rscript -e 'quit(status = if (requireNamespace("Seurat", quietly = TRUE) && requireNamespace("Matrix", quietly = TRUE)) 0 else 42)' >/dev/null 2>&1; then
  echo "skip: Seurat and Matrix are required for spatial intake smoke"
  exit 0
fi

mkdir -p "${TMP_ROOT}/refs" "${TMP_ROOT}/matrix/st_syf_1" "${TMP_ROOT}/project/config" \
  "${TMP_ROOT}/project/metadata" "${TMP_ROOT}/project/reports/intake"

cat > "${TMP_ROOT}/refs/genome.fa" <<'EOF'
>chr1
ACGT
EOF

cat > "${TMP_ROOT}/refs/genes.gtf" <<'EOF'
chrM	source	gene	1	4	.	+	.	gene_id "g1"; gene_name "MT-CO1";
EOF

Rscript - <<RSCRIPT
if (!requireNamespace("Matrix", quietly = TRUE)) quit(status = 42)
root <- "${TMP_ROOT}/matrix/st_syf_1"
counts <- Matrix::Matrix(c(5, 0, 2, 3, 4, 0), nrow = 3, sparse = TRUE)
Matrix::writeMM(counts, file.path(root, "matrix.mtx"))
write.table(data.frame("g1", "MT-CO1", "Gene Expression"), file.path(root, "features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write("spot1\nspot2", file.path(root, "barcodes.tsv"))
write.csv(data.frame(barcode = c("spot1", "spot2"), row = c(1, 2), col = c(1, 2), in_tissue = c(TRUE, TRUE)), file.path(root, "coords.csv"), row.names = FALSE)
RSCRIPT

cat > "${TMP_ROOT}/project/config/project_config.sh" <<EOF
export PROJECT_ROOT="${TMP_ROOT}/project"
export PIPELINE_ROOT="${PIPELINE_ROOT}"
export DATA_DIR="${TMP_ROOT}/project/data"
export RESULTS_DIR="${TMP_ROOT}/project/results"
export REPORT_DIR="${TMP_ROOT}/project/reports"
export EDA_REPORT_DIR="${TMP_ROOT}/project/reports/eda"
export INTAKE_REPORT_DIR="${TMP_ROOT}/project/reports/intake"
export MANIFEST_DIR="${TMP_ROOT}/project/results/manifests"
export METADATA_DIR="${TMP_ROOT}/project/metadata"
export PROJECT_CONFIG_DIR="${TMP_ROOT}/project/config"
export SAMPLE_SHEET="${TMP_ROOT}/project/metadata/samples.tsv"
export CANONICAL_SAMPLE_SHEET="${TMP_ROOT}/project/metadata/samples.canonical.tsv"
export SECTION_SHEET="${TMP_ROOT}/project/metadata/sections.tsv"
export SPATIAL_INPUT_INVENTORY_FILE="${TMP_ROOT}/project/reports/intake/spatial_input_inventory.tsv"
export SPATIAL_INTAKE_CONTRACT_FILE="${TMP_ROOT}/project/config/spatial_intake_contract.tsv"
export SPATIAL_QC_THRESHOLD_FILE="${TMP_ROOT}/project/config/spatial_qc_thresholds.tsv"
export EDA_GATE_FILE="${TMP_ROOT}/project/config/eda_gates.tsv"
export SPATIAL_RESULTS_DIR="${TMP_ROOT}/project/results/spatial"
export SPATIAL_CHECKPOINT_DIR="${TMP_ROOT}/project/results/spatial/checkpoints"
export SPATIAL_TABLE_DIR="${TMP_ROOT}/project/results/spatial/tables"
export SPATIAL_FIGURE_DIR="${TMP_ROOT}/project/results/spatial/figures"
export SPATIAL_PRE_QC_REPORT_DIR="${TMP_ROOT}/project/reports/eda/spatial_pre_qc"
export CLEAN_GTF="${TMP_ROOT}/refs/genes.gtf"
EOF

cat > "${TMP_ROOT}/project/metadata/samples.canonical.tsv" <<EOF
sample_id	condition	biological_replicate	technical_replicate	batch	input_mode	input_source	source_path	platform	gene_id_type	reference_version	group_id	timepoint	tissue	chemistry	run_main	run_velocity	run_scenic	modality	section_id	chip_id	bundle_layout	image_path	run_spatial	run_deconv	run_joint	notes
st_syf_1	syf	syf_section_1	1	st_batch1	spatial_matrix	smoke	${TMP_ROOT}/matrix/st_syf_1	generic	auto	custom_reference	syf	D0	chicken_ovary	auto	no	no	no	spatial	syf_1	syf_chip_1	generic_spatial_matrix		yes	auto	auto
EOF
cp "${TMP_ROOT}/project/metadata/samples.canonical.tsv" "${TMP_ROOT}/project/metadata/samples.tsv"

cat > "${TMP_ROOT}/project/metadata/sections.tsv" <<EOF
section_id	chip_id	tissue_block	condition	platform	bundle_root	image_lowres	image_hires	tissue_positions	scalefactors	enabled	notes
syf_1	syf_chip_1	block	syf	generic	${TMP_ROOT}/matrix/st_syf_1					yes	smoke
EOF

cat > "${TMP_ROOT}/project/reports/intake/spatial_input_inventory.tsv" <<EOF
sample_id	section_id	modality	condition	platform	bundle_layout	source_path	section_bundle_root	standardized_outs	ready	ready_status	required_relpaths	missing_required_relpaths	optional_present_relpaths	optional_missing_relpaths	notes
st_syf_1	syf_1	spatial	syf	generic	generic_spatial_matrix	${TMP_ROOT}/matrix/st_syf_1	${TMP_ROOT}/matrix/st_syf_1	${TMP_ROOT}/matrix/st_syf_1	true	true					smoke
EOF

cat > "${TMP_ROOT}/project/config/spatial_intake_contract.tsv" <<'EOF'
platform	bundle_layout	required_relpaths	optional_relpaths	standardized_outs_relpath	ready_status_when_present	notes
generic	generic_spatial_matrix	matrix.mtx,features.tsv,barcodes.tsv	coords.csv	.	true	Smoke generic spatial matrix.
EOF
cp "${PIPELINE_ROOT}/config/spatial_qc_thresholds.tsv.template" "${TMP_ROOT}/project/config/spatial_qc_thresholds.tsv"
cp "${PIPELINE_ROOT}/config/eda_gates.tsv.template" "${TMP_ROOT}/project/config/eda_gates.tsv"

export SCRNA_PIPELINE_CONFIG="${TMP_ROOT}/project/config/project_config.sh"
set -a
source "${SCRNA_PIPELINE_CONFIG}"
set +a

Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/testing/smoke_spatial_loader.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/01_build_spatial_objects.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/01a_pre_spot_qc_eda.R"

[[ -s "${TMP_ROOT}/project/results/spatial/checkpoints/st_syf_1_raw.rds" ]]
[[ -s "${TMP_ROOT}/project/reports/eda/spatial_pre_qc/report.md" ]]
[[ -s "${TMP_ROOT}/project/results/spatial/tables/spatial_pre_qc/triage_signals.tsv" ]]
grep -q $'^spatial_pre_qc\tpending' "${TMP_ROOT}/project/config/eda_gates.tsv"

echo "smoke_spatial_intake_ok"
