#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/spatial_normalize_smoke.XXXXXX")"
trap 'rm -rf "${TMP_ROOT}"' EXIT

if ! command -v Rscript >/dev/null 2>&1; then
  echo "skip: Rscript is not available"
  exit 0
fi

if ! Rscript -e 'quit(status = if (requireNamespace("Seurat", quietly = TRUE) && requireNamespace("Matrix", quietly = TRUE) && requireNamespace("sctransform", quietly = TRUE)) 0 else 42)' >/dev/null 2>&1; then
  echo "skip: Seurat, Matrix, and sctransform are required for spatial normalize smoke"
  exit 0
fi

mkdir -p "${TMP_ROOT}/refs" "${TMP_ROOT}/matrix" "${TMP_ROOT}/project/config" \
  "${TMP_ROOT}/project/metadata" "${TMP_ROOT}/project/reports/intake"

cat > "${TMP_ROOT}/refs/genes.gtf" <<'EOF'
chrM	source	gene	1	4	.	+	.	gene_id "MT-CO1"; gene_name "MT-CO1";
EOF

Rscript - <<RSCRIPT
if (!requireNamespace("Matrix", quietly = TRUE)) quit(status = 42)
set.seed(11)
root <- "${TMP_ROOT}/matrix"
features <- paste0("gene", seq_len(80))
features[1] <- "MT-CO1"
for (sec in c("sec_1", "sec_2", "sec_3")) {
  dir.create(file.path(root, sec), recursive = TRUE, showWarnings = FALSE)
  counts <- matrix(rpois(80 * 50, lambda = 3), nrow = 80)
  counts[1, ] <- rpois(50, lambda = 1)
  counts[2:8, ] <- counts[2:8, ] + 2
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  Matrix::writeMM(counts, file.path(root, sec, "matrix.mtx"))
  write.table(
    data.frame(gene_id = features, gene_name = features, type = "Gene Expression"),
    file.path(root, sec, "features.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  barcodes <- paste0(sec, "_spot", seq_len(50))
  writeLines(barcodes, file.path(root, sec, "barcodes.tsv"))
  write.csv(
    data.frame(
      barcode = barcodes,
      row = rep(seq_len(5), each = 10),
      col = rep(seq_len(10), times = 5),
      in_tissue = TRUE
    ),
    file.path(root, sec, "coords.csv"),
    row.names = FALSE
  )
}
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
export SPATIAL_OBJECT_LAYER_FILE="${TMP_ROOT}/project/config/spatial_object_layers.tsv"
export EDA_GATE_FILE="${TMP_ROOT}/project/config/eda_gates.tsv"
export SPATIAL_RESULTS_DIR="${TMP_ROOT}/project/results/spatial"
export SPATIAL_CHECKPOINT_DIR="${TMP_ROOT}/project/results/spatial/checkpoints"
export SPATIAL_TABLE_DIR="${TMP_ROOT}/project/results/spatial/tables"
export SPATIAL_FIGURE_DIR="${TMP_ROOT}/project/results/spatial/figures"
export SPATIAL_PRE_QC_REPORT_DIR="${TMP_ROOT}/project/reports/eda/spatial_pre_qc"
export SPATIAL_POST_QC_REPORT_DIR="${TMP_ROOT}/project/reports/eda/spatial_post_qc"
export SPATIAL_NORMALIZATION_VARIANTS_DIR="${TMP_ROOT}/project/results/spatial/checkpoints/normalization_variants"
export SPATIAL_NORMALIZATION_COMPARE_DIR="${TMP_ROOT}/project/reports/spatial/normalization_compare"
export SPATIAL_NORMALIZATION_OVERRIDE_FILE="${TMP_ROOT}/project/config/spatial_normalization_override.tsv"
export SPATIAL_DEFAULT_NORMALIZATION_METHOD="m3_sct_v2"
export PY_SPATIAL_BIN="${TMP_ROOT}/missing_py_spatial/bin/python"
export CLEAN_GTF="${TMP_ROOT}/refs/genes.gtf"
EOF

{
  printf 'sample_id\tcondition\tbiological_replicate\ttechnical_replicate\tbatch\tinput_mode\tinput_source\tsource_path\tplatform\tgene_id_type\treference_version\tgroup_id\ttimepoint\ttissue\tchemistry\trun_main\trun_velocity\trun_scenic\tmodality\tsection_id\tchip_id\tbundle_layout\timage_path\trun_spatial\trun_deconv\trun_joint\tnotes\n'
  for sec in sec_1 sec_2 sec_3; do
    printf 'st_%s\tsmoke\t%s\t1\tst_batch\tspatial_matrix\tsmoke\t%s/matrix/%s\tgeneric\tauto\tcustom_reference\tsmoke\tD0\tchicken_ovary\tauto\tno\tno\tno\tspatial\t%s\tchip_%s\tgeneric_spatial_matrix\t\tyes\tauto\tauto\tsmoke\n' "${sec}" "${sec}" "${TMP_ROOT}" "${sec}" "${sec}" "${sec}"
  done
} > "${TMP_ROOT}/project/metadata/samples.canonical.tsv"
cp "${TMP_ROOT}/project/metadata/samples.canonical.tsv" "${TMP_ROOT}/project/metadata/samples.tsv"

{
  printf 'section_id\tchip_id\ttissue_block\tcondition\tplatform\tbundle_root\timage_lowres\timage_hires\ttissue_positions\tscalefactors\tenabled\tnotes\n'
  for sec in sec_1 sec_2 sec_3; do
    printf '%s\tchip_%s\tblock\tsmoke\tgeneric\t%s/matrix/%s\t\t\t\t\tyes\tsmoke\n' "${sec}" "${sec}" "${TMP_ROOT}" "${sec}"
  done
} > "${TMP_ROOT}/project/metadata/sections.tsv"

{
  printf 'sample_id\tsection_id\tmodality\tcondition\tplatform\tbundle_layout\tsource_path\tsection_bundle_root\tstandardized_outs\tready\tready_status\trequired_relpaths\tmissing_required_relpaths\toptional_present_relpaths\toptional_missing_relpaths\tnotes\n'
  for sec in sec_1 sec_2 sec_3; do
    printf 'st_%s\t%s\tspatial\tsmoke\tgeneric\tgeneric_spatial_matrix\t%s/matrix/%s\t%s/matrix/%s\t%s/matrix/%s\ttrue\ttrue\t\t\t\t\tsmoke\n' "${sec}" "${sec}" "${TMP_ROOT}" "${sec}" "${TMP_ROOT}" "${sec}" "${TMP_ROOT}" "${sec}"
  done
} > "${TMP_ROOT}/project/reports/intake/spatial_input_inventory.tsv"

cat > "${TMP_ROOT}/project/config/spatial_intake_contract.tsv" <<'EOF'
platform	bundle_layout	required_relpaths	optional_relpaths	standardized_outs_relpath	ready_status_when_present	notes
generic	generic_spatial_matrix	matrix.mtx,features.tsv,barcodes.tsv	coords.csv	.	true	Smoke generic spatial matrix.
EOF

cat > "${TMP_ROOT}/project/config/spatial_qc_thresholds.tsv" <<'EOF'
section_id	qc_min_nfeature	qc_max_nfeature	qc_min_ncount	qc_max_ncount	qc_max_mito_pct	mito_set_override	spatial_aware_filter	excessive_drop_threshold
__DEFAULT__	5		20		95		true	0.5
sec_1	5		20		95		true	0.5
sec_2	5		20		95		true	0.5
sec_3	5		20		95		true	0.5
EOF
cp "${PIPELINE_ROOT}/config/spatial_object_layers.tsv.template" "${TMP_ROOT}/project/config/spatial_object_layers.tsv"
cp "${PIPELINE_ROOT}/config/spatial_normalization_override.tsv.template" "${TMP_ROOT}/project/config/spatial_normalization_override.tsv"
cp "${PIPELINE_ROOT}/config/eda_gates.tsv.template" "${TMP_ROOT}/project/config/eda_gates.tsv"

export SCRNA_PIPELINE_CONFIG="${TMP_ROOT}/project/config/project_config.sh"
set -a
source "${SCRNA_PIPELINE_CONFIG}"
set +a

Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/01_build_spatial_objects.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/01a_pre_spot_qc_eda.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/01b_spot_qc_filter.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/01c_post_spot_qc_eda.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/02_normalize_spatial.R" --methods m0,m1,m3
SPATIAL_NORMALIZE_EXPECT_M4_STATUS=skipped \
  Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/testing/smoke_spatial_normalize.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/02_normalize_spatial.R" --methods all
SPATIAL_NORMALIZE_EXPECT_M4_STATUS=failed_py_bridge \
  Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/testing/smoke_spatial_normalize.R"

echo "smoke_spatial_normalize_ok"
