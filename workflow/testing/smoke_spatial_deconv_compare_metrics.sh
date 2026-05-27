#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

export PROJECT_ROOT="${TMP_DIR}"
export PIPELINE_ROOT="${ROOT}"
export RESULTS_DIR="${TMP_DIR}/results"
export REPORT_DIR="${TMP_DIR}/reports"
export METADATA_DIR="${TMP_DIR}/metadata"
export PROJECT_CONFIG_DIR="${TMP_DIR}/config"
export MANIFEST_DIR="${TMP_DIR}/results/manifests"

mkdir -p \
  "${TMP_DIR}/results/spatial/tables/spatial_07a_deconvolution_rctd/default_deconv/all" \
  "${TMP_DIR}/results/spatial/tables/spatial_07b_deconvolution_transfer/default_deconv/all"

cat >"${TMP_DIR}/results/spatial/tables/spatial_07a_deconvolution_rctd/default_deconv/all/spot_celltype_proportions.tsv" <<'TSV'
spot_id	cell_type	proportion
s1	A	0.8
s1	B	0.2
s2	A	0.1
s2	B	0.9
TSV

cat >"${TMP_DIR}/results/spatial/tables/spatial_07b_deconvolution_transfer/default_deconv/all/spot_celltype_proportions.tsv" <<'TSV'
spot_id	cell_type	proportion
s1	A	0.7
s1	B	0.3
s2	A	0.2
s2	B	0.8
TSV

cat >"${TMP_DIR}/results/spatial/tables/spatial_07a_deconvolution_rctd/rctd_manifest.tsv" <<TSV
deconv_id	section	tool	status	reason	panorama_input_rds	n_spots	n_celltypes	runtime_sec	proportion_tsv	proportion_wide_tsv	spot_metadata_tsv	summary_tsv	method_object_rds	method_version
default_deconv	all	rctd	ok			2	2	1	${TMP_DIR}/results/spatial/tables/spatial_07a_deconvolution_rctd/default_deconv/all/spot_celltype_proportions.tsv					1.0
TSV

cat >"${TMP_DIR}/results/spatial/tables/spatial_07b_deconvolution_transfer/transfer_manifest.tsv" <<TSV
deconv_id	section	tool	status	reason	panorama_input_rds	n_spots	n_celltypes	runtime_sec	proportion_tsv	proportion_wide_tsv	spot_metadata_tsv	summary_tsv	method_object_rds	method_version
default_deconv	all	transfer	ok			2	2	1	${TMP_DIR}/results/spatial/tables/spatial_07b_deconvolution_transfer/default_deconv/all/spot_celltype_proportions.tsv					1.0
TSV

Rscript "${ROOT}/workflow/05single_script/spatial/07e_deconvolution_compare.R" >/dev/null

MATRIX="${TMP_DIR}/results/spatial/tables/spatial_07e_deconvolution_compare/method_comparison_matrix.tsv"
RANKING="${TMP_DIR}/results/spatial/tables/spatial_07e_deconvolution_compare/method_ranking.tsv"

test -s "${MATRIX}"
test -s "${RANKING}"
grep -q $'rctd\ttransfer\tall\t2\t' "${MATRIX}"
awk -F '\t' 'NR > 1 && $5 == "all" { if ($7 == "" || $8 == "" || $9 == "") exit 1; found = 1 } END { exit found ? 0 : 1 }' "${MATRIX}"
grep -qx 'rctd' "${TMP_DIR}/results/spatial/tables/spatial_07e_deconvolution_compare/recommended_method.txt"

echo "smoke_spatial_deconv_compare_metrics_ok"
