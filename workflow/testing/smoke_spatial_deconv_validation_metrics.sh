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

mkdir -p "${TMP_DIR}/results/spatial/tables/spatial_07a_deconvolution_rctd/default_deconv/all"

cat >"${TMP_DIR}/results/spatial/tables/spatial_07a_deconvolution_rctd/default_deconv/all/spot_celltype_proportions.tsv" <<'TSV'
spot_id	cell_type	proportion
s1	A	0.75
s1	B	0.25
s2	A	0.20
s2	B	0.80
TSV

cat >"${TMP_DIR}/results/spatial/tables/spatial_07a_deconvolution_rctd/rctd_manifest.tsv" <<TSV
deconv_id	section	tool	status	reason	panorama_input_rds	n_spots	n_celltypes	runtime_sec	proportion_tsv	proportion_wide_tsv	spot_metadata_tsv	summary_tsv	method_object_rds	method_version
default_deconv	all	rctd	ok			2	2	1	${TMP_DIR}/results/spatial/tables/spatial_07a_deconvolution_rctd/default_deconv/all/spot_celltype_proportions.tsv					1.0
TSV

cat >"${TMP_DIR}/truth.tsv" <<'TSV'
spot_id	cell_type	true_proportion
s1	A	0.70
s1	B	0.30
s2	A	0.25
s2	B	0.75
TSV

export SPATIAL_VALIDATION_TRUTH_TSV="${TMP_DIR}/truth.tsv"
Rscript "${ROOT}/workflow/05single_script/spatial/07f_deconvolution_validation.R" >/dev/null

MANIFEST="${TMP_DIR}/results/spatial/tables/spatial_07f_deconvolution_validation/validation_manifest.tsv"
SUMMARY="${TMP_DIR}/results/spatial/tables/spatial_07f_deconvolution_validation/deconv_validation_default/method_summary.tsv"

test -s "${MANIFEST}"
test -s "${SUMMARY}"
grep -q $'\trctd\t' "${MANIFEST}"
grep -q '^rctd' "${SUMMARY}"
awk -F '\t' 'NR > 1 { if ($4 == "" || $5 == "") exit 1; found = 1 } END { exit found ? 0 : 1 }' "${SUMMARY}"

echo "smoke_spatial_deconv_validation_metrics_ok"
