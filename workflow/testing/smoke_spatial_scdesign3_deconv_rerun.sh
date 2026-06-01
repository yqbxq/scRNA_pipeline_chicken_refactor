#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

mkdir -p "${TMP_DIR}/metadata"
cp "${ROOT}/metadata/deconv_pairs.tsv" "${TMP_DIR}/metadata/deconv_pairs.tsv"
cp "${ROOT}/metadata/spatial_reference_inventory.tsv" "${TMP_DIR}/metadata/spatial_reference_inventory.tsv"

export PROJECT_ROOT="${TMP_DIR}"
export PIPELINE_ROOT="${ROOT}"
export RESULTS_DIR="${TMP_DIR}/results"
export REPORT_DIR="${TMP_DIR}/reports"
export METADATA_DIR="${TMP_DIR}/metadata"
export PROJECT_CONFIG_DIR="${TMP_DIR}/config"
export MANIFEST_DIR="${TMP_DIR}/results/manifests"
export SPATIAL_VALIDATION_MODE="scdesign3"
export SPATIAL_VALIDATION_N_SPOTS="8"

Rscript "${ROOT}/workflow/05single_script/spatial/07f_deconvolution_validation.R" >/dev/null

VALIDATION_DIR="${TMP_DIR}/results/spatial/tables/spatial_07f_deconvolution_validation"
RUN_DIR="${VALIDATION_DIR}/deconv_validation_default"

test -s "${VALIDATION_DIR}/validation_manifest.tsv"
test -s "${RUN_DIR}/synthetic_h5ad_manifest.tsv"
test -s "${RUN_DIR}/synthetic_deconv_manifest.tsv"
grep -q 'synthetic_reference' "${RUN_DIR}/synthetic_h5ad_manifest.tsv"
grep -q 'synthetic_spatial' "${RUN_DIR}/synthetic_h5ad_manifest.tsv"
grep -q 'synthetic_prediction_source' "${VALIDATION_DIR}/validation_manifest.tsv"
grep -q 'synthetic_rerun_ok_methods' "${VALIDATION_DIR}/validation_manifest.tsv"

if awk -F'\t' '$1 == "I11_deconv_validation" && $3 == "PASS" { found = 1 } END { exit found ? 0 : 1 }' "${VALIDATION_DIR}/spatial_question_gate_status.tsv"; then
  echo "scdesign3 mode must not PASS I11 without synthetic_h5ad_rerun >=2 methods" >&2
  exit 1
fi

echo "smoke_spatial_scdesign3_deconv_rerun_ok"
