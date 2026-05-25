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
STRICT_DECONV_SMOKE="${STRICT_DECONV_SMOKE:-no}"

for script in \
  07a_deconvolution_rctd.R \
  07b_deconvolution_transfer.R \
  07c_deconvolution_card.R \
  07d_deconvolution_cell2location.R \
  07e_deconvolution_compare.R \
  07f_deconvolution_validation.R; do
  Rscript "${ROOT}/workflow/05single_script/spatial/${script}" >/dev/null
done

test -s "${TMP_DIR}/results/spatial/tables/spatial_07a_deconvolution_rctd/rctd_manifest.tsv"
test -s "${TMP_DIR}/results/spatial/tables/spatial_07b_deconvolution_transfer/transfer_manifest.tsv"
test -s "${TMP_DIR}/results/spatial/tables/spatial_07c_deconvolution_card/card_manifest.tsv"
test -s "${TMP_DIR}/results/spatial/tables/spatial_07d_deconvolution_cell2location/cell2location_manifest.tsv"
test -s "${TMP_DIR}/results/spatial/tables/spatial_07e_deconvolution_compare/deconv_compare_manifest.tsv"
test -s "${TMP_DIR}/results/spatial/tables/spatial_07f_deconvolution_validation/validation_manifest.tsv"
grep -qx 'rctd' "${TMP_DIR}/results/spatial/tables/spatial_07e_deconvolution_compare/recommended_method.txt"

if [[ "${STRICT_DECONV_SMOKE}" == "yes" ]]; then
  for manifest in \
    "${TMP_DIR}/results/spatial/tables/spatial_07a_deconvolution_rctd/rctd_manifest.tsv" \
    "${TMP_DIR}/results/spatial/tables/spatial_07b_deconvolution_transfer/transfer_manifest.tsv" \
    "${TMP_DIR}/results/spatial/tables/spatial_07c_deconvolution_card/card_manifest.tsv"; do
    grep -q $'\tok\t' "${manifest}" || {
      echo "strict deconvolution smoke expected at least one ok row: ${manifest}" >&2
      exit 1
    }
  done
fi

echo "smoke_spatial_deconv_07_ok"
