#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d)"

cleanup() {
  rm -rf "${TMP_ROOT}"
}
trap cleanup EXIT

make_mock_rds() {
  local path="$1"
  local mode="$2"
  Rscript - "${path}" "${mode}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
mode <- args[[2]]
obj <- list(
  obs = data.frame(sample_id = c("S1", "S2"), condition = c("syf", "f5"), celltype_l1 = c("GC", "TC"), stringsAsFactors = FALSE),
  var = data.frame(gene_id = c("g1", "g2"), gene_symbol = c("G1", "G2"), stringsAsFactors = FALSE),
  obsm = list(X_pca = matrix(c(1, 2, 3, 4), nrow = 2), X_umap = matrix(c(1, 2, 3, 4), nrow = 2)),
  layers = list(counts = matrix(as.integer(c(1, 0, 2, 3)), nrow = 2)),
  uns = list(module = "smoke", export_timestamp = "2026-05-26T00:00:00Z")
)
if (mode == "sct") {
  obj$layers$SCT <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
}
if (mode == "bad_contract") {
  obj$obsm$X_pca <- NULL
}
if (mode == "spatial") {
  obj$obs$section_id <- c("sec1", "sec2")
  obj$obsm$spatial <- matrix(c(10, 20, 11, 21), nrow = 2)
  obj$spatial <- list(
    tissue_positions = data.frame(spot_id = c("s1", "s2"), x = c(10, 11), y = c(20, 21)),
    scale_factors = list(spot_diameter_fullres = 65)
  )
}
saveRDS(obj, path)
RS
}

make_manifest() {
  local manifest="$1"
  local rds="$2"
  python3 - "${manifest}" "${rds}" <<'PY'
import json
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
rds = Path(sys.argv[2])
manifest.parent.mkdir(parents=True, exist_ok=True)
payload = {
    "module": "mock_upstream",
    "version": "1.0",
    "timestamp": "2026-05-26T00:00:00Z",
    "base_dir": str(rds.parent),
    "inputs": {},
    "outputs": {
        "annotated_object": {
            "path": rds.name,
            "type": "rds",
            "produced_by": "smoke",
            "row_semantics": "mock object",
        }
    },
    "depends_on": [],
}
manifest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
PY
}

run_export() {
  local case_id="$1"
  local mode="$2"
  local modality="$3"
  local case_dir="${TMP_ROOT}/${case_id}"
  mkdir -p "${case_dir}"
  local rds="${case_dir}/object.rds"
  local manifest="${case_dir}/manifest/_manifest.json"
  local target="${case_dir}/h5ad"
  make_mock_rds "${rds}" "${mode}"
  make_manifest "${manifest}" "${rds}"
  H5AD_EXPORT_DIRECT_RSCRIPT=yes \
    H5AD_EXPORT_BACKEND=mock \
    H5AD_CONTRACT_FILE="${PIPELINE_ROOT}/metadata/h5ad_export_contract.tsv.template" \
    bash "${PIPELINE_ROOT}/workflow/03stages/90_export_h5ad.sh" \
      --upstream-manifest "${manifest}" \
      --output-key annotated_object \
      --modality "${modality}" \
      --target-dir "${target}" >/tmp/h5ad_export_${case_id}.log
  echo "${target}"
}

target_c1="$(run_export c1 minimal scrna)"
test -s "${target_c1}/mock_upstream_annotated_object_scrna.h5ad"
grep -q $'ok' "${target_c1}/h5ad_export_summary.tsv"

target_c2="$(run_export c2 sct scrna)"
grep -q 'SCT' "${target_c2}/mock_upstream_annotated_object_scrna.h5ad"

target_c3="$(run_export c3 bad_contract scrna)"
test -s "${target_c3}/mock_upstream_annotated_object_scrna_h5ad_skipped.json"
grep -q $'skipped_contract' "${target_c3}/h5ad_export_summary.tsv"

target_c4="$(run_export c4 spatial spatial)"
test -s "${target_c4}/mock_upstream_annotated_object_spatial_sec1.h5ad"
test -s "${target_c4}/mock_upstream_annotated_object_spatial_sec2.h5ad"

before_mtime="$(stat -c %Y "${target_c1}/_manifest.json")"
sleep 1
H5AD_EXPORT_DIRECT_RSCRIPT=yes \
  H5AD_EXPORT_BACKEND=mock \
  H5AD_CONTRACT_FILE="${PIPELINE_ROOT}/metadata/h5ad_export_contract.tsv.template" \
  bash "${PIPELINE_ROOT}/workflow/03stages/90_export_h5ad.sh" \
    --upstream-manifest "${TMP_ROOT}/c1/manifest/_manifest.json" \
    --output-key annotated_object \
    --modality scrna \
    --target-dir "${target_c1}" >/tmp/h5ad_export_c5.log
after_mtime="$(stat -c %Y "${target_c1}/_manifest.json")"
grep -q "跳过" /tmp/h5ad_export_c5.log
[[ "${before_mtime}" == "${after_mtime}" ]]

echo "smoke_h5ad_export_ok"
