#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

mkdir -p "${TMP_DIR}/metadata" "${TMP_DIR}/config"
export PROJECT_ROOT="${TMP_DIR}"
export PIPELINE_ROOT="${ROOT}"
export RESULTS_DIR="${TMP_DIR}/results"
export REPORT_DIR="${TMP_DIR}/reports"
export METADATA_DIR="${TMP_DIR}/metadata"
export PROJECT_CONFIG_DIR="${TMP_DIR}/config"
export MANIFEST_DIR="${TMP_DIR}/results/manifests"
export PY_SPATIAL_BIN="${PY_SPATIAL_BIN:-python3}"

python3 -m py_compile "${ROOT}/workflow/04python/spatialde2_svg.py"
Rscript -e 'files <- c("workflow/05single_script/spatial/09a_spatialde2_svg.R","workflow/05single_script/spatial/09b_sparkx_svg.R","workflow/05single_script/spatial/helpers/spatial_svg_utils.R"); for (f in files) parse(f)'
Rscript "${ROOT}/workflow/05single_script/spatial/09b_sparkx_svg.R" >/dev/null

test -s "${TMP_DIR}/results/spatial/tables/09_svg/sparkx/sparkx_manifest.tsv"
grep -Eq 'skipped_no_sparkx|skipped_sparkx_adapter_unimplemented|ok' "${TMP_DIR}/results/spatial/tables/09_svg/sparkx/sparkx_manifest.tsv"
grep -q 'is_ligand_candidate' "${ROOT}/workflow/04python/spatialde2_svg.py"

echo "smoke_spatial_svg_core_methods_ok"
