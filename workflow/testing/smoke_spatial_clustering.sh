#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/spatial_clustering_smoke.XXXXXX")"
trap 'rm -rf "${TMP_ROOT}"' EXIT

source "${SCRIPT_DIR}/lib_spatial_smoke_03.sh"
spatial_smoke03_require_r
spatial_smoke03_setup "${TMP_ROOT}" "${PIPELINE_ROOT}"
spatial_smoke03_run_to_normalized "${PIPELINE_ROOT}"

Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/02a_spatial_integration_eda.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/02b_finalize_spatial_clustering.R" --backends all
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/testing/smoke_spatial_clustering.R"

echo "smoke_spatial_clustering_ok"
