#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/spatial_marker_05_smoke.XXXXXX")"
trap 'rm -rf "${TMP_ROOT}"' EXIT

source "${SCRIPT_DIR}/lib_spatial_smoke_05.sh"
spatial_smoke05_prepare_region_fixture "${TMP_ROOT}" "${PIPELINE_ROOT}"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/05_region_marker_discovery.R"
spatial_smoke05_assert_manifest_status "${SPATIAL_TABLE_DIR}/spatial_05_region_marker/marker_discovery_manifest.tsv" "ok"
echo "smoke_spatial_marker_05_ok"
