#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/spatial_de_05_smoke.XXXXXX")"
trap 'rm -rf "${TMP_ROOT}"' EXIT

source "${SCRIPT_DIR}/lib_spatial_smoke_05.sh"
spatial_smoke05_prepare_region_fixture "${TMP_ROOT}" "${PIPELINE_ROOT}"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/05_region_marker_discovery.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/05a_region_pseudobulk.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/05b_region_composition.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/05_spatial_de.R"
spatial_smoke05_assert_manifest_status "${SPATIAL_TABLE_DIR}/spatial_05_spatial_de/spatial_de_manifest.tsv" "ok"
grep -q '^spatial_marker_de	' "${PIPELINE_ROOT}/config/eda_gates.tsv.template"
echo "smoke_spatial_de_05_ok"
