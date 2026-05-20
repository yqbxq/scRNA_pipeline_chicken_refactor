#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/spatial_composition_05b_pass.XXXXXX")"
trap 'rm -rf "${TMP_ROOT}"' EXIT

source "${SCRIPT_DIR}/lib_spatial_smoke_05.sh"
spatial_smoke05_prepare_passgate_fixture "${TMP_ROOT}" "${PIPELINE_ROOT}"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/05b_region_composition.R"
spatial_smoke05_assert_manifest_status "${SPATIAL_TABLE_DIR}/spatial_05b_region_composition/composition_manifest.tsv" "ok" "yes"
echo "smoke_spatial_composition_05b_pass_ok"
