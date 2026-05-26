#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"

fail() {
  echo "[FAIL] $*" >&2
  exit 1
}

require_grep() {
  local pattern="$1"
  local path="$2"
  rg -q "${pattern}" "${path}" || fail "missing pattern ${pattern} in ${path}"
}

for name in H5AD_EXPORT_GATES H5AD_EXPORT_ON_FAILURE H5AD_PYTHON_FALLBACK_TO_RDS; do
  require_grep "${name}" workflow/02lib/common.sh
  require_grep "${name}" workflow/02lib/env_registry.sh
  require_grep "${name}" docs/env_registry.md
done

require_grep "should_export_h5ad" workflow/02lib/common.sh
require_grep "handle_h5ad_export_failure" workflow/02lib/common.sh
require_grep "export_h5ad_for_gate" workflow/02lib/common.sh

for gate in \
  03d_annotation \
  04b_subcluster \
  06d_enrichment \
  07c_communication \
  spatial_03_region \
  spatial_04b_subcluster \
  spatial_05_marker \
  spatial_06_extensions; do
  require_grep "${gate}" workflow/02lib/common.sh
done

require_grep "export_h5ad_for_gate" workflow/03stages/03_panorama.sh
require_grep "export_h5ad_for_gate" workflow/03stages/04_subcluster.sh
require_grep "export_h5ad_for_gate" workflow/03stages/06_enrichment.sh
require_grep "export_h5ad_for_gate" workflow/03stages/07_communication.sh
require_grep "export_h5ad_for_gate" workflow/03stages/50_run_spatial_pipeline.sh
require_grep "export_h5ad_for_gate" workflow/03stages/60_run_spatial_extensions.sh

require_grep "load_scrna_h5ad" workflow/04python/helpers/scrna_io.py
require_grep "load_spatial_h5ad" workflow/04python/helpers/spatial_io.py
require_grep "load_all_spatial_h5ad" workflow/04python/helpers/spatial_io.py
require_grep "DEPRECATED: P-R02-E" workflow/04python/helpers/scrna_io.py
require_grep "DEPRECATED: P-R02-E" workflow/04python/helpers/spatial_io.py
require_grep "resolve_scrna_h5ad_path" workflow/04python/cell2location_pipeline.py
require_grep "resolve_spatial_h5ad_path" workflow/04python/cell2location_pipeline.py
require_grep "resolve_spatial_h5ad_path" workflow/04python/spatial_neighborhood.py
require_grep "resolve_spatial_h5ad_path" workflow/04python/spatial_decoupler.py
require_grep "find_scrna_h5ad" workflow/04python/scvelo_pipeline.py

tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT
mkdir -p "${tmp}/90a_export_h5ad/03d_panorama" "${tmp}/90a_export_h5ad/spatial_03_region"
touch "${tmp}/90a_export_h5ad/03d_panorama/mock.h5ad"
touch "${tmp}/90a_export_h5ad/spatial_03_region/section_A.h5ad"

PYTHONPATH="${ROOT}/workflow/04python" RESULTS_DIR="${tmp}" python - <<'PY'
from helpers.scrna_io import resolve_scrna_h5ad_path
from helpers.spatial_io import resolve_spatial_h5ad_path

assert resolve_scrna_h5ad_path("03d_panorama").name == "mock.h5ad"
assert resolve_spatial_h5ad_path("spatial_03_region", section_id="section_A").name == "section_A.h5ad"
print("python_h5ad_helper_resolution_ok")
PY

echo "[PASS] H5AD hookup smoke checks passed."
