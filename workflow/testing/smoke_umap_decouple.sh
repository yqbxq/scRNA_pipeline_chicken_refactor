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

reject_grep() {
  local pattern="$1"
  local path="$2"
  if rg -q "${pattern}" "${path}"; then
    fail "unexpected pattern ${pattern} in ${path}"
  fi
}

reject_grep "RunUMAP" "workflow/05single_script/03a2_reduce_integrate.R"
reject_grep "RunUMAP" "workflow/05single_script/spatial/02a_spatial_integration_eda.R"
reject_grep "spatial_run_umap_or_fallback" "workflow/05single_script/spatial/02a_spatial_integration_eda.R"

require_grep "run_umap_only" "workflow/05single_script/03a3_compute_umap.R"
require_grep "spatial_run_umap_or_fallback" "workflow/05single_script/spatial/02a2_compute_umap_spatial.R"
require_grep "03a3_compute_umap" "workflow/03stages/03_panorama.sh"
require_grep "02a2_compute_umap_spatial.R" "workflow/03stages/50_run_spatial_pipeline.sh"
require_grep "candidate_umap_index_tsv" "workflow/05single_script/03b_integration_eda.R"

reject_grep "FindNeighbors\\([^\\n]*umap" "workflow/05single_script/03c_cluster.R"
reject_grep "FindNeighbors\\([^\\n]*umap" "workflow/05single_script/spatial/02b_finalize_spatial_clustering.R"

for name in UMAP_N_NEIGHBORS UMAP_MIN_DIST UMAP_SPREAD UMAP_SEED UMAP_METRIC UMAP_LOCAL_CONNECTIVITY MODULE_03A3_VERSION MODULE_SPATIAL_02A2_VERSION; do
  require_grep "${name}" "workflow/02lib/common.sh"
  require_grep "${name}" "workflow/02lib/env_registry.sh"
  require_grep "${name}" "docs/env_registry.md"
done

echo "[PASS] UMAP decouple smoke checks passed."
