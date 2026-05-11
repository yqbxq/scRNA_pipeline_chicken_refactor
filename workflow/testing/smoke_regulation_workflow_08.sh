#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PIPELINE_ROOT="$(cd "${WORKFLOW_ROOT}/.." && pwd)"
R_BIN="${R_BIN:-Rscript}"
RESULTS_DIR="${RESULTS_DIR:-${PIPELINE_ROOT}/results}"

require_rscript_08() {
  if ! command -v "${R_BIN}" >/dev/null 2>&1; then
    echo "Missing Rscript for smoke_regulation_workflow_08: ${R_BIN}" >&2
    exit 1
  fi
}

parse_r_08() {
  local file_path="$1"
  "${R_BIN}" -e "parse(file='${file_path}')" >/dev/null
}

check_static_08() {
  bash -n "${WORKFLOW_ROOT}/01run.sh"
  bash -n "${WORKFLOW_ROOT}/03stages/08_regulation.sh"
  bash -n "${WORKFLOW_ROOT}/05single_script/08b_scenic_grn.sh"

  require_rscript_08
  parse_r_08 "${WORKFLOW_ROOT}/05single_script/08a_scenic_export.R"
  parse_r_08 "${WORKFLOW_ROOT}/05single_script/08c_scenic_regulons.R"
  parse_r_08 "${WORKFLOW_ROOT}/05single_script/08d_scenic_downstream.R"
  parse_r_08 "${WORKFLOW_ROOT}/05single_script/08e_decoupler.R"
  parse_r_08 "${WORKFLOW_ROOT}/05single_script/08f_regulation_eda.R"
  parse_r_08 "${WORKFLOW_ROOT}/05single_script/helpers/project_paths_08.R"
  parse_r_08 "${WORKFLOW_ROOT}/05single_script/helpers/scenic_utils.R"
  parse_r_08 "${WORKFLOW_ROOT}/05single_script/helpers/decoupler_utils.R"
}

check_metadata_08() {
  require_rscript_08
  "${R_BIN}" "${WORKFLOW_ROOT}/05single_script/95_build_metadata_from_questions.R"
  "${R_BIN}" "${WORKFLOW_ROOT}/05single_script/96_validate_metadata_consistency.R"
}

check_decoupler_outputs_08() {
  local layer="${REGULATION_LAYERS:-panorama}"
  local safe_layer="${layer//[^A-Za-z0-9_.-]/_}"
  test -s "${RESULTS_DIR}/manifests/08e_decoupler/_manifest.json"
  test -s "${RESULTS_DIR}/tables/regulation/decoupler/${safe_layer}/tf_activity.tsv"
  test -s "${RESULTS_DIR}/tables/regulation/decoupler/${safe_layer}/pathway_activity.tsv"
  test -s "${RESULTS_DIR}/tables/regulation/decoupler/${safe_layer}/network_mapping_summary.tsv"
}

run_server_smoke_08() {
  env SCENIC_RESOURCES_ONLY=yes bash "${WORKFLOW_ROOT}/03stages/08_regulation.sh"
  env REGULATION_LAYERS="${REGULATION_LAYERS:-panorama}" RUN_SCENIC_GRN=no RUN_DECOUPLER=yes bash "${WORKFLOW_ROOT}/03stages/08_regulation.sh"
  check_decoupler_outputs_08
  env REGULATION_LAYERS="${REGULATION_LAYERS:-panorama}" RUN_SCENIC_GRN=yes RUN_DECOUPLER=yes bash "${WORKFLOW_ROOT}/03stages/08_regulation.sh"
}

check_static_08

if [[ "${SMOKE_08_METADATA:-no}" == "yes" ]]; then
  check_metadata_08
fi

if [[ "${SMOKE_08_SERVER:-no}" == "yes" ]]; then
  run_server_smoke_08
fi

echo "smoke_regulation_workflow_08: PASS"
