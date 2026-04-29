#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

ensure_metadata_fresh
check_stage_deps "05_deg"

MODULE_03D_MANIFEST="${MANIFEST_DIR}/03d_annotate/_manifest.json"
MODULE_05A_MANIFEST="${MANIFEST_DIR}/05a_marker_discovery/_manifest.json"
MODULE_05B_MANIFEST="${MANIFEST_DIR}/05b_pseudobulk_de/_manifest.json"
MODULE_05C_MANIFEST="${MANIFEST_DIR}/05c_composition/_manifest.json"
MODULE_05D_MANIFEST="${MANIFEST_DIR}/05d_deg_eda/_manifest.json"

ensure_eda_control_files
ensure_object_layer_config_file

require_manifest_output "${MODULE_03D_MANIFEST}" "annotated_object" >/dev/null

DEG_RERAN=0
run_deg_stage_if_stale() {
  local script_path="$1"
  local output_path="$2"
  shift 2 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "运行 ${script_path}"
    run_r_main "${script_path}"
    DEG_RERAN=1
  else
    echo "已存在且未过期，跳过: ${output_path}"
  fi
}

run_deg_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/05a_marker_discovery.R" \
  "${MODULE_05A_MANIFEST}" \
  "${LAYER_STATUS_FILE}" \
  "${COMPARISON_SHEET}" \
  "${MODULE_03D_MANIFEST}"
require_manifest_output "${MODULE_05A_MANIFEST}" "marker_discovery_manifest_tsv" >/dev/null

run_deg_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/05b_pseudobulk_de.R" \
  "${MODULE_05B_MANIFEST}" \
  "${MODULE_05A_MANIFEST}" \
  "${LAYER_STATUS_FILE}" \
  "${COMPARISON_SHEET}"
require_manifest_output "${MODULE_05B_MANIFEST}" "pseudobulk_manifest_tsv" >/dev/null

run_deg_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/05c_composition.R" \
  "${MODULE_05C_MANIFEST}" \
  "${LAYER_STATUS_FILE}" \
  "${COMPARISON_SHEET}" \
  "${MODULE_03D_MANIFEST}"
require_manifest_output "${MODULE_05C_MANIFEST}" "composition_manifest_tsv" >/dev/null

run_deg_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/05d_deg_eda.R" \
  "${MODULE_05D_MANIFEST}" \
  "${MODULE_05A_MANIFEST}" \
  "${MODULE_05B_MANIFEST}" \
  "${MODULE_05C_MANIFEST}"
DEG_REPORT="$(require_manifest_output "${MODULE_05D_MANIFEST}" "report")"
require_manifest_output "${MODULE_05D_MANIFEST}" "deg_status_matrix_tsv" >/dev/null

if [[ "${DEG_RERAN}" == "1" ]]; then
  set_eda_gate_status \
    "deg" \
    "pending" \
    "" \
    "05 DEG completed; review ${DEG_REPORT}, then approve deg before 06_enrichment"
fi

sync_workflow_gate_statuses
update_workflow_status \
  "05_deg_completed" \
  "review ${DEG_REPORT} and approve deg before bash ${PIPELINE_ROOT}/workflow/03stages/06_enrichment.sh" \
  "status.05a_marker_discovery_completed=true" \
  "status.05b_pseudobulk_de_completed=true" \
  "status.05c_composition_completed=true" \
  "status.05d_deg_eda_completed=true" \
  "status.05_deg_completed=true" \
  "status.deg_gate_passed=$(eda_gate_passed deg && echo true || echo false)"
