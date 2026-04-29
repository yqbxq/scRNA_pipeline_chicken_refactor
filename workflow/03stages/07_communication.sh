#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

check_stage_deps "07_communication"
hold_for_gate annotation

MODULE_00_MANIFEST="${ORTHOLOG_MANIFEST}"
MODULE_03D_MANIFEST="${MANIFEST_DIR}/03d_annotate/_manifest.json"
MODULE_04B_MANIFEST="${MANIFEST_DIR}/04b_subcluster_annotate/_manifest.json"
MODULE_05D_MANIFEST="${MANIFEST_DIR}/05d_deg_eda/_manifest.json"
MODULE_07A_MANIFEST="${MANIFEST_DIR}/07a_cellchat/_manifest.json"
MODULE_07B_MANIFEST="${MANIFEST_DIR}/07b_nichenet/_manifest.json"
MODULE_07C_MANIFEST="${MANIFEST_DIR}/07c_communication_eda/_manifest.json"

ensure_eda_control_files
ensure_communication_pairs_sheet
require_manifest_output "${MODULE_03D_MANIFEST}" "annotated_object" >/dev/null
require_manifest_output "${MODULE_00_MANIFEST}" "human_best" >/dev/null

COMMUNICATION_RERAN=0

run_comm_stage_if_stale() {
  local script_path="$1"
  local output_path="$2"
  shift 2 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "运行 ${script_path}"
    run_r_interaction "${script_path}"
    COMMUNICATION_RERAN=1
  else
    echo "已存在且未过期，跳过: ${output_path}"
  fi
}

run_comm_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/07a_cellchat.R" \
  "${MODULE_07A_MANIFEST}" \
  "${MODULE_03D_MANIFEST}" \
  "${MODULE_04B_MANIFEST}" \
  "${MODULE_00_MANIFEST}" \
  "${COMMUNICATION_PAIRS_SHEET}"
require_manifest_output "${MODULE_07A_MANIFEST}" "cellchat_index_tsv" >/dev/null

if [[ "${RUN_NICHENET:-yes}" != "no" ]]; then
  if is_stale_output "${MODULE_07B_MANIFEST}" \
    "${MODULE_07A_MANIFEST}" \
    "${MODULE_05D_MANIFEST}" \
    "${MODULE_00_MANIFEST}" \
    "${COMMUNICATION_PAIRS_SHEET}" \
    "${NICHENET_RESOURCE_DIR}"; then
    echo "运行 ${WORKFLOW_ROOT}/05single_script/07b_nichenet.R"
    if run_r_interaction "${WORKFLOW_ROOT}/05single_script/07b_nichenet.R"; then
      COMMUNICATION_RERAN=1
    else
      warn "07b_nichenet 失败；继续生成 07c cellchat-only 报告。常见原因是 NicheNet 资源未准备。"
    fi
  else
    echo "已存在且未过期，跳过: ${MODULE_07B_MANIFEST}"
  fi
else
  warn "RUN_NICHENET=no，跳过 07b_nichenet。"
fi

run_comm_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/07c_communication_eda.R" \
  "${MODULE_07C_MANIFEST}" \
  "${MODULE_07A_MANIFEST}" \
  "${MODULE_07B_MANIFEST}" \
  "${COMMUNICATION_PAIRS_SHEET}"
COMMUNICATION_REPORT="$(require_manifest_output "${MODULE_07C_MANIFEST}" "report_md")"

if [[ "${COMMUNICATION_RERAN}" == "1" ]]; then
  set_eda_gate_status \
    "communication" \
    "pending" \
    "" \
    "07 communication completed; review ${COMMUNICATION_REPORT}, then approve communication if desired"
fi

sync_workflow_gate_statuses
update_workflow_status \
  "07_communication_completed" \
  "review ${COMMUNICATION_REPORT}; communication gate is advisory and does not block downstream modules" \
  "status.07a_cellchat_completed=true" \
  "status.07b_nichenet_completed=$([[ -f "${MODULE_07B_MANIFEST}" ]] && echo true || echo false)" \
  "status.07c_communication_eda_completed=true" \
  "status.07_communication_completed=true" \
  "status.communication_gate_passed=$(eda_gate_passed communication && echo true || echo false)"
