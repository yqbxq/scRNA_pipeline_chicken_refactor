#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

ensure_metadata_fresh
check_stage_deps "09_trajectory"
hold_for_gate annotation
hold_for_gate subcluster

MODULE_03D_MANIFEST="${MANIFEST_DIR}/03d_annotate/_manifest.json"
MODULE_04B_MANIFEST="${MANIFEST_DIR}/04b_subcluster_annotate/_manifest.json"
MODULE_09A_MANIFEST="${MANIFEST_DIR}/09a_trajectory_inputs/_manifest.json"
MODULE_09B_MANIFEST="${MANIFEST_DIR}/09b_trajectory_inputs_eda/_manifest.json"

ensure_eda_control_files
require_manifest_output "${MODULE_03D_MANIFEST}" "annotated_object" >/dev/null

if awk -F '\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      idx[$i] = i
    }
    next
  }
  NR > 1 && $0 !~ /^#/ {
    method = tolower($(idx["method"]))
    enabled = tolower($(idx["enabled"]))
    scope = $(idx["layer_scope"])
    if (method == "trajectory" && enabled != "no" && scope != "" && scope != "panorama") {
      found = 1
    }
  }
  END { exit(found ? 0 : 1) }
' "${TRAJECTORY_PAIRS_SHEET}"; then
  [[ -s "${MODULE_04B_MANIFEST}" ]] || die "trajectory_pairs.tsv 请求 subcluster layer，但缺少 04b manifest: ${MODULE_04B_MANIFEST}"
fi

TRAJECTORY_INPUT_RERAN=0

run_trajectory_stage_if_stale() {
  local script_path="$1"
  local output_path="$2"
  shift 2 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "运行 ${script_path}"
    run_r_main "${script_path}"
    TRAJECTORY_INPUT_RERAN=1
  else
    echo "已存在且未过期，跳过: ${output_path}"
  fi
}

run_trajectory_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/09a_trajectory_inputs.R" \
  "${MODULE_09A_MANIFEST}" \
  "${MODULE_03D_MANIFEST}" \
  "${MODULE_04B_MANIFEST}" \
  "${ORTHOLOG_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}"
require_manifest_output "${MODULE_09A_MANIFEST}" "trajectory_input_index" >/dev/null

run_trajectory_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/09b_trajectory_inputs_eda.R" \
  "${MODULE_09B_MANIFEST}" \
  "${MODULE_09A_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}"
TRAJECTORY_INPUT_REPORT="$(require_manifest_output "${MODULE_09B_MANIFEST}" "report_md")"

if [[ "${TRAJECTORY_INPUT_RERAN}" == "1" ]]; then
  set_eda_gate_status \
    "trajectory_inputs" \
    "pending" \
    "" \
    "09a/09b trajectory input preparation completed; review ${TRAJECTORY_INPUT_REPORT}, then approve trajectory_inputs before 09c+."
fi

sync_workflow_gate_statuses
update_workflow_status \
  "09_trajectory_inputs_completed" \
  "review ${TRAJECTORY_INPUT_REPORT}; approve trajectory_inputs before running trajectory method scripts" \
  "status.09a_trajectory_inputs_completed=true" \
  "status.09b_trajectory_inputs_eda_completed=true" \
  "status.09_trajectory_inputs_completed=true" \
  "status.trajectory_inputs_gate_passed=$(eda_gate_passed trajectory_inputs && echo true || echo false)"
