#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

check_stage_deps "70_joint_scrna_spatial"
prepare_project_state_dirs

JOINT_R_DIR="${WORKFLOW_ROOT}/05single_script/joint"

run_future_joint_r_stage() {
  local script_name="$1"
  local manifest_name="$2"
  shift 2 || true
  local script_path="${JOINT_R_DIR}/${script_name}"
  local manifest_path="${MANIFEST_DIR}/${manifest_name}/_manifest.json"
  [[ -f "${script_path}" ]] || die "R stage 尚未实现: ${script_path}。本轮只完成非 R 脚手架。"
  run_stage_if_stale_with_runner run_r_spatial "${script_path}" "${manifest_path}" "$@"
}

run_future_joint_r_stage "00_freeze_scrna_reference.R" "joint_00_freeze_scrna_reference" "${SPATIAL_REFERENCE_INVENTORY_FILE}"
run_future_joint_r_stage "01_label_transfer_to_spatial.R" "joint_01_label_transfer" "${SPATIAL_REFERENCE_INVENTORY_FILE}"
run_future_joint_r_stage "02_joint_trajectory.R" "joint_02_trajectory"
set_eda_gate_status "joint_trajectory" "pending" "" "Review joint trajectory projection outputs."
hold_for_gate joint_trajectory

run_future_joint_r_stage "03_joint_communication_evidence.R" "joint_03_communication"
set_eda_gate_status "joint_communication" "pending" "" "Review joint communication evidence outputs."
hold_for_gate joint_communication

run_future_joint_r_stage "04_joint_neighborhood_validation.R" "joint_04_neighborhood_validation"

update_workflow_status \
  "joint_completed" \
  "review joint reports under ${JOINT_RESULTS_DIR}" \
  "status.joint_completed=true"
