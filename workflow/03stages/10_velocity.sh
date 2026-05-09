#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

ensure_metadata_fresh
check_stage_deps "10_velocity"
hold_for_gate annotation
hold_for_gate subcluster

MODULE_10A_MANIFEST="${MANIFEST_DIR}/10a_run_velocyto/_manifest.json"
MODULE_10B_MANIFEST="${MANIFEST_DIR}/10b_prepare_velocity_reference/_manifest.json"

ensure_eda_control_files
ensure_dir "${VELOCITY_INPUT_DIR}" "${VELOCITY_LOOM_DIR}" "${VELOCITY_OUTPUT_DIR}"

VELOCITY_INPUT_RERAN=0

run_velocity_input_if_stale() {
  local runner_fn="$1"
  local script_path="$2"
  local output_path="$3"
  shift 3 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "Running ${script_path}"
    "${runner_fn}" "${script_path}"
    VELOCITY_INPUT_RERAN=1
  else
    echo "Existing output is current; skipping: ${output_path}"
  fi
}

run_velocity_input_if_stale \
  run_velocity \
  "${WORKFLOW_ROOT}/05single_script/10a_run_velocyto.sh" \
  "${MODULE_10A_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}" \
  "${SAMPLE_SHEET}" \
  "${CANONICAL_SAMPLE_SHEET}" \
  "${DNBC4TOOLS_OUT_DIR}" \
  "${CELLRANGER_OUT_DIR}" \
  "${VELOCITY_GTF}"
require_manifest_output "${MODULE_10A_MANIFEST}" "velocity_loom_index" >/dev/null

run_velocity_input_if_stale \
  run_r_main \
  "${WORKFLOW_ROOT}/05single_script/10b_prepare_velocity_reference.R" \
  "${MODULE_10B_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}" \
  "${MODULE_03D_MANIFEST:-${MANIFEST_DIR}/03d_annotate/_manifest.json}" \
  "${MODULE_04B_MANIFEST:-${MANIFEST_DIR}/04b_subcluster_annotate/_manifest.json}" \
  "${LAYER_STATUS_FILE}"
require_manifest_output "${MODULE_10B_MANIFEST}" "velocity_reference_index" >/dev/null

if [[ "${VELOCITY_INPUT_RERAN}" == "1" ]]; then
  set_eda_gate_status \
    "velocity_inputs" \
    "pending" \
    "" \
    "10a/10b RNA velocity input preparation completed; review loom and reference indices before approving velocity_inputs."
fi

sync_workflow_gate_statuses
update_workflow_status \
  "10_velocity_inputs_completed" \
  "review velocity loom/reference indices and approve velocity_inputs before running downstream RNA velocity methods" \
  "status.10a_velocity_loom_completed=true" \
  "status.10b_velocity_reference_completed=true" \
  "status.10_velocity_inputs_completed=true" \
  "status.velocity_inputs_gate_passed=$(eda_gate_passed velocity_inputs && echo true || echo false)"
