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
MODULE_10C_MANIFEST="${MANIFEST_DIR}/10c_scvelo_dynamical/_manifest.json"
MODULE_10D_MANIFEST="${MANIFEST_DIR}/10d_velocyto_steady_state/_manifest.json"
MODULE_10E_MANIFEST="${MANIFEST_DIR}/10e_scvelo_drivers/_manifest.json"
MODULE_10F_MANIFEST="${MANIFEST_DIR}/10f_cellrank_fate/_manifest.json"
MODULE_10H_MANIFEST="${MANIFEST_DIR}/10h_velocity_root_terminal/_manifest.json"

ensure_eda_control_files
ensure_dir "${VELOCITY_INPUT_DIR}" "${VELOCITY_LOOM_DIR}" "${VELOCITY_OUTPUT_DIR}"

VELOCITY_INPUT_RERAN=0
VELOCITY_STAGE3_RERAN=0

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

run_velocity_stage3_if_stale() {
  local runner_fn="$1"
  local script_path="$2"
  local output_path="$3"
  shift 3 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "Running ${script_path}"
    "${runner_fn}" "${script_path}"
    VELOCITY_STAGE3_RERAN=1
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

hold_for_gate velocity_inputs

run_velocity_stage3_if_stale \
  run_scvelo \
  "${WORKFLOW_ROOT}/04python/10c_scvelo_dynamical.py" \
  "${MODULE_10C_MANIFEST}" \
  "${MODULE_10A_MANIFEST}" \
  "${MODULE_10B_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}" \
  "${WORKFLOW_ROOT}/04python/10c_scvelo_dynamical.py"
require_manifest_output "${MODULE_10C_MANIFEST}" "scvelo_index_tsv" >/dev/null
require_manifest_output "${MODULE_10C_MANIFEST}" "velocity_qc_tsv" >/dev/null

run_velocity_stage3_if_stale \
  run_r_main \
  "${WORKFLOW_ROOT}/05single_script/10d_velocyto_steady_state.R" \
  "${MODULE_10D_MANIFEST}" \
  "${MODULE_10A_MANIFEST}" \
  "${MODULE_10B_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}" \
  "${WORKFLOW_ROOT}/05single_script/10d_velocyto_steady_state.R"
require_manifest_output "${MODULE_10D_MANIFEST}" "velocyto_steady_index_tsv" >/dev/null

run_velocity_stage3_if_stale \
  run_scvelo \
  "${WORKFLOW_ROOT}/04python/10e_scvelo_drivers.py" \
  "${MODULE_10E_MANIFEST}" \
  "${MODULE_10C_MANIFEST}" \
  "${WORKFLOW_ROOT}/04python/10e_scvelo_drivers.py"
require_manifest_output "${MODULE_10E_MANIFEST}" "velocity_driver_index_tsv" >/dev/null
require_manifest_output "${MODULE_10E_MANIFEST}" "velocity_driver_overlap_tsv" >/dev/null

run_velocity_stage3_if_stale \
  run_scvelo \
  "${WORKFLOW_ROOT}/04python/10f_cellrank_fate.py" \
  "${MODULE_10F_MANIFEST}" \
  "${MODULE_10C_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}" \
  "${WORKFLOW_ROOT}/04python/10f_cellrank_fate.py"
require_manifest_output "${MODULE_10F_MANIFEST}" "cellrank_index_tsv" >/dev/null

run_velocity_stage3_if_stale \
  run_r_main \
  "${WORKFLOW_ROOT}/05single_script/10h_velocity_root_terminal.R" \
  "${MODULE_10H_MANIFEST}" \
  "${MODULE_10C_MANIFEST}" \
  "${MODULE_10F_MANIFEST}" \
  "${WORKFLOW_ROOT}/05single_script/10h_velocity_root_terminal.R"
require_manifest_output "${MODULE_10H_MANIFEST}" "velocity_root_terminal_index_tsv" >/dev/null

sync_workflow_gate_statuses
update_workflow_status \
  "10_velocity_stage3_completed" \
  "run 09_trajectory.sh next, then 10_velocity_finalize.sh to generate 10g/10i velocity EDA" \
  "status.10a_velocity_loom_completed=true" \
  "status.10b_velocity_reference_completed=true" \
  "status.10_velocity_inputs_completed=true" \
  "status.10c_scvelo_dynamical_completed=true" \
  "status.10d_velocyto_steady_completed=true" \
  "status.10e_scvelo_drivers_completed=true" \
  "status.10f_cellrank_fate_completed=true" \
  "status.10h_velocity_root_terminal_completed=true" \
  "status.10_velocity_stage3_completed=true" \
  "status.velocity_inputs_gate_passed=$(eda_gate_passed velocity_inputs && echo true || echo false)"
