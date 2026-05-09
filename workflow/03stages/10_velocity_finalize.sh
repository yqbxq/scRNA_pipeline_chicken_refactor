#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

ensure_metadata_fresh
check_stage_deps "10_velocity_finalize"

MODULE_09D_MANIFEST="${MODULE_09D_MANIFEST:-${MANIFEST_DIR}/09d_trajectory_slingshot/_manifest.json}"
MODULE_10A_MANIFEST="${MANIFEST_DIR}/10a_run_velocyto/_manifest.json"
MODULE_10B_MANIFEST="${MANIFEST_DIR}/10b_prepare_velocity_reference/_manifest.json"
MODULE_10C_MANIFEST="${MANIFEST_DIR}/10c_scvelo_dynamical/_manifest.json"
MODULE_10D_MANIFEST="${MANIFEST_DIR}/10d_velocyto_steady_state/_manifest.json"
MODULE_10E_MANIFEST="${MANIFEST_DIR}/10e_scvelo_drivers/_manifest.json"
MODULE_10F_MANIFEST="${MANIFEST_DIR}/10f_cellrank_fate/_manifest.json"
MODULE_10G_MANIFEST="${MANIFEST_DIR}/10g_velocity_consistency/_manifest.json"
MODULE_10H_MANIFEST="${MANIFEST_DIR}/10h_velocity_root_terminal/_manifest.json"
MODULE_10I_MANIFEST="${MANIFEST_DIR}/10i_velocity_eda/_manifest.json"

ensure_eda_control_files
ensure_dir "${VELOCITY_REPORT_DIR}"

VELOCITY_FINAL_RERAN=0

run_velocity_final_if_stale() {
  local runner_fn="$1"
  local script_path="$2"
  local output_path="$3"
  shift 3 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "Running ${script_path}"
    "${runner_fn}" "${script_path}"
    VELOCITY_FINAL_RERAN=1
  else
    echo "Existing output is current; skipping: ${output_path}"
  fi
}

run_velocity_final_if_stale \
  run_r_main \
  "${WORKFLOW_ROOT}/05single_script/10g_velocity_consistency.R" \
  "${MODULE_10G_MANIFEST}" \
  "${MODULE_10C_MANIFEST}" \
  "${MODULE_10D_MANIFEST}" \
  "${MODULE_10E_MANIFEST}" \
  "${MODULE_10F_MANIFEST}" \
  "${MODULE_09D_MANIFEST}" \
  "${WORKFLOW_ROOT}/05single_script/10g_velocity_consistency.R"
require_manifest_output "${MODULE_10G_MANIFEST}" "velocity_consistency_index_tsv" >/dev/null
require_manifest_output "${MODULE_10G_MANIFEST}" "velocity_consistency_summary_tsv" >/dev/null
require_manifest_output "${MODULE_10G_MANIFEST}" "velocity_split_compare_tsv" >/dev/null

run_velocity_final_if_stale \
  run_r_main \
  "${WORKFLOW_ROOT}/05single_script/10i_velocity_eda.R" \
  "${MODULE_10I_MANIFEST}" \
  "${MODULE_10A_MANIFEST}" \
  "${MODULE_10B_MANIFEST}" \
  "${MODULE_10C_MANIFEST}" \
  "${MODULE_10D_MANIFEST}" \
  "${MODULE_10E_MANIFEST}" \
  "${MODULE_10F_MANIFEST}" \
  "${MODULE_10G_MANIFEST}" \
  "${MODULE_10H_MANIFEST}" \
  "${WORKFLOW_ROOT}/05single_script/10i_velocity_eda.R"
VELOCITY_FINAL_REPORT="$(require_manifest_output "${MODULE_10I_MANIFEST}" "report_md")"
require_manifest_output "${MODULE_10I_MANIFEST}" "velocity_module_status_tsv" >/dev/null
require_manifest_output "${MODULE_10I_MANIFEST}" "velocity_triage_tsv" >/dev/null

if [[ "${VELOCITY_FINAL_RERAN}" == "1" ]]; then
  set_eda_gate_status \
    "velocity_finalize" \
    "pending" \
    "" \
    "10 RNA velocity final report completed; review ${VELOCITY_FINAL_REPORT}, then approve velocity_finalize."
fi

sync_workflow_gate_statuses
update_workflow_status \
  "10_velocity_completed" \
  "review ${VELOCITY_FINAL_REPORT} and approve velocity_finalize before treating module 10 as final" \
  "status.10a_velocity_loom_completed=true" \
  "status.10b_velocity_reference_completed=true" \
  "status.10_velocity_inputs_completed=true" \
  "status.10c_scvelo_dynamical_completed=true" \
  "status.10d_velocyto_steady_completed=true" \
  "status.10e_scvelo_drivers_completed=true" \
  "status.10f_cellrank_fate_completed=true" \
  "status.10g_velocity_consistency_completed=true" \
  "status.10h_velocity_root_terminal_completed=true" \
  "status.10_velocity_stage3_completed=true" \
  "status.10i_velocity_eda_completed=true" \
  "status.10_velocity_completed=true" \
  "status.velocity_inputs_gate_passed=$(eda_gate_passed velocity_inputs && echo true || echo false)" \
  "status.velocity_finalize_gate_passed=$(eda_gate_passed velocity_finalize && echo true || echo false)"
