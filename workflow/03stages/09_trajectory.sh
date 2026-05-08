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
MODULE_09C_MANIFEST="${MANIFEST_DIR}/09c_trajectory_root/_manifest.json"
MODULE_09D_MANIFEST="${MANIFEST_DIR}/09d_trajectory_slingshot/_manifest.json"
MODULE_09E_MANIFEST="${MANIFEST_DIR}/09e_trajectory_monocle3/_manifest.json"
MODULE_09E2_MANIFEST="${MANIFEST_DIR}/09e2_trajectory_monocle2/_manifest.json"
MODULE_09F_MANIFEST="${MANIFEST_DIR}/09f_trajectory_paga_dpt/_manifest.json"
MODULE_09G_MANIFEST="${MANIFEST_DIR}/09g_trajectory_palantir/_manifest.json"
MODULE_09H_MANIFEST="${MANIFEST_DIR}/09h_trajectory_tradeseq/_manifest.json"
MODULE_09I_MANIFEST="${MANIFEST_DIR}/09i_trajectory_consensus/_manifest.json"
MODULE_09J_MANIFEST="${MANIFEST_DIR}/09j_trajectory_velocity_link/_manifest.json"

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
TRAJECTORY_METHOD_RERAN=0
TRAJECTORY_STAGE3_RERAN=0

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

run_trajectory_method_if_stale() {
  local script_path="$1"
  local output_path="$2"
  shift 2 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "运行 ${script_path}"
    run_r_main "${script_path}"
    TRAJECTORY_METHOD_RERAN=1
  else
    echo "已存在且未过期，跳过: ${output_path}"
  fi
}

run_trajectory_stage3_if_stale() {
  local script_path="$1"
  local output_path="$2"
  shift 2 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "运行 ${script_path}"
    run_r_main "${script_path}"
    TRAJECTORY_STAGE3_RERAN=1
  else
    echo "已存在且未过期，跳过: ${output_path}"
  fi
}

trajectory_has_monocle2_optin() {
  awk -F '\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) idx[$i] = i
      next
    }
    NR > 1 && $0 !~ /^#/ {
      method = tolower($(idx["method"]))
      enabled = tolower($(idx["enabled"]))
      extra = tolower($(idx["methods_extra"]))
      if (method == "trajectory" && enabled != "no" && extra ~ /(^|[,;[:space:]])[+]?monocle2([,;[:space:]]|$)/) found = 1
    }
    END { exit(found ? 0 : 1) }
  ' "${TRAJECTORY_PAIRS_SHEET}"
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

hold_for_gate trajectory_inputs

run_trajectory_method_if_stale \
  "${WORKFLOW_ROOT}/05single_script/09c_trajectory_root.R" \
  "${MODULE_09C_MANIFEST}" \
  "${MODULE_09A_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}" \
  "${VELOCITY_OUTPUT_DIR}"
require_manifest_output "${MODULE_09C_MANIFEST}" "trajectory_root_index" >/dev/null

run_trajectory_method_if_stale \
  "${WORKFLOW_ROOT}/05single_script/09d_trajectory_slingshot.R" \
  "${MODULE_09D_MANIFEST}" \
  "${MODULE_09A_MANIFEST}" \
  "${MODULE_09C_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}"
require_manifest_output "${MODULE_09D_MANIFEST}" "slingshot_index_tsv" >/dev/null

run_trajectory_method_if_stale \
  "${WORKFLOW_ROOT}/05single_script/09e_trajectory_monocle3.R" \
  "${MODULE_09E_MANIFEST}" \
  "${MODULE_09A_MANIFEST}" \
  "${MODULE_09C_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}"
require_manifest_output "${MODULE_09E_MANIFEST}" "monocle3_index_tsv" >/dev/null

run_trajectory_method_if_stale \
  "${WORKFLOW_ROOT}/05single_script/09f_trajectory_paga_dpt.R" \
  "${MODULE_09F_MANIFEST}" \
  "${MODULE_09A_MANIFEST}" \
  "${MODULE_09C_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}"
require_manifest_output "${MODULE_09F_MANIFEST}" "paga_dpt_index_tsv" >/dev/null

if [[ "${RUN_PALANTIR:-yes}" != "no" ]]; then
  run_trajectory_method_if_stale \
    "${WORKFLOW_ROOT}/05single_script/09g_trajectory_palantir.R" \
    "${MODULE_09G_MANIFEST}" \
    "${MODULE_09A_MANIFEST}" \
    "${MODULE_09C_MANIFEST}" \
    "${TRAJECTORY_PAIRS_SHEET}"
  require_manifest_output "${MODULE_09G_MANIFEST}" "palantir_index_tsv" >/dev/null
else
  warn "RUN_PALANTIR=no，跳过 09g Palantir。"
fi

if trajectory_has_monocle2_optin; then
  run_trajectory_method_if_stale \
    "${WORKFLOW_ROOT}/05single_script/09e2_trajectory_monocle2.R" \
    "${MODULE_09E2_MANIFEST}" \
    "${MODULE_09A_MANIFEST}" \
    "${MODULE_09C_MANIFEST}" \
    "${TRAJECTORY_PAIRS_SHEET}"
  require_manifest_output "${MODULE_09E2_MANIFEST}" "monocle2_index_tsv" >/dev/null
else
  echo "trajectory_pairs.tsv 未启用 methods_extra=+monocle2，跳过 09e2。"
fi

if [[ "${TRAJECTORY_METHOD_RERAN}" == "1" ]]; then
  set_eda_gate_status \
    "trajectory_methods" \
    "pending" \
    "" \
    "09c-09g trajectory methods completed; review root and method status before approving trajectory_methods."
fi

hold_for_gate trajectory_methods

run_trajectory_stage3_if_stale \
  "${WORKFLOW_ROOT}/05single_script/09h_trajectory_tradeseq.R" \
  "${MODULE_09H_MANIFEST}" \
  "${MODULE_09A_MANIFEST}" \
  "${MODULE_09D_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}"
require_manifest_output "${MODULE_09H_MANIFEST}" "tradeseq_index_tsv" >/dev/null

run_trajectory_stage3_if_stale \
  "${WORKFLOW_ROOT}/05single_script/09i_trajectory_consensus.R" \
  "${MODULE_09I_MANIFEST}" \
  "${MODULE_09D_MANIFEST}" \
  "${MODULE_09E_MANIFEST}" \
  "${MODULE_09F_MANIFEST}" \
  "${MODULE_09G_MANIFEST}" \
  "${MODULE_09E2_MANIFEST}" \
  "${TRAJECTORY_PAIRS_SHEET}"
require_manifest_output "${MODULE_09I_MANIFEST}" "consensus_index_tsv" >/dev/null

run_trajectory_stage3_if_stale \
  "${WORKFLOW_ROOT}/05single_script/09j_trajectory_velocity_link.R" \
  "${MODULE_09J_MANIFEST}" \
  "${MODULE_09A_MANIFEST}" \
  "${MODULE_09D_MANIFEST}" \
  "${VELOCITY_OUTPUT_DIR}"
require_manifest_output "${MODULE_09J_MANIFEST}" "velocity_link_index_tsv" >/dev/null

sync_workflow_gate_statuses
update_workflow_status \
  "09_trajectory_stage3_completed" \
  "review 09h tradeSeq, 09i consensus, and 09j velocity-link outputs; next implement/run 09k/09m/09l final reporting" \
  "status.09a_trajectory_inputs_completed=true" \
  "status.09b_trajectory_inputs_eda_completed=true" \
  "status.09_trajectory_inputs_completed=true" \
  "status.09c_trajectory_root_completed=true" \
  "status.09d_trajectory_slingshot_completed=true" \
  "status.09e_trajectory_monocle3_completed=true" \
  "status.09f_trajectory_paga_dpt_completed=true" \
  "status.09g_trajectory_palantir_completed=$([[ -f "${MODULE_09G_MANIFEST}" ]] && echo true || echo false)" \
  "status.09e2_trajectory_monocle2_completed=$([[ -f "${MODULE_09E2_MANIFEST}" ]] && echo true || echo false)" \
  "status.09_trajectory_methods_completed=true" \
  "status.09h_trajectory_tradeseq_completed=true" \
  "status.09i_trajectory_consensus_completed=true" \
  "status.09j_trajectory_velocity_link_completed=true" \
  "status.09_trajectory_stage3_completed=true" \
  "status.trajectory_inputs_gate_passed=$(eda_gate_passed trajectory_inputs && echo true || echo false)" \
  "status.trajectory_methods_gate_passed=$(eda_gate_passed trajectory_methods && echo true || echo false)"
