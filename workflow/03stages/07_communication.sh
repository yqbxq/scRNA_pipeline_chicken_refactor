#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

ensure_metadata_fresh
check_stage_deps "07_communication"
hold_for_gate annotation

MODULE_00_MANIFEST="${ORTHOLOG_MANIFEST}"
MODULE_03D_MANIFEST="${MANIFEST_DIR}/03d_annotate/_manifest.json"
MODULE_04B_MANIFEST="${MANIFEST_DIR}/04b_subcluster_annotate/_manifest.json"
MODULE_04C_MANIFEST="${MANIFEST_DIR}/04c_subcluster_eda/_manifest.json"
MODULE_05D_MANIFEST="${MANIFEST_DIR}/05d_deg_eda/_manifest.json"
MODULE_07A_MANIFEST="${MANIFEST_DIR}/07a_cellchat/_manifest.json"
MODULE_07B_MANIFEST="${MANIFEST_DIR}/07b_liana_consensus/_manifest.json"
MODULE_07C_MANIFEST="${MANIFEST_DIR}/07c_nichenet/_manifest.json"
MODULE_07D_MANIFEST="${MANIFEST_DIR}/07d_communication_consensus/_manifest.json"
MODULE_07E_MANIFEST="${MANIFEST_DIR}/07e_communication_eda/_manifest.json"

COMMUNICATION_REQUIRE_FULL_PIPELINE="${COMMUNICATION_REQUIRE_FULL_PIPELINE:-no}"
LIANA_CONSENSUS_ENABLED="${LIANA_CONSENSUS_ENABLED:-yes}"
MULTINICHENET_ENABLED="${MULTINICHENET_ENABLED:-auto}"
COMMUNICATION_RERAN=0
COMMUNICATION_REPORT=""

ensure_eda_control_files
require_manifest_output "${MODULE_03D_MANIFEST}" "annotated_object" >/dev/null
require_manifest_output "${MODULE_00_MANIFEST}" "human_best" >/dev/null
require_manifest_output "${MODULE_04C_MANIFEST}" "cluster_eligibility_tsv" >/dev/null

write_placeholder_manifest_07() {
  local manifest_path="$1"
  local module_name="$2"
  local reason="$3"
  local report_path="${4:-}"
  local python_bin
  python_bin="$(detect_python)"

  if [[ -n "${report_path}" ]]; then
    ensure_dir "$(dirname "${report_path}")"
    {
      printf '# %s\n\n' "${module_name}"
      printf 'Status: placeholder\n\n'
      printf '%s\n' "${reason}"
    } > "${report_path}"
  fi

  ensure_dir "$(dirname "${manifest_path}")"
  "${python_bin}" - "${manifest_path}" "${module_name}" "${reason}" "${PROJECT_ROOT}" "${report_path}" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

manifest_path, module_name, reason, base_dir, report_path = sys.argv[1:6]
outputs = {}
if report_path:
    path = Path(report_path)
    try:
        rel = path.relative_to(Path(base_dir))
    except ValueError:
        rel = path
    outputs["report_md"] = {
        "path": str(rel),
        "format": "md",
        "module": module_name,
        "description": "placeholder communication report",
    }

manifest = {
    "module": module_name,
    "status": "placeholder",
    "reason": reason,
    "base_dir": base_dir,
    "outputs": outputs,
    "inputs": {},
    "depends_on": {},
    "created_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
}
Path(manifest_path).write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY
}

run_required_comm_r_stage() {
  local script_path="$1"
  local manifest_path="$2"
  shift 2 || true

  run_stage_if_stale_interaction \
    --script "${script_path}" \
    --manifest "${manifest_path}" \
    --inputs "$@"
  COMMUNICATION_RERAN=1
}

run_optional_comm_r_stage() {
  local module_name="$1"
  local script_path="$2"
  local manifest_path="$3"
  local report_path="${4:-}"
  shift 4 || true

  if ! is_stale_output "${manifest_path}" "$@"; then
    echo "已存在且未过期，跳过: ${manifest_path}"
    return 0
  fi

  echo "运行 ${script_path}"
  if run_r_interaction "${script_path}"; then
    COMMUNICATION_RERAN=1
    return 0
  fi

  if [[ "${COMMUNICATION_REQUIRE_FULL_PIPELINE}" == "yes" ]]; then
    die "${module_name} failed and COMMUNICATION_REQUIRE_FULL_PIPELINE=yes"
  fi

  warn "${module_name} placeholder/skip; COMMUNICATION_REQUIRE_FULL_PIPELINE=no"
  write_placeholder_manifest_07 "${manifest_path}" "${module_name}" "not yet implemented or runtime unavailable" "${report_path}"
  COMMUNICATION_RERAN=1
}

run_optional_comm_python_stage() {
  local module_name="$1"
  local script_path="$2"
  local manifest_path="$3"
  shift 3 || true
  local -a command_args=()
  local -a inputs=()

  while [[ "$#" -gt 0 ]]; do
    if [[ "${1:-}" == "--stage-inputs" ]]; then
      shift || true
      inputs=("$@")
      break
    fi
    command_args+=("$1")
    shift || true
  done

  if ! is_stale_output "${manifest_path}" "${inputs[@]}"; then
    echo "已存在且未过期，跳过: ${manifest_path}"
    return 0
  fi

  echo "运行 ${script_path}"
  local python_bin
  python_bin="$(detect_python)"
  if "${python_bin}" "${script_path}" "${command_args[@]}"; then
    COMMUNICATION_RERAN=1
    return 0
  fi

  if [[ "${COMMUNICATION_REQUIRE_FULL_PIPELINE}" == "yes" ]]; then
    die "${module_name} failed and COMMUNICATION_REQUIRE_FULL_PIPELINE=yes"
  fi

  warn "${module_name} placeholder/skip; COMMUNICATION_REQUIRE_FULL_PIPELINE=no"
  write_placeholder_manifest_07 "${manifest_path}" "${module_name}" "not yet implemented or runtime unavailable"
  COMMUNICATION_RERAN=1
}

run_required_comm_r_stage \
  "${WORKFLOW_ROOT}/05single_script/07a_cellchat.R" \
  "${MODULE_07A_MANIFEST}" \
  "${MODULE_03D_MANIFEST}" \
  "${MODULE_04B_MANIFEST}" \
  "${MODULE_04C_MANIFEST}" \
  "${MODULE_00_MANIFEST}" \
  "${COMMUNICATION_PAIRS_SHEET}"
require_manifest_output "${MODULE_07A_MANIFEST}" "cellchat_index_tsv" >/dev/null
require_manifest_output "${MODULE_07A_MANIFEST}" "cellchat_inventory_gate_log_tsv" >/dev/null

if [[ "${LIANA_CONSENSUS_ENABLED}" == "no" ]]; then
  warn "LIANA_CONSENSUS_ENABLED=no，写入 07b_liana_consensus placeholder manifest。"
  write_placeholder_manifest_07 "${MODULE_07B_MANIFEST}" "07b_liana_consensus" "disabled by LIANA_CONSENSUS_ENABLED=no"
else
  run_optional_comm_python_stage \
    "07b_liana_consensus" \
    "${WORKFLOW_ROOT}/04python/07b_liana_consensus.py" \
    "${MODULE_07B_MANIFEST}" \
    --h5ad-path "${RESULTS_DIR}/90a_export_h5ad/03d_panorama" \
    --cluster-eligibility "$(require_manifest_output "${MODULE_04C_MANIFEST}" "cluster_eligibility_tsv")" \
    --cluster-col "${COMMUNICATION_CELL_TYPE_COL}" \
    --condition-col "${LIANA_CONDITION_COL:-condition}" \
    --methods "${LIANA_METHODS_LIST}" \
    --resource "${LIANA_RESOURCE_DB}" \
    --output "${TABLE_DIR}/communication/liana_consensus/liana_consensus_lr.tsv" \
    --manifest "${MODULE_07B_MANIFEST}" \
    --ortholog-lut "${ORTHOLOG_CHICKEN_HUMAN_TSV:-${METADATA_DIR}/ortholog_chicken_human.tsv}" \
    --stage-inputs \
    "${MODULE_07A_MANIFEST}" \
    "${MODULE_04B_MANIFEST}" \
    "${MODULE_04C_MANIFEST}" \
    "${RESULTS_DIR}/90a_export_h5ad/03d_panorama" \
    "${COMMUNICATION_PAIRS_SHEET}"
fi

if [[ "${MULTINICHENET_ENABLED}" == "no" || "${RUN_NICHENET:-yes}" == "no" || "${NICHENET_MODE:-auto}" == "skip" ]]; then
  warn "MULTINICHENET_ENABLED=no、RUN_NICHENET=no 或 NICHENET_MODE=skip，写入 07c_nichenet placeholder manifest。"
  write_placeholder_manifest_07 "${MODULE_07C_MANIFEST}" "07c_nichenet" "disabled by MULTINICHENET_ENABLED/RUN_NICHENET/NICHENET_MODE"
else
  run_optional_comm_r_stage \
    "07c_nichenet" \
    "${WORKFLOW_ROOT}/05single_script/07c_nichenet.R" \
    "${MODULE_07C_MANIFEST}" \
    "" \
    "${MODULE_07A_MANIFEST}" \
    "${MODULE_05D_MANIFEST}" \
    "${MODULE_04C_MANIFEST}" \
    "${MODULE_00_MANIFEST}" \
    "${COMMUNICATION_PAIRS_SHEET}" \
    "${NICHENET_RESOURCE_DIR}"
  if [[ -f "${MODULE_07C_MANIFEST}" ]] && [[ "$(manifest_has_output_key "${MODULE_07C_MANIFEST}" "nichenet_inventory_gate_log_tsv" && echo yes || echo no)" == "yes" ]]; then
    require_manifest_output "${MODULE_07C_MANIFEST}" "nichenet_inventory_gate_log_tsv" >/dev/null
  fi
fi

run_optional_comm_r_stage \
  "07d_communication_consensus" \
  "${WORKFLOW_ROOT}/05single_script/07d_communication_consensus.R" \
  "${MODULE_07D_MANIFEST}" \
  "" \
  "${MODULE_07A_MANIFEST}" \
  "${MODULE_07B_MANIFEST}" \
  "${MODULE_07C_MANIFEST}" \
  "${COMMOT_SPATIAL_SUMMARY_TSV:-}"

run_optional_comm_r_stage \
  "07e_communication_eda" \
  "${WORKFLOW_ROOT}/05single_script/07e_communication_eda.R" \
  "${MODULE_07E_MANIFEST}" \
  "${COMMUNICATION_REPORT_DIR}/report.md" \
  "${MODULE_07A_MANIFEST}" \
  "${MODULE_07B_MANIFEST}" \
  "${MODULE_07C_MANIFEST}" \
  "${MODULE_07D_MANIFEST}" \
  "${COMMUNICATION_PAIRS_SHEET}"

COMMUNICATION_REPORT="$(require_manifest_output "${MODULE_07E_MANIFEST}" "report_md")"

export_h5ad_for_gate \
  "07e_communication" \
  "${MODULE_03D_MANIFEST}" \
  "annotated_object" \
  "scrna" \
  "07e_communication"

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
  "status.07b_liana_consensus_completed=$([[ -f "${MODULE_07B_MANIFEST}" ]] && echo true || echo false)" \
  "status.07c_nichenet_completed=$([[ -f "${MODULE_07C_MANIFEST}" ]] && echo true || echo false)" \
  "status.07d_communication_consensus_completed=$([[ -f "${MODULE_07D_MANIFEST}" ]] && echo true || echo false)" \
  "status.07e_communication_eda_completed=$([[ -f "${MODULE_07E_MANIFEST}" ]] && echo true || echo false)" \
  "status.07_communication_completed=true" \
  "status.communication_gate_passed=$(eda_gate_passed communication && echo true || echo false)"
