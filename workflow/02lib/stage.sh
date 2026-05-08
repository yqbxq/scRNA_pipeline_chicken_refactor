is_stale_output() {
  local output_path="$1"
  shift || true

  if [[ ! -e "${output_path}" ]]; then
    return 0
  fi

  local dep_path
  for dep_path in "$@"; do
    [[ -n "${dep_path}" ]] || continue
    if [[ -e "${dep_path}" && "${dep_path}" -nt "${output_path}" ]]; then
      return 0
    fi
  done

  return 1
}

check_stage_deps() {
  local stage_id="$1"
  case "${stage_id}" in
    00_ortholog|alignment|register_delivery|validate_metadata|audit_inputs|standardize_inputs|input_summary|install_envs)
      return 0
      ;;
    01_build_raw)
      require_status_flag_or_warn \
        "status.main_ready" \
        "01_build_raw 需要 main_ready=true。请先完成 metadata、audit、standardize 和 input_summary。"
      require_status_flag_or_warn \
        "status.00_ortholog_completed" \
        "01_build_raw 需要先完成 00_ortholog，并通过 manifest 暴露 cc_genes。"
      ;;
    02_qc)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.01_pre_qc_eda_completed" \
        "02_qc 需要先完成 01_build_raw，并生成 pre-QC EDA 报告。"
      require_status_flag_or_warn \
        "status.pre_qc_gate_passed" \
        "02_qc 被 pre_qc gate 阻断。请先审阅 01b pre-QC 报告，并在 eda_gates.tsv 中批准 pre_qc。"
      ;;
    03_panorama)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.post_qc_gate_passed" \
        "03_panorama 被 post_qc gate 阻断。请先审阅 02 post-QC 报告，并在 eda_gates.tsv 中批准 post_qc。"
      require_status_flag_or_warn \
        "status.02_qc_completed" \
        "03_panorama 需要 02_qc 完成（含 ambient + qc + doublet）。"
      ;;
    04_subcluster)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.annotation_gate_passed" \
        "04_subcluster 被 annotation gate 阻断。请先审阅 03e panorama 注释报告，并在 eda_gates.tsv 中批准 annotation。"
      require_status_flag_or_warn \
        "status.03_panorama_completed" \
        "04_subcluster 需要 03_panorama 整链完成。"
      ;;
    04d_cluster_robustness)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.subcluster_gate_passed" \
        "04d_cluster_robustness 被 subcluster gate 阻断。请先审阅 04c subcluster EDA 报告，并在 eda_gates.tsv 中批准 subcluster。"
      require_status_flag_or_warn \
        "status.04_subcluster_completed" \
        "04d_cluster_robustness 需要 04_subcluster 完成。"
      ;;
    05_deg)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.annotation_gate_passed" \
        "05_deg 被 annotation gate 阻断。请先审阅 03e panorama 注释报告，并在 eda_gates.tsv 中批准 annotation。"
      require_status_flag_or_warn \
        "status.03_panorama_completed" \
        "05_deg 需要 03_panorama 整链完成。"
      require_status_flag_or_warn \
        "status.subcluster_gate_passed" \
        "05_deg 被 subcluster gate 阻断。请先审阅 04c 子层注释报告，并在 eda_gates.tsv 中批准 subcluster。"
      require_status_flag_or_warn \
        "status.04_subcluster_completed" \
        "05_deg 需要 04_subcluster 完成。"
      ;;
    06_enrichment)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.05_deg_completed" \
        "06_enrichment 需要先完成 05_deg。"
      require_status_flag_or_warn \
        "status.deg_gate_passed" \
        "06_enrichment 被 deg gate 阻断。请先审阅 05d DEG 报告，并在 eda_gates.tsv 中批准 deg。"
      ;;
    07_communication)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.annotation_gate_passed" \
        "07_communication 被 annotation gate 阻断。请先审阅 03e panorama 注释报告，并在 eda_gates.tsv 中批准 annotation。"
      require_status_flag_or_warn \
        "status.03_panorama_completed" \
        "07_communication 需要 03_panorama 整链完成。"
      require_status_flag_or_warn \
        "status.00_ortholog_completed" \
        "07_communication 需要先完成 00_ortholog，并通过 manifest 暴露 human_best。"
      if [[ -f "${WORKFLOW_STATUS_FILE}" ]]; then
        local subcluster_done
        subcluster_done="$(workflow_status_get "status.04_subcluster_completed" 2>/dev/null || true)"
        if [[ "${subcluster_done}" != "true" ]]; then
          warn "04_subcluster 尚未标记完成；07 将只使用当前已注释并可发现的 layer，panorama 不受影响。"
        fi
      fi
      ;;
    08_regulation)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.annotation_gate_passed" \
        "08_regulation 被 annotation gate 阻断。请先审阅 03e panorama 注释报告，并在 eda_gates.tsv 中批准 annotation。"
      require_status_flag_or_warn \
        "status.03_panorama_completed" \
        "08_regulation 需要 03_panorama 整链完成。"
      require_status_flag_or_warn \
        "status.00_ortholog_completed" \
        "08_regulation 需要先完成 00_ortholog，并通过 manifest 暴露 human_best。"
      if [[ -f "${WORKFLOW_STATUS_FILE}" ]]; then
        local subcluster_done
        subcluster_done="$(workflow_status_get "status.04_subcluster_completed" 2>/dev/null || true)"
        if [[ "${subcluster_done}" != "true" ]]; then
          warn "04_subcluster 尚未标记完成；08 默认只使用 panorama。需要子层时请在 REGULATION_LAYERS 中显式列出已注释 layer。"
        fi
      fi
      ;;
    09_trajectory)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.annotation_gate_passed" \
        "09_trajectory 被 annotation gate 阻断。请先审阅 03e panorama 注释报告，并在 eda_gates.tsv 中批准 annotation。"
      require_status_flag_or_warn \
        "status.03_panorama_completed" \
        "09_trajectory 需要 03_panorama 整链完成。"
      require_status_flag_or_warn \
        "status.subcluster_gate_passed" \
        "09_trajectory 被 subcluster gate 阻断。请先审阅 04c subcluster EDA 报告，并在 eda_gates.tsv 中批准 subcluster。"
      require_status_flag_or_warn \
        "status.04_subcluster_completed" \
        "09_trajectory 需要 04_subcluster 完成，以便解析 GC/TC subcluster layer。"
      ;;
    *)
      warn "未定义 ${stage_id} 的依赖规则，按无依赖继续。"
      ;;
  esac
}

hold_for_gate() {
  local gate_id="$1"
  sync_workflow_gate_statuses
  if ! eda_gate_passed "${gate_id}"; then
    die "${gate_id} gate 未批准。请审阅对应 EDA 报告，并在 ${EDA_GATE_FILE} 中将 ${gate_id} 设置为 approved。"
  fi
}

manifest_output_path() {
  local manifest_path="$1"
  local output_key="$2"
  local python_bin
  python_bin="$(detect_python)"

  "${python_bin}" - "${manifest_path}" "${output_key}" <<'PY'
import json
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1])
output_key = sys.argv[2]
if not manifest_path.exists():
    raise SystemExit(f"missing manifest: {manifest_path}")

manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
entry = (manifest.get("outputs") or {}).get(output_key) or {}
raw_path = entry.get("path")
if not raw_path:
    raise SystemExit(f"manifest missing output key: {output_key}")

path = Path(raw_path)
if not path.is_absolute():
    path = Path(manifest.get("base_dir") or manifest_path.parent) / path
print(path)
PY
}

manifest_has_output_key() {
  local manifest_path="$1"
  local output_key="$2"
  local python_bin
  python_bin="$(detect_python)"

  "${python_bin}" - "${manifest_path}" "${output_key}" <<'PY'
import json
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1])
output_key = sys.argv[2]
if not manifest_path.exists():
    raise SystemExit(1)

try:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
except Exception:
    raise SystemExit(1)

outputs = manifest.get("outputs") or {}
entry = outputs.get(output_key) or {}
raise SystemExit(0 if entry.get("path") else 1)
PY
}

require_manifest_output() {
  local manifest_path="$1"
  local output_key="$2"
  local output_path
  output_path="$(manifest_output_path "${manifest_path}" "${output_key}")"
  [[ -e "${output_path}" ]] || die "manifest 输出不存在: ${output_key} -> ${output_path}"
  echo "${output_path}"
}

run_stage_if_manifest_key_missing() {
  run_stage_if_manifest_key_missing_with_runner run_r_main "$@"
}

run_stage_if_manifest_key_missing_interaction() {
  run_stage_if_manifest_key_missing_with_runner run_r_interaction "$@"
}

run_stage_if_manifest_key_missing_with_runner() {
  local runner_fn="$1"
  local script_path="$2"
  local manifest_path="$3"
  local output_key="$4"
  shift 4 || true

  if ! manifest_has_output_key "${manifest_path}" "${output_key}"; then
    echo "manifest 缺少 ${output_key}，运行 ${script_path}"
    "${runner_fn}" "${script_path}"
    return
  fi

  local output_path
  output_path="$(manifest_output_path "${manifest_path}" "${output_key}")"
  if is_stale_output "${output_path}" "$@"; then
    echo "manifest 输出 ${output_key} 缺失或过期，运行 ${script_path}"
    "${runner_fn}" "${script_path}"
  else
    echo "manifest 输出已存在且未过期，跳过: ${output_key} -> ${output_path}"
  fi
}

run_stage_if_stale() {
  run_stage_if_stale_with_runner run_r_main "$@"
}

run_stage_if_stale_interaction() {
  run_stage_if_stale_with_runner run_r_interaction "$@"
}

run_stage_if_stale_with_runner() {
  local runner_fn="$1"
  local script_path="$2"
  local output_path="$3"
  shift 3 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "运行 ${script_path}"
    "${runner_fn}" "${script_path}"
  else
    echo "已存在且未过期，跳过: ${output_path}"
  fi
}

sync_workflow_gate_statuses() {
  local current_stage
  local next_step
  current_stage="$(workflow_status_get "current_stage" 2>/dev/null || echo gate_status_synced)"
  next_step="$(workflow_status_get "next_step" 2>/dev/null || echo "")"

  update_workflow_status \
    "${current_stage}" \
    "${next_step}" \
    "status.pre_qc_gate_passed=$(eda_gate_passed pre_qc && echo true || echo false)" \
    "status.post_qc_gate_passed=$(eda_gate_passed post_qc && echo true || echo false)" \
    "status.integration_gate_passed=$(eda_gate_passed integration && echo true || echo false)" \
    "status.annotation_gate_passed=$(eda_gate_passed annotation && echo true || echo false)" \
    "status.subcluster_gate_passed=$(eda_gate_passed subcluster && echo true || echo false)" \
    "status.deg_gate_passed=$(eda_gate_passed deg && echo true || echo false)" \
    "status.communication_gate_passed=$(eda_gate_passed communication && echo true || echo false)" \
    "status.regulation_gate_passed=$(eda_gate_passed regulation && echo true || echo false)" \
    "status.trajectory_inputs_gate_passed=$(eda_gate_passed trajectory_inputs && echo true || echo false)" \
    "status.trajectory_methods_gate_passed=$(eda_gate_passed trajectory_methods && echo true || echo false)" \
    "status.trajectory_finalize_gate_passed=$(eda_gate_passed trajectory_finalize && echo true || echo false)" \
    "status.velocity_inputs_gate_passed=$(eda_gate_passed velocity_inputs && echo true || echo false)" \
    "status.velocity_finalize_gate_passed=$(eda_gate_passed velocity_finalize && echo true || echo false)"
}

update_workflow_status() {
  local stage_name="$1"
  local next_step="$2"
  shift 2 || true

  prepare_project_state_dirs

  local python_bin
  local project_input_mode
  python_bin="$(detect_python)"
  project_input_mode="${PROJECT_INPUT_MODE:-}"
  if [[ -z "${project_input_mode}" ]] && declare -F detect_project_input_mode >/dev/null 2>&1; then
    project_input_mode="$(detect_project_input_mode)"
  fi

  env \
    WORKFLOW_STATUS_FILE="${WORKFLOW_STATUS_FILE}" \
    WF_STAGE="${stage_name}" \
    WF_NEXT_STEP="${next_step}" \
    WF_PROJECT_ROOT="${PROJECT_ROOT}" \
    WF_CONFIG_FILE="${CONFIG_FILE}" \
    WF_PROJECT_INPUT_MODE="${project_input_mode}" \
    WF_SAMPLE_SHEET="${SAMPLE_SHEET}" \
    WF_CANONICAL_SAMPLE_SHEET="${CANONICAL_SAMPLE_SHEET}" \
    WF_COMPARISON_SHEET="${COMPARISON_SHEET}" \
    WF_INPUT_INVENTORY_FILE="${INPUT_INVENTORY_FILE}" \
    WF_BRANCH_READINESS_FILE="${BRANCH_READINESS_FILE}" \
    WF_INTAKE_SUMMARY_FILE="${INTAKE_SUMMARY_FILE}" \
    WF_STATUS_DIR="${STATUS_DIR}" \
    WF_EDA_REPORT_DIR="${EDA_REPORT_DIR}" \
    WF_AMBIENT_REPORT_DIR="${AMBIENT_REPORT_DIR}" \
    WF_QC_THRESHOLD_FILE="${QC_THRESHOLD_FILE}" \
    WF_EDA_GATE_FILE="${EDA_GATE_FILE}" \
    WF_ANNOTATION_HUB_PATH="${ANNOTATION_HUB_PATH}" \
    WF_WAIVER_FILE="${WAIVER_FILE}" \
    "${python_bin}" - "$@" <<'PY'
import json
import os
from pathlib import Path
from datetime import datetime

status_path = Path(os.environ["WORKFLOW_STATUS_FILE"])
if status_path.exists():
    data = json.loads(status_path.read_text(encoding="utf-8"))
else:
    data = {}

def coerce_value(raw: str):
    lowered = raw.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if lowered == "null":
        return None
    return raw

def set_nested(obj, dotted_key, value):
    current = obj
    parts = dotted_key.split(".")
    for key in parts[:-1]:
      if key not in current or not isinstance(current[key], dict):
          current[key] = {}
      current = current[key]
    current[parts[-1]] = value

def get_nested(obj, dotted_key, default=None):
    current = obj
    for key in dotted_key.split("."):
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current

for item in os.sys.argv[1:]:
    if "=" not in item:
        continue
    key, raw_value = item.split("=", 1)
    set_nested(data, key, coerce_value(raw_value))

data["updated_at"] = datetime.now().astimezone().isoformat(timespec="seconds")
data["project_root"] = os.environ["WF_PROJECT_ROOT"]
data["config_file"] = os.environ["WF_CONFIG_FILE"]
data["current_stage"] = os.environ["WF_STAGE"]
data["next_step"] = os.environ["WF_NEXT_STEP"]
project_input_mode = os.environ.get("WF_PROJECT_INPUT_MODE", "")
if project_input_mode:
    data["project_input_mode"] = project_input_mode
else:
    data.setdefault("project_input_mode", "")

paths = data.setdefault("paths", {})
paths["sample_sheet"] = os.environ["WF_SAMPLE_SHEET"]
paths["canonical_sample_sheet"] = os.environ["WF_CANONICAL_SAMPLE_SHEET"]
paths["comparison_sheet"] = os.environ["WF_COMPARISON_SHEET"]
paths["input_inventory"] = os.environ["WF_INPUT_INVENTORY_FILE"]
paths["branch_readiness"] = os.environ["WF_BRANCH_READINESS_FILE"]
paths["intake_summary"] = os.environ["WF_INTAKE_SUMMARY_FILE"]
paths["status_dir"] = os.environ["WF_STATUS_DIR"]
paths["eda_report_dir"] = os.environ["WF_EDA_REPORT_DIR"]
paths["ambient_report_dir"] = os.environ["WF_AMBIENT_REPORT_DIR"]
paths["qc_threshold_file"] = os.environ["WF_QC_THRESHOLD_FILE"]
paths["eda_gate_file"] = os.environ["WF_EDA_GATE_FILE"]

annotation_path = os.environ["WF_ANNOTATION_HUB_PATH"]
annotation_exists = Path(annotation_path).exists()
data["annotation_hub"] = {
    "path": annotation_path,
    "exists": annotation_exists,
    "produced_by": "20_run_main_pipeline.sh / r/03_annotation.R",
}

waiver_file = Path(os.environ["WF_WAIVER_FILE"])
waivers = []
if waiver_file.exists():
    lines = waiver_file.read_text(encoding="utf-8").splitlines()
    if lines:
        headers = lines[0].split("\t")
        for line in lines[1:]:
            if not line.strip():
                continue
            values = line.split("\t")
            row = {headers[i]: values[i] if i < len(values) else "" for i in range(len(headers))}
            waivers.append(row)
data["waivers"] = waivers

status = data.setdefault("status", {})
for key in (
    "metadata_valid",
    "standardized_inputs",
    "input_eda_complete",
    "ambient_branch_complete",
    "pre_qc_eda_complete",
    "pre_qc_gate_passed",
    "post_qc_eda_complete",
    "post_qc_gate_passed",
    "integration_eda_complete",
    "integration_gate_passed",
    "annotation_eda_complete",
    "annotation_gate_passed",
    "subcluster_gate_passed",
    "deg_gate_passed",
    "communication_gate_passed",
    "regulation_gate_passed",
    "trajectory_inputs_gate_passed",
    "trajectory_methods_gate_passed",
    "trajectory_finalize_gate_passed",
    "velocity_inputs_gate_passed",
    "velocity_finalize_gate_passed",
    "05_deg_completed",
    "09a_trajectory_inputs_completed",
    "09b_trajectory_inputs_eda_completed",
    "09_trajectory_inputs_completed",
    "09c_trajectory_root_completed",
    "09d_trajectory_slingshot_completed",
    "09e_trajectory_monocle3_completed",
    "09e2_trajectory_monocle2_completed",
    "09f_trajectory_paga_dpt_completed",
    "09g_trajectory_palantir_completed",
    "09_trajectory_methods_completed",
    "09h_trajectory_tradeseq_completed",
    "09i_trajectory_consensus_completed",
    "09j_trajectory_velocity_link_completed",
    "09k_trajectory_figures_completed",
    "09m_trajectory_split_compare_completed",
    "09l_trajectory_eda_completed",
    "09_trajectory_stage3_completed",
    "09_trajectory_completed",
    "main_upstream_ready",
    "velocity_upstream_ready",
    "ambient_upstream_ready",
    "scenic_upstream_ready",
    "de_replicate_ready",
):
    status.setdefault(key, False)

status["main_ready"] = bool(status.get("metadata_valid")) and bool(status.get("standardized_inputs")) and bool(status.get("main_upstream_ready"))
status["deg_ready"] = annotation_exists and bool(status.get("annotation_gate_passed"))
status["enrichment_ready"] = bool(status.get("05_deg_completed")) and bool(status.get("deg_gate_passed"))
status["trajectory_ready"] = annotation_exists and bool(status.get("annotation_gate_passed")) and bool(status.get("subcluster_gate_passed"))
status["velocity_reference_ready"] = annotation_exists and bool(status.get("velocity_upstream_ready")) and bool(status.get("annotation_gate_passed"))
status["scenic_export_ready"] = annotation_exists and bool(status.get("scenic_upstream_ready")) and bool(status.get("annotation_gate_passed"))

status_path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
PY
}

workflow_status_get() {
  local dotted_key="$1"
  [[ -f "${WORKFLOW_STATUS_FILE}" ]] || return 1

  local python_bin
  python_bin="$(detect_python)"
  "${python_bin}" - "${WORKFLOW_STATUS_FILE}" "${dotted_key}" <<'PY'
import json
import sys
from pathlib import Path

payload = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
current = payload
for key in sys.argv[2].split("."):
    if not isinstance(current, dict) or key not in current:
        raise SystemExit(1)
    current = current[key]
if isinstance(current, bool):
    print("true" if current else "false")
elif current is None:
    print("null")
else:
    print(current)
PY
}

require_status_flag_or_warn() {
  local dotted_key="$1"
  local err_message="$2"
  if [[ ! -f "${WORKFLOW_STATUS_FILE}" ]]; then
    warn "缺少 workflow 状态文件 ${WORKFLOW_STATUS_FILE}，按兼容模式继续执行。"
    return 0
  fi

  local current_value
  current_value="$(workflow_status_get "${dotted_key}" 2>/dev/null || true)"
  if [[ "${current_value}" != "true" ]]; then
    die "${err_message}"
  fi
}
