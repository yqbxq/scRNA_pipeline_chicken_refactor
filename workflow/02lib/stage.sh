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

compute_input_fingerprint() {
  local python_bin
  python_bin="$(detect_python)"

  "${python_bin}" - "$@" <<'PY'
import hashlib
import json
import os
import sys
from pathlib import Path

LARGE_SUFFIXES = {".bam", ".cram", ".fastq", ".fq", ".gz", ".h5", ".h5ad", ".loom"}
threshold_mb = int(os.environ.get("CHECKPOINT_LARGE_FILE_THRESHOLD_MB", "100") or "100")
threshold_bytes = threshold_mb * 1024 * 1024
max_manifest_depth = 5
visited_manifests = set()


def file_sha(path):
    stat = path.stat()
    suffixes = {s.lower() for s in path.suffixes}
    if stat.st_size > threshold_bytes and suffixes.intersection(LARGE_SUFFIXES):
        h = hashlib.sha256()
        with path.open("rb") as handle:
            h.update(handle.read(1024 * 1024))
        return f"large:{stat.st_size}:{stat.st_mtime_ns}:{h.hexdigest()}"

    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def is_manifest(path):
    if not path.is_file() or path.suffix.lower() != ".json":
        return False
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return False
    return isinstance(data, dict) and isinstance(data.get("outputs"), dict)


def manifest_output_paths(path, depth):
    if depth > max_manifest_depth:
        raise RuntimeError(f"manifest nesting exceeds {max_manifest_depth}: {path}")
    resolved = path.resolve()
    if resolved in visited_manifests:
        return []
    visited_manifests.add(resolved)

    data = json.loads(path.read_text(encoding="utf-8"))
    base_dir = Path(data.get("base_dir") or path.parent)
    out = []
    for entry in (data.get("outputs") or {}).values():
        if not isinstance(entry, dict) or not entry.get("path"):
            continue
        candidate = Path(str(entry["path"]))
        if not candidate.is_absolute():
            candidate = base_dir / candidate
        out.extend(expand_path(candidate, depth + 1))
    return out


def expand_path(path, depth=0):
    if not path.exists():
        return [path]
    if path.is_dir():
        files = []
        for child in path.rglob("*"):
            if not child.is_file():
                continue
            try:
                rel_parts = child.relative_to(path).parts
            except ValueError:
                rel_parts = child.parts
            if any(part.startswith(".") for part in rel_parts):
                continue
            files.append(child)
        return files
    if is_manifest(path):
        return manifest_output_paths(path, depth)
    return [path]


items = []
for raw in sys.argv[1:]:
    if not raw:
        continue
    path = Path(raw)
    for item in expand_path(path):
        key = str(item)
        if item.exists() and item.is_file():
            digest = file_sha(item)
        elif item.exists() and item.is_dir():
            digest = "dir"
        else:
            digest = "missing"
        items.append((key, digest))

agg = hashlib.sha256()
for path, digest in sorted(items):
    agg.update(f"{path}\t{digest}\n".encode("utf-8"))
print(agg.hexdigest())
PY
}

compute_script_fingerprint() {
  local script_path="${1:-}"
  local python_bin
  python_bin="$(detect_python)"

  "${python_bin}" - "${script_path}" <<'PY'
import hashlib
import re
import sys
from pathlib import Path

script = Path(sys.argv[1])
files = [script]

if script.exists():
    content = script.read_text(encoding="utf-8", errors="ignore")
    if script.suffix == ".R":
        patterns = [
            r'(?:source|source_utf8)\s*\(\s*["\']([^"\']+)["\']\s*\)',
        ]
        for pattern in patterns:
            for match in re.finditer(pattern, content):
                helper = (script.parent / match.group(1)).resolve()
                if helper.exists() and helper.is_file():
                    files.append(helper)
    elif script.suffix == ".py":
        for match in re.finditer(r'(?:from\s+helpers\.(\w+)\s+import|import\s+helpers\.(\w+))', content):
            mod = match.group(1) or match.group(2)
            helper = (script.parent / "helpers" / f"{mod}.py").resolve()
            if helper.exists() and helper.is_file():
                files.append(helper)

agg = hashlib.sha256()
seen = set()
for path in sorted(files, key=lambda p: str(p)):
    key = str(path)
    if key in seen:
        continue
    seen.add(key)
    h = hashlib.sha256()
    if path.exists() and path.is_file():
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                h.update(chunk)
        digest = h.hexdigest()
    else:
        digest = "missing"
    agg.update(f"{key}\t{digest}\n".encode("utf-8"))
print(agg.hexdigest())
PY
}

compute_params_fingerprint() {
  local param_names_csv="${1:-}"
  if [[ -z "${param_names_csv}" ]]; then
    echo ""
    return 0
  fi

  local python_bin
  python_bin="$(detect_python)"
  "${python_bin}" - "${param_names_csv}" <<'PY'
import hashlib
import os
import sys

names = [item.strip() for item in sys.argv[1].split(",") if item.strip()]
agg = hashlib.sha256()
for name in sorted(set(names)):
    agg.update(f"{name}={os.environ.get(name, '__UNSET__')}\n".encode("utf-8"))
print(agg.hexdigest())
PY
}

manifest_has_fingerprints() {
  local manifest_path="$1"
  local python_bin
  python_bin="$(detect_python)"
  "${python_bin}" - "${manifest_path}" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.exists():
    raise SystemExit(1)
try:
    data = json.loads(path.read_text(encoding="utf-8"))
except Exception:
    raise SystemExit(1)
fp = data.get("fingerprints")
if not isinstance(fp, dict) or fp.get("algorithm") != "sha256":
    raise SystemExit(1)
for key in ("input_hash", "script_hash", "params_hash"):
    if key not in fp:
        raise SystemExit(1)
raise SystemExit(0)
PY
}

read_manifest_fingerprint() {
  local manifest_path="$1"
  local fingerprint_key="$2"
  local python_bin
  python_bin="$(detect_python)"
  "${python_bin}" - "${manifest_path}" "${fingerprint_key}" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
key = sys.argv[2]
data = json.loads(path.read_text(encoding="utf-8"))
fp = data.get("fingerprints") or {}
value = fp.get(key, "")
if key == "params_names" and isinstance(value, list):
    print(",".join(str(item) for item in value if str(item)))
elif value is None:
    print("")
else:
    print(str(value))
PY
}

write_manifest_fingerprints() {
  local manifest_path="$1"
  local script_path="$2"
  local param_names_csv="$3"
  shift 3 || true

  [[ -s "${manifest_path}" ]] || return 0

  local input_hash script_hash params_hash python_bin
  if ! input_hash="$(compute_input_fingerprint "$@")"; then
    return 1
  fi
  if ! script_hash="$(compute_script_fingerprint "${script_path}")"; then
    return 1
  fi
  if ! params_hash="$(compute_params_fingerprint "${param_names_csv}")"; then
    return 1
  fi
  python_bin="$(detect_python)"

  "${python_bin}" - "${manifest_path}" "${input_hash}" "${script_hash}" "${params_hash}" "${param_names_csv}" "${CHECKPOINT_MODE:-mtime}" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

path = Path(sys.argv[1])
data = json.loads(path.read_text(encoding="utf-8"))
param_names = [item.strip() for item in sys.argv[5].split(",") if item.strip()]
data["fingerprints"] = {
    "algorithm": "sha256",
    "computed_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    "input_hash": sys.argv[2],
    "script_hash": sys.argv[3],
    "params_hash": sys.argv[4],
    "params_names": param_names,
    "checkpoint_mode": sys.argv[6],
}
path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY
}

is_stage_stale_fingerprint() {
  local script_path="$1"
  local manifest_path="$2"
  local param_names_csv="$3"
  shift 3 || true

  if ! manifest_has_fingerprints "${manifest_path}"; then
    [[ "${CHECKPOINT_FINGERPRINT_VERBOSE:-no}" == "yes" ]] && echo "fingerprint: missing saved hashes for ${manifest_path}" >&2
    return 0
  fi

  local saved_params_names
  saved_params_names="$(read_manifest_fingerprint "${manifest_path}" "params_names")"
  if [[ -n "${saved_params_names}" ]]; then
    param_names_csv="${saved_params_names}"
  fi

  local cur_input_hash cur_script_hash cur_params_hash
  local saved_input_hash saved_script_hash saved_params_hash
  if ! cur_input_hash="$(compute_input_fingerprint "$@")"; then
    return 2
  fi
  if ! cur_script_hash="$(compute_script_fingerprint "${script_path}")"; then
    return 2
  fi
  if ! cur_params_hash="$(compute_params_fingerprint "${param_names_csv}")"; then
    return 2
  fi
  if ! saved_input_hash="$(read_manifest_fingerprint "${manifest_path}" "input_hash")"; then
    return 2
  fi
  if ! saved_script_hash="$(read_manifest_fingerprint "${manifest_path}" "script_hash")"; then
    return 2
  fi
  if ! saved_params_hash="$(read_manifest_fingerprint "${manifest_path}" "params_hash")"; then
    return 2
  fi

  if [[ "${CHECKPOINT_FINGERPRINT_VERBOSE:-no}" == "yes" ]]; then
    {
      echo "fingerprint compare: ${manifest_path}"
      echo "  input:  current=${cur_input_hash} saved=${saved_input_hash}"
      echo "  script: current=${cur_script_hash} saved=${saved_script_hash}"
      echo "  params: current=${cur_params_hash} saved=${saved_params_hash}"
    } >&2
  fi

  if [[ "${cur_input_hash}" != "${saved_input_hash}" ]] || \
     [[ "${cur_script_hash}" != "${saved_script_hash}" ]] || \
     [[ "${cur_params_hash}" != "${saved_params_hash}" ]]; then
    return 0
  fi
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
    04_robustness)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.subcluster_gate_passed" \
        "04_robustness 被 subcluster gate 阻断。请先审阅 04c subcluster EDA 报告，并在 eda_gates.tsv 中批准 subcluster。"
      require_status_flag_or_warn \
        "status.04_subcluster_completed" \
        "04_robustness 需要 04_subcluster 完成。"
      require_status_flag_or_warn \
        "status.04d_cluster_robustness_completed" \
        "04_robustness 需要先完成 04d_cluster_robustness，以生成 scDesign3 target gate 表。"
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
      require_status_flag_or_warn \
        "status.scdesign3_validated_gate_passed" \
        "07_communication 被 scdesign3_validated gate 阻断。请先运行 04_robustness，审阅 04e/04f scDesign3 指标，并在 eda_gates.tsv 中批准 scdesign3_validated。"
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
    10_velocity)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.annotation_gate_passed" \
        "10_velocity 被 annotation gate 阻断。请先审阅 03e panorama 注释报告，并在 eda_gates.tsv 中批准 annotation。"
      require_status_flag_or_warn \
        "status.03_panorama_completed" \
        "10_velocity 需要 03_panorama 整链完成。"
      require_status_flag_or_warn \
        "status.subcluster_gate_passed" \
        "10_velocity 被 subcluster gate 阻断。请先审阅 04c subcluster EDA 报告，并在 eda_gates.tsv 中批准 subcluster。"
      require_status_flag_or_warn \
        "status.04_subcluster_completed" \
        "10_velocity 需要 04_subcluster 完成，以便解析 GC/TC velocity reference layer。"
      require_status_flag_or_warn \
        "status.velocity_upstream_ready" \
        "10_velocity 需要 velocity_upstream_ready=true；请先完成 intake/audit 并确认 BAM 输入可用。"
      ;;
    10_velocity_finalize)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.10_velocity_stage3_completed" \
        "10_velocity_finalize 需要先完成 10_velocity stage3。"
      ;;
    50_spatial)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.spatial_upstream_ready" \
        "50_spatial 需要 spatial_upstream_ready=true。请先完成 validate/audit/standardize/input_summary，并确认 ST bundle 可用。"
      ;;
    60_spatial_extensions)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.spatial_main_completed" \
        "60_spatial_extensions 需要先完成 50_spatial 主流程。"
      ;;
    70_joint_scrna_spatial)
      sync_workflow_gate_statuses
      require_status_flag_or_warn \
        "status.spatial_main_completed" \
        "70_joint_scrna_spatial 需要先完成 50_spatial 主流程。"
      require_status_flag_or_warn \
        "status.03_panorama_completed" \
        "70_joint_scrna_spatial 需要 scRNA panorama/annotation 完成。"
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
  shift || true

  local script_path=""
  local output_path=""
  local param_names_csv=""
  local -a inputs=()

  if [[ "${1:-}" == "--script" ]]; then
    while [[ "$#" -gt 0 ]]; do
      case "$1" in
        --script)
          script_path="$2"
          shift 2 || true
          ;;
        --manifest|--output)
          output_path="$2"
          shift 2 || true
          ;;
        --inputs)
          shift || true
          while [[ "$#" -gt 0 && "${1:-}" != --* ]]; do
            inputs+=("$1")
            shift || true
          done
          ;;
        --params)
          param_names_csv="$2"
          shift 2 || true
          ;;
        *)
          inputs+=("$1")
          shift || true
          ;;
      esac
    done
  else
    script_path="${1:-}"
    output_path="${2:-}"
    shift 2 || true
    inputs=("$@")
  fi

  [[ -n "${script_path}" ]] || die "run_stage_if_stale 缺少 script_path"
  [[ -n "${output_path}" ]] || die "run_stage_if_stale 缺少 output_path"

  local checkpoint_mode="${CHECKPOINT_MODE:-mtime}"
  local stale=1

  case "${checkpoint_mode}" in
    fingerprint)
      local fingerprint_status=0
      set +e
      is_stage_stale_fingerprint "${script_path}" "${output_path}" "${param_names_csv}" "${inputs[@]}"
      fingerprint_status="$?"
      set -e
      case "${fingerprint_status}" in
        0)
          stale=0
          ;;
        1)
          stale=1
          ;;
        *)
          case "${CHECKPOINT_FINGERPRINT_FALLBACK_ON_ERROR:-mtime}" in
            error)
              die "manifest fingerprint stale check failed: ${output_path}"
              ;;
            *)
              warn "manifest fingerprint stale check failed, falling back to mtime: ${output_path}"
              if is_stale_output "${output_path}" "${inputs[@]}"; then
                stale=0
              fi
              ;;
          esac
          ;;
      esac
      ;;
    off)
      stale=0
      ;;
    mtime|*)
      if is_stale_output "${output_path}" "${inputs[@]}"; then
        stale=0
      fi
      ;;
  esac

  if [[ "${stale}" -eq 0 ]]; then
    echo "运行 ${script_path}"
    "${runner_fn}" "${script_path}"
    if [[ "${checkpoint_mode}" == "fingerprint" ]]; then
      if ! write_manifest_fingerprints "${output_path}" "${script_path}" "${param_names_csv}" "${inputs[@]}"; then
        case "${CHECKPOINT_FINGERPRINT_FALLBACK_ON_ERROR:-mtime}" in
          error)
            die "写入 manifest fingerprint 失败: ${output_path}"
            ;;
          *)
            warn "写入 manifest fingerprint 失败，保留 mtime checkpoint: ${output_path}"
            ;;
        esac
      fi
    fi
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
    "status.velocity_finalize_gate_passed=$(eda_gate_passed velocity_finalize && echo true || echo false)" \
    "status.scdesign3_targets_gate_passed=$(eda_gate_passed scdesign3_targets && echo true || echo false)" \
    "status.scdesign3_validated_gate_passed=$(eda_gate_passed scdesign3_validated && echo true || echo false)" \
    "status.spatial_pre_qc_gate_passed=$(eda_gate_passed spatial_pre_qc && echo true || echo false)" \
    "status.spatial_post_qc_gate_passed=$(eda_gate_passed spatial_post_qc && echo true || echo false)" \
    "status.spatial_integration_gate_passed=$(eda_gate_passed spatial_integration && echo true || echo false)" \
    "status.spatial_region_annotation_gate_passed=$(eda_gate_passed spatial_region_annotation && echo true || echo false)" \
    "status.spatial_deconv_gate_passed=$(eda_gate_passed spatial_deconv && echo true || echo false)" \
    "status.joint_trajectory_gate_passed=$(eda_gate_passed joint_trajectory && echo true || echo false)" \
    "status.joint_communication_gate_passed=$(eda_gate_passed joint_communication && echo true || echo false)"
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
    "scdesign3_targets_gate_passed",
    "scdesign3_validated_gate_passed",
    "deg_gate_passed",
    "communication_gate_passed",
    "regulation_gate_passed",
    "trajectory_inputs_gate_passed",
    "trajectory_methods_gate_passed",
    "trajectory_finalize_gate_passed",
    "velocity_inputs_gate_passed",
    "velocity_finalize_gate_passed",
    "05_deg_completed",
    "04d_cluster_robustness_completed",
    "04e_scdesign3_engine_completed",
    "04f_scdesign3_finalize_completed",
    "04_robustness_completed",
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
    "10a_velocity_loom_completed",
    "10b_velocity_reference_completed",
    "10_velocity_inputs_completed",
    "10c_scvelo_dynamical_completed",
    "10d_velocyto_steady_completed",
    "10e_scvelo_drivers_completed",
    "10f_cellrank_fate_completed",
    "10g_velocity_consistency_completed",
    "10h_velocity_root_terminal_completed",
    "10_velocity_stage3_completed",
    "10i_velocity_eda_completed",
    "10_velocity_completed",
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
