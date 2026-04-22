#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DEFAULT_CONFIG_FILE="${PIPELINE_ROOT}/config/server_config.sh"
CONFIG_FILE="${SCRNA_PIPELINE_CONFIG:-${DEFAULT_CONFIG_FILE}}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "缺少配置文件: ${CONFIG_FILE}" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${CONFIG_FILE}"

# Ortholog cache is a fixed project subdirectory in standard projects.
# Keep a runtime default here so older project configs generated before this
# variable existed can still run SCENIC/ortholog steps without manual edits.
export ORTHOLOG_CACHE_DIR="${ORTHOLOG_CACHE_DIR:-${RESULTS_DIR}/ortholog_cache}"
export PROJECT_CONFIG_DIR="${PROJECT_CONFIG_DIR:-${PROJECT_ROOT}/config}"
export RAW_DIR="${RAW_DIR:-${PROJECT_ROOT}/raw}"
export RAW_RECEIVED_DIR="${RAW_RECEIVED_DIR:-${RAW_DIR}/received}"
export METADATA_DIR="${METADATA_DIR:-${PROJECT_ROOT}/metadata}"
export REPORT_DIR="${REPORT_DIR:-${PROJECT_ROOT}/reports}"
export INTAKE_REPORT_DIR="${INTAKE_REPORT_DIR:-${REPORT_DIR}/intake}"
export EDA_REPORT_DIR="${EDA_REPORT_DIR:-${REPORT_DIR}/eda}"
export STATUS_DIR="${STATUS_DIR:-${PROJECT_ROOT}/status}"
export WORKFLOW_STATUS_FILE="${WORKFLOW_STATUS_FILE:-${STATUS_DIR}/workflow_status.json}"
export SAMPLE_SHEET="${SAMPLE_SHEET:-${METADATA_DIR}/samples.tsv}"
export CANONICAL_SAMPLE_SHEET="${CANONICAL_SAMPLE_SHEET:-${METADATA_DIR}/samples.canonical.tsv}"
export COMPARISON_SHEET="${COMPARISON_SHEET:-${METADATA_DIR}/comparisons.tsv}"
export DELIVERY_MANIFEST="${DELIVERY_MANIFEST:-${METADATA_DIR}/delivery_manifest.tsv}"
export RECEIVED_FILES_MANIFEST="${RECEIVED_FILES_MANIFEST:-${METADATA_DIR}/received_files_manifest.tsv}"
export INPUT_INVENTORY_FILE="${INPUT_INVENTORY_FILE:-${INTAKE_REPORT_DIR}/input_inventory.tsv}"
export BRANCH_READINESS_FILE="${BRANCH_READINESS_FILE:-${INTAKE_REPORT_DIR}/branch_readiness.tsv}"
export INTAKE_SUMMARY_FILE="${INTAKE_SUMMARY_FILE:-${INTAKE_REPORT_DIR}/intake_summary.md}"
export WAIVER_FILE="${WAIVER_FILE:-${PROJECT_CONFIG_DIR}/waivers.tsv}"
export QC_THRESHOLD_FILE="${QC_THRESHOLD_FILE:-${PROJECT_CONFIG_DIR}/qc_thresholds.tsv}"
export EDA_GATE_FILE="${EDA_GATE_FILE:-${PROJECT_CONFIG_DIR}/eda_gates.tsv}"
export OBJECT_LAYER_CONFIG_FILE="${OBJECT_LAYER_CONFIG_FILE:-${PROJECT_CONFIG_DIR}/object_layers.tsv}"
export MARKER_PANEL_DIR="${MARKER_PANEL_DIR:-${PROJECT_CONFIG_DIR}/marker_panels}"
export MITO_GENE_LIST_FILE="${MITO_GENE_LIST_FILE:-${PROJECT_CONFIG_DIR}/mito_gene_list.txt}"
export AMBIENT_REPORT_DIR="${AMBIENT_REPORT_DIR:-${EDA_REPORT_DIR}/ambient}"
export PRE_QC_REPORT_DIR="${PRE_QC_REPORT_DIR:-${EDA_REPORT_DIR}/pre_qc}"
export POST_QC_REPORT_DIR="${POST_QC_REPORT_DIR:-${EDA_REPORT_DIR}/post_qc}"
export INTEGRATION_REPORT_DIR="${INTEGRATION_REPORT_DIR:-${EDA_REPORT_DIR}/integration}"
export ANNOTATION_REPORT_DIR="${ANNOTATION_REPORT_DIR:-${EDA_REPORT_DIR}/annotation}"
export INPUT_STANDARDIZE_MODE="${INPUT_STANDARDIZE_MODE:-symlink}"
export ANNOTATION_HUB_PATH="${ANNOTATION_HUB_PATH:-${CHECKPOINT_DIR}/03_after_annotation.rds}"
export MIN_BIOLOGICAL_REPLICATES="${MIN_BIOLOGICAL_REPLICATES:-2}"
export AMBIENT_PRIMARY_METHOD="${AMBIENT_PRIMARY_METHOD:-soupx}"
export AMBIENT_FALLBACK_METHOD="${AMBIENT_FALLBACK_METHOD:-decontx}"
export AMBIENT_APPLY_POLICY="${AMBIENT_APPLY_POLICY:-manual}"
export AMBIENT_MIN_CELLS="${AMBIENT_MIN_CELLS:-50}"
export AMBIENT_CLUSTER_DIMS="${AMBIENT_CLUSTER_DIMS:-1:20}"
export AMBIENT_CLUSTER_RESOLUTION="${AMBIENT_CLUSTER_RESOLUTION:-0.4}"
export AMBIENT_MARKER_TOP_N="${AMBIENT_MARKER_TOP_N:-3}"
export AMBIENT_RECOMMEND_MIN_CONTAMINATION="${AMBIENT_RECOMMEND_MIN_CONTAMINATION:-0.05}"
export CELLBENDER_MODE="${CELLBENDER_MODE:-stub}"
export CELLBENDER_FPR="${CELLBENDER_FPR:-0.01}"
export CELLBENDER_CUDA="${CELLBENDER_CUDA:-yes}"
export CELLBENDER_EXTRA_ARGS="${CELLBENDER_EXTRA_ARGS:-}"
export DOUBLET_RATE_PER_1K="${DOUBLET_RATE_PER_1K:-${DOUBLET_RATE:-0.008}}"
export DOUBLET_PRIMARY_CALLER="${DOUBLET_PRIMARY_CALLER:-scDblFinder}"
export DOUBLET_SECONDARY_CALLER="${DOUBLET_SECONDARY_CALLER:-DoubletFinder}"
export DOUBLET_SECONDARY_ENABLED="${DOUBLET_SECONDARY_ENABLED:-yes}"
export DOUBLET_MIN_CELLS="${DOUBLET_MIN_CELLS:-50}"

export PATH="${SRA_TOOLS_DIR}:${PATH}"

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

warn() {
  echo "[WARN] $*" >&2
}

ensure_dir() {
  mkdir -p "$@"
}

ensure_eda_control_files() {
  ensure_dir \
    "${PROJECT_CONFIG_DIR}" \
    "${MARKER_PANEL_DIR}" \
    "${AMBIENT_REPORT_DIR}" \
    "${EDA_REPORT_DIR}" \
    "${PRE_QC_REPORT_DIR}" \
    "${POST_QC_REPORT_DIR}" \
    "${INTEGRATION_REPORT_DIR}" \
    "${ANNOTATION_REPORT_DIR}"

  if [[ ! -s "${EDA_GATE_FILE}" ]]; then
    cat > "${EDA_GATE_FILE}" <<'EOF'
gate_id	status	approved_by	notes
pre_qc	pending		
post_qc	pending		
integration	pending		
annotation	pending		
EOF
  fi

  local gate_id
  for gate_id in pre_qc post_qc integration annotation; do
    if ! awk -F '\t' -v gate="${gate_id}" 'NR > 1 && $1 == gate { found = 1 } END { exit(found ? 0 : 1) }' "${EDA_GATE_FILE}" >/dev/null 2>&1; then
      printf '%s\tpending\t\t\n' "${gate_id}" >> "${EDA_GATE_FILE}"
    fi
  done

  if [[ ! -s "${QC_THRESHOLD_FILE}" ]]; then
    cat > "${QC_THRESHOLD_FILE}" <<'EOF'
sample_id	qc_min_nfeature	qc_min_ncount	qc_min_log10umi	qc_max_mito_pct
EOF
  fi
}

ensure_object_layer_config_file() {
  ensure_dir "${PROJECT_CONFIG_DIR}"
  [[ -s "${OBJECT_LAYER_CONFIG_FILE}" ]] && return 0

  {
    printf '# object_layers.tsv\n'
    printf '# - panorama 行是根层，默认从全部 post-QC 细胞建模\n'
    printf '# - subcluster 行从 parent_layer 的 metadata 中按 selection_column/selection_values 取细胞\n'
    printf '# - 默认示例保留 subcluster_1 / subcluster_2，但 selection_values 需要在审阅 panorama 后填写\n'
    printf 'layer_id\tlayer_role\tenabled\tparent_layer\tsample_include\tsample_exclude\tselection_column\tselection_values\trebuild_normalization\thvg_nfeatures\tpca_dims\ttarget_clusters\tres_range\tres_fine_step\tintegration_mode\tdescription\n'
    printf 'panorama\tpanorama\tyes\t\t\t\t\t\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${HVG_NFEATURES:-2000}" \
      "${PCA_DIMS:-1:30}" \
      "${TARGET_CLUSTERS:-15}" \
      "${RES_RANGE:-0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60}" \
      "${RES_FINE_STEP:-0.005}" \
      "${INTEGRATION_MODE:-harmony}" \
      'Root panorama object built from all post-QC cells.'
    printf 'subcluster_1\tsubcluster\tyes\tpanorama\t\t\tpanorama_cluster\t\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${HVG_NFEATURES:-2000}" \
      "1:20" \
      "8" \
      "0.10,0.15,0.20,0.25,0.30,0.35,0.40" \
      "${RES_FINE_STEP:-0.005}" \
      "${INTEGRATION_MODE:-harmony}" \
      'Edit selection_column and selection_values to define this subcluster from panorama results.'
    printf 'subcluster_2\tsubcluster\tyes\tpanorama\t\t\tpanorama_cluster\t\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${HVG_NFEATURES:-2000}" \
      "1:20" \
      "6" \
      "0.10,0.15,0.20,0.25,0.30,0.35,0.40" \
      "${RES_FINE_STEP:-0.005}" \
      "${INTEGRATION_MODE:-harmony}" \
      'Edit selection_column and selection_values to define this subcluster from panorama results.'
  } > "${OBJECT_LAYER_CONFIG_FILE}"
}

ensure_marker_panel_dir() {
  ensure_dir "${MARKER_PANEL_DIR}"

  if [[ -f "${PIPELINE_ROOT}/config/marker_panels/README.md" && ! -f "${MARKER_PANEL_DIR}/README.md" ]]; then
    cp -f "${PIPELINE_ROOT}/config/marker_panels/README.md" "${MARKER_PANEL_DIR}/README.md"
  fi
  if [[ -f "${PIPELINE_ROOT}/config/marker_panels/example_panel.tsv.template" && ! -f "${MARKER_PANEL_DIR}/example_panel.tsv.template" ]]; then
    cp -f "${PIPELINE_ROOT}/config/marker_panels/example_panel.tsv.template" "${MARKER_PANEL_DIR}/example_panel.tsv.template"
  fi
}

ensure_mito_gene_list_file() {
  ensure_dir "${PROJECT_CONFIG_DIR}"

  if [[ -s "${MITO_GENE_LIST_FILE}" ]]; then
    return 0
  fi

  if [[ -f "${PIPELINE_ROOT}/config/mito_gene_list.txt.template" ]]; then
    cp -f "${PIPELINE_ROOT}/config/mito_gene_list.txt.template" "${MITO_GENE_LIST_FILE}"
    return 0
  fi

  cat > "${MITO_GENE_LIST_FILE}" <<'EOF'
# One mito gene identifier per line.
# Use gene symbols or feature IDs that exactly match your matrix rownames.
# Leave this file unchanged when GTF-based mito detection is sufficient.
EOF
}

now_iso() {
  date --iso-8601=seconds
}

detect_python() {
  if command -v python3 >/dev/null 2>&1; then
    echo "python3"
    return
  fi

  local candidate=""
  for candidate in \
    "${R_MAIN_ENV_PREFIX:-}/bin/python" \
    "${R_SCENIC_ENV_PREFIX:-}/bin/python" \
    "${VELOCITY_ENV_PREFIX:-}/bin/python"; do
    if [[ -n "${candidate}" && -x "${candidate}" ]]; then
      echo "${candidate}"
      return
    fi
  done

  if command -v python >/dev/null 2>&1 && python - <<'PY' >/dev/null 2>&1
import sys
sys.exit(0 if sys.version_info[0] >= 3 else 1)
PY
  then
    echo "python"
    return
  fi

  die "系统中缺少可用的 Python 3 解释器，无法写入 workflow 状态文件。"
}

bool_string() {
  case "${1:-false}" in
    true|TRUE|1|yes|YES|on|ON)
      echo "true"
      ;;
    *)
      echo "false"
      ;;
  esac
}

tsv_get_col_index() {
  local file_path="$1"
  local column_name="$2"
  awk -F '\t' -v target="${column_name}" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == target) {
          print i
          found = 1
          exit 0
        }
      }
      exit(found ? 0 : 1)
    }
  ' "${file_path}"
}

tsv_require_columns() {
  local file_path="$1"
  shift
  local column_name
  [[ -f "${file_path}" ]] || die "缺少 TSV 文件: ${file_path}"
  for column_name in "$@"; do
    if ! tsv_get_col_index "${file_path}" "${column_name}" >/dev/null; then
      die "TSV 文件 ${file_path} 缺少必需列: ${column_name}"
    fi
  done
}

detect_project_input_mode() {
  local file_path="${CANONICAL_SAMPLE_SHEET}"
  if [[ ! -f "${file_path}" ]]; then
    file_path="${SAMPLE_SHEET}"
  fi
  if [[ ! -f "${file_path}" ]]; then
    echo ""
    return
  fi

  local col_idx
  if ! col_idx="$(tsv_get_col_index "${file_path}" "input_mode" 2>/dev/null)"; then
    echo ""
    return
  fi

  awk -F '\t' -v col="${col_idx}" '
    NR > 1 && $col != "" {
      seen[$col] = 1
    }
    END {
      count = 0
      for (value in seen) {
        values[++count] = value
      }
      if (count == 1) {
        print values[1]
      } else if (count > 1) {
        print "mixed"
      }
    }
  ' "${file_path}"
}

has_fastq_files() {
  local dir_path="$1"
  [[ -d "${dir_path}" ]] || return 1
  find "${dir_path}" -maxdepth 1 -type f \
    \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) \
    -print -quit | grep -q .
}

is_10x_matrix_dir() {
  local dir_path="$1"
  [[ -d "${dir_path}" ]] || return 1
  [[ -f "${dir_path}/matrix.mtx.gz" || -f "${dir_path}/matrix.mtx" ]] || return 1
  [[ -f "${dir_path}/features.tsv.gz" || -f "${dir_path}/features.tsv" || -f "${dir_path}/genes.tsv.gz" || -f "${dir_path}/genes.tsv" ]] || return 1
  [[ -f "${dir_path}/barcodes.tsv.gz" || -f "${dir_path}/barcodes.tsv" ]] || return 1
}

resolve_cellranger_sample_dir() {
  local sample_id="$1"
  local input_mode="$2"
  local source_path="$3"

  case "${input_mode}" in
    fastq)
      if [[ -d "${CELLRANGER_OUT_DIR}/${sample_id}/outs" ]]; then
        echo "${CELLRANGER_OUT_DIR}/${sample_id}"
      fi
      ;;
    cellranger_out)
      if [[ -d "${source_path}/outs" ]]; then
        echo "${source_path}"
      elif [[ -d "${source_path}/${sample_id}/outs" ]]; then
        echo "${source_path}/${sample_id}"
      elif [[ -d "${CELLRANGER_OUT_DIR}/${sample_id}/outs" ]]; then
        echo "${CELLRANGER_OUT_DIR}/${sample_id}"
      fi
      ;;
    *)
      ;;
  esac
}

resolve_matrix_source_dir() {
  local sample_id="$1"
  local input_mode="$2"
  local source_path="$3"
  local sample_dir=""

  case "${input_mode}" in
    fastq|cellranger_out)
      sample_dir="$(resolve_cellranger_sample_dir "${sample_id}" "${input_mode}" "${source_path}")"
      if [[ -n "${sample_dir}" ]]; then
        echo "${sample_dir}/outs/filtered_feature_bc_matrix"
      fi
      ;;
    matrix)
      if is_10x_matrix_dir "${source_path}"; then
        echo "${source_path}"
      elif is_10x_matrix_dir "${source_path}/filtered_feature_bc_matrix"; then
        echo "${source_path}/filtered_feature_bc_matrix"
      elif is_10x_matrix_dir "${source_path}/${sample_id}"; then
        echo "${source_path}/${sample_id}"
      elif is_10x_matrix_dir "${DATA_DIR}/${sample_id}"; then
        echo "${DATA_DIR}/${sample_id}"
      fi
      ;;
    *)
      ;;
  esac
}

prepare_project_state_dirs() {
  ensure_dir \
    "${PROJECT_CONFIG_DIR}" \
    "${RAW_RECEIVED_DIR}" \
    "${METADATA_DIR}" \
    "${REPORT_DIR}" \
    "${INTAKE_REPORT_DIR}" \
    "${EDA_REPORT_DIR}" \
    "${PRE_QC_REPORT_DIR}" \
    "${POST_QC_REPORT_DIR}" \
    "${INTEGRATION_REPORT_DIR}" \
    "${ANNOTATION_REPORT_DIR}" \
    "${STATUS_DIR}" \
    "${LOG_DIR}"
  ensure_eda_control_files
  ensure_object_layer_config_file
  ensure_marker_panel_dir
  ensure_mito_gene_list_file
}

eda_gate_status() {
  local gate_id="$1"
  ensure_eda_control_files
  awk -F '\t' -v gate="${gate_id}" '
    NR > 1 && $1 == gate {
      print tolower($2)
      found = 1
      exit
    }
    END {
      if (!found) {
        print "pending"
      }
    }
  ' "${EDA_GATE_FILE}"
}

eda_gate_passed() {
  local gate_id="$1"
  local gate_status
  gate_status="$(eda_gate_status "${gate_id}")"
  [[ "${gate_status}" == "approved" || "${gate_status}" == "yes" || "${gate_status}" == "true" ]]
}

set_eda_gate_status() {
  local gate_id="$1"
  local gate_status="$2"
  local approved_by="${3:-}"
  local notes="${4:-}"

  ensure_eda_control_files

  local python_bin
  python_bin="$(detect_python)"

  env \
    EDA_GATE_FILE="${EDA_GATE_FILE}" \
    TARGET_GATE_ID="${gate_id}" \
    TARGET_GATE_STATUS="${gate_status}" \
    TARGET_APPROVED_BY="${approved_by}" \
    TARGET_NOTES="${notes}" \
    "${python_bin}" - <<'PY'
import csv
import os
from pathlib import Path

gate_file = Path(os.environ["EDA_GATE_FILE"])
target_gate_id = os.environ["TARGET_GATE_ID"]
target_gate_status = os.environ["TARGET_GATE_STATUS"]
target_approved_by = os.environ["TARGET_APPROVED_BY"]
target_notes = os.environ["TARGET_NOTES"]

rows = []
with gate_file.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = reader.fieldnames or ["gate_id", "status", "approved_by", "notes"]
    for row in reader:
        if row["gate_id"] == target_gate_id:
            row["status"] = target_gate_status
            row["approved_by"] = target_approved_by
            row["notes"] = target_notes
        rows.append(row)

if not any(row["gate_id"] == target_gate_id for row in rows):
    rows.append({
        "gate_id": target_gate_id,
        "status": target_gate_status,
        "approved_by": target_approved_by,
        "notes": target_notes,
    })

with gate_file.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
PY
}

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
    "status.annotation_gate_passed=$(eda_gate_passed annotation && echo true || echo false)"
}

update_workflow_status() {
  local stage_name="$1"
  local next_step="$2"
  shift 2 || true

  prepare_project_state_dirs

  local python_bin
  python_bin="$(detect_python)"

  env \
    WORKFLOW_STATUS_FILE="${WORKFLOW_STATUS_FILE}" \
    WF_STAGE="${stage_name}" \
    WF_NEXT_STEP="${next_step}" \
    WF_PROJECT_ROOT="${PROJECT_ROOT}" \
    WF_CONFIG_FILE="${CONFIG_FILE}" \
    WF_PROJECT_INPUT_MODE="${PROJECT_INPUT_MODE:-$(detect_project_input_mode)}" \
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
data["project_input_mode"] = os.environ.get("WF_PROJECT_INPUT_MODE", "")

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
    "main_upstream_ready",
    "velocity_upstream_ready",
    "ambient_upstream_ready",
    "scenic_upstream_ready",
    "de_replicate_ready",
):
    status.setdefault(key, False)

status["main_ready"] = bool(status.get("metadata_valid")) and bool(status.get("standardized_inputs")) and bool(status.get("main_upstream_ready"))
status["deg_ready"] = annotation_exists and bool(status.get("annotation_gate_passed"))
status["trajectory_ready"] = annotation_exists and bool(status.get("annotation_gate_passed"))
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

detect_container_runtime() {
  if command -v apptainer >/dev/null 2>&1; then
    echo "apptainer"
    return
  fi
  if command -v singularity >/dev/null 2>&1; then
    echo "singularity"
    return
  fi
  echo ""
}

detect_conda_frontend() {
  if command -v micromamba >/dev/null 2>&1; then
    echo "micromamba"
    return
  fi
  if command -v mamba >/dev/null 2>&1; then
    echo "mamba"
    return
  fi
  if command -v conda >/dev/null 2>&1; then
    echo "conda"
    return
  fi
  echo ""
}

CONTAINER_RUNTIME="${CONTAINER_RUNTIME:-$(detect_container_runtime)}"
CONDA_FRONTEND="${CONDA_FRONTEND:-$(detect_conda_frontend)}"

declare -ag CONTAINER_BIND_ARGS=()
declare -ag BIND_PATHS=(
  "${PROJECT_ROOT}"
  "${PAPER_MATRIX_SOURCE:-}"
  "${CELLRANGER_OUT_DIR:-}"
  "${REFERENCE_DIR:-}"
  "${ENV_DIR:-}"
)

if [[ -n "${EXISTING_R_SIF:-}" ]]; then
  BIND_PATHS+=("$(dirname "${EXISTING_R_SIF}")")
fi

for bind_path in "${BIND_PATHS[@]}"; do
  [[ -n "${bind_path}" ]] || continue
  if [[ -e "${bind_path}" ]]; then
    CONTAINER_BIND_ARGS+=(--bind "${bind_path}:${bind_path}")
  fi
done

run_r_main() {
  local script_path="$1"
  shift || true

  if [[ "${USE_EXISTING_SIF}" == "yes" ]]; then
    [[ -n "${CONTAINER_RUNTIME}" ]] || die "USE_EXISTING_SIF=yes，但系统中没有 apptainer/singularity。"
    [[ -f "${EXISTING_R_SIF}" ]] || die "找不到旧 R 容器: ${EXISTING_R_SIF}"

    "${CONTAINER_RUNTIME}" exec "${CONTAINER_BIND_ARGS[@]}" "${EXISTING_R_SIF}" \
      env \
      PIPELINE_ROOT="${PIPELINE_ROOT}" \
      R_LIBS_USER="${R_LIBS_MAIN}" \
      PROJECT_ROOT="${PROJECT_ROOT}" \
      DATA_DIR="${DATA_DIR}" \
      RESULTS_DIR="${RESULTS_DIR}" \
      CHECKPOINT_DIR="${CHECKPOINT_DIR}" \
      FIGURE_DIR="${FIGURE_DIR}" \
      TABLE_DIR="${TABLE_DIR}" \
      LOG_DIR="${LOG_DIR}" \
      CELLRANGER_OUT_DIR="${CELLRANGER_OUT_DIR}" \
      VELOCITY_INPUT_DIR="${VELOCITY_INPUT_DIR}" \
      VELOCITY_LOOM_DIR="${VELOCITY_LOOM_DIR}" \
      VELOCITY_OUTPUT_DIR="${VELOCITY_OUTPUT_DIR}" \
      SCENIC_INPUT_DIR="${SCENIC_INPUT_DIR}" \
      SCENIC_OUTPUT_DIR="${SCENIC_OUTPUT_DIR}" \
      SCENIC_ORTHOLOG_MAP_FILE="${SCENIC_ORTHOLOG_MAP_FILE:-}" \
      RESOURCE_DIR="${RESOURCE_DIR}" \
      RAW_GROUP_1_NAME="${RAW_GROUP_1_NAME}" \
      RAW_GROUP_1_SAMPLES="${RAW_GROUP_1_SAMPLES}" \
      RAW_GROUP_2_NAME="${RAW_GROUP_2_NAME}" \
      RAW_GROUP_2_SAMPLES="${RAW_GROUP_2_SAMPLES}" \
      RAW_SAMPLES="${RAW_SAMPLES}" \
      DEG_IDENT_1="${DEG_IDENT_1}" \
      DEG_IDENT_2="${DEG_IDENT_2}" \
      RANDOM_SEED="${RANDOM_SEED}" \
      SAMPLE_NAMES="${SAMPLE_NAMES}" \
      OBJECT_LAYER_CONFIG_FILE="${OBJECT_LAYER_CONFIG_FILE}" \
      QC_MIN_NFEATURE="${QC_MIN_NFEATURE}" \
      QC_MIN_NCOUNT="${QC_MIN_NCOUNT}" \
      QC_MIN_LOG10UMI="${QC_MIN_LOG10UMI}" \
      QC_MAX_MITO_PCT="${QC_MAX_MITO_PCT}" \
      AMBIENT_PRIMARY_METHOD="${AMBIENT_PRIMARY_METHOD}" \
      AMBIENT_FALLBACK_METHOD="${AMBIENT_FALLBACK_METHOD}" \
      AMBIENT_APPLY_POLICY="${AMBIENT_APPLY_POLICY}" \
      AMBIENT_MIN_CELLS="${AMBIENT_MIN_CELLS}" \
      AMBIENT_CLUSTER_DIMS="${AMBIENT_CLUSTER_DIMS}" \
      AMBIENT_CLUSTER_RESOLUTION="${AMBIENT_CLUSTER_RESOLUTION}" \
      AMBIENT_MARKER_TOP_N="${AMBIENT_MARKER_TOP_N}" \
      AMBIENT_RECOMMEND_MIN_CONTAMINATION="${AMBIENT_RECOMMEND_MIN_CONTAMINATION}" \
      CELLBENDER_MODE="${CELLBENDER_MODE}" \
      CELLBENDER_FPR="${CELLBENDER_FPR}" \
      CELLBENDER_CUDA="${CELLBENDER_CUDA}" \
      CELLBENDER_EXTRA_ARGS="${CELLBENDER_EXTRA_ARGS}" \
      DOUBLET_RATE="${DOUBLET_RATE}" \
      DOUBLET_RATE_PER_1K="${DOUBLET_RATE_PER_1K}" \
      DOUBLET_PRIMARY_CALLER="${DOUBLET_PRIMARY_CALLER}" \
      DOUBLET_SECONDARY_CALLER="${DOUBLET_SECONDARY_CALLER}" \
      DOUBLET_SECONDARY_ENABLED="${DOUBLET_SECONDARY_ENABLED}" \
      DOUBLET_MIN_CELLS="${DOUBLET_MIN_CELLS}" \
      DOUBLET_DIMS="${DOUBLET_DIMS}" \
      HVG_NFEATURES="${HVG_NFEATURES}" \
      PCA_DIMS="${PCA_DIMS}" \
      TARGET_CLUSTERS="${TARGET_CLUSTERS}" \
      RES_RANGE="${RES_RANGE}" \
      TRAJECTORY_START="${TRAJECTORY_START}" \
      TRADESEQ_KNOTS="${TRADESEQ_KNOTS}" \
      ENSEMBL_MIRROR="${ENSEMBL_MIRROR}" \
      MAIN_THREADS="${MAIN_THREADS}" \
      Rscript "${script_path}" "$@"
    return
  fi

  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda，无法运行主流程 R 环境。"
  [[ -d "${R_MAIN_ENV_PREFIX}" ]] || die "主流程 R 环境不存在: ${R_MAIN_ENV_PREFIX}"

  env \
    PIPELINE_ROOT="${PIPELINE_ROOT}" \
    R_LIBS_USER="${R_LIBS_MAIN}" \
    PROJECT_ROOT="${PROJECT_ROOT}" \
    DATA_DIR="${DATA_DIR}" \
    RESULTS_DIR="${RESULTS_DIR}" \
    CHECKPOINT_DIR="${CHECKPOINT_DIR}" \
    FIGURE_DIR="${FIGURE_DIR}" \
    TABLE_DIR="${TABLE_DIR}" \
    LOG_DIR="${LOG_DIR}" \
    CELLRANGER_OUT_DIR="${CELLRANGER_OUT_DIR}" \
    VELOCITY_INPUT_DIR="${VELOCITY_INPUT_DIR}" \
    VELOCITY_LOOM_DIR="${VELOCITY_LOOM_DIR}" \
    VELOCITY_OUTPUT_DIR="${VELOCITY_OUTPUT_DIR}" \
    SCENIC_INPUT_DIR="${SCENIC_INPUT_DIR}" \
    SCENIC_OUTPUT_DIR="${SCENIC_OUTPUT_DIR}" \
    SCENIC_ORTHOLOG_MAP_FILE="${SCENIC_ORTHOLOG_MAP_FILE:-}" \
    RESOURCE_DIR="${RESOURCE_DIR}" \
    RAW_GROUP_1_NAME="${RAW_GROUP_1_NAME}" \
    RAW_GROUP_1_SAMPLES="${RAW_GROUP_1_SAMPLES}" \
    RAW_GROUP_2_NAME="${RAW_GROUP_2_NAME}" \
    RAW_GROUP_2_SAMPLES="${RAW_GROUP_2_SAMPLES}" \
    RAW_SAMPLES="${RAW_SAMPLES}" \
    DEG_IDENT_1="${DEG_IDENT_1}" \
    DEG_IDENT_2="${DEG_IDENT_2}" \
    RANDOM_SEED="${RANDOM_SEED}" \
    SAMPLE_NAMES="${SAMPLE_NAMES}" \
    OBJECT_LAYER_CONFIG_FILE="${OBJECT_LAYER_CONFIG_FILE}" \
    QC_MIN_NFEATURE="${QC_MIN_NFEATURE}" \
    QC_MIN_NCOUNT="${QC_MIN_NCOUNT}" \
    QC_MIN_LOG10UMI="${QC_MIN_LOG10UMI}" \
    QC_MAX_MITO_PCT="${QC_MAX_MITO_PCT}" \
    AMBIENT_PRIMARY_METHOD="${AMBIENT_PRIMARY_METHOD}" \
    AMBIENT_FALLBACK_METHOD="${AMBIENT_FALLBACK_METHOD}" \
    AMBIENT_APPLY_POLICY="${AMBIENT_APPLY_POLICY}" \
    AMBIENT_MIN_CELLS="${AMBIENT_MIN_CELLS}" \
    AMBIENT_CLUSTER_DIMS="${AMBIENT_CLUSTER_DIMS}" \
    AMBIENT_CLUSTER_RESOLUTION="${AMBIENT_CLUSTER_RESOLUTION}" \
    AMBIENT_MARKER_TOP_N="${AMBIENT_MARKER_TOP_N}" \
    AMBIENT_RECOMMEND_MIN_CONTAMINATION="${AMBIENT_RECOMMEND_MIN_CONTAMINATION}" \
    CELLBENDER_MODE="${CELLBENDER_MODE}" \
    CELLBENDER_FPR="${CELLBENDER_FPR}" \
    CELLBENDER_CUDA="${CELLBENDER_CUDA}" \
    CELLBENDER_EXTRA_ARGS="${CELLBENDER_EXTRA_ARGS}" \
    DOUBLET_RATE="${DOUBLET_RATE}" \
    DOUBLET_RATE_PER_1K="${DOUBLET_RATE_PER_1K}" \
    DOUBLET_PRIMARY_CALLER="${DOUBLET_PRIMARY_CALLER}" \
    DOUBLET_SECONDARY_CALLER="${DOUBLET_SECONDARY_CALLER}" \
    DOUBLET_SECONDARY_ENABLED="${DOUBLET_SECONDARY_ENABLED}" \
    DOUBLET_MIN_CELLS="${DOUBLET_MIN_CELLS}" \
    DOUBLET_DIMS="${DOUBLET_DIMS}" \
    HVG_NFEATURES="${HVG_NFEATURES}" \
    PCA_DIMS="${PCA_DIMS}" \
    TARGET_CLUSTERS="${TARGET_CLUSTERS}" \
    RES_RANGE="${RES_RANGE}" \
    TRAJECTORY_START="${TRAJECTORY_START}" \
    TRADESEQ_KNOTS="${TRADESEQ_KNOTS}" \
    ENSEMBL_MIRROR="${ENSEMBL_MIRROR}" \
    MAIN_THREADS="${MAIN_THREADS}" \
    "${CONDA_FRONTEND}" run -p "${R_MAIN_ENV_PREFIX}" Rscript "${script_path}" "$@"
}

run_r_scenic() {
  local script_path="$1"
  shift || true

  if [[ "${USE_EXISTING_SIF}" == "yes" ]]; then
    [[ -n "${CONTAINER_RUNTIME}" ]] || die "USE_EXISTING_SIF=yes，但系统中没有 apptainer/singularity。"
    [[ -f "${EXISTING_R_SIF}" ]] || die "找不到旧 R 容器: ${EXISTING_R_SIF}"

    "${CONTAINER_RUNTIME}" exec "${CONTAINER_BIND_ARGS[@]}" "${EXISTING_R_SIF}" \
      env \
      PIPELINE_ROOT="${PIPELINE_ROOT}" \
      R_LIBS_USER="${R_LIBS_SCENIC}" \
      PROJECT_ROOT="${PROJECT_ROOT}" \
      DATA_DIR="${DATA_DIR}" \
      RESULTS_DIR="${RESULTS_DIR}" \
      CHECKPOINT_DIR="${CHECKPOINT_DIR}" \
      FIGURE_DIR="${FIGURE_DIR}" \
      TABLE_DIR="${TABLE_DIR}" \
      LOG_DIR="${LOG_DIR}" \
      SCENIC_INPUT_DIR="${SCENIC_INPUT_DIR}" \
      SCENIC_OUTPUT_DIR="${SCENIC_OUTPUT_DIR}" \
      SCENIC_TF_LIST="${SCENIC_TF_LIST}" \
      SCENIC_MOTIF_ANN="${SCENIC_MOTIF_ANN}" \
      SCENIC_DB_500BP="${SCENIC_DB_500BP}" \
      SCENIC_DB_10KB="${SCENIC_DB_10KB}" \
      RANDOM_SEED="${RANDOM_SEED}" \
      Rscript "${script_path}" "$@"
    return
  fi

  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda，无法运行 SCENIC R 环境。"
  [[ -d "${R_SCENIC_ENV_PREFIX}" ]] || die "SCENIC R 环境不存在: ${R_SCENIC_ENV_PREFIX}"

  env \
    PIPELINE_ROOT="${PIPELINE_ROOT}" \
    R_LIBS_USER="${R_LIBS_SCENIC}" \
    PROJECT_ROOT="${PROJECT_ROOT}" \
    DATA_DIR="${DATA_DIR}" \
    RESULTS_DIR="${RESULTS_DIR}" \
    CHECKPOINT_DIR="${CHECKPOINT_DIR}" \
    FIGURE_DIR="${FIGURE_DIR}" \
    TABLE_DIR="${TABLE_DIR}" \
    LOG_DIR="${LOG_DIR}" \
    SCENIC_INPUT_DIR="${SCENIC_INPUT_DIR}" \
    SCENIC_OUTPUT_DIR="${SCENIC_OUTPUT_DIR}" \
    SCENIC_TF_LIST="${SCENIC_TF_LIST}" \
    SCENIC_MOTIF_ANN="${SCENIC_MOTIF_ANN}" \
    SCENIC_DB_500BP="${SCENIC_DB_500BP}" \
    SCENIC_DB_10KB="${SCENIC_DB_10KB}" \
    RANDOM_SEED="${RANDOM_SEED}" \
    "${CONDA_FRONTEND}" run -p "${R_SCENIC_ENV_PREFIX}" Rscript "${script_path}" "$@"
}

run_in_conda_prefix() {
  local prefix="$1"
  shift
  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda。"
  "${CONDA_FRONTEND}" run -p "${prefix}" "$@"
}

create_or_update_conda_env() {
  local prefix="$1"
  local yaml_path="$2"
  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda。"

  if [[ -d "${prefix}" ]]; then
    "${CONDA_FRONTEND}" env update -y -p "${prefix}" -f "${yaml_path}" --prune
  else
    "${CONDA_FRONTEND}" env create -y -p "${prefix}" -f "${yaml_path}"
  fi
}

run_pyscenic() {
  [[ -d "${PYSCENIC_ENV_PREFIX}" ]] || die "pySCENIC 环境不存在: ${PYSCENIC_ENV_PREFIX}"
  run_in_conda_prefix "${PYSCENIC_ENV_PREFIX}" "$@"
}

run_r_legacy() {
  [[ -d "${R_LEGACY_ENV_PREFIX}" ]] || die "legacy R 环境不存在: ${R_LEGACY_ENV_PREFIX}"
  run_in_conda_prefix "${R_LEGACY_ENV_PREFIX}" Rscript "$@"
}

run_r_interaction() {
  [[ -d "${R_INTERACTION_ENV_PREFIX}" ]] || die "通讯分析 R 环境不存在: ${R_INTERACTION_ENV_PREFIX}"
  run_in_conda_prefix "${R_INTERACTION_ENV_PREFIX}" Rscript "$@"
}

run_r_spatial() {
  [[ -d "${R_SPATIAL_ENV_PREFIX}" ]] || die "空间分析 R 环境不存在: ${R_SPATIAL_ENV_PREFIX}"
  run_in_conda_prefix "${R_SPATIAL_ENV_PREFIX}" Rscript "$@"
}

run_py_spatial() {
  [[ -d "${PY_SPATIAL_ENV_PREFIX}" ]] || die "空间分析 Python 环境不存在: ${PY_SPATIAL_ENV_PREFIX}"
  run_in_conda_prefix "${PY_SPATIAL_ENV_PREFIX}" "$@"
}

run_py_spatial_legacy() {
  [[ -d "${PY_SPATIAL_LEGACY_ENV_PREFIX}" ]] || die "空间 legacy Python 环境不存在: ${PY_SPATIAL_LEGACY_ENV_PREFIX}"
  run_in_conda_prefix "${PY_SPATIAL_LEGACY_ENV_PREFIX}" "$@"
}

run_py_cell2location() {
  [[ -d "${PY_CELL2LOCATION_ENV_PREFIX}" ]] || die "cell2location Python 环境不存在: ${PY_CELL2LOCATION_ENV_PREFIX}"
  run_in_conda_prefix "${PY_CELL2LOCATION_ENV_PREFIX}" "$@"
}

run_r_validation() {
  [[ -d "${R_VALIDATION_ENV_PREFIX}" ]] || die "验证分析 R 环境不存在: ${R_VALIDATION_ENV_PREFIX}"
  run_in_conda_prefix "${R_VALIDATION_ENV_PREFIX}" Rscript "$@"
}
