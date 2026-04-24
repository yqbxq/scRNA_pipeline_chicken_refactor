#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")/../02lib" && pwd)/common.sh"

prepare_project_state_dirs

[[ -f "${BRANCH_READINESS_FILE}" ]] || die "缺少 readiness 文件，请先运行 03_audit_inputs.sh"

INPUT_SHEET="${CANONICAL_SAMPLE_SHEET}"
if [[ ! -f "${INPUT_SHEET}" ]]; then
  INPUT_SHEET="${SAMPLE_SHEET}"
fi
[[ -f "${INPUT_SHEET}" ]] || die "缺少样本表: ${INPUT_SHEET}"

python_bin="$(detect_python)"
SUMMARY_ENV_FILE="${STATUS_DIR}/input_summary.env"

env \
  INPUT_SHEET="${INPUT_SHEET}" \
  DATA_DIR="${DATA_DIR}" \
  INPUT_INVENTORY_FILE="${INPUT_INVENTORY_FILE}" \
  BRANCH_READINESS_FILE="${BRANCH_READINESS_FILE}" \
  INTAKE_SUMMARY_FILE="${INTAKE_SUMMARY_FILE}" \
  SUMMARY_ENV_FILE="${SUMMARY_ENV_FILE}" \
  PIPELINE_ROOT="${PIPELINE_ROOT}" \
  "${python_bin}" - <<'PY'
import csv
import os
import shlex
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ["PIPELINE_ROOT"]) / "shell" / "04python"))

from intake_contract import is_mex_matrix_dir

input_sheet = Path(os.environ["INPUT_SHEET"])
data_dir = Path(os.environ["DATA_DIR"])
inventory_path = Path(os.environ["INPUT_INVENTORY_FILE"])
readiness_path = Path(os.environ["BRANCH_READINESS_FILE"])
summary_path = Path(os.environ["INTAKE_SUMMARY_FILE"])
env_path = Path(os.environ["SUMMARY_ENV_FILE"])
pipeline_root = Path(os.environ["PIPELINE_ROOT"])

NOTE_LABELS = {
    "main_inputs_incomplete": "主流程输入不完整",
    "velocity_inputs_incomplete": "velocity 上游输入不完整",
    "ambient_inputs_incomplete": "ambient 分支上游输入不完整",
    "scenic_inputs_incomplete": "SCENIC 上游输入不完整",
    "replicates_insufficient": "replicate 信息不足",
    "reference_version_mixed": "reference_version 混用",
    "gene_id_type_unrecognized": "gene_id_type 未识别",
    "missing_filtered_matrix": "filtered matrix 缺失",
    "missing_raw_matrix": "raw matrix 缺失",
    "ambient_branch_unavailable": "ambient 分支不可用",
    "missing_bam": "BAM 缺失",
    "missing_fastq": "FASTQ 缺失",
    "missing_cellranger_out": "cellranger_out 缺失",
    "replicate_fields_missing": "replicate/batch 字段缺失",
    "velocity_requires_bam": "velocity 上游需要 BAM 文件",
}


def describe_notes(raw: str) -> list[str]:
    notes = []
    for note in (raw or "").split(";"):
        note = note.strip()
        if not note:
            continue
        notes.append(NOTE_LABELS.get(note, note))
    return notes

with input_sheet.open("r", encoding="utf-8", newline="") as handle:
    sample_rows = list(csv.DictReader(handle, delimiter="\t"))
with inventory_path.open("r", encoding="utf-8", newline="") as handle:
    inventory_rows = list(csv.DictReader(handle, delimiter="\t"))
with readiness_path.open("r", encoding="utf-8", newline="") as handle:
    readiness_rows = list(csv.DictReader(handle, delimiter="\t"))

project_row = next((row for row in readiness_rows if row["sample_id"] == "__PROJECT__"), None)
if project_row is None:
    raise SystemExit("branch_readiness.tsv 缺少 __PROJECT__ 汇总行。")

run_main_samples = [row for row in sample_rows if row.get("run_main", "yes").lower() == "yes"]
project_input_modes = {row.get("input_mode", "").strip() for row in run_main_samples if row.get("input_mode", "").strip()}
standardized_ok = all(is_mex_matrix_dir(data_dir / row["sample_id"]) for row in run_main_samples)

inventory_by_sample = {row["sample_id"]: row for row in inventory_rows}
readiness_by_sample = {row["sample_id"]: row for row in readiness_rows}
blocked_samples = []
for row in run_main_samples:
    sample_id = row["sample_id"]
    ready = inventory_by_sample.get(sample_id, {}).get("has_filtered_matrix", "false") == "true"
    standardized = is_mex_matrix_dir(data_dir / sample_id)
    if not ready or not standardized:
        blocked_samples.append(sample_id)

main_upstream_ready = project_row["main_upstream_ready"]
velocity_upstream_ready = project_row["velocity_upstream_ready"]
ambient_upstream_ready = project_row["ambient_upstream_ready"]
scenic_upstream_ready = project_row["scenic_upstream_ready"]
de_replicate_ready = project_row["de_replicate_ready"]

if main_upstream_ready == "true" and standardized_ok:
    next_step = f"bash {pipeline_root}/shell/03stages/install_envs.sh"
elif main_upstream_ready != "true":
    if project_input_modes == {"fastq"}:
        next_step = f"bash {pipeline_root}/shell/03stages/alignment.sh"
    else:
        next_step = f"bash {pipeline_root}/shell/03stages/audit_inputs.sh"
else:
    next_step = f"bash {pipeline_root}/shell/03stages/standardize_inputs.sh"

lines = []
lines.append("# Intake Summary")
lines.append("")
lines.append(f"- 样本数: `{len(sample_rows)}`")
lines.append(f"- 主流程待标准化样本数: `{len(run_main_samples)}`")
lines.append(f"- main_upstream_ready: `{main_upstream_ready}`")
lines.append(f"- standardized_inputs: `{'true' if standardized_ok else 'false'}`")
lines.append(f"- velocity_upstream_ready: `{velocity_upstream_ready}`")
lines.append(f"- ambient_upstream_ready: `{ambient_upstream_ready}`")
lines.append(f"- ambient_raw_available: `{project_row.get('raw_matrix_available', 'false')}`")
lines.append(f"- ambient_preferred_method: `{project_row.get('ambient_preferred_method', 'mixed')}`")
lines.append(f"- scenic_upstream_ready: `{scenic_upstream_ready}`")
lines.append(f"- de_replicate_ready: `{de_replicate_ready}`")
lines.append("")
lines.append("## Blockers")
if blocked_samples:
    for sample_id in blocked_samples:
        lines.append(f"- `{sample_id}` 缺少可用标准输入或标准化目录。")
else:
    lines.append("- 当前没有主流程输入 blocker。")
if project_row.get("notes"):
    lines.append("")
    lines.append("## Project Notes")
    for note in describe_notes(project_row["notes"]):
        lines.append(f"- {note}")
lines.append("")
lines.append("## Sample Contracts")
for row in run_main_samples:
    sample_id = row["sample_id"]
    inventory_row = inventory_by_sample.get(sample_id, {})
    readiness_row = readiness_by_sample.get(sample_id, {})
    lines.append(f"- `{sample_id}`: platform=`{inventory_row.get('platform_resolved', '')}`; profile=`{inventory_row.get('feature_name_profile', '')}`; gene_id_type=`{inventory_row.get('gene_id_type_resolved', '')}`; reference_version=`{inventory_row.get('reference_version', '') or 'NA'}`")
    lines.append(f"  filtered_matrix_dir: `{inventory_row.get('filtered_matrix_dir', '') or 'NA'}`")
    lines.append(f"  raw_matrix_dir: `{inventory_row.get('raw_matrix_dir', '') or 'NA'}`")
    lines.append(f"  metrics_path: `{inventory_row.get('metrics_path', '') or 'NA'}`")
    lines.append(f"  bam_path: `{inventory_row.get('bam_path', '') or 'NA'}`")
    lines.append(
        "  ambient: "
        f"preferred=`{readiness_row.get('ambient_preferred_method', 'none')}`; "
        f"fallback=`{readiness_row.get('ambient_fallback_method', 'none')}`; "
        f"soupx_ready=`{readiness_row.get('ambient_soupx_ready', 'false')}`; "
        f"decontx_ready=`{readiness_row.get('ambient_decontx_ready', 'false')}`; "
        f"cellbender_ready=`{readiness_row.get('ambient_cellbender_ready', 'false')}`"
    )
    sample_notes = describe_notes(readiness_row.get("notes", ""))
    if sample_notes:
        lines.append(f"  notes: `{'; '.join(sample_notes)}`")
    ambient_notes = [note for note in readiness_row.get("ambient_notes", "").split(";") if note.strip()]
    if ambient_notes:
        lines.append(f"  ambient_notes: `{'; '.join(ambient_notes)}`")
lines.append("")
lines.append("## Next Step")
lines.append(f"- `{next_step}`")
summary_path.parent.mkdir(parents=True, exist_ok=True)
summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

env_lines = [
    f"MAIN_UPSTREAM_READY={shlex.quote(main_upstream_ready)}",
    f"VELOCITY_UPSTREAM_READY={shlex.quote(velocity_upstream_ready)}",
    f"AMBIENT_UPSTREAM_READY={shlex.quote(ambient_upstream_ready)}",
    f"SCENIC_UPSTREAM_READY={shlex.quote(scenic_upstream_ready)}",
    f"DE_REPLICATE_READY={shlex.quote(de_replicate_ready)}",
    f"STANDARDIZED_INPUTS={shlex.quote('true' if standardized_ok else 'false')}",
    f"NEXT_STEP={shlex.quote(next_step)}",
]
env_path.write_text("\n".join(env_lines) + "\n", encoding="utf-8")
PY

# shellcheck disable=SC1090
source "${SUMMARY_ENV_FILE}"

metadata_valid="$(workflow_status_get "status.metadata_valid" 2>/dev/null || echo false)"
update_workflow_status \
  "input_summary_complete" \
  "${NEXT_STEP}" \
  "status.metadata_valid=${metadata_valid}" \
  "status.main_upstream_ready=$(bool_string "${MAIN_UPSTREAM_READY}")" \
  "status.velocity_upstream_ready=$(bool_string "${VELOCITY_UPSTREAM_READY}")" \
  "status.ambient_upstream_ready=$(bool_string "${AMBIENT_UPSTREAM_READY}")" \
  "status.scenic_upstream_ready=$(bool_string "${SCENIC_UPSTREAM_READY}")" \
  "status.de_replicate_ready=$(bool_string "${DE_REPLICATE_READY}")" \
  "status.standardized_inputs=$(bool_string "${STANDARDIZED_INPUTS}")" \
  "status.input_eda_complete=true"

echo "intake summary 已生成: ${INTAKE_SUMMARY_FILE}"
