#!/usr/bin/env bash
set -euo pipefail

LIB_DIR="$(cd "$(dirname "$0")/../02lib" && pwd)"
source "${LIB_DIR}/common.sh"
source "${LIB_DIR}/intake.sh"

prepare_project_state_dirs

INPUT_SHEET="${CANONICAL_SAMPLE_SHEET}"
if [[ ! -f "${INPUT_SHEET}" ]]; then
  die "缺少 canonical sample sheet，请先运行 02_validate_metadata.sh"
fi

ensure_dir "${DATA_DIR}" "${LOG_DIR}"

python_bin="$(detect_python)"

env \
  INPUT_SHEET="${INPUT_SHEET}" \
  DATA_DIR="${DATA_DIR}" \
  CELLRANGER_OUT_DIR="${CELLRANGER_OUT_DIR}" \
  DNBC4TOOLS_OUT_DIR="${DNBC4TOOLS_OUT_DIR}" \
  FASTQ_DIR="${FASTQ_DIR}" \
  INPUT_INVENTORY_FILE="${INPUT_INVENTORY_FILE}" \
  INPUT_STANDARDIZE_MODE="${INPUT_STANDARDIZE_MODE}" \
  PIPELINE_ROOT="${PIPELINE_ROOT}" \
  "${python_bin}" - <<'PY'
import csv
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ["PIPELINE_ROOT"]) / "shell" / "04python"))

from intake_contract import is_mex_matrix_dir, resolve_sample_contract

input_sheet = Path(os.environ["INPUT_SHEET"])
data_dir = Path(os.environ["DATA_DIR"])
cellranger_out_dir = Path(os.environ["CELLRANGER_OUT_DIR"])
dnbc4tools_out_dir = Path(os.environ["DNBC4TOOLS_OUT_DIR"])
fastq_dir = Path(os.environ["FASTQ_DIR"])
inventory_path = Path(os.environ["INPUT_INVENTORY_FILE"])
standardize_mode = os.environ["INPUT_STANDARDIZE_MODE"]


def replace_path(dst: Path, src: Path):
    if dst.exists() or dst.is_symlink():
        if dst.is_symlink() or dst.is_file():
            dst.unlink()
        else:
            shutil.rmtree(dst)
    if standardize_mode == "symlink":
        dst.symlink_to(src, target_is_directory=src.is_dir())
    elif standardize_mode == "copy":
        if src.is_dir():
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)
    else:
        raise SystemExit("INPUT_STANDARDIZE_MODE 仅支持 symlink 或 copy")

with input_sheet.open("r", encoding="utf-8", newline="") as handle:
    sample_rows = list(csv.DictReader(handle, delimiter="\t"))

inventory_by_sample = {}
if inventory_path.exists():
    with inventory_path.open("r", encoding="utf-8", newline="") as handle:
        inventory_rows = list(csv.DictReader(handle, delimiter="\t"))
    inventory_by_sample = {row["sample_id"]: row for row in inventory_rows if row.get("sample_id")}

for row in sample_rows:
    if row.get("run_main", "yes").strip().lower() != "yes":
        continue

    sample_id = row["sample_id"].strip()
    inventory_row = inventory_by_sample.get(sample_id)
    matrix_dir = None
    raw_matrix_dir = None
    if inventory_row is not None and inventory_row.get("filtered_matrix_dir"):
        matrix_dir = Path(inventory_row["filtered_matrix_dir"])
    if inventory_row is not None and inventory_row.get("raw_matrix_dir"):
        raw_matrix_dir = Path(inventory_row["raw_matrix_dir"])
    if matrix_dir is None:
        contract = resolve_sample_contract(
            row=row,
            fastq_dir=fastq_dir,
            cellranger_out_dir=cellranger_out_dir,
            data_dir=data_dir,
            dnbc4tools_out_dir=dnbc4tools_out_dir,
        )
        if contract.get("filtered_matrix_dir"):
            matrix_dir = Path(contract["filtered_matrix_dir"])
        if contract.get("raw_matrix_dir"):
            raw_matrix_dir = Path(contract["raw_matrix_dir"])
    if matrix_dir is None or not matrix_dir.exists():
        raise SystemExit(f"无法解析 {sample_id} 的标准化输入目录")
    if not is_mex_matrix_dir(matrix_dir):
        raise SystemExit(f"输入目录不是有效的 10X matrix: {matrix_dir}")

    dst_dir = data_dir / sample_id
    replace_path(dst_dir, matrix_dir)

    if raw_matrix_dir is not None and raw_matrix_dir.exists():
        raw_dst_dir = data_dir / f"{sample_id}_raw"
        replace_path(raw_dst_dir, raw_matrix_dir)
PY

update_workflow_status \
  "inputs_standardized" \
  "bash ${PIPELINE_ROOT}/workflow/03stages/input_summary.sh" \
  "status.standardized_inputs=true"

echo "主流程输入标准化完成: ${DATA_DIR}"
