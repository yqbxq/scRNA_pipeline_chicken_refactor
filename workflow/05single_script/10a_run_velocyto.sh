#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

MODULE_NAME="10a_run_velocyto"
MODULE_VERSION="${MODULE_10A_VERSION:-${MODULE_10_VERSION:-1.0}}"
PROJECT_ROOT="${PROJECT_ROOT:-${PIPELINE_ROOT}}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_ROOT}/results}"
TABLE_DIR="${TABLE_DIR:-${RESULTS_DIR}/tables}"
MANIFEST_DIR="${MANIFEST_DIR:-${RESULTS_DIR}/manifests}"
DATA_DIR="${DATA_DIR:-${PROJECT_ROOT}/data}"
METADATA_DIR="${METADATA_DIR:-${PROJECT_ROOT}/metadata}"
VELOCITY_DIR="${VELOCITY_DIR:-${RESULTS_DIR}/velocity}"
VELOCITY_INPUT_DIR="${VELOCITY_INPUT_DIR:-${VELOCITY_DIR}/input}"
VELOCITY_LOOM_DIR="${VELOCITY_LOOM_DIR:-${VELOCITY_DIR}/loom}"
VELOCITY_OUTPUT_DIR="${VELOCITY_OUTPUT_DIR:-${VELOCITY_DIR}/output}"
DNBC4TOOLS_OUT_DIR="${DNBC4TOOLS_OUT_DIR:-${DATA_DIR}/dnbc4tools_out}"
CELLRANGER_OUT_DIR="${CELLRANGER_OUT_DIR:-}"
SAMPLE_SHEET="${SAMPLE_SHEET:-${METADATA_DIR}/samples.tsv}"
CANONICAL_SAMPLE_SHEET="${CANONICAL_SAMPLE_SHEET:-${METADATA_DIR}/samples.canonical.tsv}"
REFERENCE_GTF="${REFERENCE_GTF:-${PROJECT_ROOT}/reference/genes.gtf}"
CLEAN_GTF="${CLEAN_GTF:-${PROJECT_ROOT}/reference/genes.clean.gtf}"
VELOCITY_GTF="${VELOCITY_GTF:-${CLEAN_GTF}}"
VELOCITY_BAM_PATTERN="${VELOCITY_BAM_PATTERN:-anno_decon_sorted.bam}"
VELOCITY_H5_PATTERN="${VELOCITY_H5_PATTERN:-filtered_feature_bc_matrix.h5}"
VELOCYTO_THREADS="${VELOCYTO_THREADS:-${MAIN_THREADS:-8}}"
VELOCYTO_REPEAT_MASK_GTF="${VELOCYTO_REPEAT_MASK_GTF:-}"

VELOCITY_TABLE_DIR="${TABLE_DIR}/velocity"
VELOCITY_INPUT_TABLE_DIR="${VELOCITY_TABLE_DIR}/inputs"
VELOCITY_BARCODE_DIR="${VELOCITY_INPUT_DIR}/barcodes"
MODULE_MANIFEST="${MANIFEST_DIR}/${MODULE_NAME}/_manifest.json"
LOOM_INDEX_TSV="${VELOCITY_INPUT_TABLE_DIR}/velocity_loom_index.tsv"
LOOM_JOBS_TSV="${VELOCITY_INPUT_TABLE_DIR}/velocity_loom_jobs.tsv"

mkdir -p \
  "${VELOCITY_INPUT_DIR}" \
  "${VELOCITY_LOOM_DIR}" \
  "${VELOCITY_OUTPUT_DIR}" \
  "${VELOCITY_INPUT_TABLE_DIR}" \
  "${VELOCITY_BARCODE_DIR}" \
  "$(dirname "${MODULE_MANIFEST}")"

choose_python_10a() {
  if command -v python >/dev/null 2>&1; then
    echo "python"
    return 0
  fi
  if command -v python3 >/dev/null 2>&1; then
    echo "python3"
    return 0
  fi
  return 1
}

PYTHON_BIN="$(choose_python_10a || true)"

normalize_sample_tokens_10a() {
  awk '{
    gsub(/[,;]/, " ")
    for (i = 1; i <= NF; i++) {
      if ($i != "") print $i
    }
  }'
}

sample_sheet_path_10a() {
  if [[ -s "${CANONICAL_SAMPLE_SHEET}" ]]; then
    echo "${CANONICAL_SAMPLE_SHEET}"
  else
    echo "${SAMPLE_SHEET}"
  fi
}

detect_velocity_samples_10a() {
  if [[ -n "${VELOCITY_SAMPLE_IDS:-}" ]]; then
    printf '%s\n' "${VELOCITY_SAMPLE_IDS}" | normalize_sample_tokens_10a
    return 0
  fi

  local sheet_path
  sheet_path="$(sample_sheet_path_10a)"
  if [[ -s "${sheet_path}" && -n "${PYTHON_BIN}" ]]; then
    "${PYTHON_BIN}" - "${sheet_path}" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
with path.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader((line for line in handle if line.strip() and not line.startswith("#")), delimiter="\t")
    for row in reader:
        sample_id = (row.get("sample_id") or "").strip()
        if not sample_id:
            continue
        run_velocity = (row.get("run_velocity") or "auto").strip().lower()
        if run_velocity != "no":
            print(sample_id)
PY
    return 0
  fi

  if [[ -n "${RAW_SAMPLES:-}" ]]; then
    printf '%s\n' "${RAW_SAMPLES}" | normalize_sample_tokens_10a
    return 0
  fi

  printf '%s\n' "${SAMPLE_NAMES:-}" | normalize_sample_tokens_10a
}

first_existing_file_10a() {
  local candidate
  for candidate in "$@"; do
    [[ -n "${candidate}" ]] || continue
    if [[ -s "${candidate}" ]]; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

find_first_matching_10a() {
  local root="$1"
  local sample_id="$2"
  local pattern="$3"
  local base_name
  base_name="$(basename "${pattern}")"
  [[ -d "${root}" ]] || return 1
  find "${root}" -type f -name "${base_name}" -path "*/${sample_id}/*" -print -quit 2>/dev/null
}

resolve_bam_10a() {
  local sample_id="$1"
  local path=""
  path="$(first_existing_file_10a \
    "${DNBC4TOOLS_OUT_DIR}/${sample_id}/${VELOCITY_BAM_PATTERN}" \
    "${DNBC4TOOLS_OUT_DIR}/${sample_id}/outs/${VELOCITY_BAM_PATTERN}" \
    "${DATA_DIR}/${sample_id}/${VELOCITY_BAM_PATTERN}" \
    "${DATA_DIR}/${sample_id}/outs/${VELOCITY_BAM_PATTERN}" || true)"
  if [[ -z "${path}" ]]; then
    path="$(find_first_matching_10a "${DNBC4TOOLS_OUT_DIR}" "${sample_id}" "${VELOCITY_BAM_PATTERN}" || true)"
  fi
  echo "${path}"
}

resolve_matrix_h5_10a() {
  local sample_id="$1"
  local path=""
  path="$(first_existing_file_10a \
    "${DNBC4TOOLS_OUT_DIR}/${sample_id}/${VELOCITY_H5_PATTERN}" \
    "${DNBC4TOOLS_OUT_DIR}/${sample_id}/outs/${VELOCITY_H5_PATTERN}" \
    "${DATA_DIR}/${sample_id}/${VELOCITY_H5_PATTERN}" \
    "${DATA_DIR}/${sample_id}/outs/${VELOCITY_H5_PATTERN}" \
    "${CELLRANGER_OUT_DIR}/${sample_id}/outs/${VELOCITY_H5_PATTERN}" || true)"
  if [[ -z "${path}" ]]; then
    path="$(find_first_matching_10a "${DNBC4TOOLS_OUT_DIR}" "${sample_id}" "${VELOCITY_H5_PATTERN}" || true)"
  fi
  if [[ -z "${path}" && -n "${CELLRANGER_OUT_DIR}" ]]; then
    path="$(find_first_matching_10a "${CELLRANGER_OUT_DIR}" "${sample_id}" "${VELOCITY_H5_PATTERN}" || true)"
  fi
  echo "${path}"
}

extract_whitelist_10a() {
  local matrix_h5="$1"
  local whitelist_path="$2"
  [[ -n "${PYTHON_BIN}" ]] || return 1
  mkdir -p "$(dirname "${whitelist_path}")"
  "${PYTHON_BIN}" - "${matrix_h5}" "${whitelist_path}" <<'PY'
import sys
from pathlib import Path

try:
    import h5py
except Exception as exc:
    raise SystemExit(f"h5py unavailable: {exc}")

matrix_h5 = Path(sys.argv[1])
out_path = Path(sys.argv[2])

with h5py.File(matrix_h5, "r") as handle:
    dataset = None
    for key in ("matrix/barcodes", "barcodes"):
        if key in handle:
            dataset = handle[key]
            break
    if dataset is None:
        raise SystemExit(f"no barcode dataset found in {matrix_h5}")
    barcodes = []
    for value in dataset[:]:
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        barcodes.append(str(value).strip())

barcodes = [value for value in barcodes if value]
if not barcodes:
    raise SystemExit(f"no barcodes found in {matrix_h5}")

out_path.write_text("\n".join(barcodes) + "\n", encoding="utf-8")
PY
}

write_manifest_10a() {
  [[ -n "${PYTHON_BIN}" ]] || return 1
  "${PYTHON_BIN}" - \
    "${MODULE_MANIFEST}" \
    "${PROJECT_ROOT}" \
    "${MODULE_NAME}" \
    "${MODULE_VERSION}" \
    "${LOOM_INDEX_TSV}" \
    "${LOOM_JOBS_TSV}" \
    "${SAMPLE_SHEET}" \
    "${CANONICAL_SAMPLE_SHEET}" \
    "${DNBC4TOOLS_OUT_DIR}" \
    "${CELLRANGER_OUT_DIR}" \
    "${VELOCITY_GTF}" <<'PY'
import csv
import json
import re
import sys
from datetime import datetime
from pathlib import Path

manifest_path = Path(sys.argv[1])
project_root = Path(sys.argv[2])
module_name = sys.argv[3]
version = sys.argv[4]
index_tsv = Path(sys.argv[5])
jobs_tsv = Path(sys.argv[6])
sample_sheet = sys.argv[7]
canonical_sample_sheet = sys.argv[8]
dnbc4tools_out_dir = sys.argv[9]
cellranger_out_dir = sys.argv[10]
velocity_gtf = sys.argv[11]

def rel(path):
    path = Path(path)
    try:
        return str(path.resolve().relative_to(project_root.resolve()))
    except Exception:
        return str(path)

def entry(path, typ, semantics, schema=None):
    out = {
        "path": rel(path),
        "type": typ,
        "produced_by": module_name,
        "row_semantics": semantics,
    }
    if schema is not None:
        out["schema"] = schema
    return out

def safe_id(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    return value.strip("_") or "NA"

index_schema = {
    "sample_id": "character",
    "bam_path": "character",
    "matrix_h5": "character",
    "whitelist_path": "character",
    "loom_path": "character",
    "status": "character",
    "reason": "character",
    "runtime_s": "numeric",
    "command": "character",
}
jobs_schema = {
    "sample_id": "character",
    "bam_path": "character",
    "matrix_h5": "character",
    "whitelist_path": "character",
    "loom_path": "character",
    "command": "character",
}

outputs = {
    "velocity_loom_index": entry(index_tsv, "tsv", "one row per requested velocity loom sample", index_schema),
    "velocity_loom_jobs": entry(jobs_tsv, "tsv", "velocyto command table by sample", jobs_schema),
}

if index_tsv.exists():
    with index_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            status = (row.get("status") or "").strip()
            loom_path = row.get("loom_path") or ""
            sample_id = row.get("sample_id") or ""
            if status in {"ok", "ok_existing"} and loom_path and Path(loom_path).exists():
                outputs[f"velocity_loom__{safe_id(sample_id)}"] = entry(
                    loom_path,
                    "loom",
                    "velocyto loom file for one sample",
                )

manifest = {
    "module": module_name,
    "version": version,
    "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
    "base_dir": str(project_root),
    "inputs": {
        "sample_sheet": sample_sheet,
        "canonical_sample_sheet": canonical_sample_sheet,
        "dnbc4tools_out_dir": dnbc4tools_out_dir,
        "cellranger_out_dir": cellranger_out_dir,
        "velocity_gtf": velocity_gtf,
    },
    "outputs": outputs,
    "depends_on": {},
}

manifest_path.parent.mkdir(parents=True, exist_ok=True)
manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
PY
}

append_index_row_10a() {
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$@" >> "${LOOM_INDEX_TSV}"
}

append_job_row_10a() {
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$@" >> "${LOOM_JOBS_TSV}"
}

if [[ ! -s "${VELOCITY_GTF}" && -s "${REFERENCE_GTF}" ]]; then
  VELOCITY_GTF="${REFERENCE_GTF}"
fi

{
  printf 'sample_id\tbam_path\tmatrix_h5\twhitelist_path\tloom_path\tstatus\treason\truntime_s\tcommand\n'
} > "${LOOM_INDEX_TSV}"
{
  printf 'sample_id\tbam_path\tmatrix_h5\twhitelist_path\tloom_path\tcommand\n'
} > "${LOOM_JOBS_TSV}"

sample_ids=()
while IFS= read -r sample_id; do
  [[ -n "${sample_id}" ]] || continue
  sample_ids+=("${sample_id}")
done < <(detect_velocity_samples_10a | awk '!seen[$0]++')

if [[ "${#sample_ids[@]}" -eq 0 ]]; then
  echo "10a: no samples selected for RNA velocity loom generation."
fi

for sample_id in "${sample_ids[@]}"; do
  bam_path="$(resolve_bam_10a "${sample_id}")"
  matrix_h5="$(resolve_matrix_h5_10a "${sample_id}")"
  whitelist_path="${VELOCITY_BARCODE_DIR}/${sample_id}_barcodes.tsv"
  loom_path="${VELOCITY_LOOM_DIR}/${sample_id}.loom"
  runtime_s=0
  command_text=""

  if [[ -z "${bam_path}" || -z "${matrix_h5}" ]]; then
    reason="missing BAM or filtered_feature_bc_matrix.h5"
    append_index_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "missing_input" "${reason}" "${runtime_s}" "${command_text}"
    append_job_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "${command_text}"
    continue
  fi

  if [[ ! -s "${VELOCITY_GTF}" ]]; then
    reason="missing velocity GTF"
    append_index_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "failed_missing_gtf" "${reason}" "${runtime_s}" "${command_text}"
    append_job_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "${command_text}"
    continue
  fi

  if [[ ! -s "${whitelist_path}" || "${matrix_h5}" -nt "${whitelist_path}" ]]; then
    if ! extract_whitelist_10a "${matrix_h5}" "${whitelist_path}"; then
      reason="failed to extract whitelist barcodes from filtered_feature_bc_matrix.h5"
      append_index_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "failed_whitelist" "${reason}" "${runtime_s}" "${command_text}"
      append_job_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "${command_text}"
      continue
    fi
  fi

  cmd=(velocyto run -b "${whitelist_path}" -o "${VELOCITY_LOOM_DIR}" --sampleid "${sample_id}")
  if [[ -n "${VELOCYTO_REPEAT_MASK_GTF}" && -s "${VELOCYTO_REPEAT_MASK_GTF}" ]]; then
    cmd+=(-m "${VELOCYTO_REPEAT_MASK_GTF}")
  fi
  if [[ "${VELOCYTO_THREADS}" =~ ^[0-9]+$ && "${VELOCYTO_THREADS}" -gt 0 ]]; then
    cmd+=(--samtools-threads "${VELOCYTO_THREADS}")
  fi
  cmd+=("${bam_path}" "${VELOCITY_GTF}")
  command_text="${cmd[*]}"
  append_job_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "${command_text}"

  if [[ -s "${loom_path}" && "${loom_path}" -nt "${bam_path}" && "${loom_path}" -nt "${matrix_h5}" ]]; then
    append_index_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "ok_existing" "existing loom is current" "${runtime_s}" "${command_text}"
    continue
  fi

  if ! command -v velocyto >/dev/null 2>&1; then
    reason="velocyto executable is unavailable in VELOCITY_ENV_PREFIX"
    append_index_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "failed_no_velocyto" "${reason}" "${runtime_s}" "${command_text}"
    continue
  fi

  start_ts="$(date +%s)"
  if "${cmd[@]}"; then
    end_ts="$(date +%s)"
    runtime_s=$((end_ts - start_ts))
    if [[ -s "${loom_path}" ]]; then
      append_index_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "ok" "" "${runtime_s}" "${command_text}"
    else
      append_index_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "failed_missing_output" "velocyto completed but loom output is missing" "${runtime_s}" "${command_text}"
    fi
  else
    end_ts="$(date +%s)"
    runtime_s=$((end_ts - start_ts))
    append_index_row_10a "${sample_id}" "${bam_path}" "${matrix_h5}" "${whitelist_path}" "${loom_path}" "failed" "velocyto run failed" "${runtime_s}" "${command_text}"
  fi
done

write_manifest_10a
echo "10a completed. velocity loom index: ${LOOM_INDEX_TSV}"
