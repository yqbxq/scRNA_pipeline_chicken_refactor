#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

MODULE_00_MANIFEST="${ORTHOLOG_MANIFEST}"
MODULE_03D_MANIFEST="${MANIFEST_DIR}/03d_annotate/_manifest.json"
MODULE_04B_MANIFEST="${MANIFEST_DIR}/04b_subcluster_annotate/_manifest.json"
MODULE_08A_MANIFEST="${MANIFEST_DIR}/08a_scenic_export/_manifest.json"
MODULE_08B_MANIFEST="${MANIFEST_DIR}/08b_scenic_grn/_manifest.json"
MODULE_08C_MANIFEST="${MANIFEST_DIR}/08c_scenic_regulons/_manifest.json"
MODULE_08D_MANIFEST="${MANIFEST_DIR}/08d_scenic_downstream/_manifest.json"
MODULE_08E_MANIFEST="${MANIFEST_DIR}/08e_decoupler/_manifest.json"
MODULE_08F_MANIFEST="${MANIFEST_DIR}/08f_regulation_eda/_manifest.json"

RUN_SCENIC="${RUN_SCENIC:-yes}"
RUN_SCENIC_GRN="${RUN_SCENIC_GRN:-yes}"
RUN_DECOUPLER="${RUN_DECOUPLER:-yes}"
RUN_REGULATION_EDA="${RUN_REGULATION_EDA:-yes}"

download_if_missing_08() {
  local url="$1"
  local out="$2"
  ensure_dir "$(dirname "${out}")"
  if [[ -s "${out}" ]]; then
    echo "已存在 ${out}，跳过下载"
    return 0
  fi

  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --output "${out}" "${url}"
    return 0
  fi

  if command -v wget >/dev/null 2>&1; then
    wget -O "${out}" "${url}"
    return 0
  fi

  die "系统中既没有 curl，也没有 wget，无法下载 SCENIC 资源。"
}

write_scenic_resource_manifest_08() {
  local python_bin
  python_bin="$(detect_python)"
  env \
    SCENIC_RESOURCE_MANIFEST="${SCENIC_RESOURCE_MANIFEST}" \
    SCENIC_TF_LIST="${SCENIC_TF_LIST}" \
    SCENIC_MOTIF_ANN="${SCENIC_MOTIF_ANN}" \
    SCENIC_DB_500BP="${SCENIC_DB_500BP}" \
    SCENIC_DB_10KB="${SCENIC_DB_10KB}" \
    "${python_bin}" - <<'PY'
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

items = {
    "tf_list": os.environ["SCENIC_TF_LIST"],
    "motif_annotations": os.environ["SCENIC_MOTIF_ANN"],
    "rankings_500bp": os.environ["SCENIC_DB_500BP"],
    "rankings_10kb": os.environ["SCENIC_DB_10KB"],
}

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

resources = {}
missing = []
for key, raw in items.items():
    path = Path(raw)
    if not path.exists() or path.stat().st_size == 0:
        missing.append(str(path))
        continue
    resources[key] = {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }

if missing:
    raise SystemExit("missing SCENIC resource files: " + ", ".join(missing))

manifest = {
    "module": "scenic_resources",
    "timestamp": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
    "resources": resources,
}
out = Path(os.environ["SCENIC_RESOURCE_MANIFEST"])
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY
}

download_scenic_resources_inline() {
  ensure_dir "${SCENIC_DB_DIR}" "${LOG_DIR}" "${SCENIC_INPUT_DIR}" "${SCENIC_OUTPUT_DIR}"

  download_if_missing_08 \
    "https://raw.githubusercontent.com/aertslab/SCENIC/master/inst/extdata/hs_hgnc_tfs.txt" \
    "${SCENIC_TF_LIST}"

  download_if_missing_08 \
    "https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
    "${SCENIC_MOTIF_ANN}"

  download_if_missing_08 \
    "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather" \
    "${SCENIC_DB_500BP}"

  download_if_missing_08 \
    "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather" \
    "${SCENIC_DB_10KB}"

  write_scenic_resource_manifest_08
  update_workflow_status \
    "scenic_resources_ready" \
    "bash ${PIPELINE_ROOT}/workflow/03stages/08_regulation.sh" \
    "status.scenic_resources_ready=true"
  echo "SCENIC 资源检查/下载完成: ${SCENIC_RESOURCE_MANIFEST}"
}

if [[ "${SCENIC_RESOURCES_ONLY:-no}" == "yes" ]]; then
  download_scenic_resources_inline
  exit 0
fi

check_stage_deps "08_regulation"
hold_for_gate annotation

ensure_eda_control_files
require_manifest_output "${MODULE_00_MANIFEST}" "human_best" >/dev/null
require_manifest_output "${MODULE_03D_MANIFEST}" "annotated_object" >/dev/null

if [[ "${RUN_SCENIC}" == "yes" && "${RUN_SCENIC_GRN}" == "yes" ]]; then
  download_scenic_resources_inline
fi

REGULATION_RERAN=0
REGULATION_SCENIC_INDEX=""
REGULATION_DECOUPLER_INDEX=""
REGULATION_REPORT=""
STATUS_08A="skipped_by_config"
STATUS_08B="skipped_by_config"
STATUS_08C="skipped_by_config"
STATUS_08D="skipped_by_config"
STATUS_08E="skipped_by_config"
STATUS_08F="skipped_by_config"

run_reg_r_main_if_stale() {
  local script_path="$1"
  local output_path="$2"
  shift 2 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "运行 ${script_path}"
    run_r_main "${script_path}"
    REGULATION_RERAN=1
  else
    echo "已存在且未过期，跳过: ${output_path}"
  fi
}

run_reg_r_scenic_if_stale() {
  local script_path="$1"
  local output_path="$2"
  shift 2 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "运行 ${script_path}"
    run_r_scenic "${script_path}"
    REGULATION_RERAN=1
  else
    echo "已存在且未过期，跳过: ${output_path}"
  fi
}

run_reg_r_decoupler_if_stale() {
  local script_path="$1"
  local output_path="$2"
  shift 2 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "运行 ${script_path}"
    run_r_decoupler "${script_path}"
    REGULATION_RERAN=1
  else
    echo "已存在且未过期，跳过: ${output_path}"
  fi
}

run_reg_shell_if_stale() {
  local script_path="$1"
  local output_path="$2"
  shift 2 || true

  if is_stale_output "${output_path}" "$@"; then
    echo "运行 ${script_path}"
    bash "${script_path}"
    REGULATION_RERAN=1
  else
    echo "已存在且未过期，跳过: ${output_path}"
  fi
}

if [[ "${RUN_SCENIC}" == "yes" ]]; then
  run_reg_r_main_if_stale \
    "${WORKFLOW_ROOT}/05single_script/08a_scenic_export.R" \
    "${MODULE_08A_MANIFEST}" \
    "${MODULE_03D_MANIFEST}" \
    "${MODULE_04B_MANIFEST}" \
    "${MODULE_00_MANIFEST}"
  require_manifest_output "${MODULE_08A_MANIFEST}" "scenic_export_index_tsv" >/dev/null
  STATUS_08A="completed"

  if [[ "${RUN_SCENIC_GRN}" == "yes" ]]; then
    run_reg_shell_if_stale \
      "${WORKFLOW_ROOT}/05single_script/08b_scenic_grn.sh" \
      "${MODULE_08B_MANIFEST}" \
      "${MODULE_08A_MANIFEST}" \
      "${SCENIC_RESOURCE_MANIFEST}" \
      "${SCENIC_TF_LIST}"
    require_manifest_output "${MODULE_08B_MANIFEST}" "scenic_grn_index_tsv" >/dev/null
    STATUS_08B="completed"

    run_reg_r_scenic_if_stale \
      "${WORKFLOW_ROOT}/05single_script/08c_scenic_regulons.R" \
      "${MODULE_08C_MANIFEST}" \
      "${MODULE_08A_MANIFEST}" \
      "${MODULE_08B_MANIFEST}" \
      "${SCENIC_RESOURCE_MANIFEST}" \
      "${SCENIC_MOTIF_ANN}" \
      "${SCENIC_DB_500BP}" \
      "${SCENIC_DB_10KB}"
    require_manifest_output "${MODULE_08C_MANIFEST}" "scenic_regulons_index_tsv" >/dev/null
    STATUS_08C="completed"

    run_reg_r_scenic_if_stale \
      "${WORKFLOW_ROOT}/05single_script/08d_scenic_downstream.R" \
      "${MODULE_08D_MANIFEST}" \
      "${MODULE_08C_MANIFEST}" \
      "${MODULE_03D_MANIFEST}" \
      "${MODULE_04B_MANIFEST}"
    REGULATION_SCENIC_INDEX="$(require_manifest_output "${MODULE_08D_MANIFEST}" "scenic_downstream_index_tsv")"
    STATUS_08D="completed"
  else
    STATUS_08B="skipped_by_config"
    STATUS_08C="skipped_due_to_08b"
    STATUS_08D="skipped_due_to_08c"
  fi
fi

if [[ "${RUN_DECOUPLER}" == "yes" ]]; then
  run_reg_r_decoupler_if_stale \
    "${WORKFLOW_ROOT}/05single_script/08e_decoupler.R" \
    "${MODULE_08E_MANIFEST}" \
    "${MODULE_03D_MANIFEST}" \
    "${MODULE_04B_MANIFEST}" \
    "${MODULE_00_MANIFEST}"
  REGULATION_DECOUPLER_INDEX="$(require_manifest_output "${MODULE_08E_MANIFEST}" "decoupler_index_tsv")"
  STATUS_08E="completed"
fi

if [[ "${RUN_REGULATION_EDA}" == "yes" ]]; then
  REGULATION_EDA_DEPS=()
  for dep_path in \
    "${MODULE_08A_MANIFEST}" \
    "${MODULE_08B_MANIFEST}" \
    "${MODULE_08C_MANIFEST}" \
    "${MODULE_08D_MANIFEST}" \
    "${MODULE_08E_MANIFEST}" \
    "${MANIFEST_DIR}/05d_deg_eda/_manifest.json" \
    "${MANIFEST_DIR}/06c_enrichment_eda/_manifest.json" \
    "${MANIFEST_DIR}/07c_communication_eda/_manifest.json"; do
    if [[ -e "${dep_path}" ]]; then
      REGULATION_EDA_DEPS+=("${dep_path}")
    fi
  done
  run_reg_r_main_if_stale \
    "${WORKFLOW_ROOT}/05single_script/08f_regulation_eda.R" \
    "${MODULE_08F_MANIFEST}" \
    "${REGULATION_EDA_DEPS[@]}"
  REGULATION_REPORT="$(require_manifest_output "${MODULE_08F_MANIFEST}" "report_md")"
  if [[ "${STATUS_08D}" == "completed" && "${STATUS_08E}" == "completed" ]]; then
    STATUS_08F="completed"
  else
    STATUS_08F="completed_partial"
  fi
fi

if [[ "${REGULATION_RERAN}" == "1" ]]; then
  set_eda_gate_status \
    "regulation" \
    "pending" \
    "" \
    "08 regulation modules updated; review report ${REGULATION_REPORT:-NA}, SCENIC index ${REGULATION_SCENIC_INDEX:-NA}, and decoupleR index ${REGULATION_DECOUPLER_INDEX:-NA}."
fi

sync_workflow_gate_statuses
update_workflow_status \
  "08_regulation_completed" \
  "review ${REGULATION_REPORT:-NA}" \
  "status.08a_scenic_export_status=${STATUS_08A}" \
  "status.08b_scenic_grn_status=${STATUS_08B}" \
  "status.08c_scenic_regulons_status=${STATUS_08C}" \
  "status.08d_scenic_downstream_status=${STATUS_08D}" \
  "status.08e_decoupler_status=${STATUS_08E}" \
  "status.08f_regulation_eda_status=${STATUS_08F}" \
  "status.08_regulation_scenic_completed=$([[ "${STATUS_08D}" == "completed" ]] && echo true || echo false)" \
  "status.08_regulation_decoupler_completed=$([[ "${STATUS_08E}" == "completed" ]] && echo true || echo false)" \
  "status.08_regulation_completed=$([[ "${STATUS_08F}" == completed* ]] && echo true || echo false)" \
  "status.regulation_gate_passed=$(eda_gate_passed regulation && echo true || echo false)"
