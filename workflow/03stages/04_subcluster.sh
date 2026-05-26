#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

read_candidate_layers_count() {
  local path="$1"
  if [[ ! -s "${path}" ]]; then
    echo 0
    return
  fi
  awk 'NR == 1 { value = $1 + 0 } END { print value + 0 }' "${path}"
}

list_candidate_review_layers() {
  local review_summary="$1"
  local python_bin
  python_bin="$(detect_python)"
  "${python_bin}" - "${review_summary}" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.exists() or path.stat().st_size == 0:
    raise SystemExit(0)

with path.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row.get("mode") == "candidate":
            layer_id = row.get("layer_id", "")
            selected = row.get("selected_integration_file", "")
            if layer_id:
                print(f"{layer_id}\t{selected}")
PY
}

check_stage_deps "04_subcluster"
hold_for_gate annotation

MODULE_03D_MANIFEST="${MANIFEST_DIR}/03d_annotate/_manifest.json"
PANORAMA_ANNOTATED_OBJECT="$(require_manifest_output "${MODULE_03D_MANIFEST}" "annotated_object")"

MODULE_04A_MANIFEST="${MANIFEST_DIR}/04a_subcluster_build/_manifest.json"
MODULE_04A_REVIEW_MANIFEST="${MANIFEST_DIR}/04a_review/_manifest.json"
MODULE_04B_MANIFEST="${MANIFEST_DIR}/04b_subcluster_annotate/_manifest.json"
MODULE_04C_MANIFEST="${MANIFEST_DIR}/04c_subcluster_eda/_manifest.json"

ensure_eda_control_files
ensure_object_layer_config_file
ensure_marker_panel_dir

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/04a_subcluster_build.R" \
  "${MODULE_04A_MANIFEST}" \
  "${MODULE_03D_MANIFEST}" \
  "${PANORAMA_ANNOTATED_OBJECT}" \
  "${OBJECT_LAYER_CONFIG_FILE}" \
  "${SELECTED_INTEGRATION_FILE}"

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/04a_review.R" \
  "${MODULE_04A_REVIEW_MANIFEST}" \
  "${MODULE_04A_MANIFEST}"

SUBCLUSTER_REVIEW_SUMMARY="$(require_manifest_output "${MODULE_04A_REVIEW_MANIFEST}" "subcluster_review_summary_tsv")"
CANDIDATE_LAYERS_COUNT_FILE="$(require_manifest_output "${MODULE_04A_REVIEW_MANIFEST}" "candidate_layers_count_txt")"
CANDIDATE_LAYERS_COUNT="$(read_candidate_layers_count "${CANDIDATE_LAYERS_COUNT_FILE}")"

if [[ "${CANDIDATE_LAYERS_COUNT}" -gt 0 ]]; then
  if ! eda_gate_passed subcluster; then
    set_eda_gate_status \
      "subcluster" \
      "pending" \
      "" \
      "04a subcluster integration candidates completed; review ${SUBCLUSTER_REVIEW_SUMMARY}, then approve subcluster"
  fi
  hold_for_gate subcluster

  while IFS=$'\t' read -r layer_id selected_file; do
    [[ -n "${layer_id}" ]] || continue
    run_stage_if_manifest_key_missing \
      "${WORKFLOW_ROOT}/05single_script/04a_subcluster_build.R" \
      "${MODULE_04A_MANIFEST}" \
      "clustered_${layer_id}" \
      "${selected_file}" \
      "${OBJECT_LAYER_CONFIG_FILE}" \
      "${MODULE_03D_MANIFEST}"
  done < <(list_candidate_review_layers "${SUBCLUSTER_REVIEW_SUMMARY}")
fi

require_manifest_output "${MODULE_04A_MANIFEST}" "subcluster_index_tsv" >/dev/null
require_manifest_output "${MODULE_04A_MANIFEST}" "candidate_index_tsv" >/dev/null

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/04b_subcluster_annotate.R" \
  "${MODULE_04B_MANIFEST}" \
  "${MODULE_04A_MANIFEST}" \
  "${MARKER_PANEL_DIR}"

require_manifest_output "${MODULE_04B_MANIFEST}" "summary_tsv" >/dev/null

export_h5ad_for_gate \
  "04b_subcluster" \
  "${MODULE_04B_MANIFEST}" \
  "panorama_cell_subtype_rds" \
  "scrna" \
  "04b_subcluster"

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/04c_subcluster_eda.R" \
  "${MODULE_04C_MANIFEST}" \
  "${MODULE_04B_MANIFEST}" \
  "${COMPARISON_SHEET}"

require_manifest_output "${MODULE_04C_MANIFEST}" "subcluster_summary_tsv" >/dev/null

update_workflow_status \
  "04_subcluster_completed" \
  "next: review 04c subcluster EDA outputs; optional next bash ${PIPELINE_ROOT}/workflow/03stages/04d_cluster_robustness.sh" \
  "status.04a_subcluster_build_completed=true" \
  "status.04b_subcluster_annotated=true" \
  "status.04c_subcluster_eda_completed=true" \
  "status.04_subcluster_completed=true"
