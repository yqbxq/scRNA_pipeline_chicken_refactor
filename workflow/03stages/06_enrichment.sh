#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

check_stage_deps "06_enrichment"
hold_for_gate deg

MODULE_05D_MANIFEST="${MANIFEST_DIR}/05d_deg_eda/_manifest.json"
MODULE_06A_MANIFEST="${MANIFEST_DIR}/06a_go_enrichment/_manifest.json"
MODULE_06B_MANIFEST="${MANIFEST_DIR}/06b_kegg_enrichment/_manifest.json"
MODULE_06C_MANIFEST="${MANIFEST_DIR}/06c_enrichment_eda/_manifest.json"

ensure_eda_control_files
require_manifest_output "${MODULE_05D_MANIFEST}" "deg_status_matrix_tsv" >/dev/null

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/06a_go_enrichment.R" \
  "${MODULE_06A_MANIFEST}" \
  "${MODULE_05D_MANIFEST}" \
  "${LAYER_STATUS_FILE}" \
  "${COMPARISON_SHEET}"
require_manifest_output "${MODULE_06A_MANIFEST}" "go_enrichment_manifest_tsv" >/dev/null

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/06b_kegg_enrichment.R" \
  "${MODULE_06B_MANIFEST}" \
  "${MODULE_05D_MANIFEST}" \
  "${LAYER_STATUS_FILE}" \
  "${COMPARISON_SHEET}"
require_manifest_output "${MODULE_06B_MANIFEST}" "kegg_enrichment_manifest_tsv" >/dev/null

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/06c_enrichment_eda.R" \
  "${MODULE_06C_MANIFEST}" \
  "${MODULE_06A_MANIFEST}" \
  "${MODULE_06B_MANIFEST}"
ENRICHMENT_REPORT="$(require_manifest_output "${MODULE_06C_MANIFEST}" "report")"

update_workflow_status \
  "06_enrichment_completed" \
  "review ${ENRICHMENT_REPORT}; next: 07_communication or 08_regulation" \
  "status.06a_go_enrichment_completed=true" \
  "status.06b_kegg_enrichment_completed=true" \
  "status.06c_enrichment_eda_completed=true" \
  "status.06_enrichment_completed=true"
