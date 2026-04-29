#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

check_stage_deps "02_qc"
hold_for_gate pre_qc

MODULE_01A_MANIFEST="${MANIFEST_DIR}/01a_build_raw_objects/_manifest.json"
MODULE_01B_MANIFEST="${MANIFEST_DIR}/01b_pre_qc_eda/_manifest.json"
MODULE_02A_MANIFEST="${MANIFEST_DIR}/02a_ambient_branch/_manifest.json"
MODULE_02B1_MANIFEST="${MANIFEST_DIR}/02b1_qc_filter/_manifest.json"
MODULE_02B2_MANIFEST="${MANIFEST_DIR}/02b2_doublet/_manifest.json"
MODULE_02C_MANIFEST="${MANIFEST_DIR}/02c_post_qc_eda/_manifest.json"

[[ -f "${INPUT_INVENTORY_FILE}" ]] || die "缺少 input inventory: ${INPUT_INVENTORY_FILE}"
[[ -f "${BRANCH_READINESS_FILE}" ]] || die "缺少 branch readiness: ${BRANCH_READINESS_FILE}"
[[ -f "${QC_THRESHOLD_FILE}" ]] || die "缺少 QC threshold file: ${QC_THRESHOLD_FILE}"

RAW_OBJECT="$(require_manifest_output "${MODULE_01A_MANIFEST}" "raw_object")"
require_manifest_output "${MODULE_01B_MANIFEST}" "sample_qc_summary" >/dev/null
require_manifest_output "${MODULE_01B_MANIFEST}" "report" >/dev/null

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/02a_ambient_branch.R" \
  "${MODULE_02A_MANIFEST}" \
  "${MODULE_01A_MANIFEST}" \
  "${RAW_OBJECT}" \
  "${INPUT_INVENTORY_FILE}" \
  "${BRANCH_READINESS_FILE}"

AMBIENT_OBJECT="$(require_manifest_output "${MODULE_02A_MANIFEST}" "ambient_object")"
require_manifest_output "${MODULE_02A_MANIFEST}" "ambient_readiness" >/dev/null
require_manifest_output "${MODULE_02A_MANIFEST}" "ambient_sample_summary" >/dev/null
require_manifest_output "${MODULE_02A_MANIFEST}" "ambient_method_decisions" >/dev/null
require_manifest_output "${MODULE_02A_MANIFEST}" "ambient_marker_leakage" >/dev/null
require_manifest_output "${MODULE_02A_MANIFEST}" "cellbender_stub" >/dev/null
require_manifest_output "${MODULE_02A_MANIFEST}" "report" >/dev/null

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/02b1_qc_filter.R" \
  "${MODULE_02B1_MANIFEST}" \
  "${MODULE_02A_MANIFEST}" \
  "${AMBIENT_OBJECT}" \
  "${QC_THRESHOLD_FILE}"

QC_FILTER_OBJECT="$(require_manifest_output "${MODULE_02B1_MANIFEST}" "qc_filter_object")"
require_manifest_output "${MODULE_02B1_MANIFEST}" "qc_filter_sample_summary" >/dev/null
require_manifest_output "${MODULE_02B1_MANIFEST}" "qc_filter_diagnostics" >/dev/null

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/02b2_doublet.R" \
  "${MODULE_02B2_MANIFEST}" \
  "${MODULE_02B1_MANIFEST}" \
  "${QC_FILTER_OBJECT}"

POST_QC_OBJECT="$(require_manifest_output "${MODULE_02B2_MANIFEST}" "post_qc_object")"
require_manifest_output "${MODULE_02B2_MANIFEST}" "doublet_sample_summary" >/dev/null
require_manifest_output "${MODULE_02B2_MANIFEST}" "doublet_cell_diagnostics" >/dev/null
require_manifest_output "${MODULE_02B2_MANIFEST}" "doublet_concordance_summary" >/dev/null
require_manifest_output "${MODULE_02B2_MANIFEST}" "doublet_cluster_risk" >/dev/null

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/02c_post_qc_eda.R" \
  "${MODULE_02C_MANIFEST}" \
  "${MODULE_01A_MANIFEST}" \
  "${MODULE_01B_MANIFEST}" \
  "${MODULE_02A_MANIFEST}" \
  "${MODULE_02B1_MANIFEST}" \
  "${MODULE_02B2_MANIFEST}" \
  "${POST_QC_OBJECT}"

require_manifest_output "${MODULE_02C_MANIFEST}" "post_qc_sample_summary" >/dev/null
require_manifest_output "${MODULE_02C_MANIFEST}" "post_qc_triage" >/dev/null
require_manifest_output "${MODULE_02C_MANIFEST}" "qc_retention_comparison" >/dev/null
require_manifest_output "${MODULE_02C_MANIFEST}" "post_qc_report" >/dev/null
require_manifest_output "${MODULE_02C_MANIFEST}" "pre_vs_post_qc_violin_png" >/dev/null
require_manifest_output "${MODULE_02C_MANIFEST}" "pre_vs_post_qc_violin_pdf" >/dev/null
require_manifest_output "${MODULE_02C_MANIFEST}" "doublet_score_distribution_png" >/dev/null
require_manifest_output "${MODULE_02C_MANIFEST}" "doublet_score_distribution_pdf" >/dev/null
require_manifest_output "${MODULE_02C_MANIFEST}" "cluster_doublet_heatmap_png" >/dev/null
require_manifest_output "${MODULE_02C_MANIFEST}" "cluster_doublet_heatmap_pdf" >/dev/null

set_eda_gate_status \
  "post_qc" \
  "pending" \
  "" \
  "02_qc completed; manual review required before integration"

update_workflow_status \
  "02_qc_completed" \
  "review ${EDA_GATE_FILE} and approve post_qc before integration" \
  "status.02a_ambient_branch_completed=true" \
  "status.02b1_qc_filter_completed=true" \
  "status.02b2_doublet_completed=true" \
  "status.02c_post_qc_eda_completed=true" \
  "status.02_qc_completed=true" \
  "status.post_qc_eda_complete=true" \
  "status.post_qc_gate_passed=false"
