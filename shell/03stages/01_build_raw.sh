#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
SHELL_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${SHELL_ROOT}/.." && pwd)}"

source "${SHELL_ROOT}/02lib/common.sh"

check_stage_deps "01_build_raw"

ORTHOLOG_MANIFEST="${ORTHOLOG_MANIFEST:-${ORTHOLOG_CACHE_DIR}/_manifest.json}"
CC_GENES_OUTPUT="$(require_manifest_output "${ORTHOLOG_MANIFEST}" "cc_genes")"

MODULE_01A_MANIFEST="${MANIFEST_DIR}/01a_build_raw_objects/_manifest.json"
MODULE_01B_MANIFEST="${MANIFEST_DIR}/01b_pre_qc_eda/_manifest.json"
PRE_QC_REQUIRED_OUTPUTS=(
  "${MODULE_01A_MANIFEST}"
  "${MODULE_01B_MANIFEST}"
  "${PRE_QC_REPORT_DIR}/sample_qc_summary.tsv"
  "${PRE_QC_REPORT_DIR}/triage.tsv"
  "${PRE_QC_REPORT_DIR}/report.md"
  "${PRE_QC_REPORT_DIR}/pre_qc_violin.png"
  "${PRE_QC_REPORT_DIR}/pre_qc_violin.pdf"
  "${PRE_QC_REPORT_DIR}/pre_qc_scatter.png"
  "${PRE_QC_REPORT_DIR}/pre_qc_scatter.pdf"
  "${PRE_QC_REPORT_DIR}/pre_qc_cutoff_burden.png"
  "${PRE_QC_REPORT_DIR}/pre_qc_cutoff_burden.pdf"
)

[[ -f "${CANONICAL_SAMPLE_SHEET}" ]] || die "缺少 canonical sample sheet: ${CANONICAL_SAMPLE_SHEET}"
[[ -f "${INPUT_INVENTORY_FILE}" ]] || die "缺少 input inventory: ${INPUT_INVENTORY_FILE}"
[[ -f "${BRANCH_READINESS_FILE}" ]] || die "缺少 branch readiness: ${BRANCH_READINESS_FILE}"

run_stage_if_stale \
  "${SHELL_ROOT}/05single_script/01a_build_raw_objects.R" \
  "${MODULE_01A_MANIFEST}" \
  "${CANONICAL_SAMPLE_SHEET}" \
  "${INPUT_INVENTORY_FILE}" \
  "${CC_GENES_OUTPUT}"

RAW_OBJECT_OUTPUT="$(require_manifest_output "${MODULE_01A_MANIFEST}" "raw_object")"

run_stage_if_stale \
  "${SHELL_ROOT}/05single_script/01b_pre_qc_eda.R" \
  "${MODULE_01B_MANIFEST}" \
  "${MODULE_01A_MANIFEST}" \
  "${RAW_OBJECT_OUTPUT}" \
  "${INPUT_INVENTORY_FILE}" \
  "${BRANCH_READINESS_FILE}" \
  "${QC_THRESHOLD_FILE}"

require_manifest_output "${MODULE_01B_MANIFEST}" "sample_qc_summary" >/dev/null
require_manifest_output "${MODULE_01B_MANIFEST}" "triage" >/dev/null
require_manifest_output "${MODULE_01B_MANIFEST}" "report" >/dev/null
require_manifest_output "${MODULE_01B_MANIFEST}" "pre_qc_violin_png" >/dev/null
require_manifest_output "${MODULE_01B_MANIFEST}" "pre_qc_violin_pdf" >/dev/null
require_manifest_output "${MODULE_01B_MANIFEST}" "pre_qc_scatter_png" >/dev/null
require_manifest_output "${MODULE_01B_MANIFEST}" "pre_qc_scatter_pdf" >/dev/null
require_manifest_output "${MODULE_01B_MANIFEST}" "pre_qc_cutoff_burden_png" >/dev/null
require_manifest_output "${MODULE_01B_MANIFEST}" "pre_qc_cutoff_burden_pdf" >/dev/null

for required_output in "${PRE_QC_REQUIRED_OUTPUTS[@]}"; do
  [[ -e "${required_output}" ]] || die "01_build_raw 缺少产物: ${required_output}"
done

set_eda_gate_status \
  "pre_qc" \
  "pending" \
  "" \
  "01b_pre_qc_eda completed; manual review required before 02_qc"

update_workflow_status \
  "01_pre_qc_eda_completed" \
  "review ${EDA_GATE_FILE} and approve pre_qc before bash ${PIPELINE_ROOT}/shell/03stages/02_qc.sh" \
  "status.01a_build_raw_objects_completed=true" \
  "status.01b_pre_qc_eda_completed=true" \
  "status.01_pre_qc_eda_completed=true" \
  "status.pre_qc_eda_complete=true" \
  "status.pre_qc_gate_passed=false"
