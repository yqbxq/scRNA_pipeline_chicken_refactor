#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

ensure_metadata_fresh
check_stage_deps "04_robustness"

MODULE_04D_MANIFEST="${MANIFEST_DIR}/04d_cluster_robustness/_manifest.json"
MODULE_04E_MANIFEST="${MANIFEST_DIR}/04e_scdesign3_engine/_manifest.json"
MODULE_04F_MANIFEST="${MANIFEST_DIR}/04f_scdesign3_finalize/_manifest.json"

ensure_eda_control_files

require_manifest_output "${MODULE_04D_MANIFEST}" "target_gate_status_tsv" >/dev/null
require_manifest_output "${MODULE_04D_MANIFEST}" "scdesign3_all_questions_status_tsv" >/dev/null

hold_for_gate scdesign3_targets

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/04e_scdesign3_engine.R" \
  "${MODULE_04E_MANIFEST}" \
  "${SCDESIGN3_TARGETS_SHEET}" \
  "${SCDESIGN3_THRESHOLDS_SHEET}" \
  "${MODULE_04D_MANIFEST}"

require_manifest_output "${MODULE_04E_MANIFEST}" "target_metrics_tsv" >/dev/null
require_manifest_output "${MODULE_04E_MANIFEST}" "engine_status_tsv" >/dev/null

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/04f_scdesign3_finalize.R" \
  "${MODULE_04F_MANIFEST}" \
  "${MODULE_04E_MANIFEST}" \
  "${MODULE_04D_MANIFEST}"

require_manifest_output "${MODULE_04F_MANIFEST}" "scdesign3_all_questions_status_tsv" >/dev/null

set_eda_gate_status \
  "scdesign3_validated" \
  "pending" \
  "" \
  "review 04e ARI/NMI metrics and figures before allowing downstream core interpretation"
sync_workflow_gate_statuses

update_workflow_status \
  "04_robustness_completed" \
  "review scdesign3_validated gate; approve only after 04e/04f scDesign3 metrics are accepted" \
  "status.04e_scdesign3_engine_completed=true" \
  "status.04f_scdesign3_finalize_completed=true" \
  "status.04_robustness_completed=true" \
  "status.scdesign3_validated_gate_passed=$(eda_gate_passed scdesign3_validated && echo true || echo false)"
