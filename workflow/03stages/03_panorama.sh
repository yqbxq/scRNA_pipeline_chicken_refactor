#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

check_stage_deps "03_panorama"
hold_for_gate post_qc

MODULE_02B2_MANIFEST="${MANIFEST_DIR}/02b2_doublet/_manifest.json"
MODULE_02C_MANIFEST="${MANIFEST_DIR}/02c_post_qc_eda/_manifest.json"
POST_QC_OBJECT="$(require_manifest_output "${MODULE_02B2_MANIFEST}" "post_qc_object")"
require_manifest_output "${MODULE_02C_MANIFEST}" "report" >/dev/null

MODULE_03A1_MANIFEST="${MANIFEST_DIR}/03a1_normalize_hvg/_manifest.json"
MODULE_03A2_MANIFEST="${MANIFEST_DIR}/03a2_reduce_integrate/_manifest.json"
MODULE_03B_MANIFEST="${MANIFEST_DIR}/03b_integration_eda/_manifest.json"
MODULE_03C_MANIFEST="${MANIFEST_DIR}/03c_cluster/_manifest.json"
MODULE_03C2_MANIFEST="${MANIFEST_DIR}/03c2_cluster_marker_audit/_manifest.json"
MODULE_03C3_MANIFEST="${MANIFEST_DIR}/03c3_cluster_selection_report/_manifest.json"
MODULE_03D_MANIFEST="${MANIFEST_DIR}/03d_marker_risk/_manifest.json"
MODULE_03E_MANIFEST="${MANIFEST_DIR}/03e_panel_evidence/_manifest.json"
MODULE_03F_MANIFEST="${MANIFEST_DIR}/03f_apply_manual_annotation/_manifest.json"
MODULE_03G_MANIFEST="${MANIFEST_DIR}/03g_annotation_eda/_manifest.json"

ensure_eda_control_files
ensure_object_layer_config_file
ensure_marker_panel_dir

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/03a1_normalize_hvg.R" \
  "${MODULE_03A1_MANIFEST}" \
  "${MODULE_02B2_MANIFEST}" \
  "${MODULE_02C_MANIFEST}" \
  "${POST_QC_OBJECT}" \
  "${OBJECT_LAYER_CONFIG_FILE}"

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/03a2_reduce_integrate.R" \
  "${MODULE_03A2_MANIFEST}" \
  "${MODULE_03A1_MANIFEST}" \
  "${OBJECT_LAYER_CONFIG_FILE}"

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/03b_integration_eda.R" \
  "${MODULE_03B_MANIFEST}" \
  "${MODULE_03A2_MANIFEST}"

require_manifest_output "${MODULE_03B_MANIFEST}" "selected_integration_txt" >/dev/null
if ! eda_gate_passed integration; then
  set_eda_gate_status \
    "integration" \
    "pending" \
    "" \
    "03b integration EDA completed; review ${SELECTED_INTEGRATION_FILE}, then approve integration"
fi

update_workflow_status \
  "03b_integration_eda_completed" \
  "review ${SELECTED_INTEGRATION_FILE}, approve integration in ${EDA_GATE_FILE}, then rerun 03_panorama" \
  "status.03a1_normalize_hvg_completed=true" \
  "status.03a2_reduce_integrate_completed=true" \
  "status.03b_integration_eda_completed=true" \
  "status.integration_eda_complete=true" \
  "status.integration_gate_passed=$(eda_gate_passed integration && echo true || echo false)"

hold_for_gate integration

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/03c_cluster.R" \
  "${MODULE_03C_MANIFEST}" \
  "${MODULE_03A2_MANIFEST}" \
  "${MODULE_03B_MANIFEST}" \
  "${SELECTED_INTEGRATION_FILE}" \
  "${OBJECT_LAYER_CONFIG_FILE}"

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/03c2_cluster_marker_audit.R" \
  "${MODULE_03C2_MANIFEST}" \
  "${MODULE_03C_MANIFEST}"

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/03c3_cluster_selection_report.R" \
  "${MODULE_03C3_MANIFEST}" \
  "${MODULE_03C_MANIFEST}" \
  "${MODULE_03C2_MANIFEST}" \
  "${CLUSTER_SELECTION_FILE}"

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/03d_marker_risk.R" \
  "${MODULE_03D_MANIFEST}" \
  "${MODULE_03C3_MANIFEST}"

require_manifest_output "${MODULE_03C3_MANIFEST}" "report" >/dev/null
require_manifest_output "${MODULE_03D_MANIFEST}" "report" >/dev/null
if ! eda_gate_passed clustering; then
  set_eda_gate_status \
    "clustering" \
    "pending" \
    "" \
    "03c3 clustering selection and 03d marker risk audit completed; review clustering and marker risk reports, then approve clustering"
fi

update_workflow_status \
  "03d_marker_risk_completed" \
  "review clustering and marker risk reports, approve clustering in ${EDA_GATE_FILE}, then rerun 03_panorama" \
  "status.03a1_normalize_hvg_completed=true" \
  "status.03a2_reduce_integrate_completed=true" \
  "status.03b_integration_eda_completed=true" \
  "status.03c_cluster_completed=true" \
  "status.03c2_cluster_marker_audit_completed=true" \
  "status.03c3_cluster_selection_completed=true" \
  "status.03d_marker_risk_completed=true" \
  "status.integration_eda_complete=true" \
  "status.integration_gate_passed=$(eda_gate_passed integration && echo true || echo false)" \
  "status.clustering_gate_passed=$(eda_gate_passed clustering && echo true || echo false)"

hold_for_gate clustering

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/03e_panel_evidence.R" \
  "${MODULE_03E_MANIFEST}" \
  "${MODULE_03C3_MANIFEST}" \
  "${MODULE_03D_MANIFEST}" \
  "${MARKER_PANEL_DIR}" \
  "${MANUAL_ANNOTATION_FILE}"

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/03f_apply_manual_annotation.R" \
  "${MODULE_03F_MANIFEST}" \
  "${MODULE_03C3_MANIFEST}" \
  "${MODULE_03D_MANIFEST}" \
  "${MODULE_03E_MANIFEST}" \
  "${MANUAL_ANNOTATION_FILE}"

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/03g_annotation_eda.R" \
  "${MODULE_03G_MANIFEST}" \
  "${MODULE_03F_MANIFEST}"

require_manifest_output "${MODULE_03C3_MANIFEST}" "clustered_object" >/dev/null
require_manifest_output "${MODULE_03F_MANIFEST}" "annotated_object" >/dev/null
require_manifest_output "${MODULE_03G_MANIFEST}" "report" >/dev/null

if ! eda_gate_passed annotation; then
  set_eda_gate_status \
    "annotation" \
    "pending" \
    "" \
    "03g annotation EDA completed; manual review required before 04_subcluster"
fi

update_workflow_status \
  "03_panorama_completed" \
  "review ${EDA_GATE_FILE} and approve annotation before bash ${PIPELINE_ROOT}/workflow/03stages/04_subcluster.sh" \
  "status.03a1_normalize_hvg_completed=true" \
  "status.03a2_reduce_integrate_completed=true" \
  "status.03b_integration_eda_completed=true" \
  "status.03c_cluster_completed=true" \
  "status.03c2_cluster_marker_audit_completed=true" \
  "status.03c3_cluster_selection_completed=true" \
  "status.03d_marker_risk_completed=true" \
  "status.03e_panel_evidence_completed=true" \
  "status.03f_manual_annotation_applied=true" \
  "status.03g_annotation_eda_completed=true" \
  "status.03d_annotate_completed=true" \
  "status.03e_annotation_eda_completed=true" \
  "status.03_panorama_completed=true" \
  "status.clustering_gate_passed=$(eda_gate_passed clustering && echo true || echo false)" \
  "status.annotation_gate_passed=$(eda_gate_passed annotation && echo true || echo false)"
