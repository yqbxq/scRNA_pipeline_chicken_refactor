#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

prepare_project_state_dirs
sync_workflow_gate_statuses

require_status_flag_or_warn \
  "status.main_ready" \
  "主流程尚未 ready。请先完成 metadata/intake 审计、输入标准化和 input summary。"

ensure_dir "${LOG_DIR}"

run_stage_if_stale() {
  local script_path="$1"
  local output_path="$2"
  shift 2 || true

  STAGE_RERAN="0"
  if is_stale_output "${output_path}" "$@"; then
    echo "Run stage: $(basename "${script_path}")"
    run_r_main "${script_path}"
    STAGE_RERAN="1"
  else
    echo "Skip up-to-date stage: $(basename "${script_path}")"
  fi
}

reset_gates_if_reran() {
  local note="$1"
  shift || true
  [[ "${STAGE_RERAN}" == "1" ]] || return 0

  local gate_id
  for gate_id in "$@"; do
    set_eda_gate_status "${gate_id}" "pending" "" "${note}"
  done
  sync_workflow_gate_statuses
}

hold_for_gate() {
  local gate_id="$1"
  local report_dir="$2"
  local stage_name="$3"

  sync_workflow_gate_statuses
  if eda_gate_passed "${gate_id}"; then
    echo "Gate approved: ${gate_id}"
    return 0
  fi

  local next_step
  next_step="审阅 ${report_dir}/report.md，并在 ${EDA_GATE_FILE} 中将 ${gate_id} 的 status 改为 approved，然后重跑 bash ${PIPELINE_ROOT}/workflow/20_run_main_pipeline.sh"
  update_workflow_status "${stage_name}" "${next_step}"
  echo "Gate hold: ${gate_id}"
  echo "Review report: ${report_dir}/report.md"
  echo "Then update gate file: ${EDA_GATE_FILE}"
  exit 0
}

RAW_OUTPUT="${CHECKPOINT_DIR}/00_raw_objects.rds"
PRE_QC_REPORT="${PRE_QC_REPORT_DIR}/report.md"
AMBIENT_OUTPUT="${CHECKPOINT_DIR}/01_after_ambient_branch.rds"
AMBIENT_REPORT="${AMBIENT_REPORT_DIR}/report.md"
QC_OUTPUT="${CHECKPOINT_DIR}/01_after_qc_doublet.rds"
POST_QC_REPORT="${POST_QC_REPORT_DIR}/report.md"
REDUCTION_OUTPUT="${CHECKPOINT_DIR}/02_reduction_candidates.rds"
INTEGRATION_REPORT="${INTEGRATION_REPORT_DIR}/report.md"
CLUSTER_OUTPUT="${CHECKPOINT_DIR}/02_after_clustering.rds"
ANNOTATION_OUTPUT="${ANNOTATION_HUB_PATH}"
ANNOTATION_REPORT="${ANNOTATION_REPORT_DIR}/report.md"
MARKER_DISCOVERY_OUTPUT="${TABLE_DIR}/marker_discovery/marker_discovery_manifest.tsv"
PSEUDOBULK_OUTPUT="${TABLE_DIR}/pseudobulk_ds/pseudobulk_manifest.tsv"
COMPOSITION_OUTPUT="${TABLE_DIR}/composition/composition_manifest.tsv"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/01_build_raw_objects.R" \
  "${RAW_OUTPUT}" \
  "${CONFIG_FILE}" \
  "${SAMPLE_SHEET}" \
  "${CANONICAL_SAMPLE_SHEET}" \
  "${INPUT_INVENTORY_FILE}"
reset_gates_if_reran "raw_objects_reran" pre_qc post_qc integration annotation
update_workflow_status \
  "raw_objects_built" \
  "bash ${PIPELINE_ROOT}/workflow/20_run_main_pipeline.sh"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/01a_pre_qc_eda.R" \
  "${PRE_QC_REPORT}" \
  "${RAW_OUTPUT}" \
  "${CONFIG_FILE}" \
  "${INPUT_INVENTORY_FILE}" \
  "${BRANCH_READINESS_FILE}" \
  "${QC_THRESHOLD_FILE}"
reset_gates_if_reran "pre_qc_eda_reran" pre_qc post_qc integration annotation
update_workflow_status \
  "pre_qc_eda_completed" \
  "审阅 pre-QC EDA 报告并更新 gate 后重跑 20_run_main_pipeline.sh" \
  "status.pre_qc_eda_complete=true"
hold_for_gate "pre_qc" "${PRE_QC_REPORT_DIR}" "waiting_for_pre_qc_gate"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/01ab_ambient_branch.R" \
  "${AMBIENT_OUTPUT}" \
  "${RAW_OUTPUT}" \
  "${CONFIG_FILE}" \
  "${INPUT_INVENTORY_FILE}" \
  "${BRANCH_READINESS_FILE}" \
  "${QC_THRESHOLD_FILE}"
reset_gates_if_reran "ambient_branch_reran" post_qc integration annotation
update_workflow_status \
  "ambient_branch_completed" \
  "bash ${PIPELINE_ROOT}/workflow/20_run_main_pipeline.sh" \
  "status.ambient_branch_complete=true"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/01b_qc_doublet.R" \
  "${QC_OUTPUT}" \
  "${AMBIENT_OUTPUT}" \
  "${CONFIG_FILE}" \
  "${QC_THRESHOLD_FILE}"
reset_gates_if_reran "qc_doublet_reran" post_qc integration annotation
update_workflow_status \
  "qc_doublet_completed" \
  "bash ${PIPELINE_ROOT}/workflow/20_run_main_pipeline.sh"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/01c_post_qc_eda.R" \
  "${POST_QC_REPORT}" \
  "${RAW_OUTPUT}" \
  "${QC_OUTPUT}" \
  "${CONFIG_FILE}" \
  "${TABLE_DIR}/qc_doublet_summary.csv" \
  "${TABLE_DIR}/qc_doublet_cell_level.csv.gz"
reset_gates_if_reran "post_qc_eda_reran" post_qc integration annotation
update_workflow_status \
  "post_qc_eda_completed" \
  "审阅 post-QC EDA 报告并更新 gate 后重跑 20_run_main_pipeline.sh" \
  "status.post_qc_eda_complete=true"
hold_for_gate "post_qc" "${POST_QC_REPORT_DIR}" "waiting_for_post_qc_gate"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/02_build_reductions.R" \
  "${REDUCTION_OUTPUT}" \
  "${QC_OUTPUT}" \
  "${CONFIG_FILE}" \
  "${OBJECT_LAYER_CONFIG_FILE}"
reset_gates_if_reran "reduction_candidates_reran" integration annotation
update_workflow_status \
  "reduction_candidates_built" \
  "bash ${PIPELINE_ROOT}/workflow/20_run_main_pipeline.sh"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/02a_integration_eda.R" \
  "${INTEGRATION_REPORT}" \
  "${REDUCTION_OUTPUT}" \
  "${CONFIG_FILE}" \
  "${OBJECT_LAYER_CONFIG_FILE}"
reset_gates_if_reran "integration_eda_reran" integration annotation
update_workflow_status \
  "integration_eda_completed" \
  "审阅 integration EDA 报告并更新 gate 后重跑 20_run_main_pipeline.sh" \
  "status.integration_eda_complete=true"
hold_for_gate "integration" "${INTEGRATION_REPORT_DIR}" "waiting_for_integration_gate"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/02b_finalize_clustering.R" \
  "${CLUSTER_OUTPUT}" \
  "${REDUCTION_OUTPUT}" \
  "${CONFIG_FILE}" \
  "${INTEGRATION_REPORT}" \
  "${OBJECT_LAYER_CONFIG_FILE}"
reset_gates_if_reran "clustering_reran" annotation
update_workflow_status \
  "clustering_finalized" \
  "bash ${PIPELINE_ROOT}/workflow/20_run_main_pipeline.sh"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/03_annotation.R" \
  "${ANNOTATION_OUTPUT}" \
  "${CLUSTER_OUTPUT}" \
  "${CONFIG_FILE}"
reset_gates_if_reran "annotation_reran" annotation
update_workflow_status \
  "annotation_completed" \
  "bash ${PIPELINE_ROOT}/workflow/20_run_main_pipeline.sh"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/03a_annotation_eda.R" \
  "${ANNOTATION_REPORT}" \
  "${ANNOTATION_OUTPUT}" \
  "${CONFIG_FILE}" \
  "${BRANCH_READINESS_FILE}"
reset_gates_if_reran "annotation_eda_reran" annotation
update_workflow_status \
  "annotation_eda_completed" \
  "审阅 annotation EDA 报告并更新 gate 后重跑 20_run_main_pipeline.sh" \
  "status.annotation_eda_complete=true"
hold_for_gate "annotation" "${ANNOTATION_REPORT_DIR}" "waiting_for_annotation_gate"

sync_workflow_gate_statuses

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/04_marker_discovery.R" \
  "${MARKER_DISCOVERY_OUTPUT}" \
  "${ANNOTATION_OUTPUT}" \
  "${CONFIG_FILE}"
update_workflow_status \
  "marker_discovery_completed" \
  "bash ${PIPELINE_ROOT}/workflow/20_run_main_pipeline.sh"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/04a_pseudobulk_ds.R" \
  "${PSEUDOBULK_OUTPUT}" \
  "${ANNOTATION_OUTPUT}" \
  "${CONFIG_FILE}"
update_workflow_status \
  "pseudobulk_ds_completed" \
  "bash ${PIPELINE_ROOT}/workflow/20_run_main_pipeline.sh"

run_stage_if_stale \
  "${PIPELINE_ROOT}/workflow/r/04b_composition.R" \
  "${COMPOSITION_OUTPUT}" \
  "${ANNOTATION_OUTPUT}" \
  "${CONFIG_FILE}"
update_workflow_status \
  "composition_completed" \
  "bash ${PIPELINE_ROOT}/workflow/20_run_main_pipeline.sh"

[[ -f "${ANNOTATION_HUB_PATH}" ]] || die "主流程完成后缺少 annotation hub: ${ANNOTATION_HUB_PATH}"

update_workflow_status \
  "main_pipeline_completed" \
  "审阅 ${ANNOTATION_REPORT_DIR}/report.md、${EDA_REPORT_DIR}/pseudobulk_ds/report.md 和 ${EDA_REPORT_DIR}/composition/report.md；trajectory/velocity 等后续模块需单独运行。"

echo "主流程运行完成"
