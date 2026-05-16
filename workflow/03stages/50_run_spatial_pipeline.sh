#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

ensure_metadata_fresh
check_stage_deps "50_spatial"
prepare_project_state_dirs

SPATIAL_R_DIR="${WORKFLOW_ROOT}/05single_script/spatial"

require_spatial_intake_ready() {
  [[ -s "${SPATIAL_INPUT_INVENTORY_FILE}" ]] || die "缺少 ST input inventory，请先运行 audit_inputs.sh。"
  if ! awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; i++) idx[$i] = i; next } idx["ready"] && $idx["ready"] != "true" { bad = 1 } END { exit bad ? 1 : 0 }' "${SPATIAL_INPUT_INVENTORY_FILE}"; then
    die "ST 输入尚未全部 ready，请查看 ${SPATIAL_INPUT_INVENTORY_FILE}。"
  fi
}

run_future_spatial_r_stage() {
  local script_name="$1"
  local manifest_name="$2"
  shift 2 || true
  local script_path="${SPATIAL_R_DIR}/${script_name}"
  local manifest_path="${MANIFEST_DIR}/${manifest_name}/_manifest.json"
  [[ -f "${script_path}" ]] || die "R stage 尚未实现: ${script_path}。"
  run_stage_if_stale_with_runner run_r_spatial "${script_path}" "${manifest_path}" "$@"
}

require_spatial_intake_ready

run_future_spatial_r_stage "01_build_spatial_objects.R" "spatial_01_build_objects" "${SPATIAL_INPUT_INVENTORY_FILE}" "${SECTION_SHEET}"
run_future_spatial_r_stage "01a_pre_spot_qc_eda.R" "spatial_01a_pre_qc_eda" "${SPATIAL_QC_THRESHOLD_FILE}"
set_eda_gate_status "spatial_pre_qc" "pending" "" "Review spatial pre-QC report before spatial filtering."
hold_for_gate spatial_pre_qc

run_future_spatial_r_stage "01b_spot_qc_filter.R" "spatial_01b_qc_filter" "${SPATIAL_QC_THRESHOLD_FILE}"
run_future_spatial_r_stage "01c_post_spot_qc_eda.R" "spatial_01c_post_qc_eda"
set_eda_gate_status "spatial_post_qc" "pending" "" "Review spatial post-QC report before spatial normalization."
hold_for_gate spatial_post_qc

run_future_spatial_r_stage "02_normalize_spatial.R" "spatial_02_normalize" "${SPATIAL_OBJECT_LAYER_FILE}"
run_future_spatial_r_stage "02a_spatial_integration_eda.R" "spatial_02a_integration_eda"
set_eda_gate_status "spatial_integration" "pending" "" "Review spatial integration diagnostics before clustering."
hold_for_gate spatial_integration

run_future_spatial_r_stage "02b_finalize_spatial_clustering.R" "spatial_02b_finalize_clustering" "${SPATIAL_OBJECT_LAYER_FILE}"
run_future_spatial_r_stage "03_region_annotation.R" "spatial_03_region_annotation" "${MARKER_PANEL_DIR}"
run_future_spatial_r_stage "03a_region_annotation_eda.R" "spatial_03a_region_annotation_eda"
set_eda_gate_status "spatial_region_annotation" "pending" "" "Review spatial region annotation before marker/DE/composition."
hold_for_gate spatial_region_annotation

run_future_spatial_r_stage "04_region_marker_discovery.R" "spatial_04_region_marker"
run_future_spatial_r_stage "04a_region_pseudobulk.R" "spatial_04a_region_pseudobulk"
run_future_spatial_r_stage "04b_region_composition.R" "spatial_04b_region_composition"
run_future_spatial_r_stage "04_spatial_de.R" "spatial_04_spatial_de"

update_workflow_status \
  "spatial_main_completed" \
  "bash ${PIPELINE_ROOT}/workflow/03stages/60_run_spatial_extensions.sh" \
  "status.spatial_main_completed=true"
