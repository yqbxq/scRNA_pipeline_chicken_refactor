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
run_future_spatial_r_stage "02a2_compute_umap_spatial.R" "spatial_02a2_compute_umap" "${MANIFEST_DIR}/spatial_02a_integration_eda/_manifest.json"
set_eda_gate_status "spatial_integration" "pending" "" "Review spatial integration diagnostics before clustering."
hold_for_gate spatial_integration

run_future_spatial_r_stage "02b_finalize_spatial_clustering.R" "spatial_02b_finalize_clustering" "${SPATIAL_OBJECT_LAYER_FILE}"
run_future_spatial_r_stage "03_region_annotation.R" "spatial_03_region_annotation" "${MARKER_PANEL_DIR}"
export_h5ad_for_gate \
  "spatial_03_region" \
  "${MANIFEST_DIR}/spatial_03_region_annotation/_manifest.json" \
  "panorama_annotated" \
  "spatial" \
  "spatial_03_region"
run_future_spatial_r_stage "03a_region_annotation_eda.R" "spatial_03a_region_annotation_eda"
set_eda_gate_status "spatial_region_annotation" "pending" "" "Review spatial region annotation before sub-clustering."
hold_for_gate spatial_region_annotation

run_future_spatial_r_stage "04a_subcluster_build.R" "spatial_04a_subcluster_build" "${SPATIAL_OBJECT_LAYER_FILE}"
run_future_spatial_r_stage "04b_subcluster_annotate.R" "spatial_04b_subcluster_annotate" "${MARKER_PANEL_DIR}"
export_h5ad_for_gate \
  "spatial_04b_subcluster" \
  "${MANIFEST_DIR}/spatial_04b_subcluster_annotate/_manifest.json" \
  "panorama_subannotated" \
  "spatial" \
  "spatial_04b_subcluster"
run_future_spatial_r_stage "04c_subcluster_eda.R" "spatial_04c_subcluster_eda"
set_eda_gate_status "spatial_region_annotation" "pending" "" "Review sub-cluster annotation before marker/DE/composition."
hold_for_gate spatial_region_annotation

run_future_spatial_r_stage "05_region_marker_discovery.R" "spatial_05_region_marker"
export_h5ad_for_gate \
  "spatial_05_marker" \
  "${MANIFEST_DIR}/spatial_03_region_annotation/_manifest.json" \
  "panorama_annotated" \
  "spatial" \
  "spatial_05_marker"
run_future_spatial_r_stage "05a_region_pseudobulk.R" "spatial_05a_region_pseudobulk"
run_future_spatial_r_stage "05b_region_composition.R" "spatial_05b_region_composition"
run_future_spatial_r_stage "05_spatial_de.R" "spatial_05_spatial_de"
set_eda_gate_status "spatial_marker_de" "pending" "" "Review spatial region marker / DE / composition outputs before enrichment and deconvolution."
hold_for_gate spatial_marker_de

update_workflow_status \
  "spatial_main_completed" \
  "bash ${PIPELINE_ROOT}/workflow/03stages/60_run_spatial_extensions.sh" \
  "status.spatial_main_completed=true"
