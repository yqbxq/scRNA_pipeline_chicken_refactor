#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

RUN_SVG="yes"
RUN_DECONV="yes"
RUN_DECONV_EXTRA="yes"
RUN_DECONV_VALIDATION="yes"
RUN_NEIGHBORHOOD="yes"
RUN_NICHE="yes"
RUN_COMMOT="yes"
RUN_DECOUPLER="no"
RUN_PYSCENIC="no"
RUN_ENRICHMENT="yes"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-svg) RUN_SVG="yes"; shift ;;
    --without-svg) RUN_SVG="no"; shift ;;
    --with-deconv) RUN_DECONV="yes"; shift ;;
    --without-deconv) RUN_DECONV="no"; shift ;;
    --with-deconv-extra) RUN_DECONV_EXTRA="yes"; shift ;;
    --without-deconv-extra) RUN_DECONV_EXTRA="no"; shift ;;
    --with-deconv-validation) RUN_DECONV_VALIDATION="yes"; shift ;;
    --without-deconv-validation) RUN_DECONV_VALIDATION="no"; shift ;;
    --with-neighborhood) RUN_NEIGHBORHOOD="yes"; shift ;;
    --without-neighborhood) RUN_NEIGHBORHOOD="no"; shift ;;
    --with-niche) RUN_NICHE="yes"; shift ;;
    --without-niche) RUN_NICHE="no"; shift ;;
    --with-commot) RUN_COMMOT="yes"; shift ;;
    --without-commot) RUN_COMMOT="no"; shift ;;
    --with-decoupler) RUN_DECOUPLER="yes"; shift ;;
    --without-decoupler) RUN_DECOUPLER="no"; shift ;;
    --with-pyscenic) RUN_PYSCENIC="yes"; shift ;;
    --with-enrichment) RUN_ENRICHMENT="yes"; shift ;;
    --without-enrichment) RUN_ENRICHMENT="no"; shift ;;
    -h|--help)
      cat <<'USAGE'
Usage: 60_run_spatial_extensions.sh [flags]

Flags:
  --with-svg / --without-svg
  --with-deconv / --without-deconv / --with-deconv-extra / --without-deconv-extra
  --with-deconv-validation / --without-deconv-validation
  --with-neighborhood / --without-neighborhood
  --with-niche / --without-niche
  --with-commot / --without-commot
  --with-decoupler / --without-decoupler
  --with-pyscenic
  --with-enrichment / --without-enrichment
USAGE
      exit 0
      ;;
    *) die "Unknown argument: $1" ;;
  esac
done

check_stage_deps "60_spatial_extensions"
prepare_project_state_dirs

SPATIAL_R_DIR="${WORKFLOW_ROOT}/05single_script/spatial"

run_future_spatial_r_stage() {
  local script_name="$1"
  local manifest_name="$2"
  shift 2 || true
  local script_path="${SPATIAL_R_DIR}/${script_name}"
  local manifest_path="${MANIFEST_DIR}/${manifest_name}/_manifest.json"
  [[ -f "${script_path}" ]] || die "R stage 尚未实现: ${script_path}。本轮只完成非 R 脚手架。"
  run_stage_if_stale_with_runner run_r_spatial "${script_path}" "${manifest_path}" "$@"
}

run_spatial_python_script() {
  local script_path="$1"
  local python_bin="${PY_SPATIAL_BIN:-$(detect_python)}"
  [[ -x "${python_bin}" || "$(command -v "${python_bin}" 2>/dev/null || true)" ]] || python_bin="$(detect_python)"
  "${python_bin}" "${script_path}"
}

if [[ "${RUN_SVG}" == "yes" ]]; then
  run_future_spatial_r_stage "09a_spatialde2_svg.R" "spatial_09a_spatialde2_svg"
  run_future_spatial_r_stage "09b_sparkx_svg.R" "spatial_09b_sparkx_svg"
  run_future_spatial_r_stage "09c_svg_consensus_report.R" "spatial_09c_svg_consensus"
  if [[ "${SVG_REQUIRE_REVIEW:-no}" == "yes" ]]; then
    set_eda_gate_status "spatial_svg" "pending" "" "Review exploratory SVG gene evidence before downstream interpretation."
    hold_for_gate spatial_svg
  fi
fi
if [[ "${RUN_ENRICHMENT}" == "yes" ]]; then
  run_future_spatial_r_stage "06a_region_go_enrichment.R" "spatial_06a_region_go" "${COMPARISON_SHEET}"
  run_future_spatial_r_stage "06b_region_kegg_enrichment.R" "spatial_06b_region_kegg" "${COMPARISON_SHEET}"
  run_future_spatial_r_stage "06c_region_enrichment_eda.R" "spatial_06c_region_enrichment_eda"
  set_eda_gate_status "spatial_enrichment" "pending" "" "Review spatial region GO/KEGG enrichment before neighborhood analysis."
  hold_for_gate spatial_enrichment
fi
if [[ "${RUN_NEIGHBORHOOD}" == "yes" ]]; then
  run_future_spatial_r_stage "06d_spatial_neighborhood.R" "spatial_06d_neighborhood"
  set_eda_gate_status "spatial_neighborhood" "pending" "" "Review spatial neighborhood enrichment before niche derivation."
  hold_for_gate spatial_neighborhood
fi
if [[ "${RUN_DECONV}" == "yes" ]]; then
  run_future_spatial_r_stage "07a_deconvolution_rctd.R" "spatial_07a_deconvolution_rctd" \
    "${SPATIAL_REFERENCE_INVENTORY_FILE}" "${DECONV_PAIRS_SHEET}"
  if [[ "${RUN_DECONV_EXTRA}" == "yes" ]]; then
    run_future_spatial_r_stage "07b_deconvolution_transfer.R" "spatial_07b_deconvolution_transfer" \
      "${SPATIAL_REFERENCE_INVENTORY_FILE}" "${DECONV_PAIRS_SHEET}"
    run_future_spatial_r_stage "07c_deconvolution_card.R" "spatial_07c_deconvolution_card" \
      "${SPATIAL_REFERENCE_INVENTORY_FILE}" "${DECONV_PAIRS_SHEET}"
    run_future_spatial_r_stage "07d_deconvolution_cell2location.R" "spatial_07d_deconvolution_cell2location" \
      "${SPATIAL_REFERENCE_INVENTORY_FILE}" "${DECONV_PAIRS_SHEET}"
  fi
  run_future_spatial_r_stage "07e_deconvolution_compare.R" "spatial_07e_deconvolution_compare"
  if [[ "${RUN_DECONV_VALIDATION}" == "yes" ]]; then
    run_future_spatial_r_stage "07f_deconvolution_validation.R" "spatial_07f_deconvolution_validation" \
      "${SPATIAL_REFERENCE_INVENTORY_FILE}"
  fi
  set_eda_gate_status "spatial_deconv" "pending" "" "Review deconvolution outputs, multi-method comparison, and optional scDesign3 validation."
  hold_for_gate spatial_deconv
fi
if [[ "${RUN_COMMOT}" == "yes" && "${COMMOT_ENABLED:-auto}" != "no" ]]; then
  run_future_spatial_r_stage "08a_spatial_communication_io.R" "spatial_08a_spatial_communication_io" \
    "${MANIFEST_DIR}/spatial_07e_deconvolution_compare/_manifest.json"
  export COMMOT_INPUT_MANIFEST="${SPATIAL_TABLE_DIR}/08_commot/commot_input_manifest.tsv"
  export COMMOT_LR_CANDIDATES_PREPARED_TSV="${SPATIAL_TABLE_DIR}/08_commot/commot_lr_candidates.tsv"
  export COMMOT_OUTPUT_TSV="${SPATIAL_TABLE_DIR}/08_commot/commot_lr.tsv"
  export COMMOT_MANIFEST="${MANIFEST_DIR}/spatial_08f_commot_spatial/_manifest.json"
  run_stage_if_stale_with_runner \
    run_spatial_python_script \
    --script "${WORKFLOW_ROOT}/04python/07f_commot_spatial.py" \
    --manifest "${COMMOT_MANIFEST}" \
    --inputs "${COMMOT_INPUT_MANIFEST}" "${COMMOT_LR_CANDIDATES_PREPARED_TSV}" \
    --params COMMOT_DISTANCE_THRESHOLD,COMMOT_PERMUTATIONS,COMMOT_REQUIRE_RUNTIME,RANDOM_SEED
  run_future_spatial_r_stage "08b_spatial_communication_eda.R" "spatial_08b_spatial_communication_eda" \
    "${COMMOT_MANIFEST}"
  export COMMOT_SPATIAL_SUMMARY_TSV="${SPATIAL_TABLE_DIR}/08_commot/commot_spatial_summary.tsv"
  run_future_spatial_r_stage "08c_lr_colocalization.R" "spatial_08c_lr_colocalization" \
    "${MANIFEST_DIR}/spatial_08a_spatial_communication_io/_manifest.json" \
    "${MANIFEST_DIR}/spatial_07e_deconvolution_compare/_manifest.json"
  run_future_spatial_r_stage "08d_communication_neighborhood_consistency.R" "spatial_08d_communication_neighborhood" \
    "${MANIFEST_DIR}/spatial_08c_lr_colocalization/_manifest.json" \
    "${MANIFEST_DIR}/spatial_06d_neighborhood/_manifest.json"
  run_future_spatial_r_stage "08e_spatial_communication_consensus.R" "spatial_08e_spatial_communication_consensus" \
    "${MANIFEST_DIR}/spatial_08b_spatial_communication_eda/_manifest.json" \
    "${MANIFEST_DIR}/spatial_08c_lr_colocalization/_manifest.json" \
    "${MANIFEST_DIR}/spatial_08d_communication_neighborhood/_manifest.json"
  run_future_spatial_r_stage "08f_spatial_communication_report.R" "spatial_08f_spatial_communication_report" \
    "${MANIFEST_DIR}/spatial_08e_spatial_communication_consensus/_manifest.json"
fi
if [[ "${RUN_NICHE}" == "yes" ]]; then
  run_future_spatial_r_stage "06e_niche_derivation.R" "spatial_06e_niche"
  export_h5ad_for_gate \
    "spatial_06_extensions" \
    "${MANIFEST_DIR}/spatial_06e_niche/_manifest.json" \
    "panorama_niched" \
    "spatial" \
    "spatial_06_extensions"
  set_eda_gate_status "spatial_niche" "pending" "" "Review niche derivation before joint analyses."
  hold_for_gate spatial_niche
fi
[[ "${RUN_DECOUPLER}" == "yes" ]] && run_future_spatial_r_stage "08a_spatial_decoupler.R" "spatial_08a_decoupler"
[[ "${RUN_PYSCENIC}" == "yes" ]] && run_future_spatial_r_stage "08b_spatial_pyscenic.R" "spatial_08b_pyscenic"

update_workflow_status \
  "spatial_extensions_completed" \
  "bash ${PIPELINE_ROOT}/workflow/03stages/70_run_joint_scrna_spatial.sh" \
  "status.spatial_extensions_completed=true"
