#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

RUN_SVG="yes"
RUN_DECONV="yes"
RUN_DECONV_EXTRA="no"
RUN_NEIGHBORHOOD="yes"
RUN_DECOUPLER="yes"
RUN_PYSCENIC="no"
RUN_ENRICHMENT="yes"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-svg) RUN_SVG="yes"; shift ;;
    --without-svg) RUN_SVG="no"; shift ;;
    --with-deconv) RUN_DECONV="yes"; shift ;;
    --without-deconv) RUN_DECONV="no"; shift ;;
    --with-deconv-extra) RUN_DECONV_EXTRA="yes"; shift ;;
    --with-neighborhood) RUN_NEIGHBORHOOD="yes"; shift ;;
    --without-neighborhood) RUN_NEIGHBORHOOD="no"; shift ;;
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
  --with-deconv / --without-deconv / --with-deconv-extra
  --with-neighborhood / --without-neighborhood
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

[[ "${RUN_SVG}" == "yes" ]] && run_future_spatial_r_stage "05_svg.R" "spatial_05_svg"
if [[ "${RUN_DECONV}" == "yes" ]]; then
  run_future_spatial_r_stage "07_deconvolution.R" "spatial_07_deconvolution" "${SPATIAL_REFERENCE_INVENTORY_FILE}"
  run_future_spatial_r_stage "07a_deconvolution_validation.R" "spatial_07a_deconvolution_validation"
  set_eda_gate_status "spatial_deconv" "pending" "" "Review deconvolution and scDesign3 validation outputs."
  hold_for_gate spatial_deconv
fi
[[ "${RUN_NEIGHBORHOOD}" == "yes" ]] && run_future_spatial_r_stage "06_spatial_neighborhood.R" "spatial_06_neighborhood"
[[ "${RUN_DECOUPLER}" == "yes" ]] && run_future_spatial_r_stage "06b_spatial_decoupler.R" "spatial_06b_decoupler"
[[ "${RUN_PYSCENIC}" == "yes" ]] && run_future_spatial_r_stage "06c_spatial_pyscenic.R" "spatial_06c_pyscenic"
[[ "${RUN_ENRICHMENT}" == "yes" ]] && run_future_spatial_r_stage "08_spatial_enrichment.R" "spatial_08_enrichment"

update_workflow_status \
  "spatial_extensions_completed" \
  "bash ${PIPELINE_ROOT}/workflow/03stages/70_run_joint_scrna_spatial.sh" \
  "status.spatial_extensions_completed=true"
