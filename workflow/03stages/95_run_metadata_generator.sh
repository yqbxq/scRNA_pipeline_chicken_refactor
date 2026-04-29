#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

ensure_dir "${RESULTS_DIR}/00_validation" "${METADATA_DIR}"

R_BIN="${R_BIN:-}"
if [[ -z "${R_BIN}" && -n "${R_MAIN_ENV_PREFIX:-}" && -x "${R_MAIN_ENV_PREFIX}/bin/Rscript" ]]; then
  R_BIN="${R_MAIN_ENV_PREFIX}/bin/Rscript"
fi
if [[ -z "${R_BIN}" ]] && command -v Rscript >/dev/null 2>&1; then
  R_BIN="$(command -v Rscript)"
fi
[[ -n "${R_BIN}" && -x "${R_BIN}" ]] || die "找不到 Rscript，无法运行 metadata generator。"

env \
  PROJECT_ROOT="${PROJECT_ROOT}" \
  PIPELINE_ROOT="${PIPELINE_ROOT}" \
  METADATA_DIR="${METADATA_DIR}" \
  RESULTS_DIR="${RESULTS_DIR}" \
  TABLE_DIR="${TABLE_DIR}" \
  "${R_BIN}" "${WORKFLOW_ROOT}/05single_script/95_build_metadata_from_questions.R"

if [[ "${RUN_METADATA_VALIDATION_AFTER_GENERATE:-yes}" == "yes" ]]; then
  bash "${WORKFLOW_ROOT}/03stages/96_validate_metadata.sh"
fi
