#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

ensure_eda_control_files
ensure_dir "${TABLE_DIR}/trajectory_velocity_preflight" "${EDA_REPORT_DIR}/trajectory_velocity_preflight"

run_r_main "${WORKFLOW_ROOT}/05single_script/09_10_preflight.R"

echo "09/10 preflight report: ${EDA_REPORT_DIR}/trajectory_velocity_preflight/report.md"
