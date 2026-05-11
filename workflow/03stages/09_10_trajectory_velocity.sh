#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

run_stage_09_10() {
  local stage_script="$1"
  echo "Running ${stage_script}"
  PIPELINE_ROOT="${PIPELINE_ROOT}" bash "${STAGE_DIR}/${stage_script}"
}

run_stage_09_10 "10_velocity.sh"
run_stage_09_10 "09_trajectory.sh"
run_stage_09_10 "10_velocity_finalize.sh"
