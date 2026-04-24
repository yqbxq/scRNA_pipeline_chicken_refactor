#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
SHELL_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${SHELL_ROOT}/.." && pwd)}"

source "${SHELL_ROOT}/02lib/common.sh"

check_stage_deps "01_build_raw"

RAW_OBJECTS_R="${RAW_OBJECTS_R:-${SHELL_ROOT}/05single_script/01_build_raw_objects.R}"
RAW_OUTPUT="${CHECKPOINT_DIR}/00_raw_objects.rds"

run_stage_if_stale \
  "${RAW_OBJECTS_R}" \
  "${RAW_OUTPUT}" \
  "${CANONICAL_SAMPLE_SHEET}" \
  "${INPUT_INVENTORY_FILE}"

update_workflow_status \
  "01_build_raw_completed" \
  "" \
  "status.01_build_raw_completed=true"
