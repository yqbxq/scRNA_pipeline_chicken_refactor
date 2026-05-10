#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${STAGE_DIR}/../.." && pwd)}"

bash "${STAGE_DIR}/70_run_joint_scrna_spatial.sh" "$@"
