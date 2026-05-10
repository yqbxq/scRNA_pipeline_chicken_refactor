#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${STAGE_DIR}/../.." && pwd)}"

bash "${STAGE_DIR}/50_run_spatial_pipeline.sh" "$@"
