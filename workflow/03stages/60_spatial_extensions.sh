#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${STAGE_DIR}/../.." && pwd)}"

# Gate handling for spatial_enrichment, spatial_neighborhood, and spatial_niche
# lives in 60_run_spatial_extensions.sh; this wrapper keeps the legacy entrypoint.
bash "${STAGE_DIR}/60_run_spatial_extensions.sh" "$@"
