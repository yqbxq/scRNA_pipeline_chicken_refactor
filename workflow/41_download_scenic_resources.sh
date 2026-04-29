#!/usr/bin/env bash
set -euo pipefail

PIPELINE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEFAULT_CONFIG_FILE="${PIPELINE_ROOT}/config/server_config.sh"

if [[ -z "${SCRNA_PIPELINE_CONFIG:-}" && -f "${DEFAULT_CONFIG_FILE}" ]]; then
  export SCRNA_PIPELINE_CONFIG="${DEFAULT_CONFIG_FILE}"
fi

exec env SCENIC_RESOURCES_ONLY=yes bash "${PIPELINE_ROOT}/shell/03stages/08_regulation.sh"
