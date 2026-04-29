#!/usr/bin/env bash
set -euo pipefail

WORKFLOW_ROOT="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

usage() {
  cat <<'EOF'
Usage:
  bash workflow/01run.sh <stage>

Examples:
  bash workflow/01run.sh 08_regulation
  bash workflow/01run.sh workflow/03stages/08_regulation.sh
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 1 ]]; then
  usage
  exit 0
fi

stage_arg="$1"
shift || true

case "${stage_arg}" in
  */*)
    stage_path="${stage_arg}"
    ;;
  *.sh)
    stage_path="${WORKFLOW_ROOT}/03stages/${stage_arg}"
    ;;
  *)
    stage_path="${WORKFLOW_ROOT}/03stages/${stage_arg}.sh"
    ;;
esac

if [[ ! -f "${stage_path}" ]]; then
  echo "[ERROR] 未找到 stage: ${stage_arg} -> ${stage_path}" >&2
  exit 1
fi

if [[ ! -s "${stage_path}" ]]; then
  echo "[ERROR] stage 文件为空: ${stage_path}" >&2
  exit 1
fi

PIPELINE_ROOT="${PIPELINE_ROOT}" bash "${stage_path}" "$@"
