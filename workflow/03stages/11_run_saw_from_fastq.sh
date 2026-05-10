#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

prepare_project_state_dirs

if [[ "${ST_ENABLED:-yes}" == "no" ]]; then
  die "ST_ENABLED=no；跳过 SAW 上游入口。"
fi

die "11_run_saw_from_fastq.sh 当前是 SAW 上游占位入口。M2 只校验 SAW bundle 是否存在；真实 GEF 解析/转换在 M5 的 workflow/04python/saw_intake.py 中落地。"
