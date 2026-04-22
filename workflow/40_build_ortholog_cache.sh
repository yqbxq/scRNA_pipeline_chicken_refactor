#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

sync_workflow_gate_statuses

[[ -f "${ANNOTATION_HUB_PATH}" ]] || die "缺少 annotation hub: ${ANNOTATION_HUB_PATH}"

require_status_flag_or_warn \
  "status.annotation_gate_passed" \
  "同源映射缓存尚未 ready。请先完成 annotation EDA 审阅并批准 annotation gate。"

run_r_main "${PIPELINE_ROOT}/workflow/r/06_build_ortholog_cache.R"

update_workflow_status \
  "ortholog_cache_completed" \
  "bash ${PIPELINE_ROOT}/workflow/42_run_scenic.sh"

echo "同源映射缓存已生成: ${ORTHOLOG_CACHE_DIR}"
