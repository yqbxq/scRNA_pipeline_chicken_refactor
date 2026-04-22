#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

sync_workflow_gate_statuses

require_status_flag_or_warn \
  "status.velocity_reference_ready" \
  "scVelo 尚未 ready。请先完成主流程 annotation hub，并确认 velocity 上游输入已通过审计。"

run_r_main "${PIPELINE_ROOT}/workflow/r/06_prepare_velocity_reference.R"

export VELOCITY_INPUT_DIR VELOCITY_LOOM_DIR VELOCITY_OUTPUT_DIR SCVELO_THREADS
run_in_conda_prefix "${VELOCITY_ENV_PREFIX}" python "${PIPELINE_ROOT}/workflow/python/scvelo_pipeline.py"

update_workflow_status \
  "scvelo_completed" \
  "按需审阅 ${VELOCITY_OUTPUT_DIR} 或继续 SCENIC 分支"

echo "scVelo 运行完成"
