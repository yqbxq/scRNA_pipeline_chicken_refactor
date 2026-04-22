#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

sync_workflow_gate_statuses

require_status_flag_or_warn \
  "status.scenic_export_ready" \
  "SCENIC 尚未 ready。请先完成主流程 annotation hub，并确认主输入通过 intake 审计。"

"${PIPELINE_ROOT}/workflow/41_download_scenic_resources.sh"

if [[ -z "${SCENIC_ORTHOLOG_MAP_FILE:-}" ]]; then
  "${PIPELINE_ROOT}/workflow/40_build_ortholog_cache.sh"
  export SCENIC_ORTHOLOG_MAP_FILE="${ORTHOLOG_CACHE_DIR}/chicken_human_orthologs_for_pipeline.csv"
elif [[ ! -f "${SCENIC_ORTHOLOG_MAP_FILE}" ]]; then
  die "SCENIC_ORTHOLOG_MAP_FILE 指向的文件不存在: ${SCENIC_ORTHOLOG_MAP_FILE}"
fi

run_r_main "${PIPELINE_ROOT}/workflow/r/07_export_scenic_input.R"

run_pyscenic pyscenic grn \
  --num_workers "${SCENIC_THREADS}" \
  --method grnboost2 \
  -o "${SCENIC_OUTPUT_DIR}/adjacencies.tsv" \
  "${SCENIC_INPUT_DIR}/expr_mat_human.csv" \
  "${SCENIC_TF_LIST}"

run_r_scenic "${PIPELINE_ROOT}/workflow/r/07b_scenic_regulons_aucell.R"
run_r_scenic "${PIPELINE_ROOT}/workflow/r/08_scenic_downstream.R"

update_workflow_status \
  "scenic_completed" \
  "SCENIC 分支完成，按需审阅 ${SCENIC_OUTPUT_DIR} 和 ${FIGURE_DIR}"

echo "SCENIC 分支运行完成"
