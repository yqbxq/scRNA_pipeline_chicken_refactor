#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
SHELL_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${SHELL_ROOT}/.." && pwd)}"

source "${SHELL_ROOT}/02lib/common.sh"
source "${SHELL_ROOT}/02lib/env_registry.sh"
source "${SHELL_ROOT}/02lib/executor.sh"
source "${SHELL_ROOT}/02lib/stage.sh"

check_stage_deps "00_ortholog"

ORTHOLOG_OUTPUT="${ORTHOLOG_CACHE_DIR}/chicken_human_orthologs.csv"
CC_GENES_OUTPUT="${ORTHOLOG_CACHE_DIR}/chicken_cc_genes.rds"
ORTHOLOG_REPORT="${ORTHOLOG_CACHE_DIR}/ortholog_quality_report.md"

run_stage_if_stale \
  "${SHELL_ROOT}/05single_script/00a_build_ortholog_cache.R" \
  "${ORTHOLOG_OUTPUT}" \
  "${REFERENCE_GTF}"

run_stage_if_stale \
  "${SHELL_ROOT}/05single_script/00b_ortholog_report.R" \
  "${ORTHOLOG_REPORT}" \
  "${ORTHOLOG_OUTPUT}"

run_stage_if_stale \
  "${SHELL_ROOT}/05single_script/00c_cc_gene_mapping.R" \
  "${CC_GENES_OUTPUT}" \
  "${ORTHOLOG_OUTPUT}"

update_workflow_status \
  "00_ortholog_completed" \
  "" \
  "status.00_ortholog_completed=true"
