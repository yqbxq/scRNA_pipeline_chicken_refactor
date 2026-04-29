#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"
source "${WORKFLOW_ROOT}/02lib/env_registry.sh"
source "${WORKFLOW_ROOT}/02lib/executor.sh"
source "${WORKFLOW_ROOT}/02lib/stage.sh"

check_stage_deps "00_ortholog"

HUMAN_BEST_OUTPUT="${ORTHOLOG_CACHE_DIR}/chicken_human_orthologs.csv"
HUMAN_ALL_OUTPUT="${ORTHOLOG_CACHE_DIR}/chicken_human_orthologs_all_candidates.csv"
MOUSE_BEST_OUTPUT="${ORTHOLOG_CACHE_DIR}/chicken_mouse_orthologs.csv"
MOUSE_ALL_OUTPUT="${ORTHOLOG_CACHE_DIR}/chicken_mouse_orthologs_all_candidates.csv"
ORTHOLOG_MANIFEST="${ORTHOLOG_CACHE_DIR}/_manifest.json"
CC_GENES_OUTPUT="${ORTHOLOG_CACHE_DIR}/chicken_cc_genes.rds"
ORTHOLOG_REPORT="${ORTHOLOG_CACHE_DIR}/ortholog_quality_report.md"
ORTHOLOG_REQUIRED_OUTPUTS=(
  "${HUMAN_BEST_OUTPUT}"
  "${HUMAN_ALL_OUTPUT}"
  "${MOUSE_BEST_OUTPUT}"
  "${MOUSE_ALL_OUTPUT}"
  "${ORTHOLOG_MANIFEST}"
  "${ORTHOLOG_REPORT}"
  "${ORTHOLOG_CACHE_DIR}/figures/ortholog_identity_scatter.png"
  "${ORTHOLOG_CACHE_DIR}/figures/ortholog_goc_distribution.png"
  "${ORTHOLOG_CACHE_DIR}/figures/ortholog_unmapped_biotype.png"
  "${ORTHOLOG_CACHE_DIR}/figures/ortholog_species_comparison.png"
  "${CC_GENES_OUTPUT}"
  "${ORTHOLOG_CACHE_DIR}/chicken_cc_genes_summary.csv"
  "${ORTHOLOG_CACHE_DIR}/chicken_cc_genes_report.md"
)

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/00a_build_ortholog_cache.R" \
  "${ORTHOLOG_MANIFEST}" \
  "${REFERENCE_GTF}"

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/00b_ortholog_report.R" \
  "${ORTHOLOG_REPORT}" \
  "${ORTHOLOG_MANIFEST}" \
  "${HUMAN_BEST_OUTPUT}" \
  "${MOUSE_BEST_OUTPUT}"

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/00c_cc_gene_mapping.R" \
  "${CC_GENES_OUTPUT}" \
  "${ORTHOLOG_MANIFEST}" \
  "${HUMAN_ALL_OUTPUT}"

for required_output in "${ORTHOLOG_REQUIRED_OUTPUTS[@]}"; do
  [[ -e "${required_output}" ]] || die "00_ortholog 缺少产物: ${required_output}"
done

update_workflow_status \
  "00_ortholog_completed" \
  "" \
  "status.00_ortholog_completed=true"
