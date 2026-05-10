#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

ensure_metadata_fresh
check_stage_deps "04d_cluster_robustness"

MODULE_04B_MANIFEST="${MANIFEST_DIR}/04b_subcluster_annotate/_manifest.json"
MODULE_04C_MANIFEST="${MANIFEST_DIR}/04c_subcluster_eda/_manifest.json"
MODULE_04D_MANIFEST="${MANIFEST_DIR}/04d_cluster_robustness/_manifest.json"

ensure_eda_control_files
ensure_object_layer_config_file

require_manifest_output "${MODULE_04C_MANIFEST}" "subcluster_summary_tsv" >/dev/null

run_stage_if_stale \
  "${WORKFLOW_ROOT}/05single_script/04d_cluster_robustness.R" \
  "${MODULE_04D_MANIFEST}" \
  "${MODULE_04C_MANIFEST}" \
  "${MODULE_04B_MANIFEST}"

require_manifest_output "${MODULE_04D_MANIFEST}" "metrics_tsv" >/dev/null
require_manifest_output "${MODULE_04D_MANIFEST}" "question_gate_status_tsv" >/dev/null
require_manifest_output "${MODULE_04D_MANIFEST}" "scdesign3_all_questions_status_tsv" >/dev/null

update_workflow_status \
  "04d_cluster_robustness_completed" \
  "next: review scDesign3 all-question gate status; run locked scDesign3 simulation engine after formal 04 objects exist" \
  "status.04d_cluster_robustness_completed=true"
