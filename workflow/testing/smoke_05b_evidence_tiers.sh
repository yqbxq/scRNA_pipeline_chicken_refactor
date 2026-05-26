#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

Rscript - <<'RS'
source("workflow/05single_script/helpers/evidence_tier_dispatch_utils.R", encoding = "UTF-8")

eligibility <- data.frame(
  layer_id = c("panorama", "panorama", "panorama"),
  group_var = c("condition", "condition", "condition"),
  comparison_id = c("syf_vs_f5", "syf_vs_f5", "syf_vs_f5"),
  cluster_id = c("pGC", "eGC", "rgGC"),
  evidence_tier = c("primary", "module_score_only", "candidate_only"),
  recommended_action = c("run_pseudobulk_de", "run_module_score", "report_as_candidate"),
  reason = c("ok", "too few cells", "single sample dominance"),
  parent_cluster_id = c("", "", "GC"),
  stringsAsFactors = FALSE
)

decision <- evidence_decision_05(eligibility, "panorama", "syf_vs_f5", "condition")
stopifnot(identical(decision$evidence_tier, "candidate_only"))
stopifnot(!evidence_allows_formal_05(decision))

cluster_decision <- evidence_decision_05(eligibility, "panorama", "syf_vs_f5", "condition", cluster_id = "pGC")
stopifnot(identical(cluster_decision$evidence_tier, "primary"))
stopifnot(evidence_allows_formal_05(cluster_decision))
RS

! rg -q "run_exploratory_findmarkers_05\\(" workflow/05single_script/05b_pseudobulk_de.R
rg -q "cluster_eligibility_tsv" workflow/05single_script/05b_pseudobulk_de.R
rg -q "evidence_tier" workflow/05single_script/05b_pseudobulk_de.R
rg -q "require_manifest_output.*cluster_eligibility_tsv" workflow/03stages/05_deg.sh

echo "smoke_05b_evidence_tiers_ok"
