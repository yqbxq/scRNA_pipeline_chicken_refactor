#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

Rscript - <<'RS'
source("workflow/05single_script/helpers/evidence_tier_dispatch_utils.R", encoding = "UTF-8")

eligibility <- data.frame(
  layer_id = c("subregion", "subregion"),
  group_var = c("condition", "condition"),
  comparison_id = c("syf_vs_f5_spatial", "syf_vs_f5_spatial"),
  cluster_id = c("GC_core", "TC_edge"),
  evidence_tier = c("exploratory", "skip"),
  recommended_action = c("run_pseudobulk_de_with_warning", "skip"),
  reason = c("ST N=2", "group missing"),
  parent_cluster_id = c("GC", "TC"),
  stringsAsFactors = FALSE
)

decision <- evidence_decision_05(eligibility, "subregion", "syf_vs_f5_spatial", "condition", cluster_id = "GC_core")
stopifnot(identical(decision$evidence_tier, "exploratory"))
stopifnot(evidence_allows_formal_05(decision))

skip_decision <- evidence_decision_05(eligibility, "subregion", "syf_vs_f5_spatial", "condition", cluster_id = "TC_edge")
stopifnot(identical(skip_decision$evidence_tier, "skip"))
stopifnot(!evidence_allows_formal_05(skip_decision))
RS

rg -q "cluster_eligibility_tsv" workflow/05single_script/spatial/05a_region_pseudobulk.R
rg -q "evidence_tier" workflow/05single_script/spatial/05a_region_pseudobulk.R
rg -q "exploratory_fallback = \"\"" workflow/05single_script/spatial/05a_region_pseudobulk.R

echo "smoke_spatial_05a_evidence_tiers_ok"
