#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

Rscript - <<'RS'
source("workflow/05single_script/helpers/cell_count_inventory_utils.R", encoding = "UTF-8")
Sys.setenv(
  CELL_COUNT_INVENTORY_MIN_CELLS_PER_SAMPLE = "2",
  CELL_COUNT_INVENTORY_MIN_SAMPLES_PER_GROUP = "1",
  CELL_COUNT_INVENTORY_MIN_TOTAL_UMI = "10",
  CELL_COUNT_INVENTORY_MAX_SINGLE_SAMPLE_FRAC = "0.95"
)

meta <- data.frame(
  sub_region = rep(c("GC_core", "TC_edge"), each = 8),
  section_id = rep(c("sec1", "sec2"), each = 4, times = 2),
  region_label = rep(c("GC", "TC"), each = 8),
  condition = rep(c("SYF", "F5"), each = 4, times = 2),
  nCount_Spatial = 100,
  stringsAsFactors = FALSE
)
inventory <- build_cell_count_inventory(
  meta,
  cluster_var = "sub_region",
  sample_var = "section_id",
  group_var = "condition",
  parent_var = "region_label",
  count_var = "nCount_Spatial"
)
comparison <- data.frame(
  comparison_id = "syf_vs_f5_spatial",
  group_var = "condition",
  ident_1 = "SYF",
  ident_2 = "F5",
  enabled = "yes",
  stringsAsFactors = FALSE
)
eligibility <- classify_all_cluster_eligibilities(inventory, comparison)
stopifnot(nrow(inventory) == 4)
stopifnot(nrow(eligibility) == 2)
stopifnot(all(eligibility$evidence_tier %in% c("primary", "exploratory", "merge_to_parent", "candidate_only", "skip")))
RS

rg -q "cell_count_inventory_tsv" workflow/05single_script/spatial/04c_subcluster_eda.R
rg -q "cluster_eligibility_tsv" workflow/05single_script/spatial/04c_subcluster_eda.R
rg -q "inventory_summary_per_comparison_tsv" workflow/05single_script/spatial/04c_subcluster_eda.R
rg -q "Cell-count inventory \\+ evidence tier" workflow/05single_script/spatial/04c_subcluster_eda.R
rg -q "require_manifest_output.*cluster_eligibility_tsv" workflow/03stages/50_run_spatial_pipeline.sh

echo "smoke_spatial_04c_inventory_ok"
