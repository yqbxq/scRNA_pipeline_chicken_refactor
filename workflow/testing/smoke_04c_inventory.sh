#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

Rscript - <<'RS'
source("workflow/05single_script/helpers/cell_count_inventory_utils.R", encoding = "UTF-8")
Sys.setenv(
  CELL_COUNT_INVENTORY_MIN_CELLS_PER_SAMPLE = "2",
  CELL_COUNT_INVENTORY_MIN_SAMPLES_PER_GROUP = "2",
  CELL_COUNT_INVENTORY_MIN_TOTAL_UMI = "10",
  CELL_COUNT_INVENTORY_MAX_SINGLE_SAMPLE_FRAC = "0.95"
)

meta <- data.frame(
  cell_subtype = rep(c("pGC", "eGC"), each = 16),
  sample_id = rep(rep(c("s1", "s2", "s3", "s4"), each = 4), 2),
  condition = rep(rep(c("SYF", "SYF", "F5", "F5"), each = 4), 2),
  parent_cell_type = rep(c("GC", "GC"), each = 16),
  nCount_RNA = 100,
  stringsAsFactors = FALSE
)
inventory <- build_cell_count_inventory(
  meta,
  cluster_var = "cell_subtype",
  sample_var = "sample_id",
  group_var = "condition",
  parent_var = "parent_cell_type",
  count_var = "nCount_RNA"
)
comparison <- data.frame(
  comparison_id = "syf_vs_f5_gc",
  group_var = "condition",
  ident_1 = "SYF",
  ident_2 = "F5",
  enabled = "yes",
  stringsAsFactors = FALSE
)
eligibility <- classify_all_cluster_eligibilities(inventory, comparison)
stopifnot(nrow(inventory) == 8)
stopifnot(nrow(eligibility) == 2)
stopifnot(all(eligibility$evidence_tier == "primary"))
RS

rg -q "cell_count_inventory_tsv" workflow/05single_script/04c_subcluster_eda.R
rg -q "cluster_eligibility_tsv" workflow/05single_script/04c_subcluster_eda.R
rg -q "inventory_summary_per_comparison_tsv" workflow/05single_script/04c_subcluster_eda.R
rg -q "Cell-count inventory \\+ evidence tier" workflow/05single_script/04c_subcluster_eda.R
rg -q "require_manifest_output.*cluster_eligibility_tsv" workflow/03stages/04_subcluster.sh

echo "smoke_04c_inventory_ok"
