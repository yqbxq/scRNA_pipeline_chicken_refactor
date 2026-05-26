#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"

Rscript - <<'RS'
source("workflow/05single_script/helpers/cell_count_inventory_utils.R")

make_cells <- function(cluster, group, sample, n, parent = "") {
  data.frame(
    cluster = rep(cluster, n),
    group = rep(group, n),
    sample = rep(sample, n),
    parent = rep(parent, n),
    nCount_RNA = rep(1000, n),
    stringsAsFactors = FALSE
  )
}

meta <- do.call(rbind, list(
  make_cells("A", "syf", "syf1", 50), make_cells("A", "syf", "syf2", 50),
  make_cells("A", "f5", "f51", 50), make_cells("A", "f5", "f52", 50),
  make_cells("B", "syf", "syf1", 15), make_cells("B", "syf", "syf2", 15),
  make_cells("B", "f5", "f51", 15), make_cells("B", "f5", "f52", 15),
  make_cells("C", "syf", "syf1", 8), make_cells("C", "f5", "f51", 8),
  make_cells("D", "syf", "syf1", 10, "D_parent"), make_cells("D", "syf", "syf2", 10, "D_parent"),
  make_cells("D", "f5", "f51", 10, "D_parent"), make_cells("D", "f5", "f52", 10, "D_parent"),
  make_cells("E", "syf", "syf1", 100), make_cells("E", "f5", "f51", 1),
  make_cells("F", "f5", "f51", 50), make_cells("F", "f5", "f52", 50)
))

inventory <- build_cell_count_inventory(
  meta,
  cluster_var = "cluster",
  sample_var = "sample",
  group_var = "group",
  parent_var = "parent",
  count_var = "nCount_RNA"
)

cmp <- data.frame(
  comparison_id = "syf_vs_f5",
  ident_1 = "syf",
  ident_2 = "f5",
  stringsAsFactors = FALSE
)

expect_tier <- function(cluster, tier, action_prefix = NULL) {
  row <- cmp
  row$cluster_id <- cluster
  result <- classify_cluster_eligibility(inventory, row, NULL)
  if (!identical(result$evidence_tier, tier)) {
    stop(sprintf("cluster %s expected %s got %s: %s", cluster, tier, result$evidence_tier, result$reason), call. = FALSE)
  }
  if (!is.null(action_prefix) && !startsWith(result$recommended_action, action_prefix)) {
    stop(sprintf("cluster %s expected action prefix %s got %s", cluster, action_prefix, result$recommended_action), call. = FALSE)
  }
}

expect_tier("A", "primary")
expect_tier("B", "exploratory")
expect_tier("C", "module_score_only")
expect_tier("D", "merge_to_parent", "merge_to:D_parent")
expect_tier("E", "candidate_only")
expect_tier("F", "skip")

row <- cmp
row$cluster_id <- "A"
row$min_cells_override <- 60
override_result <- classify_cluster_eligibility(inventory, row, NULL)
stopifnot(override_result$threshold_resolved$min_cells_per_sample == 60)
stopifnot(override_result$evidence_tier == "exploratory")

thresholds <- data.frame(
  celltype = c("*", "A"),
  min_cells_per_sample = c(45, 55),
  min_samples_per_group = c(2, 2),
  min_total_umi = c(10000, 10000),
  max_single_sample_frac = c(0.8, 0.8),
  min_cells_soft_floor = c(5, 5),
  stringsAsFactors = FALSE
)
row <- cmp
row$cluster_id <- "A"
specific <- classify_cluster_eligibility(inventory, row, thresholds)
stopifnot(specific$threshold_resolved$min_cells_per_sample == 55)
stopifnot(specific$evidence_tier == "exploratory")

thresholds <- thresholds[thresholds$celltype == "*", , drop = FALSE]
star <- classify_cluster_eligibility(inventory, row, thresholds)
stopifnot(star$threshold_resolved$min_cells_per_sample == 45)
stopifnot(star$evidence_tier == "primary")

all_rows <- classify_all_cluster_eligibilities(inventory, cmp, thresholds)
stopifnot(nrow(all_rows) == length(unique(inventory$cluster_id)))
cat("smoke_cell_count_inventory_ok\n")
RS
