#!/usr/bin/env Rscript

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(getwd()), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/spatial/helpers/spatial_common.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/spatial/helpers/project_paths_spatial.R"), encoding = "UTF-8")
load_required_packages(c("Seurat", "jsonlite"))

cfg <- get_spatial_script_config()
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) > 0) args[[1]] else "one_enabled"
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

index <- spatial_read_tsv(file.path(cfg$spatial_subcluster_build_table_dir, "region_subset_index.tsv"))
evidence <- spatial_read_tsv(file.path(cfg$spatial_subcluster_annotate_table_dir, "evidence_summary.tsv"))
assignment <- spatial_read_tsv(file.path(cfg$spatial_subcluster_annotate_table_dir, "region_subset_assignment.tsv"))
triage <- spatial_read_tsv(file.path(cfg$spatial_subcluster_eda_table_dir, "triage.tsv"))

if (identical(mode, "default_disabled")) {
  assert(any(index$status == "skipped_no_enabled_layers"), "default-disabled 04a did not report skipped_no_enabled_layers")
  assert(any(evidence$status == "no_subsets_built"), "default-disabled 04b did not report no_subsets_built")
  assert(any(triage$status == "skipped"), "default-disabled 04c did not write skipped triage")
} else if (identical(mode, "one_enabled")) {
  assert(any(index$layer_id == "region_gc" & index$status == "ok"), "region_gc 04a did not finish ok")
  assert(any(evidence$layer_id == "region_gc" & evidence$status == "ok"), "region_gc 04b did not finish ok")
  assert(nrow(assignment) > 0 && any(assignment$layer_id == "region_gc"), "region_gc assignment is empty")
  assert(file.exists(cfg$spatial_panorama_subannotated_rds), "panorama_subannotated.rds is missing")
  panorama <- readRDS(cfg$spatial_panorama_subannotated_rds)
  assert("sub_region" %in% colnames(panorama@meta.data), "sub_region column is missing from panorama")
  assert(sum(!is.na(panorama$sub_region)) > 0, "sub_region projection is empty")
  figs <- list.files(file.path(cfg$spatial_subcluster_eda_figure_dir, "region_gc"), pattern = "\\.png$", full.names = TRUE)
  assert(length(figs) >= 4 && all(file.info(figs)$size > 0), "region_gc EDA figures are incomplete")
} else if (identical(mode, "one_enabled_no_panel")) {
  assert(any(evidence$layer_id == "region_gc" & evidence$status == "failed_panel_missing"), "missing panel path did not report failed_panel_missing")
  assert(nrow(triage) > 0, "04c did not write triage for missing panel case")
} else if (identical(mode, "two_enabled")) {
  assert(all(c("region_gc", "region_tc") %in% index$layer_id[index$status == "ok"]), "two-enabled 04a did not build both layers")
  assert(all(c("region_gc", "region_tc") %in% evidence$layer_id[evidence$status == "ok"]), "two-enabled 04b did not annotate both layers")
  assert(all(c("region_gc", "region_tc") %in% unique(assignment$layer_id)), "two-enabled assignment is missing a layer")
  assert(all(c("region_gc", "region_tc") %in% unique(triage$layer_id)), "two-enabled triage is missing a layer")
} else {
  stop(sprintf("unknown smoke mode: %s", mode), call. = FALSE)
}

message(sprintf("smoke_spatial_subcluster_assertions_ok:%s", mode))
