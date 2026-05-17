#!/usr/bin/env Rscript

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(getwd()), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/spatial/helpers/spatial_common.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/spatial/helpers/project_paths_spatial.R"), encoding = "UTF-8")
load_required_packages(c("Seurat", "jsonlite"))

cfg <- get_spatial_script_config()
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

expected_levels <- c("GC_rich", "TC_rich", "stroma", "vasculature", "mixed_or_uncertain")
assert(file.exists(cfg$spatial_panorama_annotated_rds), "missing spatial_panorama_annotated.rds")
panorama <- readRDS(cfg$spatial_panorama_annotated_rds)
assert("region" %in% colnames(panorama@meta.data), "region column is missing")
assert(is.factor(panorama$region), "region is not a factor")
assert(all(expected_levels %in% levels(panorama$region)), "region factor levels are incomplete")

assignment <- spatial_read_tsv(file.path(cfg$spatial_region_annotation_table_dir, "region_assignment.tsv"))
cluster_n <- length(unique(panorama$cluster_default))
assert(nrow(assignment) == cluster_n, "region_assignment.tsv row count does not match cluster count")
assert(any(assignment$region %in% c("GC_rich", "TC_rich")), "no cluster assigned to GC_rich or TC_rich")

triage <- spatial_read_tsv(file.path(cfg$spatial_region_annotation_dir, "region_triage.tsv"))
summary <- spatial_read_tsv(file.path(cfg$spatial_region_annotation_dir, "region_summary.tsv"))
assert(nrow(triage) == 6, "region_triage.tsv does not contain 6 signals")
assert(nrow(summary) > 0, "region_summary.tsv is empty")

figs <- file.path(
  cfg$spatial_region_annotation_figure_dir,
  c("spatial_region_assignment.png", "region_marker_heatmap.png", "cluster_region_sankey.png", "region_spot_barchart.png")
)
assert(all(file.exists(figs) & file.info(figs)$size > 0), "one or more region annotation figures are missing")

message("smoke_spatial_region_annotation_assertions_ok")
