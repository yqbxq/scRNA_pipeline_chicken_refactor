#!/usr/bin/env Rscript

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(getwd()), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/spatial/helpers/spatial_common.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/spatial/helpers/project_paths_spatial.R"), encoding = "UTF-8")
load_required_packages(c("Seurat", "jsonlite"))

cfg <- get_spatial_script_config()
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

assert(file.exists(cfg$spatial_panorama_integrated_rds), "missing spatial_panorama_integrated.rds")
panorama <- readRDS(cfg$spatial_panorama_integrated_rds)
assert("pca_none" %in% names(panorama@reductions), "missing pca_none reduction")
assert("umap_none" %in% names(panorama@reductions), "missing umap_none reduction")

lisi <- spatial_read_tsv(file.path(cfg$spatial_integration_compare_dir, "lisi_summary.tsv"))
bio <- spatial_read_tsv(file.path(cfg$spatial_integration_compare_dir, "biological_consistency.tsv"))
recommended <- spatial_read_tsv(file.path(cfg$spatial_integration_compare_dir, "recommended_mode.tsv"))
assert(nrow(lisi) >= 1 && any(lisi$status == "ok"), "lisi_summary.tsv has no ok row")
assert(nrow(bio) >= 1, "biological_consistency.tsv has no rows")
assert(identical(recommended$integration_mode_default[[1]], "none"), "default integration mode is not none")
if (!requireNamespace("harmony", quietly = TRUE) && "harmony" %in% lisi$integration_mode) {
  assert(lisi$status[lisi$integration_mode == "harmony"][[1]] == "failed_dependency", "harmony missing dependency path was not recorded")
}

message("smoke_spatial_integration_assertions_ok")
