#!/usr/bin/env Rscript

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(getwd()), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/spatial/helpers/spatial_common.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/spatial/helpers/project_paths_spatial.R"), encoding = "UTF-8")
load_required_packages(c("Seurat", "jsonlite"))

cfg <- get_spatial_script_config()
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

assert(file.exists(cfg$spatial_panorama_clustered_rds), "missing spatial_panorama_clustered.rds")
panorama <- readRDS(cfg$spatial_panorama_clustered_rds)
assert("cluster_default" %in% colnames(panorama@meta.data), "cluster_default column is missing")
cluster_n <- length(unique(panorama$cluster_default))
target <- as.integer(Sys.getenv("SPATIAL_CLUSTER_TARGET", unset = "4"))
assert(cluster_n >= 2 && cluster_n <= ceiling(target * 1.5), sprintf("cluster count out of range: %d", cluster_n))

metrics <- spatial_read_tsv(file.path(cfg$clustering_compare_dir, "cluster_metrics.tsv"))
selected <- spatial_read_tsv(cfg$selected_backend_tsv)
assert(any(metrics$backend == "b1_seurat_snn" & metrics$status == "ok"), "b1_seurat_snn did not succeed")
assert(selected$final_choice[[1]] == "b1_seurat_snn", "default selected backend is not b1_seurat_snn")
selected_variant <- spatial_resolve_path(selected$selected_variant_rds[[1]], cfg$project_root)
assert(file.exists(selected_variant), "selected variant RDS is missing")
assert(unname(tools::md5sum(selected_variant)) == unname(tools::md5sum(cfg$spatial_panorama_clustered_rds)), "canonical RDS is not bytewise identical to selected variant")
if (!requireNamespace("BayesSpace", quietly = TRUE)) {
  assert(any(metrics$backend == "b2_bayesspace" & metrics$status == "failed_dependency"), "BayesSpace missing dependency path was not recorded")
}
assert(any(metrics$backend == "b3_spagcn" & metrics$status == "failed_py_bridge"), "SpaGCN failed_py_bridge path was not recorded")
assert(any(metrics$backend == "b4_stagate" & metrics$status == "failed_py_bridge"), "STAGATE failed_py_bridge path was not recorded")

message("smoke_spatial_clustering_assertions_ok")
