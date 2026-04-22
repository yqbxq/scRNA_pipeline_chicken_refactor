source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "object_layer_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

input_rds <- file.path(cfg$checkpoint_dir, "01_after_qc_doublet.rds")
if (!file.exists(input_rds)) {
  stop("缺少输入检查点: 01_after_qc_doublet.rds", call. = FALSE)
}

layer_specs <- read_object_layer_specs(cfg)
panorama_id <- panorama_layer_id(layer_specs)
panorama_spec <- layer_specs[[panorama_id]]
if (!isTRUE(panorama_spec$enabled)) {
  stop("panorama 根层不能被禁用", call. = FALSE)
}

panorama_paths <- object_layer_paths(cfg, panorama_id)
ensure_object_layer_dirs(panorama_paths)

obj <- readRDS(input_rds)
obj <- maybe_join_layers(obj)

prepared <- prepare_parent_object_for_layer(obj, panorama_spec)
if (!identical(prepared$status, "ready")) {
  stop(sprintf("无法构建 panorama 对象: %s", prepared$reason), call. = FALSE)
}

panorama_obj <- build_layer_reduction_candidates(prepared$object, panorama_spec, cfg)

saveRDS(panorama_obj, panorama_paths$candidate_rds)
saveRDS(panorama_obj, file.path(cfg$checkpoint_dir, "02_reduction_candidates.rds"))

message("已保存 panorama reduction candidates: ", panorama_paths$candidate_rds)
message("兼容输出已更新: ", file.path(cfg$checkpoint_dir, "02_reduction_candidates.rds"))
