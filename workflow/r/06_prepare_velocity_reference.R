source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

ref_rds <- file.path(cfg$checkpoint_dir, "03_after_annotation.rds")
if (!file.exists(ref_rds)) {
  stop("缺少主流程注释结果: 03_after_annotation.rds", call. = FALSE)
}

reference_obj <- readRDS(ref_rds)
reference_obj <- maybe_join_layers(reference_obj)
DefaultAssay(reference_obj) <- "RNA"

if (length(cfg$raw_samples) == 0) {
  stop("server_config.sh 中 RAW_SAMPLES 为空", call. = FALSE)
}

raw_list <- list()
for (sample_id in cfg$raw_samples) {
  matrix_dir <- file.path(cfg$cellranger_out_dir, sample_id, "outs", "filtered_feature_bc_matrix")
  if (!dir.exists(matrix_dir)) {
    stop(sprintf("缺少 raw matrix 目录: %s", matrix_dir), call. = FALSE)
  }

  obj <- CreateSeuratObject(
    counts = Read10X(data.dir = matrix_dir),
    project = sample_id,
    min.features = 0
  )
  obj <- RenameCells(obj, new.names = paste0(sample_id, ":", colnames(obj)))
  obj$sample_id <- sample_id
  obj$group <- if (sample_id %in% cfg$raw_group_1_samples) {
    cfg$raw_group_1_name
  } else if (sample_id %in% cfg$raw_group_2_samples) {
    cfg$raw_group_2_name
  } else {
    "Unknown"
  }
  raw_list[[sample_id]] <- obj
}

query_obj <- raw_list[[1]]
if (length(raw_list) > 1) {
  query_obj <- merge(x = raw_list[[1]], y = raw_list[2:length(raw_list)])
}
query_obj <- maybe_join_layers(query_obj)
DefaultAssay(query_obj) <- "RNA"

query_obj <- NormalizeData(query_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
query_obj <- FindVariableFeatures(query_obj, selection.method = "vst", nfeatures = cfg$hvg_nfeatures, verbose = FALSE)
query_obj <- ScaleData(query_obj, features = rownames(query_obj), verbose = FALSE)
query_obj <- RunPCA(query_obj, features = VariableFeatures(query_obj), verbose = FALSE)
query_obj <- RunHarmony(query_obj, group.by.vars = "sample_id", theta = 3, lambda = 1, plot_convergence = FALSE)
query_obj <- RunUMAP(query_obj, reduction = "harmony", dims = cfg$pca_dims, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = reference_obj,
  query = query_obj,
  normalization.method = "LogNormalize",
  dims = cfg$pca_dims
)

celltype_predictions <- TransferData(
  anchorset = anchors,
  refdata = reference_obj$cell_type,
  dims = cfg$pca_dims
)
cluster_predictions <- TransferData(
  anchorset = anchors,
  refdata = reference_obj$seurat_clusters,
  dims = cfg$pca_dims
)

query_obj <- AddMetaData(query_obj, metadata = celltype_predictions)
query_obj$cell_type <- factor(
  query_obj$predicted.id,
  levels = levels(reference_obj$cell_type)
)

query_obj <- AddMetaData(
  query_obj,
  metadata = setNames(cluster_predictions["predicted.id"], "predicted_cluster")
)
query_obj$seurat_clusters <- factor(
  as.character(query_obj$predicted_cluster),
  levels = levels(reference_obj$seurat_clusters)
)

umap_df <- as.data.frame(Embeddings(query_obj, reduction = "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

write.csv(umap_df, file.path(cfg$velocity_input_dir, "velocity_umap.csv"))
write.csv(query_obj@meta.data, file.path(cfg$velocity_input_dir, "velocity_metadata.csv"))
saveRDS(query_obj, file.path(cfg$checkpoint_dir, "velocity_reference_query.rds"))

message("已导出 RNA velocity 所需的参考对象、UMAP 与 metadata")
