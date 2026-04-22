source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(DoubletFinder)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

inventory_df <- read_input_inventory(cfg)
sample_names <- main_sample_names(cfg)
sample_dirs <- file.path(cfg$data_dir, sample_names)
if (!all(dir.exists(sample_dirs))) {
  missing_dirs <- sample_dirs[!dir.exists(sample_dirs)]
  stop(sprintf("缺少主流程输入矩阵目录: %s", paste(missing_dirs, collapse = ", ")), call. = FALSE)
}

message("开始读取主流程输入矩阵...")
seurat_list <- lapply(sample_names, function(sample_name) {
  inventory_row <- resolve_inventory_row(cfg, sample_name, inventory_df = inventory_df)
  gene_id_type_resolved <- if (!is.null(inventory_row) && "gene_id_type_resolved" %in% colnames(inventory_row)) {
    normalize_scalar_value(inventory_row$gene_id_type_resolved[1], "unknown")
  } else {
    "unknown"
  }
  obj <- CreateSeuratObject(
    counts = Read10X(data.dir = file.path(cfg$data_dir, sample_name)),
    project = sample_name,
    min.features = 0
  )
  obj$group <- sample_name
  obj$analysis_group <- if (
    nzchar(cfg$analysis_group_1_name) &&
      sample_name %in% cfg$analysis_group_1_samples
  ) {
    cfg$analysis_group_1_name
  } else if (
    nzchar(cfg$analysis_group_2_name) &&
      sample_name %in% cfg$analysis_group_2_samples
  ) {
    cfg$analysis_group_2_name
  } else {
    sample_name
  }
  add_basic_qc_metrics(obj, cfg, declared_gene_id_type = gene_id_type_resolved)$object
})
names(seurat_list) <- sample_names

merged <- seurat_list[[1]]
if (length(seurat_list) > 1) {
  merged <- merge(
    x = seurat_list[[1]],
    y = seurat_list[2:length(seurat_list)],
    add.cell.ids = sample_names
  )
}
merged <- maybe_join_layers(merged)

filtered <- subset(
  merged,
  subset =
    nFeature_RNA >= cfg$qc_min_nfeature &
      nCount_RNA >= cfg$qc_min_ncount &
      log10GenesPerUMI >= cfg$qc_min_log10umi &
      percent.mito <= cfg$qc_max_mito_pct
)

doublet_pkg <- asNamespace("DoubletFinder")
sweep_fun <- if (exists("paramSweep_v3", envir = doublet_pkg, mode = "function")) {
  get("paramSweep_v3", envir = doublet_pkg)
} else {
  get("paramSweep", envir = doublet_pkg)
}
df_fun <- if (exists("doubletFinder_v3", envir = doublet_pkg, mode = "function")) {
  get("doubletFinder_v3", envir = doublet_pkg)
} else {
  get("doubletFinder", envir = doublet_pkg)
}

split_objs <- SplitObject(filtered, split.by = "orig.ident")
clean_list <- list()
summary_rows <- list()

for (sample_name in names(split_objs)) {
  message("运行 DoubletFinder: ", sample_name)
  obj <- split_objs[[sample_name]]
  n_after_qc <- ncol(obj)
  obj <- maybe_join_layers(obj)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = cfg$hvg_nfeatures, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = cfg$doublet_dims, verbose = FALSE)

  sweep_res <- sweep_fun(obj, PCs = cfg$doublet_dims, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  best_pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  expected_doublets <- as.integer(round(cfg$doublet_rate * ncol(obj) / 1000 * ncol(obj)))

  obj <- df_fun(
    obj,
    PCs = cfg$doublet_dims,
    pN = 0.25,
    pK = best_pk,
    nExp = expected_doublets,
    sct = FALSE
  )

  df_cols <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  if (length(df_cols) == 0) {
    stop(sprintf("DoubletFinder 没有返回分类列: %s", sample_name), call. = FALSE)
  }
  obj$Doublet_Class <- obj@meta.data[[tail(df_cols, 1)]]
  obj <- subset(obj, subset = Doublet_Class == "Singlet")

  keep_cols <- !grepl("DF.classifications|pANN", colnames(obj@meta.data))
  obj@meta.data <- obj@meta.data[, keep_cols, drop = FALSE]
  clean_list[[sample_name]] <- obj
  summary_rows[[sample_name]] <- data.frame(
    sample = sample_name,
    cells_after_qc = n_after_qc,
    singlets_after_doublet = ncol(obj),
    stringsAsFactors = FALSE
  )
}

final_obj <- clean_list[[1]]
if (length(clean_list) > 1) {
  final_obj <- merge(x = clean_list[[1]], y = clean_list[2:length(clean_list)])
}
final_obj <- maybe_join_layers(final_obj)

summary_df <- bind_rows(summary_rows)

saveRDS(final_obj, file.path(cfg$checkpoint_dir, "01_after_qc_doublet.rds"))
write.csv(summary_df, file.path(cfg$table_dir, "qc_doublet_summary.csv"), row.names = FALSE)

message("已保存: 01_after_qc_doublet.rds")
