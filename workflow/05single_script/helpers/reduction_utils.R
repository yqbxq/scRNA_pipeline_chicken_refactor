strip_reduction_state <- function(seu) {
  seu <- maybe_join_layers(seu)
  if ("reductions" %in% slotNames(seu)) {
    seu@reductions <- list()
  }
  if ("graphs" %in% slotNames(seu)) {
    seu@graphs <- list()
  }
  if ("neighbors" %in% slotNames(seu)) {
    seu@neighbors <- list()
  }
  if ("commands" %in% slotNames(seu)) {
    seu@commands <- list()
  }
  VariableFeatures(seu) <- character(0)

  stale_cols <- c("seurat_clusters", "cluster_id", "cell_type", "cell_type_confidence", "annotation_relation")
  stale_cols <- c(stale_cols, grep("_cluster$", colnames(seu@meta.data), value = TRUE))
  keep_cols <- setdiff(colnames(seu@meta.data), unique(stale_cols))
  seu@meta.data <- seu@meta.data[, keep_cols, drop = FALSE]
  seu
}

strip_reduction_state_preserving_parent_meta <- function(seu, parent_prefix = "panorama_") {
  parent_prefix <- normalize_scalar_value(parent_prefix, "panorama_")
  meta_df <- seu@meta.data
  preserve_cols <- unique(c(
    "seurat_clusters",
    "cluster_id",
    "cell_type",
    "cell_type_confidence",
    "annotation_relation",
    grep("_cluster$", colnames(meta_df), value = TRUE)
  ))
  preserve_cols <- intersect(preserve_cols, colnames(meta_df))

  preserved <- list()
  for (col in preserve_cols) {
    target_col <- if (startsWith(col, parent_prefix)) col else paste0(parent_prefix, col)
    preserved[[target_col]] <- meta_df[[col]]
  }

  seu <- strip_reduction_state(seu)
  for (target_col in names(preserved)) {
    seu[[target_col]] <- preserved[[target_col]]
  }
  seu
}

validate_regression_variables <- function(seu, vars_to_regress) {
  vars_to_regress <- vars_to_regress[nzchar(vars_to_regress)]
  if (length(vars_to_regress) == 0) {
    return(invisible(TRUE))
  }
  missing_vars <- setdiff(vars_to_regress, colnames(seu@meta.data))
  if (length(missing_vars) > 0) {
    stop(sprintf("vars_to_regress 字段在 metadata 中不存在: %s", paste(missing_vars, collapse = ",")), call. = FALSE)
  }
  invisible(TRUE)
}

normalize_layer <- function(seu, method, layer_spec, cfg) {
  method <- tolower(normalize_scalar_value(method))
  if (!method %in% LAYER_NORMALIZATION_METHODS) {
    stop(sprintf("不支持的 normalization method: %s", method), call. = FALSE)
  }
  seu <- maybe_join_layers(seu)
  vars_to_regress <- layer_spec$vars_to_regress
  validate_regression_variables(seu, vars_to_regress)
  vars_arg <- if (length(vars_to_regress) > 0) vars_to_regress else NULL

  if (identical(method, "lognorm")) {
    DefaultAssay(seu) <- "RNA"
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    seu <- FindVariableFeatures(
      seu,
      selection.method = "vst",
      nfeatures = layer_spec$hvg_nfeatures,
      verbose = FALSE
    )
    seu <- ScaleData(
      seu,
      features = rownames(seu),
      vars.to.regress = vars_arg,
      verbose = FALSE
    )
    return(seu)
  }

  DefaultAssay(seu) <- "RNA"
  seu <- SCTransform(
    seu,
    variable.features.n = layer_spec$hvg_nfeatures,
    vars.to.regress = vars_arg,
    vst.flavor = "v2",
    verbose = FALSE
  )
  DefaultAssay(seu) <- "SCT"
  seu
}

usable_reduction_dims <- function(seu, requested_dims, reduction_name) {
  if (!reduction_name %in% Reductions(seu)) {
    stop(sprintf("对象缺少 reduction: %s", reduction_name), call. = FALSE)
  }
  embed <- Embeddings(seu, reduction = reduction_name)
  max_dims <- ncol(embed)
  dims <- requested_dims[requested_dims <= max_dims]
  dims <- dims[is.finite(dims)]
  if (length(dims) == 0) {
    dims <- seq_len(max_dims)
  }
  dims
}

reduce_pca <- function(seu, layer_spec, assay, reduction_key_prefix, reduction_name = "pca") {
  DefaultAssay(seu) <- assay
  features <- VariableFeatures(seu)
  if (length(features) == 0) {
    features <- rownames(seu)
  }
  npcs <- max(layer_spec$pca_dims)
  npcs <- min(npcs, max(2L, ncol(seu) - 1L), max(2L, length(features) - 1L))
  if (!is.finite(npcs) || npcs < 2L) {
    stop("可用于 PCA 的维度不足。", call. = FALSE)
  }
  seu <- RunPCA(
    seu,
    features = features,
    npcs = npcs,
    reduction.name = reduction_name,
    reduction.key = reduction_key_prefix,
    verbose = FALSE
  )
  seu
}

run_umap_only <- function(seu, layer_spec, reduction_name, umap_name, seed = NULL) {
  dims <- usable_reduction_dims(seu, layer_spec$pca_dims, reduction_name)
  if (length(dims) >= 2) {
    if (is.null(umap_name) || !nzchar(umap_name)) {
      umap_name <- paste0("umap_", reduction_name)
    }
    if (!is.null(seed) && is.finite(seed)) {
      set.seed(as.integer(seed))
    }
    seu <- RunUMAP(
      seu,
      reduction = reduction_name,
      dims = dims,
      reduction.name = umap_name,
      reduction.key = paste0(gsub("[^A-Za-z0-9]", "", toupper(umap_name)), "_"),
      n.neighbors = as.integer(Sys.getenv("UMAP_N_NEIGHBORS", "30")),
      min.dist = as.numeric(Sys.getenv("UMAP_MIN_DIST", "0.3")),
      spread = as.numeric(Sys.getenv("UMAP_SPREAD", "1.0")),
      metric = Sys.getenv("UMAP_METRIC", "cosine"),
      local.connectivity = as.numeric(Sys.getenv("UMAP_LOCAL_CONNECTIVITY", "1")),
      verbose = FALSE
    )
  }
  seu
}

reduce_pca_umap <- function(seu, layer_spec, assay, reduction_key_prefix, umap_name = NULL) {
  .Deprecated("reduce_pca + run_umap_only")
  seu <- reduce_pca(seu, layer_spec, assay, reduction_key_prefix, reduction_name = "pca")
  run_umap_only(seu, layer_spec, "pca", umap_name %||% paste0("umap_rna_", tolower(assay)), seed = as.integer(Sys.getenv("UMAP_SEED", "42")))
}
