integration_runtime <- function(expr) {
  start_time <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, runtime_sec = proc.time()[["elapsed"]] - start_time)
}

integration_result <- function(object, method, reduction_name, umap_name, runtime_sec, downgrade_reason = "", extra = list()) {
  list(
    object = object,
    reduction_name = reduction_name,
    umap_name = umap_name,
    diagnostics = c(
      list(
        method = method,
        runtime_sec = runtime_sec,
        downgrade_reason = downgrade_reason
      ),
      extra
    )
  )
}

integrate_none <- function(seu, method, layer_spec, normalization_method, runtime_sec = 0, downgrade_reason = "") {
  integration_result(
    object = seu,
    method = method,
    reduction_name = "pca",
    umap_name = sprintf("umap_%s_%s", method, normalization_method),
    runtime_sec = runtime_sec,
    downgrade_reason = downgrade_reason
  )
}

integrate_harmony <- function(seu, group_var, layer_spec, normalization_method) {
  if (!group_var %in% colnames(seu@meta.data)) {
    return(integrate_none(seu, "harmony", layer_spec, normalization_method, downgrade_reason = "missing_group_var"))
  }
  if (length(unique(as.character(seu@meta.data[[group_var]]))) <= 1) {
    return(integrate_none(seu, "harmony", layer_spec, normalization_method, downgrade_reason = "single_group"))
  }
  if (!requireNamespace("harmony", quietly = TRUE)) {
    return(integrate_none(seu, "harmony", layer_spec, normalization_method, downgrade_reason = "harmony_unavailable"))
  }
  timed <- tryCatch(
    integration_runtime({
      harmony::RunHarmony(
        seu,
        group.by.vars = group_var,
        reduction.use = "pca",
        dims.use = usable_reduction_dims(seu, layer_spec$pca_dims, "pca"),
        theta = 2,
        lambda = 1,
        max.iter.harmony = 20,
        plot_convergence = FALSE,
        verbose = FALSE
      )
    }),
    error = function(e) list(value = seu, runtime_sec = 0, error = conditionMessage(e))
  )
  if (!is.null(timed$error) || !"harmony" %in% Reductions(timed$value)) {
    return(integrate_none(timed$value, "harmony", layer_spec, normalization_method, timed$runtime_sec, paste0("harmony_failed:", timed$error %||% "no_reduction")))
  }
  integration_result(timed$value, "harmony", "harmony", sprintf("umap_harmony_%s", normalization_method), timed$runtime_sec)
}

integrate_seurat_anchor <- function(seu, group_var, layer_spec, normalization_method, reduction_method) {
  method <- if (identical(reduction_method, "cca")) "seurat_cca" else "seurat_rpca"
  if (!group_var %in% colnames(seu@meta.data)) {
    return(integrate_none(seu, method, layer_spec, normalization_method, downgrade_reason = "missing_group_var"))
  }
  if (length(unique(as.character(seu@meta.data[[group_var]]))) <= 1) {
    return(integrate_none(seu, method, layer_spec, normalization_method, downgrade_reason = "single_group"))
  }
  timed <- tryCatch(
    integration_runtime({
      assay_name <- if (identical(normalization_method, "sct") && "SCT" %in% Assays(seu)) "SCT" else "RNA"
      DefaultAssay(seu) <- assay_name
      obj_list <- SplitObject(seu, split.by = group_var)
      obj_list <- obj_list[vapply(obj_list, ncol, integer(1)) > 0]
      if (length(obj_list) <= 1) {
        stop("single_nonempty_group", call. = FALSE)
      }
      features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = layer_spec$hvg_nfeatures)
      normalization_method_seurat <- if (identical(normalization_method, "sct")) "SCT" else "LogNormalize"
      if (identical(normalization_method_seurat, "SCT")) {
        obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features, verbose = FALSE)
      } else {
        obj_list <- lapply(obj_list, function(obj) {
          DefaultAssay(obj) <- "RNA"
          obj <- NormalizeData(obj, verbose = FALSE)
          obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = layer_spec$hvg_nfeatures, verbose = FALSE)
          obj
        })
      }
      if (identical(reduction_method, "rpca")) {
        obj_list <- lapply(obj_list, function(obj) {
          obj <- ScaleData(obj, features = features, verbose = FALSE)
          RunPCA(obj, features = features, npcs = max(layer_spec$pca_dims), verbose = FALSE)
        })
      }
      anchors <- FindIntegrationAnchors(
        object.list = obj_list,
        anchor.features = features,
        normalization.method = normalization_method_seurat,
        reduction = reduction_method,
        dims = usable_reduction_dims(seu, layer_spec$pca_dims, "pca"),
        verbose = FALSE
      )
      integrated <- IntegrateData(
        anchorset = anchors,
        normalization.method = normalization_method_seurat,
        dims = usable_reduction_dims(seu, layer_spec$pca_dims, "pca"),
        verbose = FALSE
      )
      DefaultAssay(integrated) <- "integrated"
      if (!identical(normalization_method_seurat, "SCT")) {
        integrated <- ScaleData(integrated, verbose = FALSE)
      }
      reduction_name <- if (identical(reduction_method, "cca")) "integrated_cca" else "integrated_rpca"
      RunPCA(
        integrated,
        npcs = max(layer_spec$pca_dims),
        reduction.name = reduction_name,
        reduction.key = paste0(toupper(reduction_name), "_"),
        verbose = FALSE
      )
    }),
    error = function(e) list(value = seu, runtime_sec = 0, error = conditionMessage(e))
  )
  if (!is.null(timed$error)) {
    return(integrate_none(timed$value, method, layer_spec, normalization_method, timed$runtime_sec, paste0(method, "_failed:", timed$error)))
  }
  reduction_name <- if (identical(reduction_method, "cca")) "integrated_cca" else "integrated_rpca"
  integration_result(timed$value, method, reduction_name, sprintf("umap_%s_%s", reduction_name, normalization_method), timed$runtime_sec)
}

integrate_scanorama <- function(seu, group_var, layer_spec, normalization_method) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(integrate_none(seu, "scanorama", layer_spec, normalization_method, downgrade_reason = "reticulate_unavailable"))
  }
  scanorama <- tryCatch(reticulate::import("scanorama", delay_load = FALSE), error = function(e) NULL)
  if (is.null(scanorama)) {
    return(integrate_none(seu, "scanorama", layer_spec, normalization_method, downgrade_reason = "scanorama_unavailable"))
  }
  timed <- tryCatch(
    integration_runtime({
      assay_name <- if (identical(normalization_method, "sct") && "SCT" %in% Assays(seu)) "SCT" else "RNA"
      DefaultAssay(seu) <- assay_name
      features <- VariableFeatures(seu)
      features <- intersect(features, rownames(seu))
      if (length(features) < 2) {
        stop("too_few_features", call. = FALSE)
      }
      split_cells <- split(colnames(seu), as.character(seu@meta.data[[group_var]]))
      split_cells <- split_cells[vapply(split_cells, length, integer(1)) > 0]
      data_mat <- as.matrix(GetAssayData(seu, assay = assay_name, slot = "data")[features, , drop = FALSE])
      datasets <- lapply(split_cells, function(cells) t(data_mat[, cells, drop = FALSE]))
      gene_lists <- lapply(datasets, function(x) features)
      integrated <- scanorama$integrate(datasets, gene_lists, dimred = as.integer(max(50L, max(layer_spec$pca_dims))))
      embeddings <- matrix(NA_real_, nrow = ncol(seu), ncol = ncol(integrated[[1]][[1]]))
      rownames(embeddings) <- colnames(seu)
      for (idx in seq_along(split_cells)) {
        cells <- split_cells[[idx]]
        embeddings[cells, ] <- integrated[[1]][[idx]]
      }
      colnames(embeddings) <- paste0("SCANORAMA_", seq_len(ncol(embeddings)))
      seu[["scanorama"]] <- CreateDimReducObject(embeddings = embeddings, key = "SCANORAMA_", assay = assay_name)
      seu
    }),
    error = function(e) list(value = seu, runtime_sec = 0, error = conditionMessage(e))
  )
  if (!is.null(timed$error) || !"scanorama" %in% Reductions(timed$value)) {
    return(integrate_none(timed$value, "scanorama", layer_spec, normalization_method, timed$runtime_sec, paste0("scanorama_failed:", timed$error %||% "no_reduction")))
  }
  integration_result(timed$value, "scanorama", "scanorama", sprintf("umap_scanorama_%s", normalization_method), timed$runtime_sec)
}

INTEGRATION_METHOD_REGISTRY <- list(
  none = list(label = "Unintegrated", run = function(seu, group_var, layer_spec, normalization_method) {
    integrate_none(seu, "none", layer_spec, normalization_method)
  }),
  harmony = list(label = "Harmony", run = integrate_harmony),
  seurat_cca = list(label = "Seurat CCA", run = function(seu, group_var, layer_spec, normalization_method) {
    integrate_seurat_anchor(seu, group_var, layer_spec, normalization_method, "cca")
  }),
  seurat_rpca = list(label = "Seurat RPCA", run = function(seu, group_var, layer_spec, normalization_method) {
    integrate_seurat_anchor(seu, group_var, layer_spec, normalization_method, "rpca")
  }),
  scanorama = list(label = "Scanorama", run = integrate_scanorama)
)

dispatch_integration <- function(seu, method, group_var = "orig.ident", layer_spec, normalization_method) {
  method <- tolower(normalize_scalar_value(method))
  if (!method %in% names(INTEGRATION_METHOD_REGISTRY)) {
    return(integrate_none(seu, method, layer_spec, normalization_method, downgrade_reason = "unknown_method"))
  }
  INTEGRATION_METHOD_REGISTRY[[method]]$run(seu, group_var, layer_spec, normalization_method)
}
