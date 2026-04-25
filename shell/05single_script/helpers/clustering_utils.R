evaluate_resolution_local <- function(seu, res, seed) {
  tmp <- FindClusters(
    seu,
    resolution = res,
    random.seed = seed,
    verbose = FALSE
  )
  list(object = tmp, n_clusters = length(levels(factor(tmp$seurat_clusters))))
}

run_resolution_search <- function(seu, reduction_name, layer_spec, seed, log_path = NULL) {
  dims <- usable_reduction_dims(seu, layer_spec$pca_dims, reduction_name)
  seu <- FindNeighbors(seu, reduction = reduction_name, dims = dims, verbose = FALSE)

  search_records <- list()
  best_obj <- NULL
  exact_obj <- NULL
  exact_res <- NA_real_
  exact_n <- NA_integer_

  for (res in sort(unique(layer_spec$res_range))) {
    eval_res <- evaluate_resolution_local(seu, res, seed)
    search_records[[length(search_records) + 1]] <- data.frame(
      stage = "coarse",
      resolution = res,
      n_clusters = eval_res$n_clusters,
      stringsAsFactors = FALSE
    )
    if (eval_res$n_clusters == layer_spec$target_clusters) {
      exact_obj <- eval_res$object
      exact_res <- res
      exact_n <- eval_res$n_clusters
      break
    }
  }

  if (is.na(exact_res)) {
    coarse_df <- do.call(rbind, search_records)
    best_idx <- which.min(abs(coarse_df$n_clusters - layer_spec$target_clusters))
    best_res <- coarse_df$resolution[best_idx]
    fine_grid <- seq(
      max(0, best_res - 2 * layer_spec$res_fine_step),
      best_res + 2 * layer_spec$res_fine_step,
      by = layer_spec$res_fine_step
    )
    fine_grid <- sort(unique(round(fine_grid, 6)))
    fine_grid <- setdiff(fine_grid, coarse_df$resolution)
    for (res in fine_grid) {
      eval_res <- evaluate_resolution_local(seu, res, seed)
      search_records[[length(search_records) + 1]] <- data.frame(
        stage = "fine",
        resolution = res,
        n_clusters = eval_res$n_clusters,
        stringsAsFactors = FALSE
      )
      if (eval_res$n_clusters == layer_spec$target_clusters) {
        exact_obj <- eval_res$object
        exact_res <- res
        exact_n <- eval_res$n_clusters
        break
      }
    }
  }

  search_df <- do.call(rbind, search_records)
  search_df <- dplyr::distinct(search_df, stage, resolution, .keep_all = TRUE)
  search_df <- dplyr::arrange(search_df, resolution, stage)

  fallback_used <- FALSE
  if (!is.na(exact_res)) {
    selected_resolution <- exact_res
    selected_cluster_count <- exact_n
    selected_obj <- exact_obj
  } else {
    fallback_used <- TRUE
    best_idx <- which.min(abs(search_df$n_clusters - layer_spec$target_clusters))
    selected_resolution <- search_df$resolution[best_idx]
    fallback_eval <- evaluate_resolution_local(seu, selected_resolution, seed)
    selected_obj <- fallback_eval$object
    selected_cluster_count <- fallback_eval$n_clusters
    search_df <- rbind(
      search_df,
      data.frame(stage = "fallback", resolution = selected_resolution, n_clusters = selected_cluster_count, stringsAsFactors = FALSE)
    )
  }

  search_df$selected <- search_df$resolution == selected_resolution & search_df$n_clusters == selected_cluster_count
  search_df$fallback_used <- fallback_used
  if (!is.null(log_path) && nzchar(log_path)) {
    write_tsv_local(search_df, log_path)
  }

  list(
    seu = selected_obj,
    selected_resolution = selected_resolution,
    selected_cluster_count = selected_cluster_count,
    search_table = search_df,
    fallback_used = fallback_used,
    dims = dims
  )
}

finalize_layer_object <- function(seu, reduction_name, umap_name, cluster_col) {
  if (umap_name %in% Reductions(seu)) {
    seu[["umap"]] <- seu[[umap_name]]
  }
  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    stop("对象缺少 seurat_clusters，无法 finalize。", call. = FALSE)
  }
  Idents(seu) <- "seurat_clusters"
  current_ids <- levels(factor(as.character(seu$seurat_clusters)))
  new_ids <- setNames(as.character(seq_along(current_ids)), current_ids)
  seu <- RenameIdents(seu, new_ids)
  renumbered <- as.character(Idents(seu))
  seu$seurat_clusters <- factor(renumbered, levels = as.character(seq_along(current_ids)))
  seu[[cluster_col]] <- as.character(seu$seurat_clusters)
  Idents(seu) <- cluster_col
  seu@misc$selected_reduction <- reduction_name
  seu@misc$selected_umap <- umap_name
  seu
}

build_cluster_summary_local <- function(seu, cluster_col, group_var = NULL) {
  meta_df <- tibble::rownames_to_column(seu@meta.data, "cell_id")
  if (is.null(group_var) || !nzchar(group_var) || !group_var %in% colnames(meta_df)) {
    group_var <- if ("sample_id" %in% colnames(meta_df)) "sample_id" else "orig.ident"
  }
  meta_df$plot_group <- as.character(meta_df[[group_var]])
  meta_df$cluster_value <- as.character(meta_df[[cluster_col]])
  out <- dplyr::count(meta_df, plot_group, cluster_value, name = "cell_number")
  out <- dplyr::group_by(out, plot_group)
  out <- dplyr::mutate(out, proportion = cell_number / sum(cell_number))
  out <- dplyr::ungroup(out)
  colnames(out)[colnames(out) == "cluster_value"] <- cluster_col
  out
}

upsert_layer_status <- function(path, row) {
  existing <- read_tsv_optional(path)
  if (nrow(existing) == 0) {
    out <- row
  } else {
    for (col in setdiff(colnames(row), colnames(existing))) {
      existing[[col]] <- ""
    }
    for (col in setdiff(colnames(existing), colnames(row))) {
      row[[col]] <- ""
    }
    row <- row[, colnames(existing), drop = FALSE]
    existing <- existing[existing$layer_id != row$layer_id[[1]], , drop = FALSE]
    out <- rbind(existing, row)
  }
  write_tsv_local(out, path)
  out
}
