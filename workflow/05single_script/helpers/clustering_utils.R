evaluate_resolution_local <- function(seu, res, seed, graph_name = NULL) {
  args <- list(
    object = seu,
    resolution = res,
    random.seed = seed,
    verbose = FALSE
  )
  if (!is.null(graph_name) && nzchar(graph_name) && "graphs" %in% slotNames(seu) && graph_name %in% names(seu@graphs)) {
    args$graph.name <- graph_name
  }
  tmp <- do.call(FindClusters, args)
  list(object = tmp, n_clusters = length(levels(factor(tmp$seurat_clusters))))
}

run_legacy_target_cluster_search <- function(seu, layer_spec, seed, graph_name = NULL, log_path = NULL) {
  search_records <- list()
  exact_obj <- NULL
  exact_res <- NA_real_
  exact_n <- NA_integer_

  for (res in sort(unique(layer_spec$res_range))) {
    eval_res <- evaluate_resolution_local(seu, res, seed, graph_name = graph_name)
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
      eval_res <- evaluate_resolution_local(seu, res, seed, graph_name = graph_name)
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
    fallback_eval <- evaluate_resolution_local(seu, selected_resolution, seed, graph_name = graph_name)
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
    selection_mode = "target_clusters",
    output_paths = list(resolution_search_tsv = log_path)
  )
}

run_fixed_resolution_cluster <- function(seu, layer_spec, seed, graph_name = NULL, log_path = NULL, decision_path = NULL) {
  selected_resolution <- approved_resolution_from_decision(decision_path %||% "", layer_spec)
  if (!is.finite(selected_resolution)) {
    stop(sprintf("layer `%s` fixed_resolution requires approved_resolution or approved decision file row", layer_spec$layer_id), call. = FALSE)
  }
  eval_res <- evaluate_resolution_local(seu, selected_resolution, seed, graph_name = graph_name)
  search_df <- data.frame(
    stage = "fixed_resolution",
    resolution = selected_resolution,
    n_clusters = eval_res$n_clusters,
    selected = TRUE,
    fallback_used = FALSE,
    stringsAsFactors = FALSE
  )
  if (!is.null(log_path) && nzchar(log_path)) {
    write_tsv_local(search_df, log_path)
  }
  list(
    seu = eval_res$object,
    selected_resolution = selected_resolution,
    selected_cluster_count = eval_res$n_clusters,
    search_table = search_df,
    fallback_used = FALSE,
    selection_mode = "fixed_resolution",
    output_paths = list(resolution_search_tsv = log_path)
  )
}

run_auto_resolution_recommendation <- function(seu, reduction_name, dims, layer_spec, seed, graph_name = NULL, log_path = NULL, decision_path = NULL) {
  out_dir <- if (!is.null(log_path) && nzchar(log_path)) dirname(log_path) else tempdir()
  output_paths <- list(
    resolution_search_tsv = log_path,
    resolution_metrics_tsv = file.path(out_dir, "resolution_metrics.tsv"),
    adjacent_ari_tsv = file.path(out_dir, "adjacent_ari.tsv"),
    seed_stability_tsv = file.path(out_dir, "seed_stability.tsv"),
    seed_ari_pairs_tsv = file.path(out_dir, "seed_ari_pairs.tsv"),
    cluster_size_metrics_tsv = file.path(out_dir, "cluster_size_metrics.tsv"),
    silhouette_metrics_tsv = file.path(out_dir, "silhouette_metrics.tsv"),
    ch_metrics_tsv = file.path(out_dir, "calinski_harabasz.tsv"),
    stability_plateaus_tsv = file.path(out_dir, "stability_plateaus.tsv"),
    resolution_ranking_tsv = file.path(out_dir, "resolution_ranking.tsv"),
    cluster_resolution_recommendation_tsv = file.path(out_dir, "cluster_resolution_recommendation.tsv"),
    cluster_resolution_decision_tsv = decision_path
  )

  grid <- run_resolution_grid(seu, layer_spec$res_range, seed, graph_name = graph_name)
  search_df <- grid$metrics
  search_df$selected <- FALSE
  search_df$fallback_used <- FALSE
  size_df <- compute_cluster_size_metrics(grid$labels, layer_spec$cluster_min_cells_abs, layer_spec$cluster_min_cell_fraction)
  adjacent_df <- compute_adjacent_ari(search_df, grid$labels)
  metrics_df <- attach_local_ari_metrics(search_df, adjacent_df, layer_spec$stable_local_ari)
  metrics_df <- merge(metrics_df, size_df, by = "resolution_label", all.x = TRUE, sort = FALSE)

  embedding <- Embeddings(seu, reduction = reduction_name)[, dims, drop = FALSE]
  silhouette_df <- compute_silhouette_metrics(embedding, grid$labels, layer_spec$cluster_silhouette_max_cells, seed)
  ch_df <- compute_ch_metrics(embedding, grid$labels)
  metrics_df <- merge(metrics_df, silhouette_df, by = "resolution_label", all.x = TRUE, sort = FALSE)
  metrics_df <- merge(metrics_df, ch_df, by = "resolution_label", all.x = TRUE, sort = FALSE)
  metrics_df <- metrics_df[order(metrics_df$resolution), , drop = FALSE]

  prelim <- metrics_df
  prelim$prelim_score <- 0.35 * pmin(1, ifelse(is.finite(prelim$local_ari_mean), prelim$local_ari_mean / layer_spec$stable_local_ari, 0)) +
    0.25 * cluster_rank01(prelim$silhouette_mean, high_is_good = TRUE) +
    0.20 * cluster_rank01(prelim$calinski_harabasz, high_is_good = TRUE) +
    0.20 * ifelse(is.finite(prelim$size_score), prelim$size_score, 0)
  prelim <- prelim[order(-prelim$prelim_score, prelim$resolution), , drop = FALSE]
  shortlist_n <- max(1L, as.integer(layer_spec$cluster_marker_shortlist_n))
  shortlist <- head(prelim$resolution, min(shortlist_n, nrow(prelim)))
  seeds <- unique(as.integer(layer_spec$cluster_seeds))
  seeds <- seeds[is.finite(seeds)]
  if (length(seeds) < 2L) {
    seeds <- unique(c(seed, seed + 1L))
  }
  seed_result <- compute_seed_stability(seu, shortlist, seeds, graph_name = graph_name)
  ranking_df <- rank_resolution_candidates(metrics_df, seed_result$summary, layer_spec)
  plateau_df <- detect_stability_plateaus(metrics_df, layer_spec$plateau_min_points)
  recommendation <- choose_resolution_recommendation(ranking_df, plateau_df, layer_spec$layer_id, layer_spec$score_tie_delta)

  if (!is.null(log_path) && nzchar(log_path)) {
    write_tsv_local(search_df, log_path)
  }
  write_tsv_local(metrics_df, output_paths$resolution_metrics_tsv)
  write_tsv_local(adjacent_df, output_paths$adjacent_ari_tsv)
  write_tsv_local(seed_result$summary, output_paths$seed_stability_tsv)
  write_tsv_local(seed_result$pairs, output_paths$seed_ari_pairs_tsv)
  write_tsv_local(size_df, output_paths$cluster_size_metrics_tsv)
  write_tsv_local(silhouette_df, output_paths$silhouette_metrics_tsv)
  write_tsv_local(ch_df, output_paths$ch_metrics_tsv)
  write_tsv_local(plateau_df, output_paths$stability_plateaus_tsv)
  write_tsv_local(ranking_df, output_paths$resolution_ranking_tsv)
  write_tsv_local(recommendation, output_paths$cluster_resolution_recommendation_tsv)
  if (!is.null(decision_path) && nzchar(decision_path)) {
    write_cluster_resolution_decision(decision_path, recommendation)
  }

  approved_resolution <- approved_resolution_from_decision(decision_path %||% "", layer_spec)
  if (!is.finite(approved_resolution)) {
    if (cluster_bool_env("CLUSTER_REQUIRE_APPROVAL", TRUE)) {
      stop(
        sprintf(
          "layer `%s` cluster resolution recommendation is pending approval. Review `%s`, then set status=approved and approved_resolution in `%s` before rerunning.",
          layer_spec$layer_id,
          output_paths$cluster_resolution_recommendation_tsv,
          decision_path
        ),
        call. = FALSE
      )
    }
    approved_resolution <- recommendation$recommended_resolution[[1]]
  }

  final_eval <- evaluate_resolution_local(seu, approved_resolution, seed, graph_name = graph_name)
  search_df$selected <- search_df$resolution == approved_resolution
  if (!any(search_df$selected)) {
    search_df <- rbind(
      search_df,
      data.frame(
        stage = "approved_resolution",
        resolution = approved_resolution,
        resolution_label = cluster_resolution_label(approved_resolution),
        n_clusters = final_eval$n_clusters,
        selected = TRUE,
        fallback_used = FALSE,
        stringsAsFactors = FALSE
      )
    )
  }
  if (!is.null(log_path) && nzchar(log_path)) {
    write_tsv_local(search_df, log_path)
  }
  selected_cluster_count <- final_eval$n_clusters
  list(
    seu = final_eval$object,
    selected_resolution = approved_resolution,
    selected_cluster_count = selected_cluster_count,
    search_table = search_df,
    fallback_used = FALSE,
    selection_mode = "auto_recommend",
    recommendation = recommendation,
    ranking_table = ranking_df,
    output_paths = output_paths
  )
}

run_resolution_search <- function(seu, reduction_name, layer_spec, seed, log_path = NULL, decision_path = NULL) {
  dims <- usable_reduction_dims(seu, layer_spec$pca_dims, reduction_name)
  graph_name <- "cluster_snn"
  seu <- FindNeighbors(
    seu,
    reduction = reduction_name,
    dims = dims,
    graph.name = c("cluster_nn", graph_name),
    verbose = FALSE
  )
  mode <- tolower(normalize_scalar_value(layer_spec$cluster_selection_mode, "target_clusters"))

  if (identical(mode, "auto_recommend")) {
    result <- run_auto_resolution_recommendation(
      seu,
      reduction_name = reduction_name,
      dims = dims,
      layer_spec = layer_spec,
      seed = seed,
      graph_name = graph_name,
      log_path = log_path,
      decision_path = decision_path
    )
  } else if (identical(mode, "fixed_resolution")) {
    result <- run_fixed_resolution_cluster(
      seu,
      layer_spec = layer_spec,
      seed = seed,
      graph_name = graph_name,
      log_path = log_path,
      decision_path = decision_path
    )
  } else {
    result <- run_legacy_target_cluster_search(
      seu,
      layer_spec = layer_spec,
      seed = seed,
      graph_name = graph_name,
      log_path = log_path
    )
  }

  result$dims <- dims
  result$graph_name <- graph_name
  result
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
