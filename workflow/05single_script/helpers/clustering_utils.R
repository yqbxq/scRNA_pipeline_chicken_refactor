evaluate_resolution_local <- function(seu, res, seed) {
  tmp <- FindClusters(
    seu,
    resolution = res,
    random.seed = seed,
    verbose = FALSE
  )
  list(object = tmp, n_clusters = length(levels(factor(tmp$seurat_clusters))))
}

adjusted_rand_index_local <- function(x, y) {
  x <- as.character(x)
  y <- as.character(y)
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (n < 2L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(NA_real_)
  }
  tab <- table(x, y)
  choose2 <- function(v) v * (v - 1) / 2
  sum_nij <- sum(choose2(as.numeric(tab)))
  sum_ai <- sum(choose2(rowSums(tab)))
  sum_bj <- sum(choose2(colSums(tab)))
  total <- choose2(n)
  expected <- sum_ai * sum_bj / total
  maximum <- 0.5 * (sum_ai + sum_bj)
  denom <- maximum - expected
  if (!is.finite(denom) || denom == 0) {
    return(NA_real_)
  }
  (sum_nij - expected) / denom
}

res_label_local <- function(res) {
  paste0("res_", gsub("\\.", "p", format(res, trim = TRUE, scientific = FALSE)))
}

candidate_cluster_col_local <- function(layer_id, res) {
  paste0(layer_id, "_cluster_", res_label_local(res))
}

safe_minmax_score_local <- function(x, higher_is_better = TRUE, default = 0.5) {
  out <- rep(default, length(x))
  finite <- is.finite(x)
  if (!any(finite)) {
    return(out)
  }
  rng <- range(x[finite], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    out[finite] <- default
    return(out)
  }
  out[finite] <- (x[finite] - rng[1]) / (rng[2] - rng[1])
  if (!higher_is_better) {
    out[finite] <- 1 - out[finite]
  }
  pmax(0, pmin(1, out))
}

score_near_target_local <- function(n_clusters, target_clusters) {
  if (!is.finite(target_clusters) || target_clusters <= 0) {
    return(rep(0.5, length(n_clusters)))
  }
  score <- 1 - abs(n_clusters - target_clusters) / pmax(target_clusters, n_clusters, 1)
  pmax(0, pmin(1, score))
}

sample_cells_for_metric_local <- function(cells, max_cells, seed) {
  cells <- as.character(cells)
  if (length(cells) <= max_cells) {
    return(cells)
  }
  set.seed(seed)
  sample(cells, max_cells)
}

calinski_harabasz_local <- function(embedding, cluster) {
  cluster <- as.factor(cluster)
  keep <- stats::complete.cases(embedding) & !is.na(cluster)
  embedding <- embedding[keep, , drop = FALSE]
  cluster <- droplevels(cluster[keep])
  n <- nrow(embedding)
  k <- nlevels(cluster)
  if (n <= k || k < 2L) {
    return(NA_real_)
  }
  overall <- colMeans(embedding)
  within <- 0
  between <- 0
  for (lev in levels(cluster)) {
    idx <- cluster == lev
    center <- colMeans(embedding[idx, , drop = FALSE])
    centered <- embedding[idx, , drop = FALSE] - matrix(center, nrow = sum(idx), ncol = ncol(embedding), byrow = TRUE)
    within <- within + sum(rowSums(centered^2))
    between <- between + sum(idx) * sum((center - overall)^2)
  }
  if (!is.finite(within) || within <= 0) {
    return(NA_real_)
  }
  (between / (k - 1)) / (within / (n - k))
}

mean_silhouette_local <- function(embedding, cluster) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    return(NA_real_)
  }
  cluster <- as.factor(cluster)
  keep <- stats::complete.cases(embedding) & !is.na(cluster)
  embedding <- embedding[keep, , drop = FALSE]
  cluster <- droplevels(cluster[keep])
  if (nrow(embedding) < 3L || nlevels(cluster) < 2L || nrow(embedding) <= nlevels(cluster)) {
    return(NA_real_)
  }
  out <- tryCatch(
    cluster::silhouette(as.integer(cluster), stats::dist(embedding)),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(NA_real_)
  }
  mean(out[, "sil_width"], na.rm = TRUE)
}

cluster_qc_eta2_local <- function(meta_df, cluster) {
  candidates <- intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent_mito", "percent.ribo", "percent_ribo", "S.Score", "G2M.Score"),
    colnames(meta_df)
  )
  if (length(candidates) == 0) {
    return(NA_real_)
  }
  cluster <- as.factor(cluster)
  eta_values <- vapply(candidates, function(col) {
    x <- suppressWarnings(as.numeric(meta_df[[col]]))
    keep <- is.finite(x) & !is.na(cluster)
    if (sum(keep) < 3L || length(unique(cluster[keep])) < 2L) {
      return(NA_real_)
    }
    grand <- mean(x[keep])
    groups <- split(x[keep], cluster[keep])
    ss_between <- sum(vapply(groups, function(v) length(v) * (mean(v) - grand)^2, numeric(1)))
    ss_total <- sum((x[keep] - grand)^2)
    if (!is.finite(ss_total) || ss_total <= 0) {
      return(NA_real_)
    }
    ss_between / ss_total
  }, numeric(1))
  eta_values <- eta_values[is.finite(eta_values)]
  if (length(eta_values) == 0) {
    return(NA_real_)
  }
  max(eta_values)
}

run_quality_resolution_search <- function(
  seu,
  reduction_name,
  layer_spec,
  seed,
  max_metric_cells = 2000L,
  seed_repeats = 3L,
  subsample_repeats = 2L,
  subsample_fraction = 0.8,
  candidate_top_n = 5L
) {
  dims <- usable_reduction_dims(seu, layer_spec$pca_dims, reduction_name)
  seu <- FindNeighbors(seu, reduction = reduction_name, dims = dims, verbose = FALSE)

  cells_all <- colnames(seu)
  metric_cells <- sample_cells_for_metric_local(cells_all, max_metric_cells, seed)
  embedding_all <- Embeddings(seu, reduction = reduction_name)[, dims, drop = FALSE]
  embedding_metric <- embedding_all[metric_cells, , drop = FALSE]
  meta_df <- seu@meta.data

  resolution_values <- sort(unique(layer_spec$res_range))
  if (length(resolution_values) == 0) {
    stop("resolution grid is empty.", call. = FALSE)
  }

  assignment_records <- list()
  metric_rows <- list()
  primary_assignments <- list()
  for (res in resolution_values) {
    eval_res <- evaluate_resolution_local(seu, res, seed)
    labels <- as.character(eval_res$object$seurat_clusters)
    names(labels) <- colnames(eval_res$object)
    col_name <- candidate_cluster_col_local(layer_spec$layer_id, res)
    seu[[col_name]] <- labels[colnames(seu)]
    primary_assignments[[as.character(res)]] <- labels

    sizes <- table(labels)
    small_cutoff <- max(20L, ceiling(0.005 * length(labels)))
    small_clusters <- sum(sizes < small_cutoff)
    metric_labels <- labels[metric_cells]

    seed_ari <- numeric(0)
    if (seed_repeats > 1L) {
      repeat_seeds <- seed + seq_len(seed_repeats)
      for (seed_i in repeat_seeds) {
        repeat_eval <- evaluate_resolution_local(seu, res, seed_i)
        seed_ari <- c(seed_ari, adjusted_rand_index_local(labels, as.character(repeat_eval$object$seurat_clusters)))
      }
    }

    subsample_ari <- numeric(0)
    if (subsample_repeats > 0L && subsample_fraction > 0 && subsample_fraction < 1 && length(cells_all) >= 50L) {
      n_sub <- max(20L, floor(length(cells_all) * subsample_fraction))
      n_sub <- min(length(cells_all), n_sub)
      for (rep_i in seq_len(subsample_repeats)) {
        set.seed(seed + 1000L + rep_i)
        sub_cells <- sample(cells_all, n_sub)
        sub_obj <- subset(seu, cells = sub_cells)
        sub_dims <- usable_reduction_dims(sub_obj, dims, reduction_name)
        sub_obj <- FindNeighbors(sub_obj, reduction = reduction_name, dims = sub_dims, verbose = FALSE)
        sub_eval <- evaluate_resolution_local(sub_obj, res, seed + 2000L + rep_i)
        subsample_ari <- c(
          subsample_ari,
          adjusted_rand_index_local(labels[sub_cells], as.character(sub_eval$object$seurat_clusters))
        )
      }
    }

    metric_rows[[length(metric_rows) + 1]] <- data.frame(
      resolution = res,
      candidate_id = res_label_local(res),
      candidate_cluster_col = col_name,
      n_clusters = length(sizes),
      min_cluster_size = min(as.integer(sizes)),
      median_cluster_size = stats::median(as.integer(sizes)),
      max_cluster_size = max(as.integer(sizes)),
      small_cluster_cutoff = small_cutoff,
      small_cluster_n = small_clusters,
      small_cluster_fraction = small_clusters / length(sizes),
      seed_stability_ari = if (length(seed_ari) > 0) mean(seed_ari, na.rm = TRUE) else NA_real_,
      subsample_stability_ari = if (length(subsample_ari) > 0) mean(subsample_ari, na.rm = TRUE) else NA_real_,
      silhouette_mean = mean_silhouette_local(embedding_metric, metric_labels),
      ch_index = calinski_harabasz_local(embedding_metric, metric_labels),
      max_qc_eta2 = cluster_qc_eta2_local(meta_df, labels),
      target_clusters = layer_spec$target_clusters,
      stringsAsFactors = FALSE
    )

    assignment_records[[length(assignment_records) + 1]] <- data.frame(
      cell_id = cells_all,
      resolution = res,
      candidate_id = res_label_local(res),
      candidate_cluster_col = col_name,
      cluster = labels[cells_all],
      stringsAsFactors = FALSE
    )
  }

  metrics <- dplyr::bind_rows(metric_rows)
  assignments <- dplyr::bind_rows(assignment_records)

  adjacent_rows <- list()
  if (length(resolution_values) >= 2L) {
    for (i in seq_len(length(resolution_values) - 1L)) {
      left <- as.character(resolution_values[[i]])
      right <- as.character(resolution_values[[i + 1L]])
      adjacent_rows[[length(adjacent_rows) + 1]] <- data.frame(
        resolution_left = resolution_values[[i]],
        resolution_right = resolution_values[[i + 1L]],
        ari = adjusted_rand_index_local(primary_assignments[[left]], primary_assignments[[right]]),
        stringsAsFactors = FALSE
      )
    }
  }
  adjacent_ari <- if (length(adjacent_rows) > 0) dplyr::bind_rows(adjacent_rows) else data.frame(
    resolution_left = numeric(0),
    resolution_right = numeric(0),
    ari = numeric(0),
    stringsAsFactors = FALSE
  )

  metrics$adjacent_ari_prev <- NA_real_
  metrics$adjacent_ari_next <- NA_real_
  for (i in seq_len(nrow(adjacent_ari))) {
    left <- adjacent_ari$resolution_left[[i]]
    right <- adjacent_ari$resolution_right[[i]]
    ari <- adjacent_ari$ari[[i]]
    metrics$adjacent_ari_next[metrics$resolution == left] <- ari
    metrics$adjacent_ari_prev[metrics$resolution == right] <- ari
  }
  metrics$adjacent_ari_mean <- rowMeans(metrics[, c("adjacent_ari_prev", "adjacent_ari_next"), drop = FALSE], na.rm = TRUE)
  metrics$adjacent_ari_mean[!is.finite(metrics$adjacent_ari_mean)] <- NA_real_

  metrics$seed_stability_score <- ifelse(is.finite(metrics$seed_stability_ari), pmax(0, pmin(1, metrics$seed_stability_ari)), 0.5)
  metrics$subsample_stability_score <- ifelse(is.finite(metrics$subsample_stability_ari), pmax(0, pmin(1, metrics$subsample_stability_ari)), 0.5)
  metrics$adjacent_stability_score <- ifelse(is.finite(metrics$adjacent_ari_mean), pmax(0, pmin(1, metrics$adjacent_ari_mean)), 0.5)
  metrics$stability_score <- rowMeans(metrics[, c("seed_stability_score", "subsample_stability_score", "adjacent_stability_score"), drop = FALSE], na.rm = TRUE)
  metrics$separation_score <- rowMeans(data.frame(
    silhouette = safe_minmax_score_local(metrics$silhouette_mean, higher_is_better = TRUE),
    ch = safe_minmax_score_local(metrics$ch_index, higher_is_better = TRUE)
  ), na.rm = TRUE)
  metrics$fragmentation_score <- rowMeans(data.frame(
    min_size = safe_minmax_score_local(metrics$min_cluster_size, higher_is_better = TRUE),
    small_fraction = safe_minmax_score_local(metrics$small_cluster_fraction, higher_is_better = FALSE)
  ), na.rm = TRUE)
  metrics$nuisance_score <- safe_minmax_score_local(metrics$max_qc_eta2, higher_is_better = FALSE)
  metrics$target_score <- score_near_target_local(metrics$n_clusters, metrics$target_clusters)
  for (score_col in c("stability_score", "separation_score", "fragmentation_score", "nuisance_score", "target_score")) {
    metrics[[score_col]][!is.finite(metrics[[score_col]])] <- 0.5
  }
  metrics$clustering_score <- 0.45 * metrics$stability_score +
    0.25 * metrics$separation_score +
    0.20 * metrics$fragmentation_score +
    0.05 * metrics$nuisance_score +
    0.05 * metrics$target_score
  metrics <- dplyr::arrange(metrics, dplyr::desc(clustering_score), dplyr::desc(stability_score), resolution)
  metrics$rank <- seq_len(nrow(metrics))
  metrics$selected_by_clustering_metrics <- metrics$rank == 1L
  metrics$marker_audit_candidate <- metrics$rank <= max(1L, candidate_top_n)

  selected <- metrics[metrics$rank == 1L, , drop = FALSE]
  list(
    seu = seu,
    selected_resolution = selected$resolution[[1]],
    selected_cluster_count = selected$n_clusters[[1]],
    selected_candidate_id = selected$candidate_id[[1]],
    selected_cluster_col = selected$candidate_cluster_col[[1]],
    resolution_metrics = metrics,
    resolution_ranking = metrics,
    adjacent_ari = adjacent_ari,
    cluster_assignments = assignments,
    dims = dims
  )
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

finalize_layer_object <- function(seu, reduction_name, umap_name, cluster_col, source_cluster_col = "seurat_clusters") {
  if (umap_name %in% Reductions(seu)) {
    seu[["umap"]] <- seu[[umap_name]]
  }
  if (!source_cluster_col %in% colnames(seu@meta.data)) {
    stop(sprintf("对象缺少 cluster 列 `%s`，无法 finalize。", source_cluster_col), call. = FALSE)
  }
  Idents(seu) <- source_cluster_col
  current_ids <- levels(factor(as.character(seu@meta.data[[source_cluster_col]])))
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
