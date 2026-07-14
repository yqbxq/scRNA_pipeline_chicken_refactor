cluster_resolution_label <- function(resolution) {
  gsub("\\.", "p", sprintf("%.3f", as.numeric(resolution)))
}

cluster_numeric_or <- function(value, default = NA_real_) {
  out <- suppressWarnings(as.numeric(normalize_scalar_value(value)))
  if (length(out) == 0 || is.na(out)) {
    return(default)
  }
  out
}

cluster_integer_or <- function(value, default = NA_integer_) {
  out <- suppressWarnings(as.integer(normalize_scalar_value(value)))
  if (length(out) == 0 || is.na(out)) {
    return(default)
  }
  out
}

cluster_bool_env <- function(name, default = TRUE) {
  value <- tolower(Sys.getenv(name, unset = ifelse(default, "yes", "no")))
  value %in% c("1", "true", "yes", "y")
}

cluster_require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("missing R package `%s`; install/update envs/environment_r_main.yml", pkg), call. = FALSE)
  }
}

cluster_ari <- function(x, y) {
  cluster_require_package("mclust")
  mclust::adjustedRandIndex(as.character(x), as.character(y))
}

cluster_rank01 <- function(x, high_is_good = TRUE) {
  x <- suppressWarnings(as.numeric(x))
  ok <- is.finite(x)
  out <- rep(0, length(x))
  if (sum(ok) == 0) {
    return(out)
  }
  values <- x[ok]
  if (length(unique(values)) == 1) {
    out[ok] <- 1
    return(out)
  }
  ranks <- rank(values, ties.method = "average", na.last = "keep")
  if (!high_is_good) {
    ranks <- max(ranks, na.rm = TRUE) + 1 - ranks
  }
  out[ok] <- (ranks - min(ranks, na.rm = TRUE)) / (max(ranks, na.rm = TRUE) - min(ranks, na.rm = TRUE))
  out
}

cluster_labels_from_object <- function(seu) {
  as.character(seu$seurat_clusters)
}

run_resolution_grid <- function(seu, resolutions, seed, graph_name = NULL) {
  rows <- list()
  labels <- list()
  resolutions <- sort(unique(as.numeric(resolutions[is.finite(resolutions)])))
  if (length(resolutions) == 0) {
    stop("resolution grid is empty", call. = FALSE)
  }

  for (res in resolutions) {
    eval_res <- evaluate_resolution_local(seu, res, seed, graph_name = graph_name)
    label_key <- cluster_resolution_label(res)
    current_labels <- cluster_labels_from_object(eval_res$object)
    labels[[label_key]] <- current_labels
    rows[[length(rows) + 1L]] <- data.frame(
      stage = "grid",
      resolution = res,
      resolution_label = label_key,
      n_clusters = eval_res$n_clusters,
      stringsAsFactors = FALSE
    )
  }

  list(
    metrics = dplyr::bind_rows(rows),
    labels = labels
  )
}

compute_cluster_size_metrics <- function(label_map, min_cells_abs, min_cell_fraction) {
  rows <- list()
  for (label_key in names(label_map)) {
    labels <- as.character(label_map[[label_key]])
    counts <- as.integer(table(labels))
    n_cells <- length(labels)
    min_cluster_cells <- max(as.integer(min_cells_abs), ceiling(n_cells * as.numeric(min_cell_fraction)))
    n_small <- sum(counts < min_cluster_cells)
    n_tiny <- sum(counts < 10L)
    cells_small <- sum(counts[counts < min_cluster_cells])
    rows[[length(rows) + 1L]] <- data.frame(
      resolution_label = label_key,
      min_cluster_cells_threshold = min_cluster_cells,
      min_cluster_cells = min(counts),
      n_small_clusters = n_small,
      n_tiny_clusters = n_tiny,
      fraction_small_clusters = n_small / length(counts),
      fraction_cells_in_small_clusters = cells_small / n_cells,
      size_score = ifelse(n_tiny > 0, 0, max(0, 1 - min(1, (cells_small / n_cells) / 0.10))),
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

compute_adjacent_ari <- function(metrics_df, label_map) {
  if (nrow(metrics_df) < 2) {
    return(data.frame(
      res_left = numeric(0),
      res_right = numeric(0),
      ncluster_left = integer(0),
      ncluster_right = integer(0),
      cluster_delta = integer(0),
      ari = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  metrics_df <- metrics_df[order(metrics_df$resolution), , drop = FALSE]
  rows <- list()
  for (idx in seq_len(nrow(metrics_df) - 1L)) {
    left <- metrics_df[idx, , drop = FALSE]
    right <- metrics_df[idx + 1L, , drop = FALSE]
    left_labels <- label_map[[left$resolution_label[[1]]]]
    right_labels <- label_map[[right$resolution_label[[1]]]]
    rows[[length(rows) + 1L]] <- data.frame(
      res_left = left$resolution[[1]],
      res_right = right$resolution[[1]],
      ncluster_left = left$n_clusters[[1]],
      ncluster_right = right$n_clusters[[1]],
      cluster_delta = abs(right$n_clusters[[1]] - left$n_clusters[[1]]),
      ari = cluster_ari(left_labels, right_labels),
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

attach_local_ari_metrics <- function(metrics_df, adjacent_df, stable_local_ari) {
  metrics_df$adjacent_ari_left <- NA_real_
  metrics_df$adjacent_ari_right <- NA_real_
  metrics_df$local_ari_mean <- NA_real_
  metrics_df$local_ari_min <- NA_real_
  metrics_df$neighbor_cluster_delta_max <- NA_integer_
  if (nrow(adjacent_df) == 0) {
    metrics_df$stable_point <- FALSE
    return(metrics_df)
  }

  for (idx in seq_len(nrow(metrics_df))) {
    res <- metrics_df$resolution[[idx]]
    left <- adjacent_df[adjacent_df$res_right == res, , drop = FALSE]
    right <- adjacent_df[adjacent_df$res_left == res, , drop = FALSE]
    ari_values <- c(left$ari, right$ari)
    delta_values <- c(left$cluster_delta, right$cluster_delta)
    ari_values <- ari_values[is.finite(ari_values)]
    delta_values <- delta_values[is.finite(delta_values)]
    if (nrow(left) > 0) {
      metrics_df$adjacent_ari_left[[idx]] <- left$ari[[1]]
    }
    if (nrow(right) > 0) {
      metrics_df$adjacent_ari_right[[idx]] <- right$ari[[1]]
    }
    if (length(ari_values) > 0) {
      metrics_df$local_ari_mean[[idx]] <- mean(ari_values)
      metrics_df$local_ari_min[[idx]] <- min(ari_values)
    }
    if (length(delta_values) > 0) {
      metrics_df$neighbor_cluster_delta_max[[idx]] <- max(delta_values)
    }
  }

  metrics_df$stable_point <- is.finite(metrics_df$local_ari_mean) &
    metrics_df$local_ari_mean >= stable_local_ari &
    is.finite(metrics_df$local_ari_min) &
    metrics_df$local_ari_min >= 0.85 &
    is.finite(metrics_df$neighbor_cluster_delta_max) &
    metrics_df$neighbor_cluster_delta_max <= 1
  metrics_df
}

cluster_sample_for_silhouette <- function(labels, max_cells) {
  labels <- as.character(labels)
  n <- length(labels)
  if (n <= max_cells) {
    return(seq_len(n))
  }
  groups <- split(seq_len(n), labels)
  k <- length(groups)
  base_take <- if (k * 50L <= max_cells) 50L else max(1L, floor(max_cells / k))
  selected <- unlist(lapply(groups, function(idx) sample(idx, min(length(idx), base_take))), use.names = FALSE)
  remaining_budget <- max_cells - length(selected)
  if (remaining_budget > 0) {
    pool <- setdiff(seq_len(n), selected)
    if (length(pool) > 0) {
      selected <- c(selected, sample(pool, min(length(pool), remaining_budget)))
    }
  }
  sort(unique(selected))
}

compute_silhouette_metrics <- function(embedding, label_map, max_cells, seed) {
  cluster_require_package("cluster")
  set.seed(seed)
  rows <- list()
  for (label_key in names(label_map)) {
    labels <- as.character(label_map[[label_key]])
    if (length(unique(labels)) < 2 || nrow(embedding) != length(labels)) {
      rows[[length(rows) + 1L]] <- data.frame(
        resolution_label = label_key,
        silhouette_n_cells = 0L,
        silhouette_mean = NA_real_,
        silhouette_negative_fraction = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    sampled <- cluster_sample_for_silhouette(labels, max_cells)
    sil <- cluster::silhouette(
      as.integer(factor(labels[sampled])),
      stats::dist(embedding[sampled, , drop = FALSE])
    )
    rows[[length(rows) + 1L]] <- data.frame(
      resolution_label = label_key,
      silhouette_n_cells = length(sampled),
      silhouette_mean = mean(sil[, "sil_width"], na.rm = TRUE),
      silhouette_negative_fraction = mean(sil[, "sil_width"] < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

compute_ch_index_one <- function(embedding, labels) {
  labels <- as.character(labels)
  n <- nrow(embedding)
  groups <- split(seq_len(n), labels)
  k <- length(groups)
  if (k < 2 || k >= n) {
    return(NA_real_)
  }
  overall <- colMeans(embedding)
  ss_between <- 0
  ss_within <- 0
  for (idx in groups) {
    group_embed <- embedding[idx, , drop = FALSE]
    center <- colMeans(group_embed)
    ss_between <- ss_between + length(idx) * sum((center - overall)^2)
    ss_within <- ss_within + sum(rowSums(sweep(group_embed, 2, center, "-")^2))
  }
  if (!is.finite(ss_within) || ss_within <= 0) {
    return(NA_real_)
  }
  (ss_between / (k - 1)) / (ss_within / (n - k))
}

compute_ch_metrics <- function(embedding, label_map) {
  rows <- list()
  for (label_key in names(label_map)) {
    rows[[length(rows) + 1L]] <- data.frame(
      resolution_label = label_key,
      calinski_harabasz = compute_ch_index_one(embedding, label_map[[label_key]]),
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

compute_seed_stability <- function(seu, resolutions, seeds, graph_name = NULL) {
  rows <- list()
  pair_rows <- list()
  for (res in resolutions) {
    seed_labels <- list()
    for (seed in seeds) {
      eval_res <- evaluate_resolution_local(seu, res, seed, graph_name = graph_name)
      seed_labels[[as.character(seed)]] <- cluster_labels_from_object(eval_res$object)
    }
    pairs <- utils::combn(names(seed_labels), 2, simplify = FALSE)
    ari_values <- vapply(pairs, function(pair) cluster_ari(seed_labels[[pair[[1]]]], seed_labels[[pair[[2]]]]), numeric(1))
    for (idx in seq_along(pairs)) {
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        resolution = res,
        seed_left = pairs[[idx]][[1]],
        seed_right = pairs[[idx]][[2]],
        ari = ari_values[[idx]],
        stringsAsFactors = FALSE
      )
    }
    rows[[length(rows) + 1L]] <- data.frame(
      resolution = res,
      resolution_label = cluster_resolution_label(res),
      seed_ari_median = median(ari_values, na.rm = TRUE),
      seed_ari_min = min(ari_values, na.rm = TRUE),
      seed_ari_p10 = as.numeric(stats::quantile(ari_values, 0.10, na.rm = TRUE, names = FALSE)),
      n_seed_pairs = length(ari_values),
      stringsAsFactors = FALSE
    )
  }
  list(summary = dplyr::bind_rows(rows), pairs = dplyr::bind_rows(pair_rows))
}

detect_stability_plateaus <- function(metrics_df, min_points) {
  stable <- metrics_df[isTRUE(metrics_df$stable_point) | metrics_df$stable_point, , drop = FALSE]
  if (nrow(stable) == 0) {
    return(data.frame(
      plateau_id = integer(0),
      plateau_start = numeric(0),
      plateau_end = numeric(0),
      n_points = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  stable <- stable[order(stable$resolution), , drop = FALSE]
  run_id <- cumsum(c(TRUE, diff(match(stable$resolution, metrics_df$resolution)) != 1))
  rows <- list()
  for (id in unique(run_id)) {
    block <- stable[run_id == id, , drop = FALSE]
    if (nrow(block) >= min_points) {
      rows[[length(rows) + 1L]] <- data.frame(
        plateau_id = length(rows) + 1L,
        plateau_start = min(block$resolution),
        plateau_end = max(block$resolution),
        n_points = nrow(block),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) {
    return(data.frame(
      plateau_id = integer(0),
      plateau_start = numeric(0),
      plateau_end = numeric(0),
      n_points = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(rows)
}

rank_resolution_candidates <- function(metrics_df, seed_df, layer_spec) {
  df <- merge(metrics_df, seed_df, by = c("resolution", "resolution_label"), all.x = TRUE, sort = FALSE)
  df$local_ari_score <- pmin(1, ifelse(is.finite(df$local_ari_mean), df$local_ari_mean / layer_spec$stable_local_ari, 0))
  df$seed_ari_score <- 0.7 * pmin(1, ifelse(is.finite(df$seed_ari_median), df$seed_ari_median / 0.90, 0)) +
    0.3 * pmin(1, ifelse(is.finite(df$seed_ari_min), df$seed_ari_min / 0.80, 0))
  df$silhouette_score <- cluster_rank01(df$silhouette_mean, high_is_good = TRUE)
  df$ch_score <- cluster_rank01(df$calinski_harabasz, high_is_good = TRUE)
  df$cluster_size_score <- ifelse(is.finite(df$size_score), df$size_score, 0)
  df$final_score <- 0.30 * df$seed_ari_score +
    0.20 * df$local_ari_score +
    0.20 * df$silhouette_score +
    0.15 * df$ch_score +
    0.15 * df$cluster_size_score
  df$passes_basic_filters <- df$n_clusters >= layer_spec$min_expected_clusters &
    df$n_clusters <= layer_spec$max_expected_clusters &
    df$n_tiny_clusters == 0 &
    df$fraction_cells_in_small_clusters <= 0.10 &
    ifelse(is.finite(df$seed_ari_median), df$seed_ari_median >= layer_spec$min_seed_ari, FALSE)
  df[order(-df$final_score, df$resolution), , drop = FALSE]
}

choose_resolution_recommendation <- function(ranking_df, plateau_df, layer_id, tie_delta) {
  eligible <- ranking_df[ranking_df$passes_basic_filters, , drop = FALSE]
  status <- "recommended"
  reason <- ""
  if (nrow(eligible) == 0) {
    eligible <- ranking_df
    status <- "recommended_with_warnings"
    reason <- "no candidate passed all basic filters"
  }
  if (nrow(plateau_df) > 0) {
    in_plateau <- Reduce(`|`, lapply(seq_len(nrow(plateau_df)), function(idx) {
      eligible$resolution >= plateau_df$plateau_start[[idx]] & eligible$resolution <= plateau_df$plateau_end[[idx]]
    }))
    plateau_eligible <- eligible[in_plateau, , drop = FALSE]
    if (nrow(plateau_eligible) > 0) {
      eligible <- plateau_eligible
    }
  }
  best_score <- max(eligible$final_score, na.rm = TRUE)
  tied <- eligible[eligible$final_score >= best_score - tie_delta, , drop = FALSE]
  chosen <- tied[order(tied$resolution), , drop = FALSE][1, , drop = FALSE]
  containing <- plateau_df[chosen$resolution[[1]] >= plateau_df$plateau_start & chosen$resolution[[1]] <= plateau_df$plateau_end, , drop = FALSE]
  data.frame(
    layer_id = layer_id,
    recommended_resolution = chosen$resolution[[1]],
    recommended_clusters = chosen$n_clusters[[1]],
    score = chosen$final_score[[1]],
    plateau_start = ifelse(nrow(containing) > 0, containing$plateau_start[[1]], NA_real_),
    plateau_end = ifelse(nrow(containing) > 0, containing$plateau_end[[1]], NA_real_),
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

read_cluster_resolution_decision <- function(path, layer_id) {
  df <- read_tsv_optional(path)
  if (nrow(df) == 0 || !"layer_id" %in% colnames(df)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  df[df$layer_id == layer_id, , drop = FALSE]
}

write_cluster_resolution_decision <- function(path, recommendation) {
  existing <- read_tsv_optional(path)
  base <- data.frame(
    layer_id = recommendation$layer_id,
    recommended_resolution = as.character(recommendation$recommended_resolution),
    approved_resolution = "",
    status = "pending",
    note = "",
    recommended_clusters = as.character(recommendation$recommended_clusters),
    score = as.character(signif(recommendation$score, 6)),
    plateau_start = as.character(recommendation$plateau_start),
    plateau_end = as.character(recommendation$plateau_end),
    updated_at = timestamp_now(),
    stringsAsFactors = FALSE
  )
  if (nrow(existing) > 0) {
    for (col in setdiff(colnames(base), colnames(existing))) {
      existing[[col]] <- ""
    }
    hit <- existing$layer_id == base$layer_id[[1]]
    if (any(hit)) {
      prior <- existing[hit, , drop = FALSE][1, , drop = FALSE]
      if (tolower(normalize_scalar_value(prior$status[[1]])) == "approved" && nzchar(normalize_scalar_value(prior$approved_resolution[[1]]))) {
        base$approved_resolution <- prior$approved_resolution[[1]]
        base$status <- prior$status[[1]]
        base$note <- prior$note[[1]]
      }
      existing <- existing[!hit, , drop = FALSE]
    }
    for (col in setdiff(colnames(existing), colnames(base))) {
      base[[col]] <- ""
    }
    base <- base[, colnames(existing), drop = FALSE]
    out <- rbind(existing, base)
  } else {
    out <- base
  }
  write_tsv_local(out, path)
  invisible(out)
}

approved_resolution_from_decision <- function(path, layer_spec) {
  decision <- read_cluster_resolution_decision(path, layer_spec$layer_id)
  if (nrow(decision) > 0) {
    status <- tolower(normalize_scalar_value(decision$status[[1]]))
    approved <- cluster_numeric_or(decision$approved_resolution[[1]], NA_real_)
    if (identical(status, "approved") && is.finite(approved)) {
      return(approved)
    }
  }
  if (is.finite(layer_spec$approved_resolution)) {
    return(layer_spec$approved_resolution)
  }
  NA_real_
}
