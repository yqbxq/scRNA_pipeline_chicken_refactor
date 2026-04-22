format_index_spec <- function(x) {
  x <- as.integer(x[is.finite(x)])
  if (length(x) == 0) {
    return("")
  }
  x <- sort(unique(x))
  if (length(x) >= 2 && identical(x, seq.int(min(x), max(x)))) {
    return(sprintf("%s:%s", min(x), max(x)))
  }
  paste(x, collapse = ",")
}

normalize_yes_no <- function(x, default = "yes") {
  value <- tolower(normalize_scalar_value(x, default))
  if (!value %in% c("yes", "no")) {
    value <- default
  }
  value
}

parse_integer_scalar <- function(x, default) {
  value <- suppressWarnings(as.integer(trimws(as.character(x))))
  if (length(value) == 0 || is.na(value)) {
    return(as.integer(default))
  }
  value
}

parse_numeric_scalar <- function(x, default) {
  value <- suppressWarnings(as.numeric(trimws(as.character(x))))
  if (length(value) == 0 || is.na(value)) {
    return(as.numeric(default))
  }
  value
}

parse_numeric_csv <- function(x, default) {
  values <- suppressWarnings(as.numeric(split_csv(as.character(x))))
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(as.numeric(default))
  }
  values
}

sanitize_layer_id <- function(x) {
  value <- trimws(as.character(x))
  if (!nzchar(value)) {
    return("")
  }
  gsub("[^A-Za-z0-9._-]+", "_", value)
}

normalize_integration_mode <- function(x, default = "harmony") {
  value <- tolower(normalize_scalar_value(x, default))
  if (!value %in% c("harmony", "none")) {
    value <- default
  }
  value
}

default_object_layer_table <- function(cfg) {
  data.frame(
    layer_id = c("panorama", "subcluster_1", "subcluster_2"),
    layer_role = c("panorama", "subcluster", "subcluster"),
    enabled = c("yes", "yes", "yes"),
    parent_layer = c("", "panorama", "panorama"),
    sample_include = c("", "", ""),
    sample_exclude = c("", "", ""),
    selection_column = c("", "panorama_cluster", "panorama_cluster"),
    selection_values = c("", "", ""),
    rebuild_normalization = c("yes", "yes", "yes"),
    hvg_nfeatures = c(
      cfg$hvg_nfeatures,
      cfg$hvg_nfeatures,
      cfg$hvg_nfeatures
    ),
    pca_dims = c(
      format_index_spec(cfg$pca_dims),
      "1:20",
      "1:20"
    ),
    target_clusters = c(
      cfg$target_clusters,
      8L,
      6L
    ),
    res_range = c(
      paste(cfg$res_range, collapse = ","),
      "0.10,0.15,0.20,0.25,0.30,0.35,0.40",
      "0.10,0.15,0.20,0.25,0.30,0.35,0.40"
    ),
    res_fine_step = c(
      cfg$res_fine_step,
      cfg$res_fine_step,
      cfg$res_fine_step
    ),
    integration_mode = c(
      cfg$integration_mode,
      cfg$integration_mode,
      cfg$integration_mode
    ),
    description = c(
      "Root panorama object built from all post-QC cells.",
      "Edit selection_column and selection_values to define this subcluster from panorama results.",
      "Edit selection_column and selection_values to define this subcluster from panorama results."
    ),
    stringsAsFactors = FALSE
  )
}

read_object_layer_specs <- function(cfg) {
  df <- read_tsv_optional(cfg$object_layer_config_file)
  if (nrow(df) == 0) {
    df <- default_object_layer_table(cfg)
  }

  expected_cols <- c(
    "layer_id", "layer_role", "enabled", "parent_layer",
    "sample_include", "sample_exclude",
    "selection_column", "selection_values",
    "rebuild_normalization",
    "hvg_nfeatures", "pca_dims",
    "target_clusters", "res_range", "res_fine_step",
    "integration_mode", "description"
  )
  for (col in expected_cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }

  df$layer_id <- vapply(df$layer_id, sanitize_layer_id, character(1))
  df <- df[nzchar(df$layer_id), expected_cols, drop = FALSE]
  if (nrow(df) == 0) {
    stop("对象层级配置为空，至少需要一行 panorama 配置", call. = FALSE)
  }

  panorama_hits <- which(tolower(df$layer_role) == "panorama")
  if (length(panorama_hits) == 0) {
    panorama_hits <- which(!nzchar(trimws(df$parent_layer)))
  }
  if (length(panorama_hits) == 0) {
    default_panorama <- default_object_layer_table(cfg)[1, , drop = FALSE]
    df <- rbind(default_panorama, df)
    panorama_hits <- 1L
  }
  if (length(panorama_hits) != 1) {
    stop("对象层级配置必须且只能定义一个 panorama 根层", call. = FALSE)
  }

  panorama_id <- df$layer_id[panorama_hits[1]]

  specs <- lapply(seq_len(nrow(df)), function(idx) {
    row <- df[idx, , drop = FALSE]
    layer_id <- row$layer_id[[1]]
    layer_role <- tolower(normalize_scalar_value(row$layer_role[[1]], if (identical(layer_id, panorama_id)) "panorama" else "subcluster"))
    parent_layer <- sanitize_layer_id(row$parent_layer[[1]])
    if (identical(layer_id, panorama_id)) {
      layer_role <- "panorama"
      parent_layer <- ""
    } else if (!nzchar(parent_layer)) {
      parent_layer <- panorama_id
    }

    pca_dims_raw <- normalize_scalar_value(row$pca_dims[[1]], format_index_spec(cfg$pca_dims))
    selection_column <- normalize_scalar_value(row$selection_column[[1]], "")
    selection_values <- split_csv(row$selection_values[[1]])

    list(
      layer_id = layer_id,
      layer_role = layer_role,
      enabled = !identical(normalize_yes_no(row$enabled[[1]], "yes"), "no"),
      parent_layer = parent_layer,
      sample_include = split_csv(row$sample_include[[1]]),
      sample_exclude = split_csv(row$sample_exclude[[1]]),
      selection_column = selection_column,
      selection_values = selection_values,
      rebuild_normalization = !identical(normalize_yes_no(row$rebuild_normalization[[1]], "yes"), "no"),
      hvg_nfeatures = parse_integer_scalar(row$hvg_nfeatures[[1]], cfg$hvg_nfeatures),
      pca_dims = parse_index_spec(pca_dims_raw),
      target_clusters = parse_integer_scalar(row$target_clusters[[1]], cfg$target_clusters),
      res_range = parse_numeric_csv(row$res_range[[1]], cfg$res_range),
      res_fine_step = parse_numeric_scalar(row$res_fine_step[[1]], cfg$res_fine_step),
      integration_mode = normalize_integration_mode(row$integration_mode[[1]], cfg$integration_mode),
      description = normalize_scalar_value(row$description[[1]], "")
    )
  })

  names(specs) <- vapply(specs, `[[`, character(1), "layer_id")
  specs
}

panorama_layer_id <- function(layer_specs) {
  hits <- names(layer_specs)[vapply(layer_specs, function(spec) identical(spec$layer_role, "panorama"), logical(1))]
  if (length(hits) != 1) {
    stop("对象层级配置必须且只能定义一个 panorama 根层", call. = FALSE)
  }
  hits[[1]]
}

layer_cluster_column <- function(layer_id) {
  paste0(sanitize_layer_id(layer_id), "_cluster")
}

object_layer_paths <- function(cfg, layer_id) {
  safe_id <- sanitize_layer_id(layer_id)
  checkpoint_dir <- file.path(cfg$checkpoint_dir, "layers", safe_id)
  figure_dir <- file.path(cfg$figure_dir, "layers", safe_id)
  table_dir <- file.path(cfg$table_dir, "layers", safe_id)
  report_dir <- file.path(cfg$eda_report_dir, "integration", safe_id)

  list(
    layer_id = safe_id,
    checkpoint_dir = checkpoint_dir,
    figure_dir = figure_dir,
    table_dir = table_dir,
    report_dir = report_dir,
    candidate_rds = file.path(checkpoint_dir, sprintf("%s_reduction_candidates.rds", safe_id)),
    clustered_rds = file.path(checkpoint_dir, sprintf("%s_after_clustering.rds", safe_id)),
    cluster_summary_csv = file.path(table_dir, "cluster_summary.csv"),
    resolution_search_csv = file.path(table_dir, "resolution_search.csv"),
    selected_resolution_txt = file.path(table_dir, "selected_resolution.txt"),
    pca_driver_summary_tsv = file.path(report_dir, "pca_driver_summary.tsv"),
    pc_association_tsv = file.path(report_dir, "pc_association.tsv"),
    resolution_candidates_tsv = file.path(report_dir, "resolution_candidates.tsv"),
    triage_tsv = file.path(report_dir, "integration_triage.tsv"),
    report_md = file.path(report_dir, "report.md"),
    resolution_plot_png = file.path(report_dir, "resolution_candidates.png"),
    umap_rna_by_sample_png = file.path(report_dir, "umap_rna_by_sample.png"),
    umap_harmony_by_sample_png = file.path(report_dir, "umap_harmony_by_sample.png"),
    pca_elbow_png = file.path(report_dir, "pca_elbow.png"),
    cluster_umap_png = file.path(figure_dir, sprintf("%s_cluster_umap.png", safe_id)),
    cluster_proportion_png = file.path(figure_dir, sprintf("%s_cluster_proportion.png", safe_id)),
    cluster_overview_png = file.path(figure_dir, sprintf("%s_cluster_overview.png", safe_id))
  )
}

ensure_object_layer_dirs <- function(paths) {
  dirs <- unique(unname(paths[c("checkpoint_dir", "figure_dir", "table_dir", "report_dir")]))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

apply_layer_sample_filters <- function(seu, spec) {
  meta <- seu@meta.data
  if (!"orig.ident" %in% colnames(meta)) {
    return(list(status = "invalid", reason = "对象缺少 orig.ident，无法执行 sample_include/sample_exclude 过滤", object = NULL))
  }

  keep <- rep(TRUE, nrow(meta))
  if (length(spec$sample_include) > 0) {
    keep <- keep & as.character(meta$orig.ident) %in% spec$sample_include
  }
  if (length(spec$sample_exclude) > 0) {
    keep <- keep & !as.character(meta$orig.ident) %in% spec$sample_exclude
  }

  selected_cells <- rownames(meta)[keep]
  if (length(selected_cells) == 0) {
    return(list(status = "skipped", reason = "sample_include/sample_exclude 过滤后没有剩余细胞", object = NULL))
  }

  list(status = "ready", reason = "", object = subset(seu, cells = selected_cells))
}

prepare_parent_object_for_layer <- function(parent_obj, spec) {
  filtered <- apply_layer_sample_filters(parent_obj, spec)
  if (!identical(filtered$status, "ready")) {
    return(filtered)
  }

  obj <- filtered$object
  if (identical(spec$layer_role, "panorama")) {
    return(list(status = "ready", reason = "", object = obj))
  }

  if (!nzchar(spec$selection_column)) {
    return(list(status = "skipped", reason = "selection_column 为空，无法从上游对象提取子对象", object = NULL))
  }
  if (length(spec$selection_values) == 0) {
    return(list(status = "skipped", reason = "selection_values 为空，请先在 object_layers.tsv 中填写上游选择规则", object = NULL))
  }
  if (!spec$selection_column %in% colnames(obj@meta.data)) {
    return(list(status = "skipped", reason = sprintf("上游对象缺少选择列: %s", spec$selection_column), object = NULL))
  }

  column_values <- as.character(obj@meta.data[[spec$selection_column]])
  keep <- column_values %in% spec$selection_values
  selected_cells <- rownames(obj@meta.data)[keep]
  if (length(selected_cells) == 0) {
    return(list(
      status = "skipped",
      reason = sprintf(
        "列 %s 中没有命中 selection_values=%s 的细胞",
        spec$selection_column,
        paste(spec$selection_values, collapse = ",")
      ),
      object = NULL
    ))
  }

  list(status = "ready", reason = "", object = subset(obj, cells = selected_cells))
}

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
  seu
}

usable_reduction_dims <- function(seu, requested_dims, reduction_name) {
  embed <- Embeddings(seu, reduction = reduction_name)
  max_dims <- ncol(embed)
  dims <- requested_dims[requested_dims <= max_dims]
  if (length(dims) == 0) {
    dims <- seq_len(max_dims)
  }
  dims
}

build_layer_reduction_candidates <- function(seu, spec, cfg) {
  seu <- strip_reduction_state(seu)
  DefaultAssay(seu) <- "RNA"

  if (spec$rebuild_normalization) {
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  }
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = spec$hvg_nfeatures, verbose = FALSE)
  seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
  seu <- RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE)

  pca_dims <- usable_reduction_dims(seu, spec$pca_dims, "pca")
  seu <- RunUMAP(
    seu,
    reduction = "pca",
    dims = pca_dims,
    reduction.name = "umap_rna",
    reduction.key = "UMAPRNA_",
    verbose = FALSE
  )

  if (identical(spec$integration_mode, "harmony") && length(unique(seu$orig.ident)) > 1) {
    seu <- harmony::RunHarmony(
      seu,
      group.by.vars = "orig.ident",
      theta = 3,
      lambda = 1,
      plot_convergence = FALSE
    )
    harmony_dims <- usable_reduction_dims(seu, spec$pca_dims, "harmony")
    seu <- RunUMAP(
      seu,
      reduction = "harmony",
      dims = harmony_dims,
      reduction.name = "umap_harmony",
      reduction.key = "UMAPHARM_",
      verbose = FALSE
    )
  }

  seu@misc$layer_workflow <- list(
    layer_id = spec$layer_id,
    layer_role = spec$layer_role,
    parent_layer = spec$parent_layer,
    cluster_column = layer_cluster_column(spec$layer_id),
    integration_mode = spec$integration_mode,
    hvg_nfeatures = spec$hvg_nfeatures,
    pca_dims = spec$pca_dims,
    target_clusters = spec$target_clusters,
    res_range = spec$res_range,
    res_fine_step = spec$res_fine_step
  )

  seu
}

selected_layer_reduction <- function(seu, spec) {
  if (identical(spec$integration_mode, "none") || !"harmony" %in% Reductions(seu)) {
    return("pca")
  }
  "harmony"
}

selected_layer_umap <- function(seu, reduction_name) {
  if (identical(reduction_name, "harmony") && "umap_harmony" %in% Reductions(seu)) {
    return("umap_harmony")
  }
  "umap_rna"
}

evaluate_layer_resolution <- function(seu, reduction_name, dims, res, seed) {
  tmp <- FindNeighbors(
    seu,
    reduction = reduction_name,
    dims = dims,
    verbose = FALSE
  )
  tmp <- FindClusters(
    tmp,
    resolution = res,
    random.seed = seed,
    verbose = FALSE
  )
  list(object = tmp, n_clusters = length(levels(tmp$seurat_clusters)))
}

run_layer_resolution_search <- function(seu, spec, cfg) {
  reduction_name <- selected_layer_reduction(seu, spec)
  dims <- usable_reduction_dims(seu, spec$pca_dims, reduction_name)

  search_records <- list()
  evaluated_res <- numeric(0)
  found_res <- NA_real_
  selected_obj <- NULL
  selected_cluster_count <- NA_integer_

  for (res in spec$res_range) {
    eval_res <- evaluate_layer_resolution(seu, reduction_name, dims, res, cfg$random_seed)
    search_records[[length(search_records) + 1]] <- data.frame(
      stage = "coarse",
      resolution = res,
      n_clusters = eval_res$n_clusters,
      stringsAsFactors = FALSE
    )
    evaluated_res <- c(evaluated_res, res)

    if (eval_res$n_clusters == spec$target_clusters) {
      found_res <- res
      selected_obj <- eval_res$object
      selected_cluster_count <- eval_res$n_clusters
      break
    }
  }

  if (is.na(found_res)) {
    coarse_df <- do.call(rbind, search_records)
    interval_ids <- integer(0)

    if (nrow(coarse_df) >= 2) {
      for (i in seq_len(nrow(coarse_df) - 1)) {
        left_n <- coarse_df$n_clusters[i]
        right_n <- coarse_df$n_clusters[i + 1]
        lower_n <- min(left_n, right_n)
        upper_n <- max(left_n, right_n)
        if (spec$target_clusters >= lower_n && spec$target_clusters <= upper_n) {
          interval_ids <- c(interval_ids, i)
        }
      }
    }

    if (length(interval_ids) == 0 && nrow(coarse_df) >= 2) {
      closest_idx <- which.min(abs(coarse_df$n_clusters - spec$target_clusters))
      interval_ids <- if (closest_idx == nrow(coarse_df)) closest_idx - 1L else closest_idx
    }

    for (interval_idx in unique(interval_ids)) {
      left_res <- coarse_df$resolution[interval_idx]
      right_res <- coarse_df$resolution[interval_idx + 1]
      fine_grid <- seq(left_res, right_res, by = spec$res_fine_step)
      fine_grid <- sort(unique(round(fine_grid, 6)))
      fine_grid <- fine_grid[!fine_grid %in% evaluated_res]

      for (res in fine_grid) {
        eval_res <- evaluate_layer_resolution(seu, reduction_name, dims, res, cfg$random_seed)
        search_records[[length(search_records) + 1]] <- data.frame(
          stage = "fine",
          resolution = res,
          n_clusters = eval_res$n_clusters,
          stringsAsFactors = FALSE
        )
        evaluated_res <- c(evaluated_res, res)

        if (eval_res$n_clusters == spec$target_clusters) {
          found_res <- res
          selected_obj <- eval_res$object
          selected_cluster_count <- eval_res$n_clusters
          break
        }
      }

      if (!is.na(found_res)) {
        break
      }
    }
  }

  if (is.na(found_res)) {
    search_df <- do.call(rbind, search_records)
    best_idx <- which.min(abs(search_df$n_clusters - spec$target_clusters))
    found_res <- search_df$resolution[best_idx]
    fallback_eval <- evaluate_layer_resolution(seu, reduction_name, dims, found_res, cfg$random_seed)
    selected_obj <- fallback_eval$object
    selected_cluster_count <- fallback_eval$n_clusters
    search_records[[length(search_records) + 1]] <- data.frame(
      stage = "fallback",
      resolution = found_res,
      n_clusters = selected_cluster_count,
      stringsAsFactors = FALSE
    )
  }

  search_df <- do.call(rbind, search_records)
  search_df <- dplyr::distinct(search_df, stage, resolution, .keep_all = TRUE)
  search_df <- dplyr::arrange(search_df, resolution, stage)

  list(
    object = selected_obj,
    search_df = search_df,
    selected_reduction = reduction_name,
    selected_umap = selected_layer_umap(selected_obj, reduction_name),
    selected_resolution = found_res,
    selected_cluster_count = selected_cluster_count,
    dims = dims
  )
}

finalize_layer_object <- function(search_result, spec) {
  obj <- search_result$object
  obj@misc$selected_reduction <- search_result$selected_reduction
  obj@misc$selected_umap <- search_result$selected_umap

  if (search_result$selected_umap %in% Reductions(obj)) {
    obj[["umap"]] <- obj[[search_result$selected_umap]]
  }

  Idents(obj) <- "seurat_clusters"
  current_ids <- levels(factor(obj$seurat_clusters))
  new_ids <- setNames(as.character(seq_along(current_ids)), current_ids)
  obj <- RenameIdents(obj, new_ids)
  obj$seurat_clusters <- factor(Idents(obj), levels = as.character(seq_along(current_ids)))

  cluster_col <- layer_cluster_column(spec$layer_id)
  obj[[cluster_col]] <- obj$seurat_clusters
  obj$layer_id <- spec$layer_id
  obj$layer_role <- spec$layer_role
  obj$parent_layer_id <- if (nzchar(spec$parent_layer)) spec$parent_layer else NA_character_

  obj@misc$layer_workflow <- c(
    obj@misc$layer_workflow,
    list(
      layer_id = spec$layer_id,
      layer_role = spec$layer_role,
      parent_layer = spec$parent_layer,
      cluster_column = cluster_col,
      selected_reduction = search_result$selected_reduction,
      selected_umap = search_result$selected_umap,
      selected_resolution = search_result$selected_resolution,
      selected_cluster_count = search_result$selected_cluster_count
    )
  )

  search_result$object <- obj
  search_result$cluster_column <- cluster_col
  search_result
}

build_layer_cluster_summary <- function(seu, group_var, cluster_var = "seurat_clusters") {
  summary_df <- tibble::rownames_to_column(seu@meta.data, "cell_id")
  summary_df$plot_group <- summary_df[[group_var]]
  summary_df$cluster_value <- summary_df[[cluster_var]]
  summary_df <- dplyr::count(summary_df, plot_group, cluster_value, name = "cell_number")
  summary_df <- dplyr::group_by(summary_df, plot_group)
  summary_df <- dplyr::mutate(summary_df, proportion = round(cell_number / sum(cell_number) * 100, 5))
  summary_df <- dplyr::ungroup(summary_df)
  colnames(summary_df)[colnames(summary_df) == "cluster_value"] <- cluster_var
  summary_df
}

build_layer_status_row <- function(spec, paths, status, reason = "", selected_resolution = NA_real_, selected_cluster_count = NA_integer_, selected_reduction = "", cluster_column = "") {
  data.frame(
    layer_id = spec$layer_id,
    layer_role = spec$layer_role,
    enabled = ifelse(spec$enabled, "yes", "no"),
    parent_layer = spec$parent_layer,
    selection_column = spec$selection_column,
    selection_values = paste(spec$selection_values, collapse = ","),
    sample_include = paste(spec$sample_include, collapse = ","),
    sample_exclude = paste(spec$sample_exclude, collapse = ","),
    status = status,
    reason = reason,
    selected_reduction = selected_reduction,
    selected_resolution = ifelse(is.na(selected_resolution), "", as.character(selected_resolution)),
    selected_cluster_count = ifelse(is.na(selected_cluster_count), "", as.character(selected_cluster_count)),
    cluster_column = cluster_column,
    candidate_rds = paths$candidate_rds,
    clustered_rds = paths$clustered_rds,
    report_md = paths$report_md,
    stringsAsFactors = FALSE
  )
}
