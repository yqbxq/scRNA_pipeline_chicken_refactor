source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "object_layer_helpers.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "plotting_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

compat_stage_dir <- ensure_eda_stage_dir(cfg, "integration")
layer_specs <- read_object_layer_specs(cfg)
panorama_id <- panorama_layer_id(layer_specs)
panorama_spec <- layer_specs[[panorama_id]]
panorama_paths <- object_layer_paths(cfg, panorama_id)
ensure_object_layer_dirs(panorama_paths)

input_rds <- panorama_paths$candidate_rds
if (!file.exists(input_rds)) {
  input_rds <- file.path(cfg$checkpoint_dir, "02_reduction_candidates.rds")
}
if (!file.exists(input_rds)) {
  stop("缺少输入检查点: panorama reduction candidates", call. = FALSE)
}

panorama_cache_signature_path <- function(paths) {
  file.path(paths$table_dir, "build_signature.tsv")
}

format_signature_numeric <- function(x) {
  paste(format(as.numeric(x), trim = TRUE, scientific = FALSE), collapse = ",")
}

panorama_cache_signature <- function(spec) {
  c(
    target_clusters = as.character(spec$target_clusters),
    res_range = format_signature_numeric(spec$res_range),
    integration_mode = normalize_scalar_value(spec$integration_mode),
    pca_dims = format_index_spec(spec$pca_dims)
  )
}

write_panorama_cache_signature <- function(spec, paths) {
  signature <- panorama_cache_signature(spec)
  write_tsv(
    data.frame(
      key = names(signature),
      value = unname(signature),
      stringsAsFactors = FALSE
    ),
    panorama_cache_signature_path(paths)
  )
}

read_panorama_cache_signature <- function(paths) {
  path <- panorama_cache_signature_path(paths)
  df <- read_tsv_optional(path)
  if (nrow(df) == 0 || !"key" %in% colnames(df) || !"value" %in% colnames(df)) {
    return(setNames(character(0), character(0)))
  }
  signature <- setNames(as.character(df$value), as.character(df$key))
  signature[nzchar(names(signature))]
}

describe_signature <- function(signature) {
  paste(sprintf("%s=%s", names(signature), unname(signature)), collapse = "; ")
}

compare_panorama_cache_signature <- function(spec, paths) {
  current <- panorama_cache_signature(spec)
  cached <- read_panorama_cache_signature(paths)
  signature_path <- panorama_cache_signature_path(paths)

  if (length(cached) == 0) {
    return(list(
      match = FALSE,
      reason = sprintf("缺少 panorama cache signature: %s", signature_path),
      current = current,
      cached = cached
    ))
  }

  missing_keys <- setdiff(names(current), names(cached))
  if (length(missing_keys) > 0) {
    return(list(
      match = FALSE,
      reason = sprintf("panorama cache signature 缺少字段: %s", paste(missing_keys, collapse = ", ")),
      current = current,
      cached = cached
    ))
  }

  cached <- cached[names(current)]
  if (identical(unname(current), unname(cached))) {
    return(list(
      match = TRUE,
      reason = describe_signature(current),
      current = current,
      cached = cached
    ))
  }

  diff_keys <- names(current)[current != cached]
  diff_text <- paste(
    sprintf("%s[%s -> %s]", diff_keys, unname(cached[diff_keys]), unname(current[diff_keys])),
    collapse = "; "
  )
  list(
    match = FALSE,
    reason = sprintf("panorama 参数已变更: %s", diff_text),
    current = current,
    cached = cached
  )
}

read_selected_resolution_metadata <- function(paths) {
  if (!file.exists(paths$selected_resolution_txt)) {
    return(list(
      selected_reduction = "",
      selected_resolution = NA_real_,
      selected_cluster_count = NA_integer_,
      cluster_column = ""
    ))
  }

  lines <- readLines(paths$selected_resolution_txt, warn = FALSE)
  lines <- lines[grepl("=", lines, fixed = TRUE)]
  if (length(lines) == 0) {
    return(list(
      selected_reduction = "",
      selected_resolution = NA_real_,
      selected_cluster_count = NA_integer_,
      cluster_column = ""
    ))
  }

  kv <- strsplit(lines, "=", fixed = TRUE)
  values <- setNames(
    vapply(kv, function(parts) {
      if (length(parts) < 2) {
        return("")
      }
      paste(parts[-1], collapse = "=")
    }, character(1)),
    vapply(kv, `[`, character(1), 1)
  )

  list(
    selected_reduction = normalize_scalar_value(values[["selected_reduction"]]),
    selected_resolution = suppressWarnings(as.numeric(normalize_scalar_value(values[["selected_resolution"]], ""))),
    selected_cluster_count = suppressWarnings(as.integer(normalize_scalar_value(values[["selected_cluster_count"]], ""))),
    cluster_column = normalize_scalar_value(values[["cluster_column"]])
  )
}

read_existing_cluster_summary <- function(paths) {
  if (!file.exists(paths$cluster_summary_csv)) {
    return(NULL)
  }
  tryCatch(
    read.csv(paths$cluster_summary_csv, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
}

refresh_panorama_compat_outputs <- function(obj, paths, result_meta = NULL, cluster_summary = NULL) {
  saveRDS(obj, file.path(cfg$checkpoint_dir, "02_after_clustering.rds"))

  if (!is.null(cluster_summary)) {
    write.csv(cluster_summary, file.path(cfg$table_dir, "cluster_summary.csv"), row.names = FALSE)
  }
  if (file.exists(paths$resolution_search_csv)) {
    file.copy(paths$resolution_search_csv, file.path(cfg$table_dir, "resolution_search.csv"), overwrite = TRUE)
  }
  if (file.exists(paths$selected_resolution_txt)) {
    file.copy(paths$selected_resolution_txt, file.path(cfg$table_dir, "selected_resolution.txt"), overwrite = TRUE)
  } else if (!is.null(result_meta)) {
    writeLines(
      c(
        sprintf("selected_reduction=%s", normalize_scalar_value(result_meta$selected_reduction, "")),
        sprintf("selected_resolution=%s", ifelse(is.finite(result_meta$selected_resolution), result_meta$selected_resolution, "")),
        sprintf("selected_cluster_count=%s", ifelse(is.finite(result_meta$selected_cluster_count), result_meta$selected_cluster_count, "")),
        sprintf("target_clusters=%s", panorama_spec$target_clusters),
        sprintf("cluster_column=%s", normalize_scalar_value(result_meta$cluster_column, ""))
      ),
      file.path(cfg$table_dir, "selected_resolution.txt")
    )
  }
}

reuse_panorama_checkpoint <- function(spec, paths) {
  finalized_rds <- first_existing_path(c(paths$clustered_rds, file.path(cfg$checkpoint_dir, "02_after_clustering.rds")))
  if (!nzchar(finalized_rds)) {
    return(list(reuse = FALSE, reason = "缺少 panorama finalized checkpoint"))
  }

  signature_check <- compare_panorama_cache_signature(spec, paths)
  if (!isTRUE(signature_check$match)) {
    return(list(reuse = FALSE, reason = signature_check$reason))
  }

  cached_obj <- readRDS(finalized_rds)
  cached_obj <- maybe_join_layers(cached_obj)
  result_meta <- read_selected_resolution_metadata(paths)
  cluster_summary <- read_existing_cluster_summary(paths)
  reason <- sprintf("复用 panorama finalized checkpoint: %s (`%s`)", signature_check$reason, finalized_rds)

  refresh_panorama_compat_outputs(
    obj = cached_obj,
    paths = paths,
    result_meta = result_meta,
    cluster_summary = cluster_summary
  )

  write_layer_report(
    spec,
    paths,
    status = "reused",
    reason = reason,
    result = result_meta,
    cluster_summary = cluster_summary
  )

  list(
    reuse = TRUE,
    object = cached_obj,
    reason = reason,
    status_row = build_layer_status_row(
      spec,
      paths,
      status = "reused",
      reason = reason,
      selected_resolution = result_meta$selected_resolution,
      selected_cluster_count = result_meta$selected_cluster_count,
      selected_reduction = result_meta$selected_reduction,
      cluster_column = result_meta$cluster_column
    )
  )
}

build_resolution_plot <- function(search_df, target_clusters, title_text) {
  ggplot(search_df, aes(x = resolution, y = n_clusters, color = stage)) +
    geom_line() +
    geom_point(size = 2) +
    geom_hline(yintercept = target_clusters, linetype = "dashed", color = "#D62728") +
    scale_color_manual(values = c(coarse = "#2F6CB3", fine = "#1A9850", fallback = "#FF8C00"), drop = FALSE) +
    theme_classic(base_size = 11) +
    labs(
      title = title_text,
      x = "Resolution",
      y = "Cluster count",
      color = "search stage"
    )
}

write_layer_report <- function(spec, paths, status, reason = "", result = NULL, cluster_summary = NULL) {
  lines <- c(
    sprintf("# %s Layer Report", spec$layer_id),
    "",
    sprintf("- status: `%s`", status),
    sprintf("- layer_role: `%s`", spec$layer_role),
    sprintf("- parent_layer: `%s`", ifelse(nzchar(spec$parent_layer), spec$parent_layer, "none")),
    sprintf("- object layer config: `%s`", cfg$object_layer_config_file),
    sprintf("- enabled: `%s`", ifelse(spec$enabled, "yes", "no"))
  )

  if (nzchar(spec$description)) {
    lines <- c(lines, sprintf("- description: %s", spec$description))
  }
  if (length(spec$sample_include) > 0) {
    lines <- c(lines, sprintf("- sample_include: `%s`", paste(spec$sample_include, collapse = ",")))
  }
  if (length(spec$sample_exclude) > 0) {
    lines <- c(lines, sprintf("- sample_exclude: `%s`", paste(spec$sample_exclude, collapse = ",")))
  }
  if (!identical(spec$layer_role, "panorama")) {
    lines <- c(
      lines,
      sprintf("- selection_column: `%s`", spec$selection_column),
      sprintf("- selection_values: `%s`", paste(spec$selection_values, collapse = ","))
    )
  }
  if (nzchar(reason)) {
    lines <- c(lines, sprintf("- note: %s", reason))
  }

  lines <- c(
    lines,
    sprintf("- hvg_nfeatures: `%s`", spec$hvg_nfeatures),
    sprintf("- pca_dims: `%s`", format_index_spec(spec$pca_dims)),
    sprintf("- target_clusters: `%s`", spec$target_clusters),
    sprintf("- res_range: `%s`", paste(spec$res_range, collapse = ",")),
    sprintf("- res_fine_step: `%s`", spec$res_fine_step),
    sprintf("- integration_mode: `%s`", spec$integration_mode),
    ""
  )

  if (!is.null(result)) {
    lines <- c(
      lines,
      "## Finalized Model",
      sprintf("- selected_reduction: `%s`", result$selected_reduction),
      sprintf("- selected_resolution: `%s`", result$selected_resolution),
      sprintf("- selected_cluster_count: `%s`", result$selected_cluster_count),
      sprintf("- cluster_column: `%s`", result$cluster_column),
      sprintf("- candidate_rds: `%s`", paths$candidate_rds),
      sprintf("- clustered_rds: `%s`", paths$clustered_rds),
      sprintf("- resolution_search: `%s`", paths$resolution_search_csv),
      sprintf("- cluster_summary: `%s`", paths$cluster_summary_csv),
      sprintf("- resolution_plot: `%s`", paths$resolution_plot_png)
    )

    if (!is.null(cluster_summary)) {
      lines <- c(lines, sprintf("- cluster_summary_rows: `%s`", nrow(cluster_summary)))
    }
  }

  write_markdown(lines, paths$report_md)
}

save_layer_cluster_outputs <- function(obj, spec, paths, result, compat_panorama = FALSE) {
  ensure_object_layer_dirs(paths)

  group_var <- resolve_plot_group_var(obj)
  cluster_summary <- build_layer_cluster_summary(obj, group_var = group_var)
  write.csv(cluster_summary, paths$cluster_summary_csv, row.names = FALSE)
  write.csv(result$search_df, paths$resolution_search_csv, row.names = FALSE)
  writeLines(
    c(
      sprintf("selected_reduction=%s", result$selected_reduction),
      sprintf("selected_resolution=%s", result$selected_resolution),
      sprintf("selected_cluster_count=%s", result$selected_cluster_count),
      sprintf("target_clusters=%s", spec$target_clusters),
      sprintf("cluster_column=%s", result$cluster_column)
    ),
    paths$selected_resolution_txt
  )
  saveRDS(obj, paths$clustered_rds)
  if (identical(spec$layer_id, panorama_id)) {
    write_panorama_cache_signature(spec, paths)
  }

  resolution_plot <- build_resolution_plot(
    result$search_df,
    spec$target_clusters,
    sprintf("%s resolution search", spec$layer_id)
  )
  save_plot_dual(resolution_plot, paths$resolution_plot_png, width = 7, height = 5)

  tryCatch({
    cluster_umap <- plot_cluster_split_umap_paper(obj, group_var = group_var, cluster_var = "seurat_clusters")
    cluster_prop <- plot_cluster_proportion_paper(obj, group_var = group_var, cluster_var = "seurat_clusters")
    cluster_overview <- (cluster_umap | cluster_prop) +
      plot_layout(widths = c(2.4, 1)) +
      plot_annotation(tag_levels = "A")

    save_plot_dual(cluster_umap, paths$cluster_umap_png, width = 14, height = 6)
    save_plot_dual(cluster_prop, paths$cluster_proportion_png, width = 6.5, height = 6)
    save_plot_dual(cluster_overview, paths$cluster_overview_png, width = 16, height = 6.5)

    if (isTRUE(compat_panorama)) {
      save_plot_dual(cluster_umap, file.path(cfg$figure_dir, "Figure_2A.png"), width = 14, height = 6)
      save_plot_dual(cluster_prop, file.path(cfg$figure_dir, "Figure_2B.png"), width = 6.5, height = 6)
      save_plot_dual(cluster_overview, file.path(cfg$figure_dir, "Figure_2.png"), width = 16, height = 6.5)
    }
  }, error = function(e) {
    warning(sprintf("%s 聚类绘图失败，但检查点与表格已保存: %s", spec$layer_id, conditionMessage(e)))
  })

  if (isTRUE(compat_panorama)) {
    write.csv(cluster_summary, file.path(cfg$table_dir, "cluster_summary.csv"), row.names = FALSE)
    write.csv(result$search_df, file.path(cfg$table_dir, "resolution_search.csv"), row.names = FALSE)
    writeLines(
      c(
        sprintf("selected_reduction=%s", result$selected_reduction),
        sprintf("selected_resolution=%s", result$selected_resolution),
        sprintf("selected_cluster_count=%s", result$selected_cluster_count),
        sprintf("target_clusters=%s", spec$target_clusters),
        sprintf("cluster_column=%s", result$cluster_column)
      ),
      file.path(cfg$table_dir, "selected_resolution.txt")
    )
    saveRDS(obj, file.path(cfg$checkpoint_dir, "02_after_clustering.rds"))
  }

  cluster_summary
}

build_layer <- function(base_obj, spec, save_candidate = TRUE) {
  paths <- object_layer_paths(cfg, spec$layer_id)
  ensure_object_layer_dirs(paths)

  if (!isTRUE(spec$enabled)) {
    write_layer_report(spec, paths, status = "disabled", reason = "layer.enabled=no")
    return(list(
      status_row = build_layer_status_row(spec, paths, status = "disabled", reason = "layer.enabled=no")
    ))
  }

  prepared <- prepare_parent_object_for_layer(base_obj, spec)
  if (!identical(prepared$status, "ready")) {
    write_layer_report(spec, paths, status = prepared$status, reason = prepared$reason)
    return(list(
      status_row = build_layer_status_row(spec, paths, status = prepared$status, reason = prepared$reason)
    ))
  }

  candidate_obj <- build_layer_reduction_candidates(prepared$object, spec, cfg)
  if (isTRUE(save_candidate)) {
    saveRDS(candidate_obj, paths$candidate_rds)
  }

  result <- run_layer_resolution_search(candidate_obj, spec, cfg)
  result <- finalize_layer_object(result, spec)

  cluster_summary <- save_layer_cluster_outputs(
    obj = result$object,
    spec = spec,
    paths = paths,
    result = result,
    compat_panorama = identical(spec$layer_id, panorama_id)
  )
  write_layer_report(spec, paths, status = "built", result = result, cluster_summary = cluster_summary)

  list(
    object = result$object,
    paths = paths,
    result = result,
    cluster_summary = cluster_summary,
    status_row = build_layer_status_row(
      spec,
      paths,
      status = "built",
      selected_resolution = result$selected_resolution,
      selected_cluster_count = result$selected_cluster_count,
      selected_reduction = result$selected_reduction,
      cluster_column = result$cluster_column
    )
  )
}

panorama_build <- reuse_panorama_checkpoint(panorama_spec, panorama_paths)
if (isTRUE(panorama_build$reuse)) {
  message(panorama_build$reason)
} else {
  message("panorama finalized checkpoint 不复用，执行重建: ", panorama_build$reason)
  panorama_input <- readRDS(input_rds)
  panorama_input <- maybe_join_layers(panorama_input)
  panorama_build <- build_layer(panorama_input, panorama_spec, save_candidate = FALSE)
}
if (!"object" %in% names(panorama_build)) {
  stop("panorama finalize 失败，无法继续构建子对象层", call. = FALSE)
}

built_objects <- list()
built_objects[[panorama_id]] <- panorama_build$object
layer_status_rows <- list(panorama_build$status_row)

pending_layer_ids <- setdiff(names(layer_specs), panorama_id)
processed_ids <- character(0)
remaining_ids <- pending_layer_ids

while (length(remaining_ids) > 0) {
  progressed <- FALSE

  for (layer_id in remaining_ids) {
    spec <- layer_specs[[layer_id]]
    if (!spec$parent_layer %in% names(built_objects)) {
      next
    }

    child_build <- build_layer(built_objects[[spec$parent_layer]], spec, save_candidate = TRUE)
    layer_status_rows[[length(layer_status_rows) + 1]] <- child_build$status_row
    if ("object" %in% names(child_build)) {
      built_objects[[layer_id]] <- child_build$object
    }
    processed_ids <- c(processed_ids, layer_id)
    progressed <- TRUE
  }

  remaining_ids <- setdiff(pending_layer_ids, processed_ids)
  if (!progressed) {
    for (layer_id in remaining_ids) {
      spec <- layer_specs[[layer_id]]
      paths <- object_layer_paths(cfg, spec$layer_id)
      ensure_object_layer_dirs(paths)
      reason <- sprintf("等待 parent_layer=%s 但未找到可用上游对象，请检查 object_layers.tsv 是否存在循环或上游层被禁用", spec$parent_layer)
      write_layer_report(spec, paths, status = "blocked", reason = reason)
      layer_status_rows[[length(layer_status_rows) + 1]] <- build_layer_status_row(spec, paths, status = "blocked", reason = reason)
    }
    break
  }
}

layer_status_df <- bind_rows(layer_status_rows)
write_tsv(layer_status_df, file.path(cfg$table_dir, "layer_status.tsv"))

selection_skips <- layer_status_df %>%
  filter(status == "skipped", grepl("selection_values", reason, fixed = TRUE))

manifest_lines <- c(
  "# Object Layer Manifest",
  "",
  sprintf("- object layer config: `%s`", cfg$object_layer_config_file),
  sprintf("- panorama compatibility checkpoint: `%s`", file.path(cfg$checkpoint_dir, "02_after_clustering.rds")),
  sprintf("- layer status table: `%s`", file.path(cfg$table_dir, "layer_status.tsv")),
  "",
  "## Layers"
)
for (i in seq_len(nrow(layer_status_df))) {
  row <- layer_status_df[i, , drop = FALSE]
  manifest_lines <- c(
    manifest_lines,
    sprintf("- `%s` [%s]: %s", row$layer_id, row$status, row$report_md)
  )
}
if (nrow(selection_skips) > 0) {
  manifest_lines <- c(
    manifest_lines,
    "",
    "## Follow-up Required",
    sprintf(
      "- 子对象未构建，请编辑 `%s` 中这些 layer 的 `selection_values` 后重跑 `02b_finalize_clustering.R`: `%s`",
      cfg$object_layer_config_file,
      paste(selection_skips$layer_id, collapse = "`, `")
    )
  )
}
write_markdown(manifest_lines, file.path(compat_stage_dir, "layer_manifest.md"))

message("已保存 panorama 兼容输出: ", file.path(cfg$checkpoint_dir, "02_after_clustering.rds"))
message("layer status 清单: ", file.path(cfg$table_dir, "layer_status.tsv"))
if (nrow(selection_skips) > 0) {
  message(
    "子对象未构建，请编辑 ",
    cfg$object_layer_config_file,
    " 填写 selection_values 后重跑 02b_finalize_clustering.R：",
    paste(selection_skips$layer_id, collapse = ", ")
  )
}
