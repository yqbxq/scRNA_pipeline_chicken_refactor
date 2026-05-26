#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
)

source_utf8 <- function(path) {
  source(path, encoding = "UTF-8")
}

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "comparison_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "cell_count_inventory_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble", "tidyr", "jsonlite"))

cfg <- get_single_script_config_04()
module_name <- "04c_subcluster_eda"
prepare_dirs_04(cfg)
set.seed(cfg$random_seed)

empty_subcluster_summary_04c <- function() {
  data.frame(
    layer_id = character(0),
    cluster_count = integer(0),
    target_clusters = integer(0),
    determined_n = integer(0),
    tentative_n = integer(0),
    undetermined_n = integer(0),
    determined_fraction = numeric(0),
    subset_pass_count = integer(0),
    condition_split_path = character(0),
    panorama_cross_path = character(0),
    comparison_subset_summary_path = character(0),
    summary_plot_png = character(0),
    annotated_rds = character(0),
    report_md = character(0),
    stringsAsFactors = FALSE
  )
}

empty_condition_split_04c <- function() {
  data.frame(
    layer_id = character(0),
    condition_var = character(0),
    condition = character(0),
    cluster = character(0),
    cell_type = character(0),
    n_cells = integer(0),
    fraction_of_layer = numeric(0),
    stringsAsFactors = FALSE
  )
}

empty_panorama_cross_04c <- function() {
  data.frame(
    layer_id = character(0),
    panorama_cell_type = character(0),
    subcluster_cell_type = character(0),
    n_cells = integer(0),
    fraction_of_layer = numeric(0),
    stringsAsFactors = FALSE
  )
}

empty_comparison_summary_04c <- function() {
  data.frame(
    layer_id = character(0),
    comparison_id = character(0),
    group_var = character(0),
    ident_1 = character(0),
    ident_2 = character(0),
    ident_1_n = integer(0),
    ident_2_n = integer(0),
    subset_n_cells = integer(0),
    subset_column = character(0),
    subset_value = character(0),
    status = character(0),
    stringsAsFactors = FALSE
  )
}

empty_comparison_split_04c <- function() {
  data.frame(
    layer_id = character(0),
    comparison_id = character(0),
    group_var = character(0),
    group_value = character(0),
    cluster = character(0),
    cell_type = character(0),
    n_cells = integer(0),
    stringsAsFactors = FALSE
  )
}

empty_cell_count_inventory_04c <- function() {
  data.frame(
    cluster_id = character(0),
    sample_id = character(0),
    group_id = character(0),
    layer_id = character(0),
    group_var = character(0),
    parent_cluster_id = character(0),
    n_cells = integer(0),
    total_umi = numeric(0),
    median_umi_per_cell = numeric(0),
    fraction_of_sample = numeric(0),
    fraction_of_group = numeric(0),
    fraction_of_cluster = numeric(0),
    stringsAsFactors = FALSE
  )
}

empty_cluster_eligibility_04c <- function() {
  data.frame(
    layer_id = character(0),
    group_var = character(0),
    comparison_id = character(0),
    cluster_id = character(0),
    group_a = character(0),
    group_b = character(0),
    evidence_tier = character(0),
    recommended_action = character(0),
    failed_rules = character(0),
    reason = character(0),
    n_samples_a = integer(0),
    n_samples_b = integer(0),
    min_cells_a = integer(0),
    min_cells_b = integer(0),
    total_umi_a = numeric(0),
    total_umi_b = numeric(0),
    max_single_sample_frac = numeric(0),
    parent_cluster_id = character(0),
    min_cells_per_sample = numeric(0),
    min_samples_per_group = numeric(0),
    min_total_umi = numeric(0),
    max_single_sample_frac_threshold = numeric(0),
    min_cells_soft_floor = numeric(0),
    inventory_version = character(0),
    stringsAsFactors = FALSE
  )
}

empty_inventory_summary_04c <- function() {
  data.frame(
    layer_id = character(0),
    group_var = character(0),
    comparison_id = character(0),
    evidence_tier = character(0),
    n_cluster_comparisons = integer(0),
    stringsAsFactors = FALSE
  )
}

manifest_output_or_empty <- function(manifest, key) {
  entry <- manifest$outputs[[key]]
  if (is.null(entry) || is.null(entry$path)) {
    return("")
  }
  resolve_output_local(manifest, key)
}

manifest_annotated_keys_04c <- function(manifest) {
  output_names <- names(manifest$outputs %||% list())
  output_names[startsWith(output_names, "annotated_")]
}

layer_id_from_annotated_key_04c <- function(key) {
  sub("^annotated_", "", key)
}

first_present_col_04c <- function(meta_df, candidates) {
  hit <- candidates[candidates %in% colnames(meta_df)]
  if (length(hit) == 0) "" else hit[[1]]
}

count_col_for_inventory_04c <- function(meta_df) {
  first_present_col_04c(meta_df, c("nCount_RNA", "nCount_SCT", "nCount_Spatial", grep("^nCount_", colnames(meta_df), value = TRUE)))
}

sample_col_for_inventory_04c <- function(meta_df) {
  first_present_col_04c(meta_df, c("sample_id", "orig.ident", "sample", "sample_name", "section_id"))
}

parent_col_for_inventory_04c <- function(meta_df, cell_type_col) {
  first_present_col_04c(meta_df, unique(c("parent_cell_type", "panorama_cell_type", "cell_type_parent")))
}

inventory_cluster_col_04c <- function(meta_df, cluster_col, cell_type_col) {
  first_present_col_04c(meta_df, unique(c(cell_type_col, cluster_col, "cell_subtype", "cell_type", "annotation_label", "seurat_clusters")))
}

build_layer_inventory_04c <- function(seu, layer_id, comparisons, cluster_col, cell_type_col) {
  meta <- seu@meta.data
  cluster_var <- inventory_cluster_col_04c(meta, cluster_col, cell_type_col)
  sample_var <- sample_col_for_inventory_04c(meta)
  count_var <- count_col_for_inventory_04c(meta)
  if (!nzchar(cluster_var) || !nzchar(sample_var)) {
    return(empty_cell_count_inventory_04c())
  }
  if (!nzchar(count_var)) {
    count_var <- ".inventory_count_value"
    meta[[count_var]] <- 1
  }
  group_vars <- unique(as.character(comparisons$group_var %||% character(0)))
  group_vars <- group_vars[nzchar(group_vars) & group_vars %in% colnames(meta)]
  if (length(group_vars) == 0) {
    fallback_group <- first_present_col_04c(meta, c("condition", "group_id", "stage", "sample_id", "orig.ident"))
    if (nzchar(fallback_group)) {
      group_vars <- fallback_group
    }
  }
  if (length(group_vars) == 0) {
    return(empty_cell_count_inventory_04c())
  }

  parent_var <- parent_col_for_inventory_04c(meta, cell_type_col)
  rows <- lapply(group_vars, function(group_var) {
    inv <- build_cell_count_inventory(
      meta,
      cluster_var = cluster_var,
      sample_var = sample_var,
      group_var = group_var,
      layer_var = NULL,
      parent_var = if (nzchar(parent_var)) parent_var else NULL,
      count_var = count_var
    )
    inv$layer_id <- layer_id
    inv$group_var <- group_var
    inv[, c(
      "cluster_id", "sample_id", "group_id", "layer_id", "group_var", "parent_cluster_id",
      "n_cells", "total_umi", "median_umi_per_cell",
      "fraction_of_sample", "fraction_of_group", "fraction_of_cluster"
    ), drop = FALSE]
  })
  if (length(rows) == 0) empty_cell_count_inventory_04c() else dplyr::bind_rows(rows)
}

build_cluster_eligibility_04c <- function(inventory_df, comparisons) {
  if (nrow(inventory_df) == 0 || nrow(comparisons) == 0) {
    return(empty_cluster_eligibility_04c())
  }
  comparisons <- comparisons[comparisons$enabled != "no", , drop = FALSE]
  if (nrow(comparisons) == 0) {
    return(empty_cluster_eligibility_04c())
  }
  rows <- list()
  layer_ids <- unique(as.character(inventory_df$layer_id))
  for (layer_id in layer_ids) {
    layer_comparisons <- comparison_rows_for_layer_04c(comparisons, layer_id)
    if (nrow(layer_comparisons) == 0) {
      next
    }
    group_vars <- unique(as.character(inventory_df$group_var[inventory_df$layer_id == layer_id]))
    group_vars <- group_vars[nzchar(group_vars)]
    for (group_var in group_vars) {
      comp_subset <- layer_comparisons[layer_comparisons$group_var == group_var, , drop = FALSE]
      inv_subset <- inventory_df[inventory_df$layer_id == layer_id & inventory_df$group_var == group_var, , drop = FALSE]
      if (nrow(comp_subset) == 0 || nrow(inv_subset) == 0) {
        next
      }
      eligible <- classify_all_cluster_eligibilities(inv_subset, comp_subset, cfg$cell_count_threshold_file %||% NULL)
      if (nrow(eligible) == 0) {
        next
      }
      eligible$layer_id <- layer_id
      eligible$group_var <- group_var
      rows[[length(rows) + 1]] <- eligible[, c(
        "layer_id", "group_var", setdiff(colnames(eligible), c("layer_id", "group_var"))
      ), drop = FALSE]
    }
  }
  if (length(rows) == 0) empty_cluster_eligibility_04c() else dplyr::bind_rows(rows)
}

summarize_inventory_eligibility_04c <- function(eligibility_df) {
  if (nrow(eligibility_df) == 0) {
    return(empty_inventory_summary_04c())
  }
  eligibility_df %>%
    dplyr::count(layer_id, group_var, comparison_id, evidence_tier, name = "n_cluster_comparisons") %>%
    dplyr::arrange(layer_id, comparison_id, evidence_tier)
}

format_fraction_04c <- function(x) {
  if (is.na(x)) "NA" else sprintf("%.3f", x)
}

confidence_counts_04c <- function(annotation_table) {
  if (nrow(annotation_table) == 0 || !"confidence" %in% colnames(annotation_table)) {
    return(list(determined_n = 0L, tentative_n = 0L, undetermined_n = 0L, determined_fraction = 0))
  }
  confidence <- as.character(annotation_table$confidence)
  determined_n <- sum(confidence == "确定", na.rm = TRUE)
  tentative_n <- sum(confidence == "暂定", na.rm = TRUE)
  undetermined_n <- sum(confidence == "未定", na.rm = TRUE)
  list(
    determined_n = determined_n,
    tentative_n = tentative_n,
    undetermined_n = undetermined_n,
    determined_fraction = safe_rate(determined_n, length(confidence))
  )
}

target_clusters_for_layer_04c <- function(layer_df, layer_id) {
  hit <- layer_df[layer_df$layer_id == layer_id, , drop = FALSE]
  if (nrow(hit) == 0) {
    return(NA_integer_)
  }
  as.integer(hit$target_clusters[[1]])
}

build_condition_split_04c <- function(seu, layer_id, cluster_col, cell_type_col) {
  meta <- seu@meta.data
  condition_col <- first_present_col_04c(meta, c("condition", "group_id", "sample_id", "orig.ident"))
  if (!nzchar(condition_col) || nrow(meta) == 0) {
    return(empty_condition_split_04c())
  }
  cluster <- if (nzchar(cluster_col) && cluster_col %in% colnames(meta)) as.character(meta[[cluster_col]]) else ""
  cell_type <- if (nzchar(cell_type_col) && cell_type_col %in% colnames(meta)) as.character(meta[[cell_type_col]]) else ""
  out <- data.frame(
    layer_id = layer_id,
    condition_var = condition_col,
    condition = as.character(meta[[condition_col]]),
    cluster = cluster,
    cell_type = cell_type,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::count(layer_id, condition_var, condition, cluster, cell_type, name = "n_cells") %>%
    dplyr::mutate(fraction_of_layer = safe_rate(n_cells, sum(n_cells)))
  out
}

build_panorama_cross_04c <- function(seu, layer_id, cell_type_col) {
  meta <- seu@meta.data
  panorama_col <- first_present_col_04c(meta, c("panorama_cell_type", "parent_cell_type"))
  if (!nzchar(panorama_col) || !nzchar(cell_type_col) || !cell_type_col %in% colnames(meta) || nrow(meta) == 0) {
    return(empty_panorama_cross_04c())
  }
  data.frame(
    layer_id = layer_id,
    panorama_cell_type = as.character(meta[[panorama_col]]),
    subcluster_cell_type = as.character(meta[[cell_type_col]]),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::count(layer_id, panorama_cell_type, subcluster_cell_type, name = "n_cells") %>%
    dplyr::mutate(fraction_of_layer = safe_rate(n_cells, sum(n_cells)))
}

comparison_rows_for_layer_04c <- function(comparisons, layer_id) {
  if (nrow(comparisons) == 0) {
    return(comparisons)
  }
  comparisons <- comparisons[comparisons$enabled != "no", , drop = FALSE]
  if (nrow(comparisons) == 0) {
    return(comparisons)
  }
  keep <- vapply(seq_len(nrow(comparisons)), function(i) {
    comparison_applies_to_layer(comparisons[i, , drop = FALSE], layer_id)
  }, logical(1))
  comparisons[keep, , drop = FALSE]
}

build_comparison_tables_04c <- function(seu, layer_id, comparisons, cluster_col, cell_type_col) {
  if (nrow(comparisons) == 0) {
    return(list(summary = empty_comparison_summary_04c(), split = empty_comparison_split_04c(), triage = empty_triage_df(include_sample = TRUE)))
  }

  summary_rows <- list()
  split_rows <- list()
  triage_rows <- list()
  for (idx in seq_len(nrow(comparisons))) {
    comparison_row <- comparisons[idx, , drop = FALSE]
    comparison_id <- normalize_scalar_value(comparison_row$comparison_id[[1]], sprintf("comparison_%s", idx))
    group_var <- normalize_scalar_value(comparison_row$group_var[[1]])
    ident_1 <- normalize_scalar_value(comparison_row$ident_1[[1]])
    ident_2 <- normalize_scalar_value(comparison_row$ident_2[[1]])
    subset_column <- normalize_scalar_value(comparison_row$subset_column[[1]])
    subset_value <- normalize_scalar_value(comparison_row$subset_value[[1]])
    subset_requested <- nzchar(subset_column) || nzchar(subset_value)
    subset_seu <- seu
    subset_n <- ncol(seu)
    status <- "ok"

    if (isTRUE(subset_requested)) {
      meta <- seu@meta.data
      if (!nzchar(subset_column) || !subset_column %in% colnames(meta)) {
        status <- "subset_column_missing"
        subset_n <- 0L
        subset_seu <- NULL
        triage_rows[[length(triage_rows) + 1]] <- make_triage_row(
          sample_id = layer_id,
          severity = "low",
          signal_id = "subcluster_subset_no_match",
          evidence = sprintf(
            "status=subset_column_missing; layer_id=%s; comparison_id=%s; subset_column=%s; subset_value=%s",
            layer_id,
            comparison_id,
            subset_column,
            subset_value
          )
        )
      } else {
        keep_values <- split_csv_local(subset_value)
        keep <- as.character(meta[[subset_column]]) %in% keep_values
        if (!any(keep)) {
          status <- "subset_no_match"
          subset_n <- 0L
          subset_seu <- NULL
          triage_rows[[length(triage_rows) + 1]] <- make_triage_row(
            sample_id = layer_id,
            severity = "low",
            signal_id = "subcluster_subset_no_match",
            evidence = sprintf(
              "status=subset_no_match; user-requested subset matched 0 cells; layer_id=%s; comparison_id=%s; subset_column=%s; subset_value=%s",
              layer_id,
              comparison_id,
              subset_column,
              subset_value
            )
          )
        } else {
          subset_seu <- apply_subset_filter(seu, comparison_row)
          subset_n <- ncol(subset_seu)
        }
      }
    }

    meta <- if (is.null(subset_seu)) data.frame(stringsAsFactors = FALSE) else subset_seu@meta.data
    ident_1_n <- NA_integer_
    ident_2_n <- NA_integer_
    if (identical(status, "ok")) {
      if (subset_n == 0) {
        status <- "empty_subset"
      } else if (!nzchar(group_var) || !group_var %in% colnames(meta)) {
        status <- "missing_group_var"
      } else {
        group <- as.character(meta[[group_var]])
        ident_1_n <- sum(group == ident_1, na.rm = TRUE)
        ident_2_n <- sum(group == ident_2, na.rm = TRUE)
      }
    }

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      layer_id = layer_id,
      comparison_id = comparison_id,
      group_var = group_var,
      ident_1 = ident_1,
      ident_2 = ident_2,
      ident_1_n = ident_1_n,
      ident_2_n = ident_2_n,
      subset_n_cells = subset_n,
      subset_column = subset_column,
      subset_value = subset_value,
      status = status,
      stringsAsFactors = FALSE
    )

    if (identical(status, "ok") && subset_n > 0) {
      cluster <- if (nzchar(cluster_col) && cluster_col %in% colnames(meta)) as.character(meta[[cluster_col]]) else ""
      cell_type <- if (nzchar(cell_type_col) && cell_type_col %in% colnames(meta)) as.character(meta[[cell_type_col]]) else ""
      split_rows[[length(split_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        comparison_id = comparison_id,
        group_var = group_var,
        group_value = as.character(meta[[group_var]]),
        cluster = cluster,
        cell_type = cell_type,
        stringsAsFactors = FALSE
      ) %>%
        dplyr::count(layer_id, comparison_id, group_var, group_value, cluster, cell_type, name = "n_cells")
    }
  }

  list(
    summary = if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else empty_comparison_summary_04c(),
    split = if (length(split_rows) > 0) dplyr::bind_rows(split_rows) else empty_comparison_split_04c(),
    triage = if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)
  )
}

write_summary_plot_04c <- function(path, layer_id, condition_split) {
  ensure_dir(dirname(path))
  if (is.null(condition_split) || nrow(condition_split) == 0) {
    plot_obj <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = sprintf("Subcluster summary: %s", layer_id))
  } else {
    plot_obj <- ggplot2::ggplot(
      condition_split,
      ggplot2::aes(x = cluster, y = n_cells, fill = condition)
    ) +
      ggplot2::geom_col(position = "stack", width = 0.8) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(
        title = sprintf("Subcluster condition split: %s", layer_id),
        x = "Cluster",
        y = "Cells",
        fill = "Condition"
      )
  }
  save_plot_local(plot_obj, path, width = 7, height = 4.5)
}

write_layer_report_04c <- function(path, layer_id, annotated_rds, summary_row, condition_split, panorama_cross, comparison_summary, comparison_split) {
  lines <- c(
    sprintf("# 04c Subcluster EDA: %s", layer_id),
    "",
    sprintf("- annotated_rds: `%s`", annotated_rds),
    sprintf("- cluster_count: `%s`", summary_row$cluster_count[[1]]),
    sprintf("- target_clusters: `%s`", summary_row$target_clusters[[1]]),
    sprintf("- determined_fraction: `%s`", format_fraction_04c(summary_row$determined_fraction[[1]])),
    sprintf("- subset_pass_count: `%s`", summary_row$subset_pass_count[[1]]),
    sprintf("- summary_plot_png: `%s`", summary_row$summary_plot_png[[1]]),
    "",
    "## Annotation Confidence",
    render_markdown_table_local(summary_row[, c("determined_n", "tentative_n", "undetermined_n", "determined_fraction"), drop = FALSE]),
    "",
    "## Condition Split",
    render_markdown_table_local(head(condition_split, 80)),
    "",
    "## Panorama Vs Subcluster",
    render_markdown_table_local(head(panorama_cross, 80)),
    "",
    "## Comparison Subsets",
    render_markdown_table_local(head(comparison_summary, 80)),
    "",
    "## Comparison Condition Split",
    render_markdown_table_local(head(comparison_split, 80))
  )
  ensure_dir(dirname(path))
  write_markdown_local(lines, path)
}

manifest_04b <- read_manifest_local(cfg$module_04b_manifest_path)
annotated_keys <- manifest_annotated_keys_04c(manifest_04b)
annotation_summary <- read_tsv_optional(manifest_output_or_empty(manifest_04b, "summary_tsv"))

layer_df <- validate_layer_config(read_object_layer_config(cfg))
comparisons <- read_comparisons_with_subset(cfg)

summary_rows <- list()
condition_rows <- list()
panorama_rows <- list()
comparison_summary_rows <- list()
comparison_split_rows <- list()
triage_rows <- list()
inventory_rows <- list()
output_entries <- list()

for (annotated_key in annotated_keys) {
  layer_id <- layer_id_from_annotated_key_04c(annotated_key)
  annotated_rds <- resolve_output_local(manifest_04b, annotated_key)
  if (!file.exists(annotated_rds)) {
    warning(sprintf("Skipping missing annotated RDS: %s -> %s", annotated_key, annotated_rds), call. = FALSE)
    next
  }

  message("04c EDA layer: ", layer_id)
  seu <- readRDS(annotated_rds)
  seu <- maybe_join_layers(seu)
  meta <- seu@meta.data
  cluster_col <- paste0(layer_id, "_cluster")
  if (!cluster_col %in% colnames(meta)) {
    cluster_col <- first_present_col_04c(meta, grep("_cluster$", colnames(meta), value = TRUE))
  }
  cell_type_col <- paste0(layer_id, "_cell_type")
  if (!cell_type_col %in% colnames(meta)) {
    cell_type_col <- first_present_col_04c(meta, c("cell_type", "annotation_label"))
  }

  annotation_table_path <- manifest_output_or_empty(manifest_04b, paste0("annotation_table_tsv_", layer_id))
  annotation_table <- read_tsv_optional(annotation_table_path)
  counts <- confidence_counts_04c(annotation_table)
  cluster_count <- if (nzchar(cluster_col) && cluster_col %in% colnames(meta)) {
    length(unique(as.character(meta[[cluster_col]])))
  } else if (nrow(annotation_table) > 0) {
    nrow(annotation_table)
  } else {
    0L
  }

  condition_split <- build_condition_split_04c(seu, layer_id, cluster_col, cell_type_col)
  panorama_cross <- build_panorama_cross_04c(seu, layer_id, cell_type_col)
  layer_comparisons <- comparison_rows_for_layer_04c(comparisons, layer_id)
  comparison_tables <- build_comparison_tables_04c(seu, layer_id, layer_comparisons, cluster_col, cell_type_col)
  layer_inventory <- build_layer_inventory_04c(seu, layer_id, layer_comparisons, cluster_col, cell_type_col)

  table_dir <- layer_subcluster_table_dir_04(cfg, layer_id)
  report_dir <- layer_subcluster_report_dir_04(cfg, layer_id)
  ensure_dir(table_dir)
  ensure_dir(report_dir)

  condition_split_tsv <- file.path(table_dir, sprintf("condition_split_%s.tsv", layer_id))
  panorama_cross_tsv <- file.path(table_dir, sprintf("panorama_subcluster_crosstab_%s.tsv", layer_id))
  comparison_summary_tsv <- file.path(table_dir, sprintf("comparison_subset_summary_%s.tsv", layer_id))
  comparison_split_tsv <- file.path(table_dir, sprintf("comparison_condition_split_%s.tsv", layer_id))
  summary_plot_png <- file.path(report_dir, sprintf("subcluster_summary_%s.png", layer_id))
  report_md <- file.path(report_dir, sprintf("subcluster_eda_%s.md", layer_id))

  write_tsv_local(condition_split, condition_split_tsv)
  write_tsv_local(panorama_cross, panorama_cross_tsv)
  write_tsv_local(comparison_tables$summary, comparison_summary_tsv)
  write_tsv_local(comparison_tables$split, comparison_split_tsv)
  write_summary_plot_04c(summary_plot_png, layer_id, condition_split)

  subset_pass_count <- if (nrow(comparison_tables$summary) > 0) {
    sum(comparison_tables$summary$status == "ok", na.rm = TRUE)
  } else {
    0L
  }
  summary_row <- data.frame(
    layer_id = layer_id,
    cluster_count = as.integer(cluster_count),
    target_clusters = target_clusters_for_layer_04c(layer_df, layer_id),
    determined_n = as.integer(counts$determined_n),
    tentative_n = as.integer(counts$tentative_n),
    undetermined_n = as.integer(counts$undetermined_n),
    determined_fraction = as.numeric(counts$determined_fraction),
    subset_pass_count = as.integer(subset_pass_count),
    condition_split_path = normalizePath(condition_split_tsv, winslash = "/", mustWork = FALSE),
    panorama_cross_path = normalizePath(panorama_cross_tsv, winslash = "/", mustWork = FALSE),
    comparison_subset_summary_path = normalizePath(comparison_summary_tsv, winslash = "/", mustWork = FALSE),
    summary_plot_png = normalizePath(summary_plot_png, winslash = "/", mustWork = FALSE),
    annotated_rds = normalizePath(annotated_rds, winslash = "/", mustWork = FALSE),
    report_md = normalizePath(report_md, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  write_layer_report_04c(
    report_md,
    layer_id,
    annotated_rds,
    summary_row,
    condition_split,
    panorama_cross,
    comparison_tables$summary,
    comparison_tables$split
  )

  summary_rows[[length(summary_rows) + 1]] <- summary_row
  condition_rows[[length(condition_rows) + 1]] <- condition_split
  panorama_rows[[length(panorama_rows) + 1]] <- panorama_cross
  comparison_summary_rows[[length(comparison_summary_rows) + 1]] <- comparison_tables$summary
  comparison_split_rows[[length(comparison_split_rows) + 1]] <- comparison_tables$split
  triage_rows[[length(triage_rows) + 1]] <- comparison_tables$triage
  inventory_rows[[length(inventory_rows) + 1]] <- layer_inventory

  output_entries[[paste0("report_md_", layer_id)]] <- build_output_entry(report_md, "md", module_name, sprintf("04c subcluster review report for %s", layer_id), base_dir = cfg$project_root)
  output_entries[[paste0("summary_plot_png_", layer_id)]] <- build_output_entry(summary_plot_png, "png", module_name, sprintf("04c summary plot for %s", layer_id), base_dir = cfg$project_root)
  output_entries[[paste0("condition_split_tsv_", layer_id)]] <- build_output_entry(condition_split_tsv, "tsv", module_name, sprintf("condition split for %s", layer_id), base_dir = cfg$project_root, schema = infer_schema_from_df(condition_split))
  output_entries[[paste0("panorama_subcluster_crosstab_tsv_", layer_id)]] <- build_output_entry(panorama_cross_tsv, "tsv", module_name, sprintf("panorama vs subcluster crosstab for %s", layer_id), base_dir = cfg$project_root, schema = infer_schema_from_df(panorama_cross))
  output_entries[[paste0("comparison_subset_summary_tsv_", layer_id)]] <- build_output_entry(comparison_summary_tsv, "tsv", module_name, sprintf("comparison subset summary for %s", layer_id), base_dir = cfg$project_root, schema = infer_schema_from_df(comparison_tables$summary))
  output_entries[[paste0("comparison_condition_split_tsv_", layer_id)]] <- build_output_entry(comparison_split_tsv, "tsv", module_name, sprintf("comparison condition split for %s", layer_id), base_dir = cfg$project_root, schema = infer_schema_from_df(comparison_tables$split))
}

summary_df <- if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else empty_subcluster_summary_04c()
condition_df <- if (length(condition_rows) > 0) dplyr::bind_rows(condition_rows) else empty_condition_split_04c()
panorama_df <- if (length(panorama_rows) > 0) dplyr::bind_rows(panorama_rows) else empty_panorama_cross_04c()
comparison_summary_df <- if (length(comparison_summary_rows) > 0) dplyr::bind_rows(comparison_summary_rows) else empty_comparison_summary_04c()
comparison_split_df <- if (length(comparison_split_rows) > 0) dplyr::bind_rows(comparison_split_rows) else empty_comparison_split_04c()
triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)
inventory_df <- if (length(inventory_rows) > 0) dplyr::bind_rows(inventory_rows) else empty_cell_count_inventory_04c()
eligibility_df <- build_cluster_eligibility_04c(inventory_df, comparisons)
inventory_summary_df <- summarize_inventory_eligibility_04c(eligibility_df)

subcluster_summary_tsv <- file.path(cfg$subcluster_table_dir, "subcluster_summary.tsv")
condition_split_tsv <- file.path(cfg$subcluster_table_dir, "subcluster_condition_split.tsv")
panorama_cross_tsv <- file.path(cfg$subcluster_table_dir, "subcluster_panorama_crosstab.tsv")
comparison_summary_tsv <- file.path(cfg$subcluster_table_dir, "comparison_subset_summary.tsv")
comparison_split_tsv <- file.path(cfg$subcluster_table_dir, "comparison_condition_split.tsv")
triage_tsv <- file.path(cfg$subcluster_table_dir, "subcluster_eda_triage.tsv")
cell_count_inventory_tsv <- file.path(cfg$subcluster_table_dir, "cell_count_inventory.tsv")
cluster_eligibility_tsv <- file.path(cfg$subcluster_table_dir, "cluster_eligibility.tsv")
inventory_summary_tsv <- file.path(cfg$subcluster_table_dir, "inventory_summary_per_comparison.tsv")

write_tsv_local(summary_df, subcluster_summary_tsv)
write_tsv_local(condition_df, condition_split_tsv)
write_tsv_local(panorama_df, panorama_cross_tsv)
write_tsv_local(comparison_summary_df, comparison_summary_tsv)
write_tsv_local(comparison_split_df, comparison_split_tsv)
write_tsv_local(triage_df, triage_tsv)
write_tsv_local(inventory_df, cell_count_inventory_tsv)
write_tsv_local(eligibility_df, cluster_eligibility_tsv)
write_tsv_local(inventory_summary_df, inventory_summary_tsv)

report_path <- file.path(cfg$subcluster_report_dir, "04c_subcluster_eda.md")
report_lines <- c(
  "# 04c Subcluster EDA",
  "",
  sprintf("- annotated_layers: `%s`", length(summary_rows)),
  sprintf("- subcluster_summary: `%s`", subcluster_summary_tsv),
  sprintf("- comparison_subset_summary: `%s`", comparison_summary_tsv),
  "",
  "## Layer Summary",
  render_markdown_table_local(summary_df),
  "",
  "## Comparison Subsets",
  render_markdown_table_local(head(comparison_summary_df, 100)),
  "",
  "## Triage",
  render_markdown_table_local(head(triage_df, 100)),
  "",
  "## Cell-count inventory + evidence tier 建议",
  sprintf("- cell_count_inventory_tsv: `%s`", cell_count_inventory_tsv),
  sprintf("- cluster_eligibility_tsv: `%s`", cluster_eligibility_tsv),
  sprintf("- inventory_summary_per_comparison_tsv: `%s`", inventory_summary_tsv),
  "",
  "### Evidence Tier Summary",
  render_markdown_table_local(inventory_summary_df),
  "",
  "### Eligibility Preview",
  render_markdown_table_local(head(eligibility_df[, intersect(c("layer_id", "group_var", "comparison_id", "cluster_id", "evidence_tier", "recommended_action", "reason"), colnames(eligibility_df)), drop = FALSE], 100)),
  "",
  "## Per-Layer Reports",
  render_markdown_table_local(summary_df[, c("layer_id", "report_md"), drop = FALSE])
)
ensure_dir(dirname(report_path))
write_markdown_local(report_lines, report_path)

output_entries$subcluster_summary_tsv <- build_output_entry(subcluster_summary_tsv, "tsv", module_name, "one row per reviewed subcluster layer", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df))
output_entries$condition_split_tsv <- build_output_entry(condition_split_tsv, "tsv", module_name, "condition split across subcluster layers", base_dir = cfg$project_root, schema = infer_schema_from_df(condition_df))
output_entries$panorama_subcluster_crosstab_tsv <- build_output_entry(panorama_cross_tsv, "tsv", module_name, "panorama vs subcluster cell type crosstab", base_dir = cfg$project_root, schema = infer_schema_from_df(panorama_df))
output_entries$comparison_subset_summary_tsv <- build_output_entry(comparison_summary_tsv, "tsv", module_name, "comparison subset summary across subcluster layers", base_dir = cfg$project_root, schema = infer_schema_from_df(comparison_summary_df))
output_entries$comparison_condition_split_tsv <- build_output_entry(comparison_split_tsv, "tsv", module_name, "comparison condition split across subcluster layers", base_dir = cfg$project_root, schema = infer_schema_from_df(comparison_split_df))
output_entries$subcluster_eda_triage_tsv <- build_output_entry(triage_tsv, "tsv", module_name, "one row per subcluster EDA triage signal", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df))
output_entries$cell_count_inventory_tsv <- build_output_entry(cell_count_inventory_tsv, "tsv", module_name, "cell count inventory by layer, cluster, sample, and group", base_dir = cfg$project_root, schema = infer_schema_from_df(inventory_df))
output_entries$cluster_eligibility_tsv <- build_output_entry(cluster_eligibility_tsv, "tsv", module_name, "cluster x comparison evidence tiers and recommended actions", base_dir = cfg$project_root, schema = infer_schema_from_df(eligibility_df))
output_entries$inventory_summary_per_comparison_tsv <- build_output_entry(inventory_summary_tsv, "tsv", module_name, "evidence tier counts per layer and comparison", base_dir = cfg$project_root, schema = infer_schema_from_df(inventory_summary_df))
output_entries$report <- build_output_entry(report_path, "md", module_name, "subcluster EDA review report", base_dir = cfg$project_root)

if (file.exists(cfg$module_04c_manifest_path)) {
  unlink(cfg$module_04c_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_04c_manifest_path,
  new_outputs = output_entries,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_04b_manifest = cfg$module_04b_manifest_path,
    comparison_sheet = cfg$comparison_sheet,
    cell_count_threshold_file = cfg$cell_count_threshold_file,
    object_layer_config = cfg$object_layer_config_file,
    annotation_summary_tsv = manifest_output_or_empty(manifest_04b, "summary_tsv")
  ),
  version = env_or_default_03("MODULE_04C_VERSION", cfg$module_version),
  depends_on = list(module_04b = cfg$module_04b_manifest_path)
)

message("04c completed. reviewed layers: ", length(summary_rows))
