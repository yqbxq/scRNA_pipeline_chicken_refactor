`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

inventory_env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

inventory_numeric_or <- function(value, default) {
  out <- suppressWarnings(as.numeric(value))
  if (length(out) == 0 || is.na(out[[1]])) default else out[[1]]
}

inventory_integer_or <- function(value, default) {
  out <- suppressWarnings(as.integer(value))
  if (length(out) == 0 || is.na(out[[1]])) default else out[[1]]
}

inventory_blank_or_dash <- function(value) {
  value <- trimws(as.character(value %||% ""))
  !nzchar(value) || identical(value, "-") || is.na(value)
}

inventory_cell <- function(row, name, default = "") {
  if (is.null(row) || !name %in% colnames(row)) {
    return(default)
  }
  value <- row[[name]][[1]]
  if (is.na(value)) default else trimws(as.character(value))
}

inventory_meta_data <- function(seurat_obj) {
  if (is.data.frame(seurat_obj)) {
    return(seurat_obj)
  }
  if (is.list(seurat_obj) && is.data.frame(seurat_obj$meta.data)) {
    return(seurat_obj$meta.data)
  }
  meta <- tryCatch(seurat_obj@meta.data, error = function(e) NULL)
  if (!is.data.frame(meta)) {
    stop("build_cell_count_inventory requires a Seurat object or metadata data.frame", call. = FALSE)
  }
  meta
}

inventory_safe_sum <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  sum(x[is.finite(x)], na.rm = TRUE)
}

inventory_safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else stats::median(x, na.rm = TRUE)
}

build_cell_count_inventory <- function(
  seurat_obj,
  cluster_var = "cell_subtype",
  sample_var = "sample_id",
  group_var = "stage",
  layer_var = NULL,
  parent_var = NULL,
  count_var = "nCount_RNA"
) {
  meta <- inventory_meta_data(seurat_obj)
  required_cols <- c(cluster_var, sample_var, group_var, count_var)
  missing <- setdiff(required_cols, colnames(meta))
  if (length(missing) > 0) {
    stop(sprintf("build_cell_count_inventory missing metadata columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  work <- data.frame(
    cluster_id = as.character(meta[[cluster_var]]),
    sample_id = as.character(meta[[sample_var]]),
    group_id = as.character(meta[[group_var]]),
    layer_id = if (!is.null(layer_var) && layer_var %in% colnames(meta)) as.character(meta[[layer_var]]) else "__all__",
    parent_cluster_id = if (!is.null(parent_var) && parent_var %in% colnames(meta)) as.character(meta[[parent_var]]) else "",
    count_value = suppressWarnings(as.numeric(meta[[count_var]])),
    stringsAsFactors = FALSE
  )
  work <- work[nzchar(work$cluster_id) & nzchar(work$sample_id) & nzchar(work$group_id), , drop = FALSE]
  work$layer_id[is.na(work$layer_id) | !nzchar(work$layer_id)] <- "__all__"
  work$parent_cluster_id[is.na(work$parent_cluster_id)] <- ""
  work$count_value[!is.finite(work$count_value)] <- 0

  if (nrow(work) == 0) {
    return(data.frame(
      cluster_id = character(),
      sample_id = character(),
      group_id = character(),
      layer_id = character(),
      parent_cluster_id = character(),
      n_cells = integer(),
      total_umi = numeric(),
      median_umi_per_cell = numeric(),
      fraction_of_sample = numeric(),
      fraction_of_group = numeric(),
      fraction_of_cluster = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  group_cols <- c("cluster_id", "sample_id", "group_id", "layer_id", "parent_cluster_id")
  n_cells <- stats::aggregate(work$count_value, work[group_cols], length)
  total_umi <- stats::aggregate(work$count_value, work[group_cols], inventory_safe_sum)
  median_umi <- stats::aggregate(work$count_value, work[group_cols], inventory_safe_median)
  colnames(n_cells)[ncol(n_cells)] <- "n_cells"
  colnames(total_umi)[ncol(total_umi)] <- "total_umi"
  colnames(median_umi)[ncol(median_umi)] <- "median_umi_per_cell"

  inventory <- merge(n_cells, total_umi, by = group_cols, all = TRUE)
  inventory <- merge(inventory, median_umi, by = group_cols, all = TRUE)
  inventory$n_cells <- as.integer(inventory$n_cells)
  inventory$total_umi <- as.numeric(inventory$total_umi)

  sample_totals <- stats::aggregate(inventory$n_cells, inventory["sample_id"], sum)
  group_totals <- stats::aggregate(inventory$n_cells, inventory["group_id"], sum)
  cluster_totals <- stats::aggregate(inventory$n_cells, inventory["cluster_id"], sum)
  names(sample_totals) <- c("sample_id", "sample_total")
  names(group_totals) <- c("group_id", "group_total")
  names(cluster_totals) <- c("cluster_id", "cluster_total")
  inventory <- merge(inventory, sample_totals, by = "sample_id", all.x = TRUE)
  inventory <- merge(inventory, group_totals, by = "group_id", all.x = TRUE)
  inventory <- merge(inventory, cluster_totals, by = "cluster_id", all.x = TRUE)
  inventory$fraction_of_sample <- inventory$n_cells / pmax(inventory$sample_total, 1)
  inventory$fraction_of_group <- inventory$n_cells / pmax(inventory$group_total, 1)
  inventory$fraction_of_cluster <- inventory$n_cells / pmax(inventory$cluster_total, 1)
  inventory$sample_total <- NULL
  inventory$group_total <- NULL
  inventory$cluster_total <- NULL
  inventory$layer_id[inventory$layer_id == "__all__"] <- NA_character_

  inventory <- inventory[
    order(inventory$cluster_id, inventory$group_id, inventory$sample_id),
    c(
      "cluster_id", "sample_id", "group_id", "layer_id", "parent_cluster_id",
      "n_cells", "total_umi", "median_umi_per_cell",
      "fraction_of_sample", "fraction_of_group", "fraction_of_cluster"
    ),
    drop = FALSE
  ]
  rownames(inventory) <- NULL
  attr(inventory, "inventory_version") <- inventory_env_or_default("MODULE_INVENTORY_VERSION", "1.0")
  attr(inventory, "built_at") <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  inventory
}

read_threshold_overrides <- function(threshold_overrides_tsv = NULL) {
  if (is.null(threshold_overrides_tsv) || identical(threshold_overrides_tsv, "")) {
    path <- inventory_env_or_default("CELL_COUNT_INVENTORY_THRESHOLD_OVERRIDE_TSV", "")
    if (!nzchar(path) || !file.exists(path)) {
      return(NULL)
    }
    threshold_overrides_tsv <- path
  }
  if (is.data.frame(threshold_overrides_tsv)) {
    return(threshold_overrides_tsv)
  }
  if (!file.exists(threshold_overrides_tsv) || file.info(threshold_overrides_tsv)$size == 0) {
    return(NULL)
  }
  utils::read.delim(threshold_overrides_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")
}

merge_threshold_row <- function(thresholds, row) {
  map <- c(
    min_cells_per_sample = "min_cells_per_sample",
    min_samples_per_group = "min_samples_per_group",
    min_total_umi = "min_total_umi",
    max_single_sample_frac = "max_single_sample_frac",
    min_cells_soft_floor = "min_cells_soft_floor"
  )
  for (name in names(map)) {
    col <- map[[name]]
    if (col %in% colnames(row) && !inventory_blank_or_dash(row[[col]][[1]])) {
      thresholds[[name]] <- inventory_numeric_or(row[[col]][[1]], thresholds[[name]])
    }
  }
  thresholds
}

resolve_thresholds <- function(comparison_row = NULL, threshold_overrides = NULL, cluster_id = NULL) {
  thresholds <- list(
    min_cells_per_sample = inventory_integer_or(inventory_env_or_default("CELL_COUNT_INVENTORY_MIN_CELLS_PER_SAMPLE", "20"), 20L),
    min_samples_per_group = inventory_integer_or(inventory_env_or_default("CELL_COUNT_INVENTORY_MIN_SAMPLES_PER_GROUP", "2"), 2L),
    min_total_umi = inventory_integer_or(inventory_env_or_default("CELL_COUNT_INVENTORY_MIN_TOTAL_UMI", "10000"), 10000L),
    max_single_sample_frac = inventory_numeric_or(inventory_env_or_default("CELL_COUNT_INVENTORY_MAX_SINGLE_SAMPLE_FRAC", "0.8"), 0.8),
    min_cells_soft_floor = inventory_integer_or(inventory_env_or_default("CELL_COUNT_INVENTORY_MIN_CELLS_SOFT_FLOOR", "5"), 5L)
  )

  overrides <- read_threshold_overrides(threshold_overrides)
  if (!is.null(overrides) && nrow(overrides) > 0 && "celltype" %in% colnames(overrides)) {
    star <- overrides[as.character(overrides$celltype) == "*", , drop = FALSE]
    if (nrow(star) > 0) {
      thresholds <- merge_threshold_row(thresholds, star[1, , drop = FALSE])
    }
    if (!is.null(cluster_id) && nzchar(as.character(cluster_id))) {
      specific <- overrides[as.character(overrides$celltype) == as.character(cluster_id), , drop = FALSE]
      if (nrow(specific) > 0) {
        thresholds <- merge_threshold_row(thresholds, specific[1, , drop = FALSE])
      }
    }
  }

  if (!is.null(comparison_row) && nrow(comparison_row) > 0) {
    if (!inventory_blank_or_dash(inventory_cell(comparison_row, "min_cells_override"))) {
      thresholds$min_cells_per_sample <- inventory_integer_or(inventory_cell(comparison_row, "min_cells_override"), thresholds$min_cells_per_sample)
    }
    if (!inventory_blank_or_dash(inventory_cell(comparison_row, "min_samples_override"))) {
      thresholds$min_samples_per_group <- inventory_integer_or(inventory_cell(comparison_row, "min_samples_override"), thresholds$min_samples_per_group)
    }
    if (!inventory_blank_or_dash(inventory_cell(comparison_row, "total_umi_override"))) {
      thresholds$min_total_umi <- inventory_integer_or(inventory_cell(comparison_row, "total_umi_override"), thresholds$min_total_umi)
    }
    if (!inventory_blank_or_dash(inventory_cell(comparison_row, "single_sample_frac_override"))) {
      thresholds$max_single_sample_frac <- inventory_numeric_or(inventory_cell(comparison_row, "single_sample_frac_override"), thresholds$max_single_sample_frac)
    }
  }
  thresholds
}

inventory_group_stats <- function(rows) {
  if (nrow(rows) == 0) {
    return(list(n_samples = 0L, min_cells = 0L, total_umi = 0, max_single_sample_frac = 0))
  }
  list(
    n_samples = length(unique(rows$sample_id)),
    min_cells = min(rows$n_cells, na.rm = TRUE),
    total_umi = sum(rows$total_umi, na.rm = TRUE),
    max_single_sample_frac = max(rows$fraction_of_cluster, na.rm = TRUE)
  )
}

inventory_parent_cluster_id <- function(inventory, cluster_id) {
  rows <- inventory[inventory$cluster_id == cluster_id, "parent_cluster_id", drop = TRUE]
  rows <- unique(rows[!is.na(rows) & nzchar(rows)])
  if (length(rows) == 0) NA_character_ else rows[[1]]
}

inventory_reason <- function(tier, failed_rules, stats_a, stats_b) {
  sprintf(
    "tier=%s; failed_rules=%s; samples=%s/%s; min_cells=%s/%s; max_single_sample_frac=%.3f",
    tier,
    if (length(failed_rules) == 0) "none" else paste(failed_rules, collapse = ","),
    stats_a$n_samples,
    stats_b$n_samples,
    stats_a$min_cells,
    stats_b$min_cells,
    max(stats_a$max_single_sample_frac, stats_b$max_single_sample_frac, na.rm = TRUE)
  )
}

classify_cluster_eligibility <- function(inventory, comparison_row, threshold_overrides = NULL) {
  cluster_id <- inventory_cell(comparison_row, "cluster_id")
  if (!nzchar(cluster_id)) {
    stop("classify_cluster_eligibility requires comparison_row$cluster_id", call. = FALSE)
  }
  group_a_id <- inventory_cell(comparison_row, "ident_1", inventory_cell(comparison_row, "group_a"))
  group_b_id <- inventory_cell(comparison_row, "ident_2", inventory_cell(comparison_row, "group_b"))
  comparison_id <- inventory_cell(comparison_row, "comparison_id", "comparison")
  thresholds <- resolve_thresholds(comparison_row, threshold_overrides, cluster_id)

  cluster_rows <- inventory[inventory$cluster_id == cluster_id, , drop = FALSE]
  group_a <- cluster_rows[cluster_rows$group_id == group_a_id, , drop = FALSE]
  group_b <- cluster_rows[cluster_rows$group_id == group_b_id, , drop = FALSE]
  stats_a <- inventory_group_stats(group_a)
  stats_b <- inventory_group_stats(group_b)
  max_single_sample_frac <- max(c(stats_a$max_single_sample_frac, stats_b$max_single_sample_frac), na.rm = TRUE)
  if (!is.finite(max_single_sample_frac)) {
    max_single_sample_frac <- 0
  }

  failed_rules <- character()
  if (stats_a$n_samples < thresholds$min_samples_per_group || stats_b$n_samples < thresholds$min_samples_per_group) {
    failed_rules <- c(failed_rules, "per_group_min_samples")
  }
  if (stats_a$min_cells < thresholds$min_cells_per_sample || stats_b$min_cells < thresholds$min_cells_per_sample) {
    failed_rules <- c(failed_rules, "per_sample_min_cells")
  }
  if (stats_a$total_umi < thresholds$min_total_umi || stats_b$total_umi < thresholds$min_total_umi) {
    failed_rules <- c(failed_rules, "total_umi_min")
  }
  if (max_single_sample_frac > thresholds$max_single_sample_frac) {
    failed_rules <- c(failed_rules, "single_sample_dominance")
  }

  parent_id <- inventory_parent_cluster_id(inventory, cluster_id)
  if (length(failed_rules) == 0) {
    tier <- "primary"
    action <- "run_pseudobulk_de"
  } else if ("single_sample_dominance" %in% failed_rules) {
    tier <- "candidate_only"
    action <- "report_as_candidate"
  } else if ("per_group_min_samples" %in% failed_rules &&
             (stats_a$n_samples == 0L || stats_b$n_samples == 0L)) {
    tier <- "skip"
    action <- "skip"
  } else if ("per_sample_min_cells" %in% failed_rules) {
    if (!is.na(parent_id) && nzchar(parent_id) && identical(tolower(inventory_env_or_default("CELL_COUNT_INVENTORY_ENABLE_MERGE_TO_PARENT", "yes")), "yes")) {
      tier <- "merge_to_parent"
      action <- paste0("merge_to:", parent_id)
    } else if (!"per_group_min_samples" %in% failed_rules &&
               stats_a$min_cells >= thresholds$min_cells_soft_floor &&
               stats_b$min_cells >= thresholds$min_cells_soft_floor) {
      tier <- "exploratory"
      action <- "run_pseudobulk_de_with_warning"
    } else if (stats_a$min_cells >= thresholds$min_cells_soft_floor && stats_b$min_cells >= thresholds$min_cells_soft_floor) {
      tier <- "module_score_only"
      action <- "run_module_score"
    } else {
      tier <- "skip"
      action <- "skip"
    }
  } else {
    tier <- "exploratory"
    action <- "run_pseudobulk_de_with_warning"
  }

  list(
    cluster_id = cluster_id,
    comparison_id = comparison_id,
    group_a = group_a_id,
    group_b = group_b_id,
    evidence_tier = tier,
    recommended_action = action,
    failed_rules = failed_rules,
    reason = inventory_reason(tier, failed_rules, stats_a, stats_b),
    evidence_table = rbind(group_a, group_b),
    threshold_resolved = thresholds,
    n_samples_a = stats_a$n_samples,
    n_samples_b = stats_b$n_samples,
    min_cells_a = stats_a$min_cells,
    min_cells_b = stats_b$min_cells,
    total_umi_a = stats_a$total_umi,
    total_umi_b = stats_b$total_umi,
    max_single_sample_frac = max_single_sample_frac,
    parent_cluster_id = parent_id
  )
}

as_eligibility_row <- function(result) {
  data.frame(
    comparison_id = result$comparison_id,
    cluster_id = result$cluster_id,
    group_a = result$group_a,
    group_b = result$group_b,
    evidence_tier = result$evidence_tier,
    recommended_action = result$recommended_action,
    failed_rules = paste(result$failed_rules, collapse = ","),
    reason = result$reason,
    n_samples_a = result$n_samples_a,
    n_samples_b = result$n_samples_b,
    min_cells_a = result$min_cells_a,
    min_cells_b = result$min_cells_b,
    total_umi_a = result$total_umi_a,
    total_umi_b = result$total_umi_b,
    max_single_sample_frac = result$max_single_sample_frac,
    parent_cluster_id = result$parent_cluster_id,
    min_cells_per_sample = result$threshold_resolved$min_cells_per_sample,
    min_samples_per_group = result$threshold_resolved$min_samples_per_group,
    min_total_umi = result$threshold_resolved$min_total_umi,
    max_single_sample_frac_threshold = result$threshold_resolved$max_single_sample_frac,
    min_cells_soft_floor = result$threshold_resolved$min_cells_soft_floor,
    inventory_version = inventory_env_or_default("MODULE_INVENTORY_VERSION", "1.0"),
    stringsAsFactors = FALSE
  )
}

classify_all_cluster_eligibilities <- function(inventory, comparisons_tsv, threshold_overrides_tsv = NULL) {
  comparisons <- if (is.data.frame(comparisons_tsv)) {
    comparisons_tsv
  } else {
    utils::read.delim(comparisons_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")
  }
  overrides <- read_threshold_overrides(threshold_overrides_tsv)
  rows <- list()
  for (cluster_id in unique(inventory$cluster_id)) {
    for (i in seq_len(nrow(comparisons))) {
      comparison_row <- comparisons[i, , drop = FALSE]
      comparison_row$cluster_id <- cluster_id
      rows[[length(rows) + 1L]] <- as_eligibility_row(classify_cluster_eligibility(inventory, comparison_row, overrides))
    }
  }
  if (length(rows) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}
