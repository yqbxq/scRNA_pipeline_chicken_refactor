empty_df_05 <- function(cols = character()) {
  out <- as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols), stringsAsFactors = FALSE)
  out
}

empty_marker_result_05 <- function() {
  empty_df_05(c(
    "gene", "comparison_id", "layer_id", "cluster_id", "annotation_label",
    "analysis_mode", "analysis_unit", "gene_program_role",
    "inference_scope", "inference_status", "subset_column", "subset_value"
  ))
}

empty_marker_summary_05 <- function() {
  empty_df_05(c(
    "comparison_id", "layer_id", "cluster_id", "annotation_label",
    "analysis_mode", "analysis_unit",
    "ident_1_n", "ident_2_n", "min_cells_per_group",
    "result_available", "status", "reason"
  ))
}

empty_gate_summary_05 <- function() {
  data.frame(
    group_id = character(0),
    sample_n = integer(0),
    biological_replicate_n = integer(0),
    pass = logical(0),
    stringsAsFactors = FALSE
  )
}

fill_blank_from_05 <- function(primary, fallback, default = "") {
  primary <- as.character(primary)
  fallback <- as.character(fallback)
  primary[is.na(primary)] <- ""
  fallback[is.na(fallback)] <- ""
  primary <- trimws(primary)
  fallback <- trimws(fallback)
  miss <- !nzchar(primary)
  primary[miss] <- fallback[miss]
  primary[!nzchar(primary)] <- default
  primary
}

standardize_design_metadata_df_05 <- function(meta_df) {
  meta_df <- as.data.frame(meta_df, stringsAsFactors = FALSE)
  required_cols <- c(
    "sample_id", "orig.ident", "group_id", "condition", "analysis_group",
    "group", "biological_replicate", "batch", "tissue"
  )
  for (col in required_cols) {
    if (!col %in% colnames(meta_df)) {
      meta_df[[col]] <- ""
    }
  }

  meta_df$sample_id <- fill_blank_from_05(meta_df$sample_id, meta_df$orig.ident, "")
  meta_df$group_id <- fill_blank_from_05(meta_df$group_id, meta_df$condition, "")
  meta_df$group_id <- fill_blank_from_05(meta_df$group_id, meta_df$analysis_group, "")
  meta_df$group_id <- fill_blank_from_05(meta_df$group_id, meta_df$group, "")
  meta_df$group_id <- fill_blank_from_05(meta_df$group_id, meta_df$sample_id, "")
  meta_df$condition <- fill_blank_from_05(meta_df$condition, meta_df$group_id, "")
  meta_df$analysis_group <- fill_blank_from_05(meta_df$analysis_group, meta_df$group_id, "")
  meta_df$group <- fill_blank_from_05(meta_df$group, meta_df$group_id, "")
  meta_df$biological_replicate <- fill_blank_from_05(meta_df$biological_replicate, meta_df$sample_id, "")
  meta_df$batch <- fill_blank_from_05(meta_df$batch, rep("default", nrow(meta_df)), "default")
  meta_df$tissue <- fill_blank_from_05(meta_df$tissue, rep("*", nrow(meta_df)), "*")
  meta_df
}

standardize_design_metadata_05 <- function(seu) {
  seu@meta.data <- standardize_design_metadata_df_05(seu@meta.data)
  seu
}

read_deg_comparison_sheet <- function(cfg) {
  expected_cols <- c(
    "comparison_id", "source_question_id", "contrast_axis", "ident_1", "ident_2",
    "enabled", "group_var", "batch_var", "layer_scope",
    "min_biological_replicates", "subset_column", "subset_value",
    "analysis_mode", "analysis_unit", "stat_level", "aggregation_group_var",
    "composition_group_var", "produces_gene_program", "gene_program_role",
    "force_exploratory", "min_cells_per_group", "logfc_threshold"
  )

  df <- read_tsv_optional(cfg$comparison_sheet)
  if (nrow(df) == 0 && nzchar(cfg$deg_ident_1 %||% "") && nzchar(cfg$deg_ident_2 %||% "")) {
    df <- data.frame(
      comparison_id = sprintf("%s_vs_%s", cfg$deg_ident_1, cfg$deg_ident_2),
      source_question_id = "",
      contrast_axis = "contrast_only",
      ident_1 = cfg$deg_ident_1,
      ident_2 = cfg$deg_ident_2,
      enabled = "yes",
      group_var = "group_id",
      batch_var = "batch",
      layer_scope = "*",
      min_biological_replicates = cfg$min_biological_replicates,
      subset_column = "",
      subset_value = "",
      analysis_mode = "condition_within_type",
      analysis_unit = "whole_layer",
      stat_level = "pseudobulk_formal_if_replicates",
      aggregation_group_var = "all_cells",
      composition_group_var = "",
      produces_gene_program = "yes",
      gene_program_role = "condition_deg",
      force_exploratory = "no",
      min_cells_per_group = cfg$deg_default_min_cells_per_group,
      logfc_threshold = cfg$deg_default_logfc_threshold,
      stringsAsFactors = FALSE
    )
  }
  if (nrow(df) == 0) {
    return(empty_df_05(expected_cols))
  }

  defaults <- list(
    enabled = "yes",
    source_question_id = "",
    contrast_axis = "",
    group_var = "group_id",
    batch_var = "batch",
    layer_scope = "*",
    min_biological_replicates = cfg$min_biological_replicates,
    subset_column = "",
    subset_value = "",
    analysis_mode = "",
    analysis_unit = "",
    stat_level = "",
    aggregation_group_var = "",
    composition_group_var = "",
    produces_gene_program = "yes",
    gene_program_role = "",
    force_exploratory = "no",
    min_cells_per_group = cfg$deg_default_min_cells_per_group,
    logfc_threshold = cfg$deg_default_logfc_threshold
  )
  for (col in expected_cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- if (col %in% names(defaults)) defaults[[col]] else ""
    }
  }
  df <- df[, expected_cols, drop = FALSE]
  for (col in setdiff(expected_cols, c("min_biological_replicates", "min_cells_per_group", "logfc_threshold"))) {
    df[[col]] <- vapply(df[[col]], normalize_scalar_value, character(1))
  }
  df$enabled <- normalize_flag(df$enabled, "yes")
  df$group_var[!nzchar(df$group_var)] <- "group_id"
  df$batch_var[!nzchar(df$batch_var)] <- "batch"
  df$layer_scope[!nzchar(df$layer_scope)] <- "*"
  df$analysis_mode[!nzchar(df$analysis_mode)] <- ifelse(df$contrast_axis %in% "composition", "composition", "subtype_pairwise")
  df$analysis_unit[!nzchar(df$analysis_unit)] <- ifelse(df$analysis_mode == "composition", "sample_level", "whole_layer")
  df$stat_level[!nzchar(df$stat_level)] <- ifelse(
    df$analysis_mode == "condition_within_type",
    "pseudobulk_formal_if_replicates",
    ifelse(df$analysis_mode == "composition", "composition_formal_if_replicates", "cell_level_exploratory")
  )
  df$aggregation_group_var[is.na(df$aggregation_group_var)] <- ""
  df$composition_group_var[is.na(df$composition_group_var)] <- ""
  df$produces_gene_program <- normalize_flag(df$produces_gene_program, "yes")
  df$gene_program_role[!nzchar(df$gene_program_role)] <- ifelse(df$analysis_mode == "composition", "none", "subtype_pairwise_deg")
  df$force_exploratory <- normalize_flag(df$force_exploratory, "no")
  df$min_biological_replicates <- suppressWarnings(as.integer(df$min_biological_replicates))
  df$min_biological_replicates[is.na(df$min_biological_replicates)] <- cfg$min_biological_replicates
  df$min_cells_per_group <- suppressWarnings(as.integer(df$min_cells_per_group))
  df$min_cells_per_group[is.na(df$min_cells_per_group)] <- cfg$deg_default_min_cells_per_group
  df$logfc_threshold <- suppressWarnings(as.numeric(df$logfc_threshold))
  df$logfc_threshold[is.na(df$logfc_threshold)] <- cfg$deg_default_logfc_threshold
  df[df$enabled != "no", , drop = FALSE]
}

deg_layer_status <- function(cfg) {
  df <- read_tsv_optional(cfg$layer_status_file)
  if (nrow(df) == 0) {
    stop(sprintf("Missing layer status table: %s", cfg$layer_status_file), call. = FALSE)
  }
  expected <- c("layer_id", "status", "layer_role", "cluster_column", "clustered_rds", "annotated_rds")
  for (col in expected) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  df <- df[nzchar(df$layer_id), expected, drop = FALSE]
  df$layer_role[!nzchar(df$layer_role) & df$layer_id == cfg$panorama_layer_id] <- "panorama"
  df$layer_role[!nzchar(df$layer_role)] <- "subcluster"

  object_path <- ifelse(nzchar(df$annotated_rds), df$annotated_rds, df$clustered_rds)
  keep <- nzchar(object_path) & file.exists(object_path)
  df[keep, , drop = FALSE]
}

resolve_cluster_var_05 <- function(seu, layer_row) {
  meta_cols <- colnames(seu@meta.data)
  layer_id <- layer_row$layer_id[[1]]
  candidates <- c(
    "cluster_id",
    normalize_scalar_value(layer_row$cluster_column[[1]]),
    paste0(layer_id, "_cluster"),
    "seurat_clusters"
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[candidates %in% meta_cols]
  if (length(hit) == 0) {
    stop(sprintf("Layer %s has no usable cluster column.", layer_id), call. = FALSE)
  }
  hit[[1]]
}

load_layer_for_deg <- function(cfg, layer_row) {
  obj_path <- normalize_scalar_value(layer_row$annotated_rds[[1]])
  if (!nzchar(obj_path) || !file.exists(obj_path)) {
    obj_path <- normalize_scalar_value(layer_row$clustered_rds[[1]])
  }
  if (!nzchar(obj_path) || !file.exists(obj_path)) {
    stop(sprintf("Layer %s has no readable object path.", layer_row$layer_id[[1]]), call. = FALSE)
  }

  seu <- readRDS(obj_path)
  seu <- maybe_join_layers(seu)
  seu <- standardize_design_metadata_05(seu)
  cluster_var <- resolve_cluster_var_05(seu, layer_row)
  seu$cluster_id <- as.character(seu@meta.data[[cluster_var]])
  if (!"annotation_label" %in% colnames(seu@meta.data)) {
    layer_cell_type <- paste0(layer_row$layer_id[[1]], "_cell_type")
    if (layer_cell_type %in% colnames(seu@meta.data)) {
      seu$annotation_label <- as.character(seu@meta.data[[layer_cell_type]])
    } else if ("cell_type" %in% colnames(seu@meta.data)) {
      seu$annotation_label <- as.character(seu$cell_type)
    } else {
      seu$annotation_label <- seu$cluster_id
    }
  }
  seu@misc$deg_object_path <- obj_path
  seu@misc$deg_cluster_var <- cluster_var
  seu
}

resolve_comparison_vars_05 <- function(seu, comparison_row) {
  group_var <- normalize_scalar_value(comparison_row$group_var[[1]], "group_id")
  batch_var <- normalize_scalar_value(comparison_row$batch_var[[1]], "batch")
  subset_column <- normalize_scalar_value(comparison_row$subset_column[[1]])
  subset_value <- normalize_scalar_value(comparison_row$subset_value[[1]])
  list(
    comparison_id = normalize_scalar_value(comparison_row$comparison_id[[1]]),
    ident_1 = normalize_scalar_value(comparison_row$ident_1[[1]]),
    ident_2 = normalize_scalar_value(comparison_row$ident_2[[1]]),
    group_var = group_var,
    batch_var = batch_var,
    subset_column = subset_column,
    subset_value = subset_value,
    subset_values = split_csv_local(subset_value),
    analysis_mode = normalize_scalar_value(comparison_row$analysis_mode[[1]], "subtype_pairwise"),
    analysis_unit = normalize_scalar_value(comparison_row$analysis_unit[[1]], "whole_layer"),
    stat_level = normalize_scalar_value(comparison_row$stat_level[[1]], "cell_level_exploratory"),
    aggregation_group_var = normalize_scalar_value(comparison_row$aggregation_group_var[[1]]),
    composition_group_var = normalize_scalar_value(comparison_row$composition_group_var[[1]]),
    produces_gene_program = normalize_flag(comparison_row$produces_gene_program[[1]], "yes"),
    gene_program_role = normalize_scalar_value(comparison_row$gene_program_role[[1]], "subtype_pairwise_deg"),
    force_exploratory = normalize_flag(comparison_row$force_exploratory[[1]], "no") %in% c("yes", "true", "1", "on"),
    min_biological_replicates = as.integer(comparison_row$min_biological_replicates[[1]]),
    min_cells_per_group = as.integer(comparison_row$min_cells_per_group[[1]]),
    logfc_threshold = as.numeric(comparison_row$logfc_threshold[[1]])
  )
}

subset_cells_for_comparison <- function(seu, vars) {
  meta <- seu@meta.data
  keep <- rep(TRUE, nrow(meta))
  names(keep) <- rownames(meta)

  if (nzchar(vars$subset_column) || nzchar(vars$subset_value)) {
    if (!nzchar(vars$subset_column) || !vars$subset_column %in% colnames(meta)) {
      return(list(object = NULL, status = "subset_column_missing", reason = sprintf("subset_column missing: %s", vars$subset_column), subset_n = 0L))
    }
    if (length(vars$subset_values) == 0) {
      return(list(object = NULL, status = "subset_value_empty", reason = "subset_value is empty", subset_n = 0L))
    }
    keep <- keep & as.character(meta[[vars$subset_column]]) %in% vars$subset_values
  }

  if (!vars$group_var %in% colnames(meta)) {
    return(list(object = NULL, status = "missing_group_var", reason = sprintf("group_var missing: %s", vars$group_var), subset_n = sum(keep)))
  }

  if (identical(vars$ident_2, "__rest__")) {
    keep <- keep & nzchar(as.character(meta[[vars$group_var]]))
  } else {
    keep <- keep & as.character(meta[[vars$group_var]]) %in% c(vars$ident_1, vars$ident_2)
  }
  if (!any(keep)) {
    return(list(object = NULL, status = "empty_subset", reason = "comparison matched 0 cells", subset_n = 0L))
  }
  obj <- subset(seu, cells = names(keep)[keep])
  obj <- maybe_join_layers(obj)
  list(object = obj, status = "ok", reason = "", subset_n = ncol(obj))
}

prepare_comparison_group_05 <- function(seu, vars) {
  group_var <- vars$group_var
  if (!identical(vars$ident_2, "__rest__")) {
    return(list(object = seu, group_var = group_var))
  }
  meta <- seu@meta.data
  rest_col <- ".deg_comparison_group"
  values <- ifelse(as.character(meta[[group_var]]) == vars$ident_1, vars$ident_1, "__rest__")
  seu[[rest_col]] <- values
  list(object = seu, group_var = rest_col)
}

resolve_aggregation_group_var_05 <- function(seu, vars) {
  aggregation_group_var <- normalize_scalar_value(vars$aggregation_group_var)
  if (!nzchar(aggregation_group_var) || identical(aggregation_group_var, "all_cells")) {
    seu$deg_aggregation_group <- "all_cells"
    return(list(object = seu, group_var = "deg_aggregation_group", status = "ok", reason = ""))
  }
  if (!aggregation_group_var %in% colnames(seu@meta.data)) {
    return(list(object = NULL, group_var = "", status = "missing_aggregation_group_var", reason = sprintf("aggregation_group_var missing: %s", aggregation_group_var)))
  }
  seu$deg_aggregation_group <- as.character(seu@meta.data[[aggregation_group_var]])
  list(object = seu, group_var = "deg_aggregation_group", status = "ok", reason = "")
}

resolve_composition_group_var_05 <- function(seu, vars) {
  composition_group_var <- normalize_scalar_value(vars$composition_group_var, "cluster_id")
  if (identical(composition_group_var, "cell_class_for_composition")) {
    source_col <- if ("cell_subtype" %in% colnames(seu@meta.data)) "cell_subtype" else if ("cell_type" %in% colnames(seu@meta.data)) "cell_type" else ""
    if (!nzchar(source_col)) {
      return(list(object = NULL, group_var = "", status = "missing_composition_group_var", reason = "cell_class_for_composition requires cell_subtype or cell_type"))
    }
    values <- as.character(seu@meta.data[[source_col]])
    gc_tokens <- c("GC", "pGC", "eGC", "rgGC", "lGC")
    values[values %in% gc_tokens] <- "GC"
    values[values == "TC"] <- "TC"
    seu$cell_class_for_composition <- values
    return(list(object = seu, group_var = "cell_class_for_composition", status = "ok", reason = ""))
  }
  if (!composition_group_var %in% colnames(seu@meta.data)) {
    return(list(object = NULL, group_var = "", status = "missing_composition_group_var", reason = sprintf("composition_group_var missing: %s", composition_group_var)))
  }
  list(object = seu, group_var = composition_group_var, status = "ok", reason = "")
}

replicate_gate_summary_05 <- function(meta_df, vars, replicate_var = "biological_replicate", sample_var = "sample_id") {
  required <- c(vars$group_var, replicate_var, sample_var)
  missing_cols <- setdiff(required, colnames(meta_df))
  if (length(missing_cols) > 0) {
    summary_df <- data.frame(
      group_id = c(vars$ident_1, vars$ident_2),
      sample_n = NA_integer_,
      biological_replicate_n = NA_integer_,
      pass = FALSE,
      stringsAsFactors = FALSE
    )
    return(list(pass = FALSE, reason = sprintf("missing metadata fields: %s", paste(missing_cols, collapse = ", ")), summary = summary_df))
  }

  group_rows <- lapply(c(vars$ident_1, vars$ident_2), function(group_id) {
    grp <- meta_df[as.character(meta_df[[vars$group_var]]) == group_id, , drop = FALSE]
    reps <- unique(normalize_flag(grp[[replicate_var]], ""))
    reps <- reps[nzchar(reps)]
    samples <- unique(normalize_flag(grp[[sample_var]], ""))
    samples <- samples[nzchar(samples)]
    data.frame(
      group_id = group_id,
      sample_n = length(samples),
      biological_replicate_n = length(reps),
      pass = length(reps) >= vars$min_biological_replicates,
      stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, group_rows)
  list(
    pass = all(summary_df$pass),
    reason = if (all(summary_df$pass)) {
      sprintf("each group has biological_replicate >= %s", vars$min_biological_replicates)
    } else {
      sprintf("at least one group has biological_replicate < %s", vars$min_biological_replicates)
    },
    summary = summary_df
  )
}

annotation_warning_banner_05 <- function(gate_summary, vars, forced = FALSE) {
  prefix <- if (forced) "EXPLICIT FORCE_EXPLORATORY OVERRIDE. " else ""
  if (nrow(gate_summary) == 0) {
    return(paste0(prefix, "WARNING: replicate metadata are missing; results are exploratory only and cannot support formal conclusions."))
  }
  rep_counts <- suppressWarnings(as.integer(gate_summary$biological_replicate_n))
  rep_counts <- rep_counts[is.finite(rep_counts)]
  if (length(rep_counts) > 0 && all(rep_counts == 1)) {
    return(paste0(prefix, "WARNING: each group has N=1; results are exploratory only and cannot support formal conclusions."))
  }
  paste0(
    prefix,
    sprintf("WARNING: at least one group has biological replicates below %s; results are exploratory only and cannot support formal conclusions.", vars$min_biological_replicates)
  )
}

cluster_label_for_05 <- function(seu, cluster_id) {
  hit <- as.character(seu$cluster_id) == cluster_id
  labels <- unique(as.character(seu$annotation_label[hit]))
  labels <- labels[nzchar(labels)]
  if (length(labels) == 0) cluster_id else labels[[1]]
}

aggregate_cluster_sample_counts_05 <- function(seu, cluster_var = "cluster_id", sample_var = "sample_id", group_var = "group_id", batch_var = "batch", replicate_var = "biological_replicate") {
  counts <- get_assay_matrix(seu, assay = "RNA", type = "counts")
  meta_df <- standardize_design_metadata_df_05(seu@meta.data)
  meta_df <- meta_df[colnames(counts), , drop = FALSE]
  if (!cluster_var %in% colnames(meta_df)) {
    stop(sprintf("Object is missing cluster column for aggregation: %s", cluster_var), call. = FALSE)
  }

  key_df <- unique(meta_df[, c(cluster_var, sample_var), drop = FALSE])
  key_df <- key_df[nzchar(as.character(key_df[[cluster_var]])) & nzchar(as.character(key_df[[sample_var]])), , drop = FALSE]
  if (nrow(key_df) == 0) {
    return(list(counts = counts[, FALSE], metadata = empty_df_05(c("pseudobulk_id", "cluster_id", "sample_id", "group_id", "batch", "biological_replicate", "annotation_label", "n_cells"))))
  }
  key_df <- key_df[order(key_df[[cluster_var]], key_df[[sample_var]]), , drop = FALSE]
  colnames(key_df) <- c("cluster_id", "sample_id")

  agg_cols <- lapply(seq_len(nrow(key_df)), function(idx) {
    hit <- meta_df[[cluster_var]] == key_df$cluster_id[idx] & meta_df[[sample_var]] == key_df$sample_id[idx]
    Matrix::rowSums(counts[, rownames(meta_df)[hit], drop = FALSE])
  })
  agg_mat <- do.call(cbind, agg_cols)
  if (is.null(dim(agg_mat))) {
    agg_mat <- matrix(agg_mat, ncol = 1)
  }
  rownames(agg_mat) <- rownames(counts)
  colnames(agg_mat) <- paste(key_df$cluster_id, key_df$sample_id, sep = "__")
  agg_mat <- Matrix::Matrix(agg_mat, sparse = TRUE)

  meta_rows <- lapply(seq_len(nrow(key_df)), function(idx) {
    hit <- meta_df[[cluster_var]] == key_df$cluster_id[idx] & meta_df[[sample_var]] == key_df$sample_id[idx]
    subset_df <- meta_df[hit, , drop = FALSE]
    data.frame(
      pseudobulk_id = paste(key_df$cluster_id[idx], key_df$sample_id[idx], sep = "__"),
      cluster_id = key_df$cluster_id[idx],
      sample_id = key_df$sample_id[idx],
      group_id = normalize_scalar_value(subset_df[[group_var]][1], ""),
      batch = if (batch_var %in% colnames(subset_df)) normalize_scalar_value(subset_df[[batch_var]][1], "") else "",
      biological_replicate = normalize_scalar_value(subset_df[[replicate_var]][1], ""),
      annotation_label = if ("annotation_label" %in% colnames(subset_df)) normalize_scalar_value(subset_df$annotation_label[1], "") else "",
      n_cells = nrow(subset_df),
      stringsAsFactors = FALSE
    )
  })
  list(counts = agg_mat, metadata = do.call(rbind, meta_rows))
}

sample_level_proportion_summary_05 <- function(seu, cluster_var = "cluster_id", sample_var = "sample_id", group_var = "group_id") {
  meta_df <- standardize_design_metadata_df_05(seu@meta.data)
  if (!cluster_var %in% colnames(meta_df)) {
    stop(sprintf("Object is missing cluster column for composition: %s", cluster_var), call. = FALSE)
  }
  summary_df <- meta_df %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::mutate(cluster_value = .data[[cluster_var]]) %>%
    dplyr::count(.data[[sample_var]], .data[[group_var]], cluster_value, name = "cell_number") %>%
    dplyr::group_by(.data[[sample_var]]) %>%
    dplyr::mutate(sample_total = sum(cell_number), proportion = cell_number / sample_total) %>%
    dplyr::ungroup()
  colnames(summary_df)[1:3] <- c("sample_id", "group_id", "cluster_id")
  summary_df
}

read_tsv_if_path_05 <- function(path) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path) || !file.exists(path)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  read_tsv_optional(path)
}

count_significant_rows_05 <- function(path, alpha = 0.05) {
  df <- read_tsv_if_path_05(path)
  if (nrow(df) == 0) {
    return(list(total_n = 0L, significant_n = 0L))
  }
  candidates <- c("p_val_adj", "p_adj", "adj.P.Val", "FDR", "fdr", "p_adj.loc", "p_adj.glb")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    return(list(total_n = nrow(df), significant_n = NA_integer_))
  }
  p <- suppressWarnings(as.numeric(df[[hit[[1]]]]))
  list(total_n = nrow(df), significant_n = sum(is.finite(p) & p < alpha))
}

render_status_for_report_05 <- function(status) {
  status <- normalize_scalar_value(status, "skipped")
  if (!nzchar(status)) "skipped" else status
}
