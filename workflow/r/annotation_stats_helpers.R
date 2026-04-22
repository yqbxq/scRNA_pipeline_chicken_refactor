source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))

if (!exists("sanitize_layer_id", mode = "function")) {
  sanitize_layer_id <- function(x) {
    value <- trimws(as.character(x))
    if (!nzchar(value)) {
      return("")
    }
    gsub("[^A-Za-z0-9._-]+", "_", value)
  }
}

trim_character <- function(x, default = "") {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  if (!missing(default)) {
    x[!nzchar(x)] <- default
  }
  x
}

fill_blank_from <- function(primary, fallback, default = "") {
  primary <- trim_character(primary, "")
  fallback <- trim_character(fallback, "")
  miss <- !nzchar(primary)
  primary[miss] <- fallback[miss]
  if (!missing(default)) {
    primary[!nzchar(primary)] <- default
  }
  primary
}

annotation_confidence_levels <- function() {
  c("未定", "暂定", "确定")
}

cap_annotation_confidence <- function(value, ceiling = "") {
  levels <- annotation_confidence_levels()
  value <- normalize_scalar_value(value, "未定")
  ceiling <- normalize_scalar_value(ceiling, "")
  if (!value %in% levels) {
    value <- "未定"
  }
  if (!nzchar(ceiling) || !ceiling %in% levels) {
    return(value)
  }
  levels[min(match(value, levels), match(ceiling, levels))]
}

standardize_design_metadata_df <- function(meta_df) {
  meta_df <- as.data.frame(meta_df, stringsAsFactors = FALSE)
  required_cols <- c(
    "sample_id", "orig.ident", "group_id", "condition", "analysis_group", "group",
    "biological_replicate", "batch", "tissue"
  )
  for (col in required_cols) {
    if (!col %in% colnames(meta_df)) {
      meta_df[[col]] <- ""
    }
  }

  meta_df$sample_id <- fill_blank_from(meta_df$sample_id, meta_df$orig.ident, "")
  meta_df$group_id <- fill_blank_from(meta_df$group_id, meta_df$condition, "")
  meta_df$group_id <- fill_blank_from(meta_df$group_id, meta_df$analysis_group, "")
  meta_df$group_id <- fill_blank_from(meta_df$group_id, meta_df$group, "")
  meta_df$group_id <- fill_blank_from(meta_df$group_id, meta_df$sample_id, "")
  meta_df$condition <- fill_blank_from(meta_df$condition, meta_df$group_id, "")
  meta_df$analysis_group <- fill_blank_from(meta_df$analysis_group, meta_df$group_id, "")
  meta_df$group <- fill_blank_from(meta_df$group, meta_df$group_id, "")
  meta_df$biological_replicate <- fill_blank_from(meta_df$biological_replicate, meta_df$sample_id, "")
  meta_df$batch <- trim_character(meta_df$batch, "default")
  meta_df$tissue <- trim_character(meta_df$tissue, "*")

  meta_df
}

standardize_design_metadata <- function(seu) {
  seu@meta.data <- standardize_design_metadata_df(seu@meta.data)
  seu
}

read_layer_status_table <- function(cfg) {
  path <- file.path(cfg$table_dir, "layer_status.tsv")
  df <- read_tsv_optional(path)
  if (nrow(df) == 0) {
    stop(sprintf("缺少 layer status 表: %s", path), call. = FALSE)
  }
  required <- c("layer_id", "status", "clustered_rds", "cluster_column")
  for (col in required) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  df
}

built_layer_status <- function(cfg) {
  df <- read_layer_status_table(cfg)
  df[df$status %in% c("built", "reused"), , drop = FALSE]
}

layer_annotation_paths <- function(cfg, layer_id) {
  safe_id <- sanitize_layer_id(layer_id)
  list(
    layer_id = safe_id,
    table_dir = file.path(cfg$table_dir, "annotation", "layers", safe_id),
    report_dir = file.path(cfg$eda_report_dir, "annotation", "layers", safe_id),
    annotation_table_tsv = file.path(cfg$table_dir, "annotation", "layers", safe_id, "annotation_table.tsv"),
    evidence_table_tsv = file.path(cfg$table_dir, "annotation", "layers", safe_id, "annotation_evidence.tsv"),
    module_score_tsv = file.path(cfg$table_dir, "annotation", "layers", safe_id, "module_score_summary.tsv"),
    report_md = file.path(cfg$eda_report_dir, "annotation", "layers", safe_id, "report.md"),
    top_markers_dotplot_png = file.path(cfg$eda_report_dir, "annotation", "layers", safe_id, "top_markers_dotplot.png"),
    panel_dotplot_png = file.path(cfg$eda_report_dir, "annotation", "layers", safe_id, "panel_validation_dotplot.png"),
    module_score_heatmap_png = file.path(cfg$eda_report_dir, "annotation", "layers", safe_id, "module_score_heatmap.png"),
    annotated_rds = file.path(cfg$checkpoint_dir, "layers", safe_id, sprintf("%s_after_annotation.rds", safe_id))
  )
}

marker_discovery_paths <- function(cfg, layer_id, comparison_id = NULL) {
  safe_layer <- sanitize_layer_id(layer_id)
  if (is.null(comparison_id) || !nzchar(comparison_id)) {
    return(list(
      table_dir = file.path(cfg$table_dir, "marker_discovery", "layers", safe_layer),
      cluster_markers_tsv = file.path(cfg$table_dir, "marker_discovery", "layers", safe_layer, "cluster_markers.tsv")
    ))
  }
  safe_cmp <- sanitize_layer_id(comparison_id)
  list(
    table_dir = file.path(cfg$table_dir, "marker_discovery", "comparisons", safe_cmp, safe_layer),
    exploratory_tsv = file.path(cfg$table_dir, "marker_discovery", "comparisons", safe_cmp, safe_layer, "exploratory_cell_level_findmarkers.tsv"),
    exploratory_summary_tsv = file.path(cfg$table_dir, "marker_discovery", "comparisons", safe_cmp, safe_layer, "exploratory_summary.tsv")
  )
}

pseudobulk_paths <- function(cfg, layer_id, comparison_id = NULL) {
  safe_layer <- sanitize_layer_id(layer_id)
  if (is.null(comparison_id) || !nzchar(comparison_id)) {
    return(list(
      table_dir = file.path(cfg$table_dir, "pseudobulk_ds", "layers", safe_layer),
      aggregation_rds = file.path(cfg$table_dir, "pseudobulk_ds", "layers", safe_layer, "pseudobulk_counts.rds"),
      aggregation_tsv = file.path(cfg$table_dir, "pseudobulk_ds", "layers", safe_layer, "pseudobulk_column_metadata.tsv"),
      aggregation_manifest_tsv = file.path(cfg$table_dir, "pseudobulk_ds", "layers", safe_layer, "aggregation_manifest.tsv")
    ))
  }
  safe_cmp <- sanitize_layer_id(comparison_id)
  list(
    table_dir = file.path(cfg$table_dir, "pseudobulk_ds", "comparisons", safe_cmp, safe_layer),
    ds_results_tsv = file.path(cfg$table_dir, "pseudobulk_ds", "comparisons", safe_cmp, safe_layer, "formal_pseudobulk_ds.tsv"),
    gate_summary_tsv = file.path(cfg$table_dir, "pseudobulk_ds", "comparisons", safe_cmp, safe_layer, "replicate_gate_summary.tsv"),
    status_tsv = file.path(cfg$table_dir, "pseudobulk_ds", "comparisons", safe_cmp, safe_layer, "ds_status.tsv"),
    result_rds = file.path(cfg$table_dir, "pseudobulk_ds", "comparisons", safe_cmp, safe_layer, "formal_pseudobulk_ds.rds")
  )
}

composition_paths <- function(cfg, layer_id, comparison_id = NULL) {
  safe_layer <- sanitize_layer_id(layer_id)
  if (is.null(comparison_id) || !nzchar(comparison_id)) {
    return(list(
      table_dir = file.path(cfg$table_dir, "composition", "layers", safe_layer),
      proportion_tsv = file.path(cfg$table_dir, "composition", "layers", safe_layer, "sample_level_proportions.tsv"),
      proportion_plot_png = file.path(cfg$eda_report_dir, "composition", "layers", safe_layer, "sample_level_proportions.png")
    ))
  }
  safe_cmp <- sanitize_layer_id(comparison_id)
  list(
    table_dir = file.path(cfg$table_dir, "composition", "comparisons", safe_cmp, safe_layer),
    formal_results_tsv = file.path(cfg$table_dir, "composition", "comparisons", safe_cmp, safe_layer, "formal_propeller.tsv"),
    gate_summary_tsv = file.path(cfg$table_dir, "composition", "comparisons", safe_cmp, safe_layer, "replicate_gate_summary.tsv"),
    status_tsv = file.path(cfg$table_dir, "composition", "comparisons", safe_cmp, safe_layer, "composition_status.tsv")
  )
}

ensure_parent_dirs <- function(paths) {
  dirs <- unique(dirname(unlist(paths, use.names = FALSE)))
  dirs <- dirs[nzchar(dirs)]
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

read_marker_panel_rows <- function(cfg, layer_id = NULL, tissue = NULL) {
  if (!dir.exists(cfg$marker_panel_dir)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  panel_files <- list.files(cfg$marker_panel_dir, pattern = "\\.tsv$", full.names = TRUE, ignore.case = TRUE)
  if (length(panel_files) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  parsed <- lapply(panel_files, function(path) {
    df <- read_tsv_optional(path)
    if (nrow(df) == 0) {
      return(NULL)
    }
    if (!"layer_id" %in% colnames(df)) {
      df$layer_id <- "*"
    }
    if (!"celltype" %in% colnames(df) && "cell_type" %in% colnames(df)) {
      df$celltype <- df$cell_type
    }
    if (!"evidence_source" %in% colnames(df)) {
      df$evidence_source <- "unspecified"
    }
    for (col in c("tissue", "panel_name", "evidence_note", "confidence_ceiling")) {
      if (!col %in% colnames(df)) {
        df[[col]] <- ""
      }
    }
    required <- c("layer_id", "celltype", "gene", "evidence_source")
    if (!all(required %in% colnames(df))) {
      warning(sprintf("跳过 marker panel（缺少必需列）: %s", path), call. = FALSE)
      return(NULL)
    }
    df <- df[, c(required, "tissue", "panel_name", "evidence_note", "confidence_ceiling"), drop = FALSE]
    for (col in colnames(df)) {
      df[[col]] <- trim_character(df[[col]], "")
    }
    df$panel_file <- basename(path)
    df <- df[nzchar(df$celltype) & nzchar(df$gene), , drop = FALSE]
    if (nrow(df) == 0) {
      return(NULL)
    }
    df
  })

  parsed <- Filter(Negate(is.null), parsed)
  if (length(parsed) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  out <- unique(do.call(rbind, parsed))
  if (!is.null(layer_id) && nzchar(layer_id)) {
    out <- out[out$layer_id %in% c("*", layer_id), , drop = FALSE]
  }
  if (!is.null(tissue) && nzchar(tissue) && "tissue" %in% colnames(out)) {
    out <- out[out$tissue %in% c("", "*", tissue), , drop = FALSE]
  }
  out
}

panel_gene_sets <- function(panel_df) {
  if (nrow(panel_df) == 0) {
    return(list())
  }
  split(unique(panel_df$gene), panel_df$celltype)
}

empty_marker_table <- function() {
  data.frame(
    p_val = numeric(0),
    avg_log2FC = numeric(0),
    pct.1 = numeric(0),
    pct.2 = numeric(0),
    p_val_adj = numeric(0),
    cluster = character(0),
    gene = character(0),
    stringsAsFactors = FALSE
  )
}

run_cluster_marker_discovery <- function(seu, cluster_var = "seurat_clusters", only_pos = TRUE, min_pct = 0.25, logfc_threshold = 0) {
  Idents(seu) <- cluster_var
  cluster_levels <- unique(as.character(Idents(seu)))
  cluster_levels <- cluster_levels[nzchar(cluster_levels)]
  if (length(cluster_levels) < 2) {
    return(empty_marker_table())
  }
  markers <- tryCatch(
    FindAllMarkers(
      seu,
      only.pos = only_pos,
      min.pct = min_pct,
      logfc.threshold = logfc_threshold,
      test.use = "wilcox",
      verbose = FALSE
    ),
    error = function(e) empty_marker_table()
  )
  if (!"gene" %in% colnames(markers)) {
    markers$gene <- rownames(markers)
  }
  if (!"cluster" %in% colnames(markers)) {
    markers$cluster <- as.character(Idents(seu))[match(rownames(markers), colnames(seu))]
  }
  markers$cluster <- as.character(markers$cluster)
  markers$gene <- as.character(markers$gene)
  markers
}

top_cluster_markers <- function(marker_df, cluster_id, top_n = 10) {
  if (nrow(marker_df) == 0) {
    return(marker_df)
  }
  df <- marker_df[marker_df$cluster == cluster_id, , drop = FALSE]
  if (nrow(df) == 0) {
    return(df)
  }
  df <- df[order(-df$avg_log2FC, df$p_val_adj), , drop = FALSE]
  df[seq_len(min(top_n, nrow(df))), , drop = FALSE]
}

format_data_evidence <- function(marker_df, top_n = 5) {
  if (nrow(marker_df) == 0) {
    return("no_positive_markers")
  }
  marker_df <- marker_df[order(-marker_df$avg_log2FC, marker_df$p_val_adj), , drop = FALSE]
  marker_df <- marker_df[seq_len(min(top_n, nrow(marker_df))), , drop = FALSE]
  paste(sprintf("%s(%.2f)", marker_df$gene, marker_df$avg_log2FC), collapse = ", ")
}

panel_overlap_summary <- function(cluster_marker_df, panel_df, top_n = 20) {
  if (nrow(panel_df) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  top_genes <- unique(cluster_marker_df$gene[seq_len(min(top_n, nrow(cluster_marker_df)))])
  top_genes <- top_genes[nzchar(top_genes)]
  split_rows <- split(panel_df, panel_df$celltype)
  out <- lapply(names(split_rows), function(label) {
    df <- split_rows[[label]]
    overlap_df <- df[df$gene %in% top_genes, , drop = FALSE]
    overlap_genes <- unique(overlap_df$gene)
    data.frame(
      celltype = label,
      overlap_n = length(overlap_genes),
      overlap_genes = paste(overlap_genes, collapse = ","),
      evidence_source = paste(unique(overlap_df$evidence_source), collapse = ","),
      confidence_ceiling = paste(unique(overlap_df$confidence_ceiling[nzchar(overlap_df$confidence_ceiling)]), collapse = ","),
      panel_files = paste(unique(df$panel_file), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out[order(-out$overlap_n, out$celltype), , drop = FALSE]
}

compute_module_score_summary <- function(seu, cluster_var, panel_df, seed = 42) {
  if (nrow(panel_df) == 0) {
    return(list(object = seu, summary = data.frame(stringsAsFactors = FALSE)))
  }
  gene_sets <- panel_gene_sets(panel_df)
  if (length(gene_sets) == 0) {
    return(list(object = seu, summary = data.frame(stringsAsFactors = FALSE)))
  }

  score_cols <- list()
  expr_mat <- get_assay_matrix(seu, assay = "RNA", type = "data")
  for (label in names(gene_sets)) {
    genes <- intersect(unique(gene_sets[[label]]), rownames(seu))
    if (length(genes) == 0) {
      next
    }
    safe_label <- gsub("[^A-Za-z0-9]+", "_", label)
    score_col <- paste0("module_check_", safe_label)
    seu[[score_col]] <- Matrix::colMeans(expr_mat[genes, , drop = FALSE])
    score_cols[[label]] <- score_col
  }

  if (length(score_cols) == 0) {
    return(list(object = seu, summary = data.frame(stringsAsFactors = FALSE)))
  }

  meta_df <- seu@meta.data
  clusters <- sort(unique(as.character(meta_df[[cluster_var]])))
  summary_rows <- lapply(clusters, function(cluster_id) {
    hit <- as.character(meta_df[[cluster_var]]) == cluster_id
    scores <- vapply(score_cols, function(col) stats::median(meta_df[[col]][hit]), numeric(1))
    score_order <- order(scores, decreasing = TRUE)
    best_label <- names(scores)[score_order[1]]
    best_score <- scores[score_order[1]]
    second_score <- if (length(scores) >= 2) scores[score_order[2]] else NA_real_
    data.frame(
      cluster_id = cluster_id,
      best_module_label = best_label,
      best_module_score = best_score,
      second_module_score = second_score,
      score_margin = best_score - second_score,
      module_score_check = sprintf("best=%s; margin=%.3f; score=%.3f", best_label, best_score - second_score, best_score),
      stringsAsFactors = FALSE
    )
  })
  list(object = seu, summary = do.call(rbind, summary_rows))
}

annotation_warning_banner <- function(summary_df, min_reps = 2) {
  if (nrow(summary_df) == 0) {
    return("**⚠️ 缺少 replicate 信息，无法进行正式统计推断，以下结果仅为探索性分析**")
  }
  rep_counts <- unique(summary_df$biological_replicate_n)
  rep_counts <- rep_counts[is.finite(rep_counts)]
  if (length(rep_counts) == 1 && rep_counts[1] == 1) {
    return("**⚠️ 每组 N=1，无法进行正式统计推断，以下结果仅为探索性分析**")
  }
  sprintf("**⚠️ 每组 biological replicate < %s，无法进行正式统计推断，以下结果仅为探索性分析**", min_reps)
}

replicate_gate_summary <- function(meta_df, ident_1, ident_2, group_var = "group_id", replicate_var = "biological_replicate", sample_var = "sample_id", min_reps = 2) {
  required <- c(group_var, replicate_var, sample_var)
  missing_cols <- setdiff(required, colnames(meta_df))
  if (length(missing_cols) > 0) {
    summary_df <- data.frame(
      group_id = c(ident_1, ident_2),
      sample_n = NA_integer_,
      biological_replicate_n = NA_integer_,
      pass = FALSE,
      stringsAsFactors = FALSE
    )
    return(list(
      pass = FALSE,
      reason = sprintf("缺少关键 metadata 字段: %s", paste(missing_cols, collapse = ", ")),
      summary = summary_df
    ))
  }

  subset_df <- meta_df[meta_df[[group_var]] %in% c(ident_1, ident_2), , drop = FALSE]
  group_rows <- lapply(c(ident_1, ident_2), function(group_id) {
    grp <- subset_df[subset_df[[group_var]] == group_id, , drop = FALSE]
    reps <- unique(trim_character(grp[[replicate_var]], ""))
    reps <- reps[nzchar(reps)]
    samples <- unique(trim_character(grp[[sample_var]], ""))
    samples <- samples[nzchar(samples)]
    data.frame(
      group_id = group_id,
      sample_n = length(samples),
      biological_replicate_n = length(reps),
      pass = length(reps) >= min_reps,
      stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, group_rows)
  list(
    pass = all(summary_df$pass),
    reason = if (all(summary_df$pass)) {
      sprintf("每组 biological replicate >= %s，可进行正式推断", min_reps)
    } else {
      sprintf("至少一组 biological replicate < %s，formal test blocked", min_reps)
    },
    summary = summary_df
  )
}

read_comparison_sheet <- function(cfg) {
  df <- read_tsv_optional(cfg$comparison_sheet)
  if (nrow(df) == 0 && nzchar(cfg$deg_ident_1) && nzchar(cfg$deg_ident_2)) {
    df <- data.frame(
      comparison_id = sprintf("%s_vs_%s", cfg$deg_ident_1, cfg$deg_ident_2),
      ident_1 = cfg$deg_ident_1,
      ident_2 = cfg$deg_ident_2,
      enabled = "yes",
      group_var = "group_id",
      batch_var = "batch",
      layer_scope = "*",
      min_biological_replicates = cfg$min_biological_replicates,
      stringsAsFactors = FALSE
    )
  }
  if (nrow(df) == 0) {
    return(df)
  }
  defaults <- list(
    enabled = "yes",
    group_var = "group_id",
    batch_var = "batch",
    layer_scope = "*",
    min_biological_replicates = cfg$min_biological_replicates
  )
  for (col in c("comparison_id", "ident_1", "ident_2", names(defaults))) {
    if (!col %in% colnames(df)) {
      df[[col]] <- if (col %in% names(defaults)) defaults[[col]] else ""
    }
  }
  df$enabled <- normalize_flag(df$enabled, "yes")
  df$group_var <- trim_character(df$group_var, "group_id")
  df$batch_var <- trim_character(df$batch_var, "batch")
  df$layer_scope <- trim_character(df$layer_scope, "*")
  df$min_biological_replicates <- suppressWarnings(as.integer(df$min_biological_replicates))
  df$min_biological_replicates[is.na(df$min_biological_replicates)] <- cfg$min_biological_replicates
  df[df$enabled != "no", , drop = FALSE]
}

comparison_applies_to_layer <- function(comparison_row, layer_id) {
  scopes <- split_csv(comparison_row$layer_scope[[1]])
  if (length(scopes) == 0 || "*" %in% scopes) {
    return(TRUE)
  }
  layer_id %in% scopes
}

aggregate_cluster_sample_counts <- function(seu, cluster_var = "cluster_id", sample_var = "sample_id", group_var = "group_id", batch_var = "batch", replicate_var = "biological_replicate") {
  counts <- get_assay_matrix(seu, assay = "RNA", type = "counts")
  meta_df <- standardize_design_metadata_df(seu@meta.data)
  if (!cluster_var %in% colnames(meta_df)) {
    stop(sprintf("对象缺少 aggregation 所需 cluster 列: %s", cluster_var), call. = FALSE)
  }

  key_df <- unique(meta_df[, c(cluster_var, sample_var), drop = FALSE])
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
      batch = normalize_scalar_value(subset_df[[batch_var]][1], ""),
      biological_replicate = normalize_scalar_value(subset_df[[replicate_var]][1], ""),
      annotation_label = if ("annotation_label" %in% colnames(subset_df)) normalize_scalar_value(subset_df$annotation_label[1], "") else "",
      n_cells = nrow(subset_df),
      stringsAsFactors = FALSE
    )
  })

  list(counts = agg_mat, metadata = do.call(rbind, meta_rows))
}

sample_level_proportion_summary <- function(seu, cluster_var = "cluster_id", sample_var = "sample_id", group_var = "group_id") {
  meta_df <- standardize_design_metadata_df(seu@meta.data)
  if (!cluster_var %in% colnames(meta_df)) {
    stop(sprintf("对象缺少 composition 所需 cluster 列: %s", cluster_var), call. = FALSE)
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
