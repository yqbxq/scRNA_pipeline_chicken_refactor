annotation_confidence_levels <- function() {
  c("未定", "暂定", "确定")
}

cap_annotation_confidence <- function(value, ceiling = "") {
  levels <- annotation_confidence_levels()
  value <- normalize_scalar_value(value, "未定")
  ceiling <- normalize_scalar_value(ceiling)
  if (!value %in% levels) {
    value <- "未定"
  }
  if (!nzchar(ceiling) || !ceiling %in% levels) {
    return(value)
  }
  levels[min(match(value, levels), match(ceiling, levels))]
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
      df[[col]] <- vapply(df[[col]], normalize_scalar_value, character(1))
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

empty_marker_table_local <- function() {
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

run_cluster_marker_discovery <- function(seu, cluster_var, only_pos = TRUE, min_pct = 0.25, logfc_threshold = 0.25) {
  if (!cluster_var %in% colnames(seu@meta.data)) {
    stop(sprintf("对象缺少 cluster 列: %s", cluster_var), call. = FALSE)
  }
  Idents(seu) <- cluster_var
  if (length(unique(as.character(Idents(seu)))) < 2) {
    return(empty_marker_table_local())
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
    error = function(e) empty_marker_table_local()
  )
  if (!"gene" %in% colnames(markers)) {
    markers$gene <- rownames(markers)
  }
  markers$cluster <- as.character(markers$cluster)
  markers$gene <- as.character(markers$gene)
  markers
}

top_cluster_markers_local <- function(marker_df, cluster_id, top_n = 20) {
  if (nrow(marker_df) == 0) {
    return(marker_df)
  }
  df <- marker_df[as.character(marker_df$cluster) == as.character(cluster_id), , drop = FALSE]
  if (nrow(df) == 0) {
    return(df)
  }
  score_col <- if ("avg_log2FC" %in% colnames(df)) "avg_log2FC" else if ("avg_logFC" %in% colnames(df)) "avg_logFC" else ""
  if (nzchar(score_col)) {
    df <- df[order(-df[[score_col]], df$p_val_adj), , drop = FALSE]
  }
  df[seq_len(min(top_n, nrow(df))), , drop = FALSE]
}

format_data_evidence_local <- function(marker_df, top_n = 5) {
  if (nrow(marker_df) == 0) {
    return("no_positive_markers")
  }
  score_col <- if ("avg_log2FC" %in% colnames(marker_df)) "avg_log2FC" else if ("avg_logFC" %in% colnames(marker_df)) "avg_logFC" else ""
  if (nzchar(score_col)) {
    marker_df <- marker_df[order(-marker_df[[score_col]], marker_df$p_val_adj), , drop = FALSE]
  }
  marker_df <- marker_df[seq_len(min(top_n, nrow(marker_df))), , drop = FALSE]
  if (nzchar(score_col)) {
    return(paste(sprintf("%s(%.2f)", marker_df$gene, marker_df[[score_col]]), collapse = ", "))
  }
  paste(marker_df$gene, collapse = ", ")
}

panel_overlap_metrics <- function(cluster_marker_df, panel_df, top_n = 20) {
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

panel_gene_sets_local <- function(panel_df) {
  if (nrow(panel_df) == 0) {
    return(list())
  }
  split(unique(panel_df$gene), panel_df$celltype)
}

compute_module_score_summary <- function(seu, cluster_var, panel_df, seed = 42) {
  if (nrow(panel_df) == 0) {
    return(list(object = seu, summary = data.frame(stringsAsFactors = FALSE)))
  }
  gene_sets <- panel_gene_sets_local(panel_df)
  gene_sets <- lapply(gene_sets, function(genes) intersect(unique(genes), rownames(seu)))
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) > 0]
  if (length(gene_sets) == 0) {
    return(list(object = seu, summary = data.frame(stringsAsFactors = FALSE)))
  }

  score_cols <- list()
  scored <- tryCatch(
    AddModuleScore(seu, features = gene_sets, name = "module_check_", seed = seed),
    error = function(e) NULL
  )
  if (!is.null(scored)) {
    new_cols <- setdiff(colnames(scored@meta.data), colnames(seu@meta.data))
    seu <- scored
    for (idx in seq_along(gene_sets)) {
      if (idx <= length(new_cols)) {
        score_cols[[names(gene_sets)[idx]]] <- new_cols[[idx]]
      }
    }
  }
  if (length(score_cols) == 0) {
    expr_mat <- GetAssayData(seu, assay = DefaultAssay(seu), slot = "data")
    for (label in names(gene_sets)) {
      safe_label <- gsub("[^A-Za-z0-9]+", "_", label)
      score_col <- paste0("module_check_", safe_label)
      seu[[score_col]] <- Matrix::colMeans(expr_mat[gene_sets[[label]], , drop = FALSE])
      score_cols[[label]] <- score_col
    }
  }

  meta_df <- seu@meta.data
  clusters <- sort(unique(as.character(meta_df[[cluster_var]])))
  summary_rows <- lapply(clusters, function(cluster_id) {
    hit <- as.character(meta_df[[cluster_var]]) == cluster_id
    scores <- vapply(score_cols, function(col) stats::median(meta_df[[col]][hit], na.rm = TRUE), numeric(1))
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

relation_from_overlap_local <- function(best_overlap, second_overlap, panel_present) {
  if (!panel_present || best_overlap <= 0) {
    return("无关")
  }
  if (second_overlap >= best_overlap && best_overlap > 0) {
    return("冲突")
  }
  "支持"
}

confidence_from_overlap_local <- function(best_overlap, second_overlap, evidence_source, confidence_ceiling = "", panel_present = TRUE) {
  if (!panel_present || best_overlap <= 0 || second_overlap >= best_overlap) {
    return("未定")
  }
  source_tokens <- tolower(split_csv_local(evidence_source))
  indirect_only <- length(source_tokens) > 0 && all(grepl("ortholog|candidate|indirect|putative", source_tokens))
  confidence <- if (best_overlap >= 2 && !indirect_only && second_overlap == 0) "确定" else "暂定"
  cap_annotation_confidence(confidence, confidence_ceiling)
}

build_annotation_decision <- function(cluster_id, cluster_markers, overlap_df, module_df, panel_present) {
  best_row <- if (nrow(overlap_df) > 0) overlap_df[1, , drop = FALSE] else NULL
  best_overlap <- if (!is.null(best_row)) best_row$overlap_n[[1]] else 0L
  second_overlap <- if (nrow(overlap_df) >= 2) overlap_df$overlap_n[[2]] else 0L
  literature_candidate <- if (!is.null(best_row) && best_overlap > 0) best_row$celltype[[1]] else ""
  confidence <- confidence_from_overlap_local(
    best_overlap,
    second_overlap,
    if (!is.null(best_row)) best_row$evidence_source[[1]] else "",
    if (!is.null(best_row)) best_row$confidence_ceiling[[1]] else "",
    panel_present = panel_present
  )
  relation <- relation_from_overlap_local(best_overlap, second_overlap, panel_present)
  final_annotation <- if (confidence == "未定" || !nzchar(literature_candidate)) {
    sprintf("Uncertain-%s", cluster_id)
  } else {
    literature_candidate
  }
  data_candidate <- if (nrow(cluster_markers) > 0) paste(head(cluster_markers$gene, 3), collapse = "+") else sprintf("Uncertain-%s", cluster_id)
  module_row <- if (nrow(module_df) > 0) module_df[module_df$cluster_id == cluster_id, , drop = FALSE] else data.frame()
  module_check <- if (nrow(module_row) > 0) {
    agreement <- if (nzchar(literature_candidate) && module_row$best_module_label[[1]] == literature_candidate) {
      "一致"
    } else if (nzchar(literature_candidate)) {
      "不一致"
    } else {
      "候选缺失"
    }
    sprintf("%s; agreement=%s", module_row$module_score_check[[1]], agreement)
  } else if (panel_present) {
    "panel_present_but_module_unavailable"
  } else {
    "panel_absent"
  }
  list(
    candidate_celltype = final_annotation,
    data_candidate = data_candidate,
    literature_candidate = literature_candidate,
    confidence = confidence,
    relation = relation,
    data_evidence = format_data_evidence_local(cluster_markers, top_n = 5),
    literature_evidence = if (!is.null(best_row) && best_overlap > 0) {
      sprintf("%s: %s [%s]", best_row$celltype[[1]], normalize_scalar_value(best_row$overlap_genes[[1]], "none"), normalize_scalar_value(best_row$evidence_source[[1]], "unspecified"))
    } else if (panel_present) {
      "panel_present_no_overlap"
    } else {
      "panel_absent"
    },
    module_score_check = module_check
  )
}

plot_annotation_top_markers <- function(seu, marker_df, cluster_var, output_png) {
  if (nrow(marker_df) == 0) {
    return(invisible(NULL))
  }
  genes <- marker_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_head(n = 3) %>%
    dplyr::pull(gene) %>%
    unique()
  genes <- intersect(genes, rownames(seu))
  if (length(genes) == 0) {
    return(invisible(NULL))
  }
  plot_obj <- DotPlot(seu, features = genes, group.by = cluster_var) +
    RotatedAxis() +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::labs(title = "Top cluster markers")
  save_plot_local(plot_obj, output_png, width = 10, height = 6)
}

plot_annotation_panel_validation <- function(seu, panel_df, cluster_var, output_png) {
  if (nrow(panel_df) == 0 || !"gene" %in% colnames(panel_df)) {
    return(invisible(NULL))
  }
  genes <- intersect(unique(panel_df$gene), rownames(seu))
  if (length(genes) == 0) {
    return(invisible(NULL))
  }
  plot_obj <- DotPlot(seu, features = genes, group.by = cluster_var) +
    RotatedAxis() +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::labs(title = "Panel validation markers")
  save_plot_local(plot_obj, output_png, width = 10, height = 6)
}

plot_module_score_heatmap_local <- function(module_df, output_png) {
  if (nrow(module_df) == 0 || !"best_module_label" %in% colnames(module_df)) {
    return(invisible(NULL))
  }
  heatmap_df <- module_df %>%
    dplyr::select(cluster_id, best_module_label, best_module_score, score_margin) %>%
    tidyr::pivot_longer(cols = c(best_module_score, score_margin), names_to = "metric", values_to = "value")
  plot_obj <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = metric, y = factor(cluster_id), fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#2F6CB3") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::labs(title = "Module score check", x = NULL, y = "Cluster", fill = "Value")
  save_plot_local(plot_obj, output_png, width = 6, height = 4.5)
}

annotate_one_layer <- function(seu, layer_id, marker_panels_df, cfg, paths_module) {
  cluster_var <- paths_module$cluster_var
  if (!cluster_var %in% colnames(seu@meta.data)) {
    stop(sprintf("对象缺少 cluster 列: %s", cluster_var), call. = FALSE)
  }
  ensure_dir(paths_module$table_dir)
  ensure_dir(paths_module$figure_dir)

  marker_df <- run_cluster_marker_discovery(seu, cluster_var = cluster_var, only_pos = TRUE, min_pct = 0.25, logfc_threshold = 0.25)
  write_tsv_local(marker_df, paths_module$cluster_markers_tsv)

  module_payload <- compute_module_score_summary(seu, cluster_var = cluster_var, panel_df = marker_panels_df, seed = cfg$random_seed)
  seu <- module_payload$object
  module_df <- module_payload$summary
  write_tsv_local(module_df, paths_module$module_score_summary_tsv)

  cluster_ids <- sort(unique(as.character(seu@meta.data[[cluster_var]])))
  cluster_sizes <- table(as.character(seu@meta.data[[cluster_var]]))
  annotation_rows <- list()
  evidence_rows <- list()
  for (cluster_id in cluster_ids) {
    cluster_markers <- top_cluster_markers_local(marker_df, cluster_id, top_n = 20)
    overlap_df <- panel_overlap_metrics(cluster_markers, marker_panels_df, top_n = 20)
    decision <- build_annotation_decision(cluster_id, cluster_markers, overlap_df, module_df, panel_present = nrow(marker_panels_df) > 0)
    annotation_rows[[length(annotation_rows) + 1]] <- data.frame(
      layer_id = layer_id,
      cluster = cluster_id,
      cluster_id = cluster_id,
      cluster_column = cluster_var,
      candidate_celltype = decision$candidate_celltype,
      final_annotation = decision$candidate_celltype,
      data_driven_candidate = decision$data_candidate,
      literature_candidate = decision$literature_candidate,
      confidence = decision$confidence,
      relation = decision$relation,
      evidence_relation = decision$relation,
      data_evidence = decision$data_evidence,
      literature_evidence = decision$literature_evidence,
      module_score_check = decision$module_score_check,
      n_cells = unname(cluster_sizes[[cluster_id]]),
      marker_panel_available = ifelse(nrow(marker_panels_df) > 0, "yes", "no"),
      stringsAsFactors = FALSE
    )
    if (nrow(overlap_df) == 0) {
      evidence_rows[[length(evidence_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        cluster_id = cluster_id,
        celltype = "",
        overlap_n = 0,
        overlap_genes = "",
        evidence_source = ifelse(nrow(marker_panels_df) > 0, "panel_present_no_overlap", "panel_absent"),
        confidence_ceiling = "",
        panel_files = ifelse(nrow(marker_panels_df) > 0, paste(sort(unique(marker_panels_df$panel_file)), collapse = ","), ""),
        stringsAsFactors = FALSE
      )
    } else {
      overlap_df$layer_id <- layer_id
      overlap_df$cluster_id <- cluster_id
      evidence_rows[[length(evidence_rows) + 1]] <- overlap_df[, c("layer_id", "cluster_id", "celltype", "overlap_n", "overlap_genes", "evidence_source", "confidence_ceiling", "panel_files")]
    }
  }

  annotation_table <- dplyr::bind_rows(annotation_rows)
  evidence_table <- dplyr::bind_rows(evidence_rows)
  write_tsv_local(annotation_table, paths_module$annotation_table_tsv)
  write_tsv_local(evidence_table, paths_module$annotation_evidence_tsv)

  cluster_key <- as.character(seu@meta.data[[cluster_var]])
  seu$annotation_label <- annotation_table$final_annotation[match(cluster_key, annotation_table$cluster_id)]
  seu$annotation_confidence <- annotation_table$confidence[match(cluster_key, annotation_table$cluster_id)]
  seu$annotation_evidence_relation <- annotation_table$relation[match(cluster_key, annotation_table$cluster_id)]
  seu$cell_type <- seu$annotation_label
  seu$cell_type_confidence <- seu$annotation_confidence
  seu$annotation_relation <- seu$annotation_evidence_relation

  tryCatch(plot_annotation_top_markers(seu, marker_df, cluster_var, paths_module$top_markers_dotplot_png), error = function(e) warning(conditionMessage(e), call. = FALSE))
  tryCatch(plot_annotation_panel_validation(seu, marker_panels_df, cluster_var, paths_module$panel_validation_dotplot_png), error = function(e) warning(conditionMessage(e), call. = FALSE))
  tryCatch(plot_module_score_heatmap_local(module_df, paths_module$module_score_heatmap_png), error = function(e) warning(conditionMessage(e), call. = FALSE))

  list(
    seu = seu,
    annotation_table = annotation_table,
    evidence_table = evidence_table,
    module_score_summary = module_df,
    marker_table = marker_df
  )
}
