source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "annotation_stats_helpers.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "plotting_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(tidyr)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

layer_status_df <- built_layer_status(cfg)
if (nrow(layer_status_df) == 0) {
  stop("layer_status.tsv 中没有 status=built 的层，无法执行 annotation。", call. = FALSE)
}

annotation_root_dir <- file.path(cfg$table_dir, "annotation")
dir.create(annotation_root_dir, recursive = TRUE, showWarnings = FALSE)
annotation_manifest_path <- file.path(annotation_root_dir, "layer_annotation_manifest.tsv")

choose_cluster_var <- function(obj, layer_row) {
  cluster_var <- normalize_scalar_value(layer_row$cluster_column[[1]], "")
  if (nzchar(cluster_var) && cluster_var %in% colnames(obj@meta.data)) {
    return(cluster_var)
  }
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    return("seurat_clusters")
  }
  stop(sprintf("layer `%s` 缺少 cluster 列。", layer_row$layer_id[[1]]), call. = FALSE)
}

relation_from_overlap <- function(best_overlap, second_overlap, panel_present) {
  if (!panel_present) {
    return("仅数据驱动")
  }
  if (best_overlap <= 0) {
    return("无文献验证")
  }
  if (second_overlap >= best_overlap && best_overlap > 0) {
    return("冲突")
  }
  if (best_overlap >= 2 && second_overlap == 0) {
    return("一致")
  }
  "部分一致"
}

confidence_from_overlap <- function(best_overlap, second_overlap, evidence_source, confidence_ceiling = "", panel_present = TRUE) {
  if (!panel_present || best_overlap <= 0) {
    return("未定")
  }
  if (second_overlap >= best_overlap && best_overlap > 0) {
    return("未定")
  }

  source_tokens <- trimws(unlist(strsplit(normalize_scalar_value(evidence_source, ""), ",", fixed = TRUE)))
  source_tokens <- tolower(source_tokens[nzchar(source_tokens)])
  indirect_only <- length(source_tokens) > 0 && all(grepl("ortholog|cross[_ -]?species|candidate|indirect|putative", source_tokens))

  confidence <- if (best_overlap >= 2 && !indirect_only && second_overlap == 0) {
    "确定"
  } else if (best_overlap >= 1) {
    "暂定"
  } else {
    "未定"
  }
  cap_annotation_confidence(confidence, confidence_ceiling)
}

format_literature_evidence <- function(best_row, panel_present) {
  if (!panel_present) {
    return("panel_absent")
  }
  if (is.null(best_row) || nrow(best_row) == 0 || best_row$overlap_n[[1]] <= 0) {
    return("panel_present_no_overlap")
  }
  sprintf(
    "%s: %s [%s]",
    best_row$celltype[[1]],
    normalize_scalar_value(best_row$overlap_genes[[1]], "none"),
    normalize_scalar_value(best_row$evidence_source[[1]], "unspecified")
  )
}

make_uncertain_label <- function(cluster_id) {
  sprintf("Uncertain-%s", cluster_id)
}

plot_layer_top_markers <- function(seu, marker_df, cluster_var, output_png) {
  genes <- marker_df %>%
    group_by(cluster) %>%
    slice_head(n = 3) %>%
    pull(gene) %>%
    unique()
  genes <- intersect(genes, rownames(seu))
  if (length(genes) == 0) {
    return(invisible(NULL))
  }
  plot_obj <- DotPlot(seu, features = genes, group.by = cluster_var) +
    RotatedAxis() +
    theme_classic(base_size = 10) +
    labs(title = "Top cluster markers")
  save_plot_dual(plot_obj, output_png, width = 10, height = 6)
}

plot_layer_panel_validation <- function(seu, panel_df, cluster_var, output_png) {
  if (nrow(panel_df) == 0 || !"gene" %in% colnames(panel_df)) {
    return(invisible(NULL))
  }
  genes <- intersect(unique(panel_df$gene), rownames(seu))
  if (length(genes) == 0) {
    return(invisible(NULL))
  }
  plot_obj <- DotPlot(seu, features = genes, group.by = cluster_var) +
    RotatedAxis() +
    theme_classic(base_size = 10) +
    labs(title = "Panel validation markers")
  save_plot_dual(plot_obj, output_png, width = 10, height = 6)
}

plot_module_score_heatmap <- function(module_df, output_png) {
  if (nrow(module_df) == 0 || !"best_module_label" %in% colnames(module_df)) {
    return(invisible(NULL))
  }
  heatmap_df <- module_df %>%
    select(cluster_id, best_module_label, best_module_score, score_margin) %>%
    pivot_longer(cols = c(best_module_score, score_margin), names_to = "metric", values_to = "value")
  plot_obj <- ggplot(heatmap_df, aes(x = metric, y = factor(cluster_id), fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = paper_feature_palette()) +
    theme_classic(base_size = 10) +
    labs(title = "Module score check", x = NULL, y = "Cluster", fill = "Value")
  save_plot_dual(plot_obj, output_png, width = 6, height = 4.5)
}

write_layer_report <- function(layer_row, ann_paths, marker_paths, annotation_df, evidence_df, panel_df) {
  relation_counts <- if (nrow(annotation_df) > 0) table(annotation_df$evidence_relation) else integer(0)
  confidence_counts <- if (nrow(annotation_df) > 0) table(annotation_df$confidence) else integer(0)
  report_lines <- c(
    sprintf("# Annotation Report: %s", layer_row$layer_id[[1]]),
    "",
    sprintf("- layer_role: `%s`", normalize_scalar_value(layer_row$layer_role[[1]], "")),
    sprintf("- clustered_rds: `%s`", layer_row$clustered_rds[[1]]),
    sprintf("- cluster markers: `%s`", marker_paths$cluster_markers_tsv),
    sprintf("- annotation table: `%s`", ann_paths$annotation_table_tsv),
    sprintf("- annotation evidence: `%s`", ann_paths$evidence_table_tsv),
    sprintf("- module score summary: `%s`", ann_paths$module_score_tsv),
    sprintf("- marker panel rows used: `%s`", nrow(panel_df)),
    sprintf("- marker panel files used: `%s`", if (nrow(panel_df) > 0) paste(sort(unique(panel_df$panel_file)), collapse = ", ") else "none"),
    "",
    "## Confidence Summary"
  )
  if (length(confidence_counts) == 0) {
    report_lines <- c(report_lines, "- no annotation rows")
  } else {
    for (name in names(confidence_counts)) {
      report_lines <- c(report_lines, sprintf("- `%s`: `%s`", name, confidence_counts[[name]]))
    }
  }
  report_lines <- c(report_lines, "", "## Evidence Relation Summary")
  if (length(relation_counts) == 0) {
    report_lines <- c(report_lines, "- no evidence rows")
  } else {
    for (name in names(relation_counts)) {
      report_lines <- c(report_lines, sprintf("- `%s`: `%s`", name, relation_counts[[name]]))
    }
  }
  report_lines <- c(report_lines, "", "## Cluster Review")
  for (i in seq_len(nrow(annotation_df))) {
    row <- annotation_df[i, , drop = FALSE]
    report_lines <- c(
      report_lines,
      sprintf(
        "- cluster `%s`: final=`%s`; confidence=`%s`; relation=`%s`",
        row$cluster_id,
        row$final_annotation,
        row$confidence,
        row$evidence_relation
      ),
      sprintf("  data_evidence: %s", row$data_evidence),
      sprintf("  literature_evidence: %s", row$literature_evidence),
      sprintf("  module_score_check: %s", row$module_score_check)
    )
  }
  write_markdown(report_lines, ann_paths$report_md)
}

layer_manifest_rows <- list()
panorama_object_saved <- FALSE

for (idx in seq_len(nrow(layer_status_df))) {
  layer_row <- layer_status_df[idx, , drop = FALSE]
  layer_id <- layer_row$layer_id[[1]]
  message("开始 annotation layer: ", layer_id)

  obj <- readRDS(layer_row$clustered_rds[[1]])
  obj <- maybe_join_layers(obj)
  obj <- standardize_design_metadata(obj)

  cluster_var <- choose_cluster_var(obj, layer_row)
  obj$cluster_id <- as.character(obj@meta.data[[cluster_var]])

  tissue_values <- unique(trim_character(obj$tissue, "*"))
  tissue_values <- tissue_values[nzchar(tissue_values) & tissue_values != "*"]
  panel_df <- read_marker_panel_rows(cfg, layer_id = layer_id, tissue = if (length(tissue_values) == 1) tissue_values[[1]] else NULL)

  marker_paths <- marker_discovery_paths(cfg, layer_id)
  ann_paths <- layer_annotation_paths(cfg, layer_id)
  ensure_parent_dirs(c(marker_paths, ann_paths))

  marker_df <- run_cluster_marker_discovery(obj, cluster_var = cluster_var, only_pos = TRUE, min_pct = 0.25, logfc_threshold = 0)
  write_tsv(marker_df, marker_paths$cluster_markers_tsv)

  module_payload <- compute_module_score_summary(obj, cluster_var = cluster_var, panel_df = panel_df, seed = cfg$random_seed)
  obj <- module_payload$object
  module_df <- module_payload$summary
  write_tsv(module_df, ann_paths$module_score_tsv)

  cluster_ids <- sort(unique(as.character(obj@meta.data[[cluster_var]])))
  cluster_sizes <- table(as.character(obj@meta.data[[cluster_var]]))
  annotation_rows <- vector("list", length(cluster_ids))
  evidence_rows <- list()

  for (j in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[[j]]
    cluster_markers <- top_cluster_markers(marker_df, cluster_id, top_n = 20)
    overlap_df <- panel_overlap_summary(cluster_markers, panel_df, top_n = 20)
    best_row <- if (nrow(overlap_df) > 0) overlap_df[1, , drop = FALSE] else NULL
    best_overlap <- if (!is.null(best_row)) best_row$overlap_n[[1]] else 0L
    second_overlap <- if (nrow(overlap_df) >= 2) overlap_df$overlap_n[[2]] else 0L
    relation <- relation_from_overlap(best_overlap, second_overlap, panel_present = nrow(panel_df) > 0)
    confidence <- confidence_from_overlap(
      best_overlap = best_overlap,
      second_overlap = second_overlap,
      evidence_source = if (!is.null(best_row)) best_row$evidence_source[[1]] else "",
      confidence_ceiling = if (!is.null(best_row)) best_row$confidence_ceiling[[1]] else "",
      panel_present = nrow(panel_df) > 0
    )
    literature_candidate <- if (!is.null(best_row) && best_overlap > 0) best_row$celltype[[1]] else ""
    final_annotation <- if (confidence == "未定" || !nzchar(literature_candidate)) {
      make_uncertain_label(cluster_id)
    } else {
      literature_candidate
    }
    data_candidate <- if (nrow(cluster_markers) > 0) paste(head(cluster_markers$gene, 3), collapse = "+") else make_uncertain_label(cluster_id)
    module_row <- if (nrow(module_df) > 0) module_df[module_df$cluster_id == cluster_id, , drop = FALSE] else data.frame()
    module_check <- if (nrow(module_row) > 0) {
      agreement <- if (nzchar(literature_candidate) && module_row$best_module_label[[1]] == literature_candidate) "一致" else if (nzchar(literature_candidate)) "不一致" else "候选缺失"
      sprintf("%s; agreement=%s", module_row$module_score_check[[1]], agreement)
    } else if (nrow(panel_df) > 0) {
      "panel_present_but_module_unavailable"
    } else {
      "panel_absent"
    }

    annotation_rows[[j]] <- data.frame(
      layer_id = layer_id,
      cluster_id = cluster_id,
      cluster_column = cluster_var,
      final_annotation = final_annotation,
      data_driven_candidate = data_candidate,
      data_evidence = format_data_evidence(cluster_markers, top_n = 5),
      literature_candidate = literature_candidate,
      literature_evidence = format_literature_evidence(best_row, panel_present = nrow(panel_df) > 0),
      module_score_check = module_check,
      confidence = confidence,
      evidence_relation = relation,
      n_cells = unname(cluster_sizes[[cluster_id]]),
      marker_panel_available = ifelse(nrow(panel_df) > 0, "yes", "no"),
      stringsAsFactors = FALSE
    )

    if (nrow(overlap_df) == 0) {
      evidence_rows[[length(evidence_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        cluster_id = cluster_id,
        celltype = "",
        overlap_n = 0,
        overlap_genes = "",
        evidence_source = if (nrow(panel_df) > 0) "panel_present_no_overlap" else "panel_absent",
        confidence_ceiling = "",
        panel_files = if (nrow(panel_df) > 0) paste(sort(unique(panel_df$panel_file)), collapse = ",") else "",
        stringsAsFactors = FALSE
      )
    } else {
      overlap_df$layer_id <- layer_id
      overlap_df$cluster_id <- cluster_id
      evidence_rows[[length(evidence_rows) + 1]] <- overlap_df[, c("layer_id", "cluster_id", "celltype", "overlap_n", "overlap_genes", "evidence_source", "confidence_ceiling", "panel_files")]
    }
  }

  annotation_df <- bind_rows(annotation_rows)
  evidence_df <- bind_rows(evidence_rows)
  write_tsv(annotation_df, ann_paths$annotation_table_tsv)
  write_tsv(evidence_df, ann_paths$evidence_table_tsv)

  obj$annotation_label <- annotation_df$final_annotation[match(obj$cluster_id, annotation_df$cluster_id)]
  obj$annotation_confidence <- annotation_df$confidence[match(obj$cluster_id, annotation_df$cluster_id)]
  obj$annotation_evidence_relation <- annotation_df$evidence_relation[match(obj$cluster_id, annotation_df$cluster_id)]
  obj$cell_type <- factor(obj$annotation_label)
  saveRDS(obj, ann_paths$annotated_rds)

  if (!panorama_object_saved && normalize_scalar_value(layer_row$layer_role[[1]], "") == "panorama") {
    saveRDS(obj, file.path(cfg$checkpoint_dir, "03_after_annotation.rds"))
    panorama_object_saved <- TRUE
  }

  tryCatch({
    plot_layer_top_markers(obj, marker_df, cluster_var = cluster_var, output_png = ann_paths$top_markers_dotplot_png)
  }, error = function(e) {
    warning(sprintf("layer `%s` top marker DotPlot 失败: %s", layer_id, conditionMessage(e)), call. = FALSE)
  })
  tryCatch({
    plot_layer_panel_validation(obj, panel_df, cluster_var = cluster_var, output_png = ann_paths$panel_dotplot_png)
  }, error = function(e) {
    warning(sprintf("layer `%s` panel DotPlot 失败: %s", layer_id, conditionMessage(e)), call. = FALSE)
  })
  tryCatch({
    plot_module_score_heatmap(module_df, ann_paths$module_score_heatmap_png)
  }, error = function(e) {
    warning(sprintf("layer `%s` module score heatmap 失败: %s", layer_id, conditionMessage(e)), call. = FALSE)
  })

  write_layer_report(layer_row, ann_paths, marker_paths, annotation_df, evidence_df, panel_df)

  layer_manifest_rows[[length(layer_manifest_rows) + 1]] <- data.frame(
    layer_id = layer_id,
    layer_role = normalize_scalar_value(layer_row$layer_role[[1]], ""),
    cluster_count = length(cluster_ids),
    determined_n = sum(annotation_df$confidence == "确定"),
    tentative_n = sum(annotation_df$confidence == "暂定"),
    undetermined_n = sum(annotation_df$confidence == "未定"),
    marker_panel_rows = nrow(panel_df),
    cluster_markers_tsv = marker_paths$cluster_markers_tsv,
    annotation_table_tsv = ann_paths$annotation_table_tsv,
    report_md = ann_paths$report_md,
    annotated_rds = ann_paths$annotated_rds,
    stringsAsFactors = FALSE
  )
}

layer_manifest_df <- bind_rows(layer_manifest_rows)
write_tsv(layer_manifest_df, annotation_manifest_path)

if (!panorama_object_saved) {
  first_layer_rds <- layer_manifest_df$annotated_rds[[1]]
  if (nzchar(first_layer_rds) && file.exists(first_layer_rds)) {
    file.copy(first_layer_rds, file.path(cfg$checkpoint_dir, "03_after_annotation.rds"), overwrite = TRUE)
  }
}

message("annotation 完成，manifest: ", annotation_manifest_path)
