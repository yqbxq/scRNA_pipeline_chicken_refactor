#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source_utf8 <- function(path) source(path, encoding = "UTF-8")

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_05.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_06.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_07.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "deg_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_mapping_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_08.R"))
source_utf8(file.path(.script_dir, "helpers", "scenic_utils.R"))

load_required_packages(c("Seurat", "SCENIC", "AUCell", "dplyr", "tidyr", "tibble", "ggplot2", "pheatmap", "jsonlite"))

cfg <- get_single_script_config_08()
module_name <- "08d_scenic_downstream"
prepare_dirs_08(cfg)
set.seed(cfg$random_seed)

if (!file.exists(cfg$module_08c_manifest_path) || !file.exists(cfg$scenic_regulons_index_tsv)) {
  stop(sprintf("缺少 08c SCENIC regulons 输出: %s", cfg$scenic_regulons_index_tsv), call. = FALSE)
}

module_colors <- c(
  "1" = "#8BC34A",
  "2" = "#F48FB1",
  "3" = "#26C6DA",
  "4" = "#CE93D8"
)

tasks <- read_tsv_optional(cfg$scenic_regulons_index_tsv)
if (nrow(tasks) == 0) {
  stop("08d 没有可处理的 SCENIC regulon layer", call. = FALSE)
}

index_rows <- list()
dynamic_outputs <- list()

for (idx in seq_len(nrow(tasks))) {
  task <- tasks[idx, , drop = FALSE]
  layer_id <- task$layer_id[[1]]
  safe_layer <- safe_id_08(layer_id)
  message("08d SCENIC downstream layer: ", layer_id)

  input_rds <- task$object_rds[[1]]
  auc_csv <- task$auc_matrix_csv[[1]]
  if (!file.exists(input_rds)) {
    stop(sprintf("Missing layer object: %s", input_rds), call. = FALSE)
  }
  if (!file.exists(auc_csv)) {
    stop(sprintf("Missing SCENIC AUC matrix: %s", auc_csv), call. = FALSE)
  }

  paths <- scenic_downstream_paths_08(cfg, layer_id)
  ensure_dir(paths$table_dir)
  ensure_dir(paths$figure_dir)
  ensure_dir(paths$checkpoint_dir)

  obj <- readRDS(input_rds)
  obj <- maybe_join_layers(obj)

  auc_raw <- read.csv(auc_csv, row.names = 1, check.names = FALSE)
  if (nrow(auc_raw) == ncol(obj)) {
    auc_cells_x_regulons <- auc_raw
  } else if (ncol(auc_raw) == ncol(obj)) {
    auc_cells_x_regulons <- t(auc_raw)
  } else {
    stop(sprintf("auc_matrix.csv dimensions do not match layer %s object", layer_id), call. = FALSE)
  }

  common_cells <- intersect(rownames(auc_cells_x_regulons), colnames(obj))
  if (length(common_cells) == 0) {
    stop(sprintf("No overlapping cells between layer %s AUC matrix and Seurat object", layer_id), call. = FALSE)
  }

  auc_cells_x_regulons <- auc_cells_x_regulons[common_cells, , drop = FALSE]
  obj <- subset(obj, cells = common_cells)

  auc_features_x_cells <- t(as.matrix(auc_cells_x_regulons))
  rownames(auc_features_x_cells) <- gsub("\\(_\\+\\)|\\(\\+\\)", "", rownames(auc_features_x_cells))
  obj[["SCENIC"]] <- CreateAssayObject(data = auc_features_x_cells)

  cell_type_col <- resolve_cell_type_col_07(obj, layer_id, task$layer_role[[1]])
  if (!nzchar(cell_type_col)) {
    stop(sprintf("layer %s 缺少 cell type/cluster metadata，无法计算 RSS", layer_id), call. = FALSE)
  }
  cell_types <- as.character(obj@meta.data[[cell_type_col]])
  names(cell_types) <- colnames(obj)
  celltype_colors <- named_discrete_palette_08(unique(cell_types))
  rss_mat <- calcRSS(AUC = auc_features_x_cells, cellAnnotation = cell_types)

  rss_df <- as.data.frame(rss_mat)
  rss_df$regulon <- rownames(rss_df)
  write.csv(rss_df, paths$rss_matrix_csv, row.names = FALSE)

  top_regulons <- dplyr::bind_rows(lapply(colnames(rss_mat), function(ct) {
    tibble::tibble(
      cell_type = ct,
      regulon = rownames(rss_mat),
      rss = rss_mat[, ct]
    ) %>%
      dplyr::arrange(dplyr::desc(rss)) %>%
      dplyr::slice_head(n = 5)
  })) %>%
    dplyr::distinct(regulon, .keep_all = TRUE)
  write.csv(top_regulons, paths$top_regulons_by_celltype_csv, row.names = FALSE)

  selected_cells <- unlist(lapply(names(celltype_colors), function(ct) {
    ct_cells <- names(cell_types)[cell_types == ct]
    if (length(ct_cells) > 250) {
      sample(ct_cells, 250)
    } else {
      ct_cells
    }
  }), use.names = FALSE)
  selected_cells <- selected_cells[selected_cells %in% colnames(obj)]
  selected_cells <- selected_cells[order(factor(cell_types[selected_cells], levels = names(celltype_colors)))]

  heatmap_regs <- top_regulons$regulon[top_regulons$regulon %in% rownames(auc_features_x_cells)]
  heatmap_mat <- auc_features_x_cells[heatmap_regs, selected_cells, drop = FALSE]
  heatmap_scaled <- t(scale(t(heatmap_mat)))
  heatmap_scaled[is.na(heatmap_scaled)] <- 0
  heatmap_scaled[heatmap_scaled > 2] <- 2
  heatmap_scaled[heatmap_scaled < -2] <- -2

  anno_col <- data.frame(CellType = factor(cell_types[selected_cells], levels = names(celltype_colors)))
  rownames(anno_col) <- selected_cells

  png(paths$figure_6a_png, width = 1600, height = 950, res = 150)
  pheatmap(
    heatmap_scaled,
    color = colorRampPalette(c("#313695", "white", "#A50026"))(100),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    annotation_col = anno_col,
    annotation_colors = list(CellType = celltype_colors),
    main = "Figure 6A: regulon activity heatmap"
  )
  dev.off()

  rss_long <- dplyr::bind_rows(lapply(colnames(rss_mat), function(ct) {
    tibble::tibble(
      cell_type = ct,
      regulon = rownames(rss_mat),
      rss = rss_mat[, ct]
    ) %>%
      dplyr::arrange(dplyr::desc(rss)) %>%
      dplyr::mutate(rank = dplyr::row_number())
  }))

  rss_long <- rss_long %>%
    dplyr::group_by(cell_type) %>%
    dplyr::mutate(label = ifelse(rank <= 3, regulon, NA_character_)) %>%
    dplyr::ungroup()

  fig6b <- ggplot(rss_long, aes(x = rank, y = rss)) +
    geom_point(color = "grey25", size = 1.5) +
    geom_point(data = subset(rss_long, !is.na(label)), color = "red3", size = 2.2) +
    geom_text(
      data = subset(rss_long, !is.na(label)),
      aes(label = label),
      hjust = -0.05,
      vjust = -0.2,
      size = 3
    ) +
    facet_wrap(~cell_type, ncol = 2, scales = "free_y") +
    labs(x = "Rank", y = "RSS") +
    theme_classic(base_size = 12)

  ggsave(paths$figure_6b_png, fig6b, width = 10, height = 8, dpi = 300, bg = "white")
  ggsave(paths$figure_6b_pdf, fig6b, width = 10, height = 8, bg = "white")

  csi_result <- calculate_csi_long_08(
    auc_features_x_cells,
    file.path(.script_dir, "vendor", "calculate_csi.R")
  )
  csi_long <- csi_result$csi_long
  if (identical(csi_result$source, "vendor")) {
    writeLines(
      "scFunctions package was unavailable, so the original calculate_csi source was used directly.",
      paths$csi_source_note_txt
    )
  }

  csi_matrix <- NULL
  if (!is.null(csi_long) && nrow(csi_long) > 0) {
    csi_matrix <- csi_long %>%
      tidyr::pivot_wider(names_from = regulon_2, values_from = CSI) %>%
      tibble::column_to_rownames("regulon_1") %>%
      as.matrix()
    write.csv(csi_matrix, paths$csi_matrix_csv)
  }

  if (!is.null(csi_matrix) && nrow(csi_matrix) > 1) {
    overlap_regs <- intersect(heatmap_regs, rownames(csi_matrix))
    csi_sub <- csi_matrix[overlap_regs, overlap_regs, drop = FALSE]
    csi_sub[is.na(csi_sub)] <- 0
    hr <- hclust(as.dist(1 - csi_sub), method = "ward.D2")
    modules <- cutree(hr, k = min(4, nrow(csi_sub)))

    anno_reg <- data.frame(csi_module = factor(modules))
    rownames(anno_reg) <- names(modules)

    png(paths$figure_6c_png, width = 1300, height = 1100, res = 150)
    pheatmap(
      csi_sub,
      color = colorRampPalette(c("#440154", "#21908C", "#FDE725"))(100),
      cluster_rows = hr,
      cluster_cols = hr,
      annotation_row = anno_reg,
      annotation_col = anno_reg,
      annotation_colors = list(csi_module = module_colors),
      main = "Figure 6C: CSI clustering"
    )
    dev.off()

    module_activity <- do.call(rbind, lapply(sort(unique(modules)), function(mod_id) {
      regs <- names(modules[modules == mod_id])
      module_mean <- sapply(names(celltype_colors), function(ct) {
        ct_cells <- names(cell_types)[cell_types == ct]
        mean(auc_features_x_cells[regs, ct_cells, drop = FALSE])
      })
      data.frame(module = as.character(mod_id), t(module_mean), check.names = FALSE)
    }))
    rownames(module_activity) <- paste0("Module_", module_activity$module)
    module_activity$module <- NULL
    module_activity <- as.matrix(module_activity)

    png(paths$figure_6d_png, width = 900, height = 700, res = 150)
    pheatmap(
      module_activity,
      color = colorRampPalette(c("#440154", "#21908C", "#FDE725"))(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      main = "Figure 6D: CSI module activity"
    )
    dev.off()
  }

  saveRDS(obj, paths$scenic_integrated_object_rds)

  if (identical(layer_id, cfg$panorama_layer_id)) {
    file.copy(paths$rss_matrix_csv, file.path(cfg$scenic_output_dir, "rss_matrix.csv"), overwrite = TRUE)
    file.copy(paths$top_regulons_by_celltype_csv, file.path(cfg$scenic_output_dir, "top_regulons_by_celltype.csv"), overwrite = TRUE)
    if (file.exists(paths$csi_matrix_csv)) {
      file.copy(paths$csi_matrix_csv, file.path(cfg$scenic_output_dir, "csi_matrix.csv"), overwrite = TRUE)
    }
    if (file.exists(paths$csi_source_note_txt)) {
      file.copy(paths$csi_source_note_txt, file.path(cfg$scenic_output_dir, "csi_original_source_used.txt"), overwrite = TRUE)
    }
    for (figure_path in c(paths$figure_6a_png, paths$figure_6b_png, paths$figure_6b_pdf, paths$figure_6c_png, paths$figure_6d_png)) {
      if (file.exists(figure_path)) {
        file.copy(figure_path, file.path(cfg$figure_dir, basename(figure_path)), overwrite = TRUE)
      }
    }
    file.copy(paths$scenic_integrated_object_rds, file.path(cfg$checkpoint_dir, "scenic_integrated_object.rds"), overwrite = TRUE)
  }

  index_rows[[length(index_rows) + 1L]] <- data.frame(
    layer_id = layer_id,
    layer_role = task$layer_role[[1]],
    object_rds = normalize_path_07(input_rds),
    auc_matrix_csv = normalize_path_07(auc_csv),
    rss_matrix_csv = normalize_path_07(paths$rss_matrix_csv),
    top_regulons_by_celltype_csv = normalize_path_07(paths$top_regulons_by_celltype_csv),
    csi_matrix_csv = if (file.exists(paths$csi_matrix_csv)) normalize_path_07(paths$csi_matrix_csv) else "",
    scenic_integrated_object_rds = normalize_path_07(paths$scenic_integrated_object_rds),
    figure_6a_png = normalize_path_07(paths$figure_6a_png),
    figure_6b_png = normalize_path_07(paths$figure_6b_png),
    figure_6b_pdf = normalize_path_07(paths$figure_6b_pdf),
    figure_6c_png = if (file.exists(paths$figure_6c_png)) normalize_path_07(paths$figure_6c_png) else "",
    figure_6d_png = if (file.exists(paths$figure_6d_png)) normalize_path_07(paths$figure_6d_png) else "",
    cell_type_col = cell_type_col,
    cell_n = ncol(obj),
    regulon_n = nrow(auc_features_x_cells),
    stringsAsFactors = FALSE
  )

  dynamic_outputs[[paste0("rss_matrix_csv_", safe_layer)]] <- build_output_entry(paths$rss_matrix_csv, "csv", module_name, "SCENIC RSS matrix", base_dir = cfg$project_root, schema = infer_schema_from_df(rss_df))
  dynamic_outputs[[paste0("top_regulons_by_celltype_csv_", safe_layer)]] <- build_output_entry(paths$top_regulons_by_celltype_csv, "csv", module_name, "top SCENIC regulons by cell type", base_dir = cfg$project_root, schema = infer_schema_from_df(top_regulons))
  if (file.exists(paths$csi_matrix_csv)) {
    dynamic_outputs[[paste0("csi_matrix_csv_", safe_layer)]] <- build_output_entry(paths$csi_matrix_csv, "csv", module_name, "SCENIC CSI matrix", base_dir = cfg$project_root)
  }
  dynamic_outputs[[paste0("scenic_integrated_object_rds_", safe_layer)]] <- build_output_entry(paths$scenic_integrated_object_rds, "rds", module_name, "Seurat object with SCENIC assay", base_dir = cfg$project_root)
  for (figure_name in c("figure_6a_png", "figure_6b_png", "figure_6b_pdf", "figure_6c_png", "figure_6d_png")) {
    if (file.exists(paths[[figure_name]])) {
      dynamic_outputs[[paste0(figure_name, "_", safe_layer)]] <- build_output_entry(paths[[figure_name]], tools::file_ext(paths[[figure_name]]), module_name, "SCENIC downstream figure", base_dir = cfg$project_root)
    }
  }
}

index_df <- dplyr::bind_rows(index_rows)
write_tsv_local(index_df, cfg$scenic_downstream_index_tsv)

if (file.exists(cfg$module_08d_manifest_path)) {
  unlink(cfg$module_08d_manifest_path)
}
fixed_outputs <- list(
  scenic_downstream_index_tsv = build_output_entry(cfg$scenic_downstream_index_tsv, "tsv", module_name, "one row per SCENIC downstream layer", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df))
)
write_manifest_local(
  manifest_path = cfg$module_08d_manifest_path,
  new_outputs = c(fixed_outputs, dynamic_outputs),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_08c = cfg$module_08c_manifest_path
  ),
  version = cfg$module_version,
  depends_on = list(
    module_08c = cfg$module_08c_manifest_path,
    module_03d = cfg$module_03d_manifest_path,
    module_04b = cfg$module_04b_manifest_path
  )
)

message("08d completed. index: ", cfg$scenic_downstream_index_tsv)
