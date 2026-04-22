source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SCENIC)
  library(AUCell)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
})

named_discrete_palette <- function(levels, base_colors = NULL) {
  levels <- unique(as.character(levels))
  levels <- levels[nzchar(levels) & !is.na(levels)]
  if (length(levels) == 0) {
    return(setNames(character(0), character(0)))
  }

  if (is.null(base_colors) || length(base_colors) == 0) {
    base_colors <- c(
      "#79BB7D", "#C4A9D6", "#FDBE6F", "#2F6CB3", "#E30073",
      "#A54F37", "#69706A", "#0E9B84", "#D95F02", "#6B63B7"
    )
  }
  if (length(levels) > length(base_colors)) {
    base_colors <- c(base_colors, scales::hue_pal()(length(levels) - length(base_colors)))
  }
  setNames(base_colors[seq_along(levels)], levels)
}

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

input_rds <- file.path(cfg$checkpoint_dir, "03_after_annotation.rds")
auc_csv <- file.path(cfg$scenic_output_dir, "auc_matrix.csv")

if (!file.exists(input_rds)) {
  stop("Missing checkpoint: 03_after_annotation.rds", call. = FALSE)
}
if (!file.exists(auc_csv)) {
  stop("Missing SCENIC AUC matrix: auc_matrix.csv", call. = FALSE)
}

module_colors <- c(
  "1" = "#8BC34A",
  "2" = "#F48FB1",
  "3" = "#26C6DA",
  "4" = "#CE93D8"
)

obj <- readRDS(input_rds)
obj <- maybe_join_layers(obj)

auc_raw <- read.csv(auc_csv, row.names = 1, check.names = FALSE)
if (nrow(auc_raw) == ncol(obj)) {
  auc_cells_x_regulons <- auc_raw
} else if (ncol(auc_raw) == ncol(obj)) {
  auc_cells_x_regulons <- t(auc_raw)
} else {
  stop("auc_matrix.csv dimensions do not match the Seurat object", call. = FALSE)
}

common_cells <- intersect(rownames(auc_cells_x_regulons), colnames(obj))
if (length(common_cells) == 0) {
  stop("No overlapping cells between auc_matrix.csv and the Seurat object", call. = FALSE)
}

auc_cells_x_regulons <- auc_cells_x_regulons[common_cells, , drop = FALSE]
obj <- subset(obj, cells = common_cells)

auc_features_x_cells <- t(as.matrix(auc_cells_x_regulons))
rownames(auc_features_x_cells) <- gsub("\\(_\\+\\)|\\(\\+\\)", "", rownames(auc_features_x_cells))
obj[["SCENIC"]] <- CreateAssayObject(data = auc_features_x_cells)

cell_types <- as.character(obj$cell_type)
names(cell_types) <- colnames(obj)
celltype_colors <- named_discrete_palette(unique(cell_types))
rss_mat <- calcRSS(AUC = auc_features_x_cells, cellAnnotation = cell_types)

rss_df <- as.data.frame(rss_mat)
rss_df$regulon <- rownames(rss_df)
write.csv(rss_df, file.path(cfg$scenic_output_dir, "rss_matrix.csv"), row.names = FALSE)

top_regulons <- bind_rows(lapply(colnames(rss_mat), function(ct) {
  tibble(
    cell_type = ct,
    regulon = rownames(rss_mat),
    rss = rss_mat[, ct]
  ) %>%
    arrange(desc(rss)) %>%
    slice_head(n = 5)
})) %>%
  distinct(regulon, .keep_all = TRUE)
write.csv(top_regulons, file.path(cfg$scenic_output_dir, "top_regulons_by_celltype.csv"), row.names = FALSE)

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

png(file.path(cfg$figure_dir, "Figure_6A_Regulon_Heatmap.png"), width = 1600, height = 950, res = 150)
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

rss_long <- bind_rows(lapply(colnames(rss_mat), function(ct) {
  tibble(
    cell_type = ct,
    regulon = rownames(rss_mat),
    rss = rss_mat[, ct]
  ) %>%
    arrange(desc(rss)) %>%
    mutate(rank = row_number())
}))

rss_long <- rss_long %>%
  group_by(cell_type) %>%
  mutate(label = ifelse(rank <= 3, regulon, NA_character_)) %>%
  ungroup()

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

ggsave(file.path(cfg$figure_dir, "Figure_6B_RSS_Ranks.png"), fig6b, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(cfg$figure_dir, "Figure_6B_RSS_Ranks.pdf"), fig6b, width = 10, height = 8, bg = "white")

csi_long <- NULL
if (requireNamespace("scFunctions", quietly = TRUE)) {
  auc_rankings <- AUCell_buildRankings(auc_features_x_cells, plotStats = FALSE, nCores = 1)
  auc_rankings@assays@data@listData$AUC <- auc_features_x_cells
  csi_long <- scFunctions::calculate_csi(auc_rankings, calc_extended = FALSE, verbose = FALSE)
} else {
  source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "vendor", "calculate_csi_scfunctions_original.R"))
  auc_rankings <- AUCell_buildRankings(auc_features_x_cells, plotStats = FALSE, nCores = 1)
  auc_rankings@assays@data@listData$AUC <- auc_features_x_cells
  csi_long <- calculate_csi(auc_rankings, calc_extended = FALSE, verbose = FALSE)
  writeLines(
    "scFunctions package was unavailable, so the original calculate_csi source was used directly.",
    file.path(cfg$scenic_output_dir, "csi_original_source_used.txt")
  )
}

csi_matrix <- NULL
if (!is.null(csi_long) && nrow(csi_long) > 0) {
  csi_matrix <- csi_long %>%
    pivot_wider(names_from = regulon_2, values_from = CSI) %>%
    column_to_rownames("regulon_1") %>%
    as.matrix()
  write.csv(csi_matrix, file.path(cfg$scenic_output_dir, "csi_matrix.csv"))
}

if (!is.null(csi_matrix) && nrow(csi_matrix) > 1) {
  overlap_regs <- intersect(heatmap_regs, rownames(csi_matrix))
  csi_sub <- csi_matrix[overlap_regs, overlap_regs, drop = FALSE]
  csi_sub[is.na(csi_sub)] <- 0
  hr <- hclust(as.dist(1 - csi_sub), method = "ward.D2")
  modules <- cutree(hr, k = min(4, nrow(csi_sub)))

  anno_reg <- data.frame(csi_module = factor(modules))
  rownames(anno_reg) <- names(modules)

  png(file.path(cfg$figure_dir, "Figure_6C_CSI_Clustering_Heatmap.png"), width = 1300, height = 1100, res = 150)
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

  png(file.path(cfg$figure_dir, "Figure_6D_CSI_Module_Activity_Heatmap.png"), width = 900, height = 700, res = 150)
  pheatmap(
    module_activity,
    color = colorRampPalette(c("#440154", "#21908C", "#FDE725"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Figure 6D: CSI module activity"
  )
  dev.off()
}

saveRDS(obj, file.path(cfg$checkpoint_dir, "scenic_integrated_object.rds"))
message("SCENIC downstream complete: Figure_6A-6D and RSS/CSI matrices are ready.")
