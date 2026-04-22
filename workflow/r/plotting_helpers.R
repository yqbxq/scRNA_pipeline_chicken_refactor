suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

paper_cluster_palette_15 <- function() {
  c(
    "1" = "#74BA73",
    "2" = "#C8C6E4",
    "3" = "#FDBE6F",
    "4" = "#2E63AF",
    "5" = "#E30073",
    "6" = "#A54F37",
    "7" = "#69706A",
    "8" = "#0E9B84",
    "9" = "#D95F02",
    "10" = "#6B63B7",
    "11" = "#D90445",
    "12" = "#1893D1",
    "13" = "#FDBA12",
    "14" = "#4C2A7A",
    "15" = "#D6457A"
  )
}

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
    extra_needed <- length(levels) - length(base_colors)
    extra_colors <- scales::hue_pal()(extra_needed)
    base_colors <- c(base_colors, extra_colors)
  }

  setNames(base_colors[seq_along(levels)], levels)
}

paper_celltype_palette <- function(levels = NULL) {
  named_discrete_palette(levels)
}

paper_feature_palette <- function() {
  c(
    "#2C7BB6", "#00A6CA", "#00CCBC", "#90EB9D",
    "#FFFF8C", "#F9D057", "#F29E2E", "#E76818", "#D7191C"
  )
}

paper_dotplot_palette <- function() {
  paper_feature_palette()
}

paper_heatmap_palette <- function() {
  c("#3B4992FF", "#FFFFFF", "#EE0000FF")
}

resolve_plot_group_var <- function(seu) {
  if ("group_id" %in% colnames(seu@meta.data)) {
    group_values <- unique(stats::na.omit(as.character(seu$group_id)))
    if (length(group_values) >= 2) {
      return("group_id")
    }
  }
  if ("condition" %in% colnames(seu@meta.data)) {
    group_values <- unique(stats::na.omit(as.character(seu$condition)))
    if (length(group_values) >= 2) {
      return("condition")
    }
  }
  if ("analysis_group" %in% colnames(seu@meta.data)) {
    group_values <- unique(stats::na.omit(as.character(seu$analysis_group)))
    if (length(group_values) >= 2) {
      return("analysis_group")
    }
  }
  "orig.ident"
}

umap_plot_df <- function(seu, vars) {
  umap_mat <- as.data.frame(Embeddings(seu, reduction = "umap"))
  if (ncol(umap_mat) < 2) {
    stop("对象缺少 UMAP 降维结果，无法绘图", call. = FALSE)
  }
  colnames(umap_mat)[1:2] <- c("UMAP_1", "UMAP_2")
  meta_df <- FetchData(seu, vars = unique(vars))
  df <- cbind(umap_mat, meta_df[rownames(umap_mat), , drop = FALSE])
  tibble::rownames_to_column(df, "cell_id")
}

cluster_label_text <- function(seu, cluster_var = "seurat_clusters") {
  counts <- table(seu[[cluster_var, drop = TRUE]])
  counts <- counts[names(counts) %in% names(paper_cluster_palette_15())]
  setNames(
    paste0(names(counts), "-", as.integer(counts), " cells"),
    names(counts)
  )
}

cluster_centers_by_group <- function(df, group_var, cluster_var) {
  df %>%
    group_by(.data[[group_var]], .data[[cluster_var]]) %>%
    summarise(
      UMAP_1 = stats::median(UMAP_1),
      UMAP_2 = stats::median(UMAP_2),
      .groups = "drop"
    )
}

celltype_centers_by_group <- function(df, group_var, celltype_var = "cell_type") {
  df %>%
    group_by(.data[[group_var]], .data[[celltype_var]]) %>%
    summarise(
      UMAP_1 = stats::median(UMAP_1),
      UMAP_2 = stats::median(UMAP_2),
      .groups = "drop"
    )
}

save_plot_dual <- function(plot_obj, png_path, width, height, dpi = 300, pdf_path = NULL) {
  if (is.null(pdf_path)) {
    pdf_path <- sub("\\.png$", ".pdf", png_path)
  }
  ggsave(png_path, plot_obj, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(pdf_path, plot_obj, width = width, height = height, bg = "white")
}

plot_cluster_split_umap_paper <- function(seu, group_var = NULL, cluster_var = "seurat_clusters", pt_size = 0.28) {
  if (is.null(group_var)) {
    group_var <- resolve_plot_group_var(seu)
  }

  palette_15 <- paper_cluster_palette_15()
  df <- umap_plot_df(seu, c(group_var, cluster_var))
  df[[cluster_var]] <- factor(as.character(df[[cluster_var]]), levels = names(palette_15))
  centers <- cluster_centers_by_group(df, group_var, cluster_var)
  label_map <- cluster_label_text(seu, cluster_var = cluster_var)

  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = .data[[cluster_var]])) +
    geom_point(size = pt_size, alpha = 0.9) +
    geom_text(
      data = centers,
      aes(x = UMAP_1, y = UMAP_2, label = .data[[cluster_var]]),
      color = "black",
      fontface = "bold",
      size = 4.4,
      inherit.aes = FALSE
    ) +
    facet_wrap(stats::as.formula(paste("~", group_var)), nrow = 1) +
    scale_color_manual(values = palette_15, labels = label_map, drop = FALSE) +
    labs(x = "UMAP_1", y = "UMAP_2", color = NULL) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 18, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
}

plot_cluster_proportion_paper <- function(seu, group_var = NULL, cluster_var = "seurat_clusters") {
  if (is.null(group_var)) {
    group_var <- resolve_plot_group_var(seu)
  }

  palette_15 <- paper_cluster_palette_15()
  df <- seu@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    mutate(
      cluster_value = factor(as.character(.data[[cluster_var]]), levels = names(palette_15)),
      group_value = .data[[group_var]]
    ) %>%
    group_by(group_value, cluster_value) %>%
    summarise(count = dplyr::n(), .groups = "drop") %>%
    group_by(group_value) %>%
    mutate(proportion = count / sum(count))

  ggplot(df, aes(x = group_value, y = proportion, fill = cluster_value)) +
    geom_col(width = 0.84, color = NA, position = "fill") +
    scale_fill_manual(values = palette_15, drop = FALSE) +
    scale_y_continuous(labels = label_percent(accuracy = 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion [%]", fill = "clusters") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold")
    )
}

plot_figure2_paper <- function(seu, group_var = NULL, cluster_var = "seurat_clusters") {
  p_umap <- plot_cluster_split_umap_paper(seu, group_var = group_var, cluster_var = cluster_var)
  p_prop <- plot_cluster_proportion_paper(seu, group_var = group_var, cluster_var = cluster_var)

  (p_umap | p_prop) +
    plot_layout(widths = c(2.4, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 26, face = "bold"))
}

plot_celltype_split_umap_paper <- function(seu, group_var = NULL, celltype_var = "cell_type", pt_size = 0.32) {
  if (is.null(group_var)) {
    group_var <- resolve_plot_group_var(seu)
  }

  palette_ct <- paper_celltype_palette(unique(as.character(df[[celltype_var]])))
  df <- umap_plot_df(seu, c(group_var, celltype_var))
  df[[celltype_var]] <- factor(as.character(df[[celltype_var]]), levels = names(palette_ct))
  centers <- celltype_centers_by_group(df, group_var, celltype_var)

  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = .data[[celltype_var]])) +
    geom_point(size = pt_size, alpha = 0.9) +
    geom_text(
      data = centers,
      aes(x = UMAP_1, y = UMAP_2, label = .data[[celltype_var]]),
      fontface = "bold",
      size = 5,
      color = "black",
      inherit.aes = FALSE
    ) +
    facet_wrap(stats::as.formula(paste("~", group_var)), nrow = 1) +
    scale_color_manual(values = palette_ct, drop = FALSE) +
    labs(x = "UMAP_1", y = "UMAP_2", color = NULL) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 18, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 12)
    )
}

plot_feature_umap_paper <- function(seu, gene, low_high = NULL, pt_size = 0.22) {
  if (is.null(low_high)) {
    low_high <- paper_feature_palette()
  }
  df <- umap_plot_df(seu, gene)
  colnames(df)[colnames(df) == gene] <- "expr_value"

  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = expr_value)) +
    geom_point(size = pt_size) +
    scale_color_gradientn(colors = low_high) +
    labs(title = gene, x = "UMAP_1", y = "UMAP_2", color = NULL) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11),
      legend.position = "right"
    )
}

plot_feature_rows_paper <- function(seu, gene_rows) {
  row_plots <- lapply(gene_rows, function(genes) {
    plots <- lapply(genes, function(gene) plot_feature_umap_paper(seu, gene))
    wrap_plots(plots, nrow = 1, guides = "collect") &
      theme(legend.position = "right")
  })
  wrap_plots(row_plots, ncol = 1)
}

plot_marker_dotplot_paper <- function(seu, genes, celltype_var = "cell_type") {
  DotPlot(
    seu,
    features = genes,
    group.by = celltype_var,
    col.min = -1,
    col.max = 1
  ) +
    scale_color_gradientn(
      colours = paper_dotplot_palette(),
      limits = c(-1, 1),
      oob = squish,
      name = "Exp.avg"
    ) +
    coord_flip() +
    labs(x = "Cell type", y = "Gene") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 11, face = "italic"),
      legend.position = "right"
    )
}

plot_celltype_proportion_paper <- function(seu, group_var = NULL, celltype_var = "cell_type") {
  if (is.null(group_var)) {
    group_var <- resolve_plot_group_var(seu)
  }

  palette_ct <- paper_celltype_palette(unique(as.character(seu[[celltype_var, drop = TRUE]])))
  df <- seu@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    mutate(
      celltype_value = factor(as.character(.data[[celltype_var]]), levels = names(palette_ct)),
      group_value = .data[[group_var]]
    ) %>%
    group_by(group_value, celltype_value) %>%
    summarise(count = dplyr::n(), .groups = "drop") %>%
    group_by(group_value) %>%
    mutate(proportion = count / sum(count))

  ggplot(df, aes(x = group_value, y = proportion, fill = celltype_value)) +
    geom_col(width = 0.84, color = NA, position = "fill") +
    scale_fill_manual(values = palette_ct, drop = FALSE) +
    scale_y_continuous(labels = label_percent(accuracy = 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion [%]", fill = "Cell type") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold")
    )
}

plot_figure3_bottom_paper <- function(seu, genes, group_var = NULL, celltype_var = "cell_type") {
  dot_plot <- plot_marker_dotplot_paper(seu, genes = genes, celltype_var = celltype_var)
  prop_plot <- plot_celltype_proportion_paper(seu, group_var = group_var, celltype_var = celltype_var)

  (dot_plot | prop_plot) +
    plot_layout(widths = c(1.9, 1))
}
