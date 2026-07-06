suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(42)

env_or <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

project_root <- env_or("PROJECT_ROOT", "/home/user_test/syf_f5/pre-sub/F6_1_panorama_project_panel_v3")
input_checkpoint <- env_or(
  "INPUT_CHECKPOINT",
  file.path(project_root, "results", "checkpoints", "05_after_manual_panorama_ovary_annotation.rds")
)
gc_panel <- env_or(
  "GC_PANEL",
  file.path(project_root, "config", "marker_panels", "GC_subcluster_panel.tsv")
)
parent_cluster_col <- env_or("PARENT_CLUSTER_COL", "seurat_clusters")
gc_compartment <- env_or("GC_COMPARTMENT", "GC-like / follicular somatic compartment")
resolution <- as.numeric(env_or("GC_RESOLUTION", "0.35"))
dims_max <- as.integer(env_or("GC_DIMS_MAX", "20"))

checkpoint_dir <- file.path(project_root, "results", "checkpoints")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures", "GC_manual_subcluster")
report_dir <- file.path(project_root, "reports", "eda", "GC_manual_subcluster")
annotation_dir <- file.path(table_dir, "annotation", "layers", "GC_manual_subcluster")
marker_dir <- file.path(table_dir, "marker_discovery", "layers", "GC_manual_subcluster")

dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

strip_ensembl_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

match_features <- function(features, genes) {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0) return(character(0))
  feature_upper <- toupper(features)
  feature_stripped <- strip_ensembl_version(features)
  gene_upper <- toupper(genes)
  gene_stripped <- strip_ensembl_version(genes)
  hit <- features %in% genes |
    feature_upper %in% gene_upper |
    feature_stripped %in% gene_stripped |
    toupper(feature_stripped) %in% gene_upper
  unique(features[hit])
}

get_norm_data <- function(obj) {
  tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
  )
}

logfc_column <- function(markers) {
  candidates <- c("avg_log2FC", "avg_logFC")
  hit <- candidates[candidates %in% colnames(markers)]
  if (length(hit) > 0) hit[[1]] else NA_character_
}

top_marker_table <- function(markers, n = 50) {
  if (nrow(markers) == 0) return(markers)
  fc_col <- logfc_column(markers)
  markers$.rank_fc <- if (!is.na(fc_col)) markers[[fc_col]] else 0
  if ("p_val_adj" %in% colnames(markers)) {
    markers <- markers %>% arrange(cluster, p_val_adj, desc(.rank_fc))
  } else {
    markers <- markers %>% arrange(cluster, desc(.rank_fc))
  }
  markers %>%
    group_by(cluster) %>%
    slice_head(n = n) %>%
    ungroup() %>%
    select(-.rank_fc)
}

run_gc_clustering <- function(sub_obj) {
  DefaultAssay(sub_obj) <- "RNA"
  log_msg("NormalizeData / FindVariableFeatures / PCA / clustering for manual GC cells")
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  sub_obj <- FindVariableFeatures(sub_obj, nfeatures = min(2000, nrow(sub_obj)), verbose = FALSE)
  var_features <- VariableFeatures(sub_obj)
  if (length(var_features) < 50) var_features <- rownames(sub_obj)
  sub_obj <- ScaleData(sub_obj, features = var_features, verbose = FALSE)
  npcs <- min(dims_max, length(var_features), max(2, ncol(sub_obj) - 1))
  sub_obj <- RunPCA(sub_obj, features = var_features, npcs = npcs, verbose = FALSE)
  dims_use <- seq_len(min(dims_max, ncol(Embeddings(sub_obj, "pca"))))
  sub_obj <- FindNeighbors(sub_obj, dims = dims_use, verbose = FALSE)
  sub_obj <- FindClusters(sub_obj, resolution = resolution, random.seed = 42, verbose = FALSE)
  if (ncol(sub_obj) >= 30) {
    sub_obj <- RunUMAP(
      sub_obj,
      dims = dims_use,
      n.neighbors = min(30, ncol(sub_obj) - 1),
      verbose = FALSE
    )
  }
  sub_obj$GC_manual_subcluster_id <- as.character(sub_obj$seurat_clusters)
  sub_obj
}

find_markers_safe <- function(obj, cluster_col) {
  Idents(obj) <- obj[[cluster_col, drop = TRUE]]
  clusters <- unique(as.character(Idents(obj)))
  if (length(clusters) < 2) return(data.frame())
  log_msg("FindAllMarkers for ", cluster_col, " (", length(clusters), " clusters)")
  markers <- tryCatch(
    FindAllMarkers(
      obj,
      only.pos = TRUE,
      assay = "RNA",
      min.pct = 0.10,
      logfc.threshold = 0.25,
      max.cells.per.ident = 1500,
      random.seed = 42,
      verbose = FALSE
    ),
    error = function(e) {
      warning("FindAllMarkers failed: ", conditionMessage(e), call. = FALSE)
      data.frame()
    }
  )
  if (nrow(markers) > 0 && !"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
  markers
}

annotate_gc_subclusters <- function(sub_obj, markers, panel_path) {
  panel <- read.delim(panel_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  panel <- panel[panel$layer_id == "GC_subcluster", , drop = FALSE]
  if (!"confidence_ceiling" %in% colnames(panel)) panel$confidence_ceiling <- ""
  matched <- lapply(panel$gene, function(gene) match_features(rownames(sub_obj), gene))
  panel$feature <- vapply(matched, function(x) if (length(x) > 0) x[[1]] else "", character(1))
  panel$matched <- nzchar(panel$feature)
  write_tsv(panel, file.path(annotation_dir, "panel_feature_match.tsv"))

  matched_panel <- panel[panel$matched, , drop = FALSE]
  clusters <- sort(unique(as.character(sub_obj$GC_manual_subcluster_id)))
  top_markers <- top_marker_table(markers, n = 100)
  data_mat <- get_norm_data(sub_obj)
  evidence <- list()

  for (cluster in clusters) {
    cells <- colnames(sub_obj)[as.character(sub_obj$GC_manual_subcluster_id) == cluster]
    cluster_markers <- character(0)
    if (nrow(top_markers) > 0) {
      cluster_markers <- unique(as.character(top_markers$gene[as.character(top_markers$cluster) == cluster]))
    }
    cluster_marker_upper <- toupper(cluster_markers)

    for (celltype in sort(unique(matched_panel$celltype))) {
      panel_ct <- matched_panel[matched_panel$celltype == celltype, , drop = FALSE]
      panel_features <- unique(panel_ct$feature)
      panel_features <- panel_features[nzchar(panel_features)]
      if (length(panel_features) == 0) next
      expr_sub <- data_mat[panel_features, cells, drop = FALSE]
      module_mean <- mean(Matrix::colMeans(expr_sub), na.rm = TRUE)
      pct_expr <- mean(Matrix::rowMeans(expr_sub > 0), na.rm = TRUE)
      overlap <- panel_features[toupper(panel_features) %in% cluster_marker_upper]
      core_features <- panel_ct$feature[panel_ct$confidence_ceiling == "确定"]
      core_overlap <- unique(overlap[overlap %in% core_features])
      score <- length(core_overlap) * 8 + length(setdiff(overlap, core_overlap)) * 2 + module_mean + pct_expr
      evidence[[length(evidence) + 1]] <- data.frame(
        layer_id = "GC_manual_subcluster",
        subcluster = cluster,
        candidate_celltype = celltype,
        n_cells = length(cells),
        matched_panel_feature_count = length(panel_features),
        marker_overlap_n = length(overlap),
        core_marker_overlap_n = length(core_overlap),
        module_score_mean = module_mean,
        pct_cells_expressing_panel_mean = pct_expr,
        annotation_score = score,
        marker_overlap_genes = paste(overlap, collapse = ","),
        core_marker_overlap_genes = paste(core_overlap, collapse = ","),
        matched_panel_features = paste(panel_features, collapse = ","),
        stringsAsFactors = FALSE
      )
    }
  }

  evidence_df <- if (length(evidence) > 0) bind_rows(evidence) else data.frame()
  write_tsv(evidence_df, file.path(annotation_dir, "annotation_evidence.tsv"))

  counts <- sub_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    count(GC_manual_subcluster_id, name = "n_cells")

  if (nrow(evidence_df) == 0) {
    annotation <- counts %>%
      transmute(
        layer_id = "GC_manual_subcluster",
        subcluster = GC_manual_subcluster_id,
        n_cells = n_cells,
        candidate_celltype = paste0("GC_manual_subcluster_", subcluster),
        annotation_label = paste0("GC_manual_subcluster_", subcluster),
        confidence = "no_panel_match",
        annotation_score = 0,
        marker_overlap_n = 0,
        core_marker_overlap_n = 0,
        module_score_mean = 0,
        marker_overlap_genes = "",
        core_marker_overlap_genes = "",
        matched_panel_features = ""
      )
  } else {
    annotation <- evidence_df %>%
      group_by(subcluster) %>%
      arrange(desc(annotation_score), desc(core_marker_overlap_n), desc(marker_overlap_n), desc(module_score_mean), .by_group = TRUE) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      mutate(
        passes_label_rule = core_marker_overlap_n >= 1 | marker_overlap_n >= 2 | module_score_mean > 0.20,
        annotation_label = ifelse(
          passes_label_rule,
          paste0(candidate_celltype, "_candidate"),
          paste0("GC_manual_subcluster_", subcluster)
        ),
        confidence = case_when(
          !passes_label_rule ~ "exploratory_no_strong_panel_label",
          core_marker_overlap_n >= 2 ~ "supported_core_markers",
          core_marker_overlap_n >= 1 & marker_overlap_n >= 2 ~ "supported",
          core_marker_overlap_n >= 1 ~ "tentative_core_marker",
          marker_overlap_n >= 2 ~ "tentative_multi_marker",
          TRUE ~ "module_score_only"
        )
      ) %>%
      select(
        layer_id, subcluster, n_cells, candidate_celltype, annotation_label,
        confidence, annotation_score, marker_overlap_n, core_marker_overlap_n,
        module_score_mean, pct_cells_expressing_panel_mean, marker_overlap_genes,
        core_marker_overlap_genes, matched_panel_features
      )
  }

  annotation <- annotation %>% arrange(as.integer(subcluster))
  write_tsv(annotation, file.path(annotation_dir, "annotation_table.tsv"))
  list(annotation = annotation, evidence = evidence_df, panel = panel)
}

plot_if_umap <- function(obj, group_col, title, file_name, width = 8, height = 6) {
  if (!"umap" %in% names(obj@reductions)) return(invisible(FALSE))
  p <- DimPlot(obj, reduction = "umap", group.by = group_col, label = TRUE, repel = TRUE) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(figure_dir, file_name), p, width = width, height = height, dpi = 180, limitsize = FALSE)
  invisible(TRUE)
}

log_msg("Project root: ", project_root)
log_msg("Input checkpoint: ", input_checkpoint)
log_msg("GC panel: ", gc_panel)

if (!file.exists(input_checkpoint)) stop("Missing input checkpoint: ", input_checkpoint, call. = FALSE)
if (!file.exists(gc_panel)) stop("Missing GC panel: ", gc_panel, call. = FALSE)

obj <- readRDS(input_checkpoint)
if (!parent_cluster_col %in% colnames(obj@meta.data)) {
  stop("Missing parent cluster column: ", parent_cluster_col, call. = FALSE)
}
if (!"panorama_ovary_compartment" %in% colnames(obj@meta.data)) {
  stop("Missing panorama_ovary_compartment; run manual panorama annotation first.", call. = FALSE)
}

obj$parent_panorama_cluster <- as.character(obj@meta.data[[parent_cluster_col]])
gc_cells <- colnames(obj)[obj$panorama_ovary_compartment == gc_compartment]
if (length(gc_cells) < 50) stop("Too few manual GC cells: ", length(gc_cells), call. = FALSE)
log_msg("Manual GC cells: ", length(gc_cells))

sub_obj <- subset(obj, cells = gc_cells)
sub_obj <- run_gc_clustering(sub_obj)
markers <- find_markers_safe(sub_obj, "GC_manual_subcluster_id")
write_tsv(markers, file.path(marker_dir, "cluster_markers.tsv"))
top_markers <- top_marker_table(markers, n = 50)
write_tsv(top_markers, file.path(marker_dir, "top50_cluster_markers.tsv"))

ann <- annotate_gc_subclusters(sub_obj, markers, gc_panel)
annotation <- ann$annotation
cluster_to_label <- setNames(annotation$annotation_label, annotation$subcluster)
cluster_to_candidate <- setNames(annotation$candidate_celltype, annotation$subcluster)
cluster_to_confidence <- setNames(annotation$confidence, annotation$subcluster)

sub_clusters <- as.character(sub_obj$GC_manual_subcluster_id)
sub_obj$GC_manual_subcluster_label <- unname(cluster_to_label[sub_clusters])
sub_obj$GC_manual_subcluster_candidate <- unname(cluster_to_candidate[sub_clusters])
sub_obj$GC_manual_subcluster_confidence <- unname(cluster_to_confidence[sub_clusters])

cell_assignments <- sub_obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  transmute(
    cell_id = cell_id,
    parent_panorama_cluster = parent_panorama_cluster,
    parent_manual_annotation = panorama_ovary_annotation,
    GC_manual_subcluster_id = GC_manual_subcluster_id,
    GC_manual_subcluster_label = GC_manual_subcluster_label,
    GC_manual_subcluster_candidate = GC_manual_subcluster_candidate,
    GC_manual_subcluster_confidence = GC_manual_subcluster_confidence
  )
write_tsv(cell_assignments, file.path(annotation_dir, "cell_assignments.tsv"))

cluster_counts <- cell_assignments %>%
  count(GC_manual_subcluster_id, GC_manual_subcluster_label, GC_manual_subcluster_candidate, GC_manual_subcluster_confidence, name = "n_cells") %>%
  mutate(percent_gc_cells = round(100 * n_cells / sum(n_cells), 4)) %>%
  arrange(as.integer(GC_manual_subcluster_id))
write_tsv(cluster_counts, file.path(annotation_dir, "cluster_cell_counts.tsv"))

parent_composition <- cell_assignments %>%
  count(GC_manual_subcluster_id, parent_panorama_cluster, parent_manual_annotation, name = "n_cells") %>%
  group_by(GC_manual_subcluster_id) %>%
  mutate(percent_subcluster = round(100 * n_cells / sum(n_cells), 4)) %>%
  ungroup() %>%
  arrange(as.integer(GC_manual_subcluster_id), desc(n_cells))
write_tsv(parent_composition, file.path(annotation_dir, "parent_panorama_cluster_composition.tsv"))

top_gene_summary <- if (nrow(top_markers) > 0) {
  top_markers %>%
    group_by(cluster) %>%
    summarise(top_genes = paste(head(unique(gene), 12), collapse = ","), .groups = "drop") %>%
    rename(GC_manual_subcluster_id = cluster)
} else {
  data.frame(GC_manual_subcluster_id = character(0), top_genes = character(0))
}

cluster_summary <- cluster_counts %>%
  left_join(annotation, by = c("GC_manual_subcluster_id" = "subcluster")) %>%
  left_join(top_gene_summary, by = "GC_manual_subcluster_id") %>%
  arrange(as.integer(GC_manual_subcluster_id))

major_parent <- parent_composition %>%
  group_by(GC_manual_subcluster_id) %>%
  arrange(desc(n_cells), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    GC_manual_subcluster_id = GC_manual_subcluster_id,
    major_parent_panorama_cluster = parent_panorama_cluster,
    major_parent_manual_annotation = parent_manual_annotation,
    major_parent_percent = percent_subcluster
  )

refine_label <- function(candidate, confidence, marker_genes, core_genes, top_genes, major_parent) {
  text <- toupper(paste(candidate, confidence, marker_genes, core_genes, top_genes, major_parent, sep = " "))
  if (grepl("TPX2|KIF4B|TOP2A|SMC2|SPAG5|CKAP2|KIF15|CKS1B|PRC1|CENPF", text, perl = TRUE)) {
    return(c("proliferating_GC_like", "增殖型颗粒细胞样细胞", "supported_by_cell_cycle_markers", "Top markers include cell-cycle/proliferation genes."))
  }
  if (grepl("\\bRLN3\\b|\\bFABP5\\b|\\bTSPAN6\\b", text, perl = TRUE) || identical(candidate, "pGC")) {
    return(c("pGC_like", "原始/早期颗粒细胞样细胞", "tentative_panel_marker", "Primitive GC panel evidence is present."))
  }
  if (grepl("\\bCYP11A1\\b|\\bFGL2\\b|\\bZP3\\b", text, perl = TRUE)) {
    return(c("steroidogenic_rgGC_like", "类固醇合成/快速生长期颗粒细胞样细胞", "tentative_panel_or_module", "rgGC/steroidogenic panel evidence is present but should be reviewed."))
  }
  if (grepl("\\bHBA1\\b|\\bHBBA\\b|\\bACTA2\\b|\\bTMSB4X\\b", text, perl = TRUE)) {
    return(c("mixed_GC_candidate_Hb_ACTA2_like", "混合型GC候选/Hb-ACTA2相关细胞", "review_mixed_signal", "Top markers include Hb/ACTA2-like mixed signals; review before treating as clean GC."))
  }
  if (grepl("\\bBMPR1B\\b|\\bESRRG\\b|\\bKCNQ1\\b|\\bNCAM2\\b", text, perl = TRUE)) {
    return(c("GC_like_BMPR1B_ESRRG_candidate", "BMPR1B/ESRRG 型颗粒细胞样候选细胞", "exploratory_marker_pattern", "Top markers define a GC-like exploratory pattern not covered by the small subtype panel."))
  }
  if (grepl("WEAK_GC_LIKE", toupper(major_parent), fixed = TRUE)) {
    return(c("weak_GC_like", "弱颗粒细胞样细胞", "parent_annotation_based", "Refined from the dominant parent panorama annotation."))
  }
  if (grepl("GC_LIKE_CANDIDATE", toupper(major_parent), fixed = TRUE)) {
    return(c("GC_like_candidate", "颗粒细胞样候选细胞", "parent_annotation_based", "Refined from the dominant parent panorama annotation."))
  }
  c(paste0(candidate, "_candidate"), paste0(candidate, " 候选细胞"), "panel_candidate", "Fallback to the best small-panel candidate.")
}

refined <- mapply(
  refine_label,
  candidate = as.character(cluster_summary$GC_manual_subcluster_candidate),
  confidence = as.character(cluster_summary$GC_manual_subcluster_confidence),
  marker_genes = as.character(cluster_summary$marker_overlap_genes),
  core_genes = as.character(cluster_summary$core_marker_overlap_genes),
  top_genes = as.character(cluster_summary$top_genes),
  major_parent = as.character(major_parent$major_parent_manual_annotation[match(cluster_summary$GC_manual_subcluster_id, major_parent$GC_manual_subcluster_id)]),
  SIMPLIFY = TRUE
)

cluster_summary <- cluster_summary %>%
  left_join(major_parent, by = "GC_manual_subcluster_id") %>%
  mutate(
    GC_manual_refined_subtype = refined[1, ],
    GC_manual_refined_subtype_zh = refined[2, ],
    GC_manual_refined_confidence = refined[3, ],
    GC_manual_refined_note = refined[4, ]
  )
write_tsv(cluster_summary, file.path(annotation_dir, "cluster_summary.tsv"))

refined_label_map <- setNames(cluster_summary$GC_manual_refined_subtype, cluster_summary$GC_manual_subcluster_id)
refined_label_zh_map <- setNames(cluster_summary$GC_manual_refined_subtype_zh, cluster_summary$GC_manual_subcluster_id)
refined_confidence_map <- setNames(cluster_summary$GC_manual_refined_confidence, cluster_summary$GC_manual_subcluster_id)
sub_obj$GC_manual_refined_subtype <- unname(refined_label_map[sub_clusters])
sub_obj$GC_manual_refined_subtype_zh <- unname(refined_label_zh_map[sub_clusters])
sub_obj$GC_manual_refined_confidence <- unname(refined_confidence_map[sub_clusters])

cell_assignments <- cell_assignments %>%
  mutate(
    GC_manual_refined_subtype = unname(refined_label_map[GC_manual_subcluster_id]),
    GC_manual_refined_subtype_zh = unname(refined_label_zh_map[GC_manual_subcluster_id]),
    GC_manual_refined_confidence = unname(refined_confidence_map[GC_manual_subcluster_id])
  )
write_tsv(cell_assignments, file.path(annotation_dir, "cell_assignments.tsv"))

cluster_counts <- cell_assignments %>%
  count(
    GC_manual_subcluster_id,
    GC_manual_subcluster_label,
    GC_manual_subcluster_candidate,
    GC_manual_subcluster_confidence,
    GC_manual_refined_subtype,
    GC_manual_refined_subtype_zh,
    GC_manual_refined_confidence,
    name = "n_cells"
  ) %>%
  mutate(percent_gc_cells = round(100 * n_cells / sum(n_cells), 4)) %>%
  arrange(as.integer(GC_manual_subcluster_id))
write_tsv(cluster_counts, file.path(annotation_dir, "cluster_cell_counts.tsv"))

obj_gc_subcluster <- rep(NA_character_, ncol(obj))
obj_gc_label <- rep(NA_character_, ncol(obj))
obj_gc_candidate <- rep(NA_character_, ncol(obj))
obj_gc_confidence <- rep(NA_character_, ncol(obj))
obj_gc_refined <- rep(NA_character_, ncol(obj))
obj_gc_refined_zh <- rep(NA_character_, ncol(obj))
obj_gc_refined_confidence <- rep(NA_character_, ncol(obj))
names(obj_gc_subcluster) <- colnames(obj)
names(obj_gc_label) <- colnames(obj)
names(obj_gc_candidate) <- colnames(obj)
names(obj_gc_confidence) <- colnames(obj)
names(obj_gc_refined) <- colnames(obj)
names(obj_gc_refined_zh) <- colnames(obj)
names(obj_gc_refined_confidence) <- colnames(obj)
obj_gc_subcluster[colnames(sub_obj)] <- as.character(sub_obj$GC_manual_subcluster_id)
obj_gc_label[colnames(sub_obj)] <- as.character(sub_obj$GC_manual_subcluster_label)
obj_gc_candidate[colnames(sub_obj)] <- as.character(sub_obj$GC_manual_subcluster_candidate)
obj_gc_confidence[colnames(sub_obj)] <- as.character(sub_obj$GC_manual_subcluster_confidence)
obj_gc_refined[colnames(sub_obj)] <- as.character(sub_obj$GC_manual_refined_subtype)
obj_gc_refined_zh[colnames(sub_obj)] <- as.character(sub_obj$GC_manual_refined_subtype_zh)
obj_gc_refined_confidence[colnames(sub_obj)] <- as.character(sub_obj$GC_manual_refined_confidence)
obj$GC_manual_subcluster_id <- obj_gc_subcluster
obj$GC_manual_subcluster_label <- obj_gc_label
obj$GC_manual_subcluster_candidate <- obj_gc_candidate
obj$GC_manual_subcluster_confidence <- obj_gc_confidence
obj$GC_manual_refined_subtype <- obj_gc_refined
obj$GC_manual_refined_subtype_zh <- obj_gc_refined_zh
obj$GC_manual_refined_confidence <- obj_gc_refined_confidence

saveRDS(sub_obj, file.path(checkpoint_dir, "06_GC_manual_subcluster_after_annotation.rds"))
saveRDS(obj, file.path(checkpoint_dir, "06_after_manual_GC_subcluster_annotation.rds"))

plot_if_umap(sub_obj, "GC_manual_subcluster_id", "Manual GC subclusters", "GC_manual_subcluster_umap.png")
plot_if_umap(sub_obj, "GC_manual_subcluster_label", "Manual GC subcluster labels", "GC_manual_subcluster_label_umap.png", width = 9, height = 6)
plot_if_umap(sub_obj, "GC_manual_refined_subtype", "Manual GC refined subtypes", "GC_manual_refined_subtype_umap.png", width = 10, height = 6)
plot_if_umap(sub_obj, "parent_panorama_cluster", "Manual GC cells by panorama cluster", "GC_manual_parent_panorama_cluster_umap.png")

panel_features <- ann$panel$feature[ann$panel$matched]
panel_features <- unique(panel_features[nzchar(panel_features)])
if (length(panel_features) > 0) {
  p <- DotPlot(sub_obj, features = panel_features, group.by = "GC_manual_subcluster_id") +
    RotatedAxis() +
    ggtitle("GC subtype panel features")
  ggsave(file.path(figure_dir, "GC_manual_subcluster_panel_dotplot.png"), p, width = max(7, length(panel_features) * 0.35), height = 4.8, dpi = 180, limitsize = FALSE)
}

report_lines <- c(
  "# Manual GC Subcluster",
  "",
  paste0("- project_root: ", project_root),
  paste0("- input_checkpoint: ", input_checkpoint),
  paste0("- parent GC compartment: ", gc_compartment),
  paste0("- GC cells: ", length(gc_cells)),
  paste0("- resolution: ", resolution),
  paste0("- dimensions: ", dims_max),
  paste0("- subclusters: ", nrow(cluster_counts)),
  paste0("- matched GC panel genes: ", sum(ann$panel$matched), " / ", nrow(ann$panel)),
  "",
  "## Cluster Counts",
  "",
  paste(capture.output(print(cluster_counts, row.names = FALSE)), collapse = "\n"),
  "",
  "## Cluster Summary",
  "",
  paste(capture.output(print(cluster_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Outputs",
  "",
  paste0("- ", file.path(annotation_dir, "annotation_table.tsv")),
  paste0("- ", file.path(annotation_dir, "cluster_summary.tsv")),
  paste0("- ", file.path(annotation_dir, "cell_assignments.tsv")),
  paste0("- ", file.path(annotation_dir, "parent_panorama_cluster_composition.tsv")),
  paste0("- ", file.path(marker_dir, "cluster_markers.tsv")),
  paste0("- ", file.path(marker_dir, "top50_cluster_markers.tsv")),
  paste0("- ", file.path(checkpoint_dir, "06_GC_manual_subcluster_after_annotation.rds")),
  paste0("- ", file.path(checkpoint_dir, "06_after_manual_GC_subcluster_annotation.rds")),
  paste0("- ", file.path(figure_dir, "GC_manual_subcluster_umap.png")),
  paste0("- ", file.path(figure_dir, "GC_manual_subcluster_label_umap.png")),
  paste0("- ", file.path(figure_dir, "GC_manual_refined_subtype_umap.png")),
  paste0("- ", file.path(figure_dir, "GC_manual_parent_panorama_cluster_umap.png")),
  paste0("- ", file.path(figure_dir, "GC_manual_subcluster_panel_dotplot.png"))
)
writeLines(report_lines, file.path(report_dir, "report.md"))

log_msg("Done. GC manual subcluster table: ", file.path(annotation_dir, "annotation_table.tsv"))
log_msg("Done. Parent checkpoint: ", file.path(checkpoint_dir, "06_after_manual_GC_subcluster_annotation.rds"))
