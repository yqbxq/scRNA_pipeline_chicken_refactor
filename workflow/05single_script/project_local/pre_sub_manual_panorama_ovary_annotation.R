suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

env_or <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

project_root <- env_or("PROJECT_ROOT", "/home/user_test/syf_f5/pre-sub/F6_1_panorama_project_panel_v3")
annotation_set <- "Panorama ovary annotation"
annotation_set_zh <- "卵巢全景注释结果"
cluster_col <- env_or("CLUSTER_COL", "seurat_clusters")

checkpoint_dir <- file.path(project_root, "results", "checkpoints")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures", "panorama_manual_annotation")
report_dir <- file.path(project_root, "reports", "eda", "panorama_manual_annotation")
config_dir <- file.path(project_root, "config")
manual_dir <- file.path(table_dir, "panorama_manual_annotation")
annotation_layer_dir <- file.path(table_dir, "annotation", "layers", "panorama")

dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annotation_layer_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

pick_checkpoint <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    stop(
      "No input checkpoint found. Tried: ",
      paste(paths, collapse = ", "),
      call. = FALSE
    )
  }
  existing[[1]]
}

manual_map <- tibble::tribble(
  ~compartment_order, ~compartment, ~compartment_zh, ~cluster, ~manual_annotation, ~manual_annotation_zh, ~expected_cells, ~expected_percent,
  1L, "Vascular / perivascular compartment", "血管 / 周血管相关细胞", "0", "Endothelial_APLN_APLNR_like", "APLN/APLNR 相关内皮样细胞", 1930L, 15.77,
  1L, "Vascular / perivascular compartment", "血管 / 周血管相关细胞", "1", "Perivascular_smooth_muscle_like", "周血管平滑肌样细胞", 1698L, 13.88,
  1L, "Vascular / perivascular compartment", "血管 / 周血管相关细胞", "8", "proliferating_vascular_like", "增殖型血管相关细胞", 536L, 4.38,
  1L, "Vascular / perivascular compartment", "血管 / 周血管相关细胞", "10", "Perivascular_smooth_muscle_like / proliferating smooth muscle-like", "增殖型周血管平滑肌样细胞", 491L, 4.01,
  1L, "Vascular / perivascular compartment", "血管 / 周血管相关细胞", "14", "Endothelial", "内皮细胞", 222L, 1.81,
  2L, "GC-like / follicular somatic compartment", "颗粒细胞样 / 卵泡体细胞相关细胞", "2", "GC_like / weak_GC_like", "颗粒细胞样 / 弱颗粒细胞样细胞", 1333L, 10.89,
  2L, "GC-like / follicular somatic compartment", "颗粒细胞样 / 卵泡体细胞相关细胞", "4", "GC_like_candidate / follicular_somatic_like_candidate", "颗粒细胞样候选 / 卵泡体细胞样候选细胞", 708L, 5.79,
  2L, "GC-like / follicular somatic compartment", "颗粒细胞样 / 卵泡体细胞相关细胞", "7", "GC_like", "颗粒细胞样细胞", 578L, 4.72,
  2L, "GC-like / follicular somatic compartment", "颗粒细胞样 / 卵泡体细胞相关细胞", "13", "proliferating_GC_like", "增殖型颗粒细胞样细胞", 255L, 2.08,
  3L, "TC / stromal / steroidogenic compartment", "膜层 / 基质 / 类固醇合成相关细胞", "5", "stromal_theca_externa_like", "基质-卵泡外膜层样细胞", 693L, 5.66,
  3L, "TC / stromal / steroidogenic compartment", "膜层 / 基质 / 类固醇合成相关细胞", "6", "steroidogenic_or_mixed_candidate", "类固醇合成或混合候选细胞", 664L, 5.43,
  3L, "TC / stromal / steroidogenic compartment", "膜层 / 基质 / 类固醇合成相关细胞", "11", "TC_broad", "膜层细胞大类", 416L, 3.40,
  3L, "TC / stromal / steroidogenic compartment", "膜层 / 基质 / 类固醇合成相关细胞", "16", "Stromal_fibroblast_like", "基质成纤维样细胞", 191L, 1.56,
  3L, "TC / stromal / steroidogenic compartment", "膜层 / 基质 / 类固醇合成相关细胞", "17", "steroidogenic_TC_like / theca_interna_like", "类固醇合成型膜层细胞样 / 卵泡内膜层样细胞", 123L, 1.01,
  4L, "Blood-derived cells", "血液来源细胞", "3", "Erythroid", "红系细胞", 1174L, 9.59,
  4L, "Blood-derived cells", "血液来源细胞", "12", "Thrombocyte_platelet_like", "血栓细胞 / 血小板样细胞", 391L, 3.20,
  5L, "Immune compartment", "免疫细胞", "9", "T_NK_like", "T / NK 样细胞", 495L, 4.05,
  5L, "Immune compartment", "免疫细胞", "15", "Myeloid_macrophage_like", "髓系 / 巨噬细胞样细胞", 215L, 1.76,
  5L, "Immune compartment", "免疫细胞", "18", "T_like / Immune_T_like", "T 细胞样 / 免疫 T 样细胞", 62L, 0.51,
  5L, "Immune compartment", "免疫细胞", "19", "B_cell_plasma_like", "B 细胞 / 浆细胞样细胞", 62L, 0.51
) %>%
  mutate(
    annotation_set = annotation_set,
    annotation_set_zh = annotation_set_zh,
    display_label = paste0(manual_annotation, " (", manual_annotation_zh, ")"),
    manual_confidence = "manual_curated",
    evidence_source = "user_manual_revision_20260706",
    note = "Manual revision supplied by user for F6_1 panorama panel v3."
  )

input_checkpoint <- pick_checkpoint(c(
  file.path(checkpoint_dir, "04_after_subcluster_annotation.rds"),
  file.path(checkpoint_dir, "03_after_annotation.rds"),
  file.path(checkpoint_dir, "02_after_clustering.rds")
))

log_msg("Project root: ", project_root)
log_msg("Input checkpoint: ", input_checkpoint)
obj <- readRDS(input_checkpoint)

if (!cluster_col %in% colnames(obj@meta.data)) {
  stop("Missing cluster column in object metadata: ", cluster_col, call. = FALSE)
}

obj_clusters <- as.character(obj@meta.data[[cluster_col]])
actual_counts <- as.data.frame(table(obj_clusters), stringsAsFactors = FALSE)
colnames(actual_counts) <- c("cluster", "actual_cells")
actual_counts$actual_cells <- as.integer(actual_counts$actual_cells)
total_cells <- ncol(obj)

manual_table <- manual_map %>%
  left_join(actual_counts, by = "cluster") %>%
  mutate(
    actual_cells = ifelse(is.na(actual_cells), 0L, actual_cells),
    actual_percent = round(100 * actual_cells / total_cells, 4),
    count_status = case_when(
      actual_cells == expected_cells ~ "match_expected",
      actual_cells == 0L ~ "missing_in_object",
      TRUE ~ "differs_from_expected"
    ),
    cell_type = manual_annotation,
    annotation_label = manual_annotation,
    confidence = manual_confidence
  ) %>%
  arrange(compartment_order, as.integer(cluster))

unknown_clusters <- setdiff(sort(unique(obj_clusters)), manual_table$cluster)
if (length(unknown_clusters) > 0) {
  extra <- actual_counts %>%
    filter(cluster %in% unknown_clusters) %>%
    mutate(
      compartment_order = 99L,
      compartment = "Unassigned",
      compartment_zh = "未分配",
      manual_annotation = "Unassigned",
      manual_annotation_zh = "未分配",
      expected_cells = NA_integer_,
      expected_percent = NA_real_,
      annotation_set = annotation_set,
      annotation_set_zh = annotation_set_zh,
      display_label = "Unassigned (未分配)",
      manual_confidence = "unassigned",
      evidence_source = "object_cluster_not_in_user_manual_revision",
      note = "Cluster exists in object but was not included in the supplied manual annotation.",
      actual_percent = round(100 * actual_cells / total_cells, 4),
      count_status = "extra_cluster_in_object",
      cell_type = manual_annotation,
      annotation_label = manual_annotation,
      confidence = manual_confidence
    )
  manual_table <- bind_rows(manual_table, extra) %>%
    arrange(compartment_order, suppressWarnings(as.integer(cluster)))
}

map_label <- setNames(manual_table$manual_annotation, manual_table$cluster)
map_label_zh <- setNames(manual_table$manual_annotation_zh, manual_table$cluster)
map_compartment <- setNames(manual_table$compartment, manual_table$cluster)
map_compartment_zh <- setNames(manual_table$compartment_zh, manual_table$cluster)
map_display <- setNames(manual_table$display_label, manual_table$cluster)

obj$panorama_ovary_annotation <- unname(map_label[obj_clusters])
obj$panorama_ovary_annotation_zh <- unname(map_label_zh[obj_clusters])
obj$panorama_ovary_compartment <- unname(map_compartment[obj_clusters])
obj$panorama_ovary_compartment_zh <- unname(map_compartment_zh[obj_clusters])
obj$panorama_ovary_display_label <- unname(map_display[obj_clusters])
obj$panorama_annotation_set <- annotation_set
obj$panorama_annotation_set_zh <- annotation_set_zh
obj$cell_type_manual <- obj$panorama_ovary_annotation
obj$cell_type_manual_zh <- obj$panorama_ovary_annotation_zh

if ("cell_type" %in% colnames(obj@meta.data)) {
  obj$cell_type_auto_before_manual <- obj$cell_type
}
if ("annotation_label" %in% colnames(obj@meta.data)) {
  obj$annotation_label_auto_before_manual <- obj$annotation_label
}
if ("panorama_annotation" %in% colnames(obj@meta.data)) {
  obj$panorama_annotation_auto_before_manual <- obj$panorama_annotation
}

obj$cell_type <- obj$panorama_ovary_annotation
obj$cell_type_zh <- obj$panorama_ovary_annotation_zh
obj$annotation_label <- obj$panorama_ovary_annotation
obj$annotation_label_zh <- obj$panorama_ovary_annotation_zh
obj$panorama_annotation <- obj$panorama_ovary_annotation
obj$panorama_annotation_zh <- obj$panorama_ovary_annotation_zh

cell_assignments <- obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  transmute(
    cell_id = cell_id,
    cluster = as.character(.data[[cluster_col]]),
    annotation_set = panorama_annotation_set,
    annotation_set_zh = panorama_annotation_set_zh,
    panorama_ovary_compartment = panorama_ovary_compartment,
    panorama_ovary_compartment_zh = panorama_ovary_compartment_zh,
    panorama_ovary_annotation = panorama_ovary_annotation,
    panorama_ovary_annotation_zh = panorama_ovary_annotation_zh,
    panorama_ovary_display_label = panorama_ovary_display_label
  )

compartment_summary <- cell_assignments %>%
  count(
    annotation_set,
    annotation_set_zh,
    panorama_ovary_compartment,
    panorama_ovary_compartment_zh,
    name = "n_cells"
  ) %>%
  mutate(percent_cells = round(100 * n_cells / sum(n_cells), 4)) %>%
  arrange(desc(n_cells))

annotation_summary <- cell_assignments %>%
  count(
    annotation_set,
    annotation_set_zh,
    panorama_ovary_compartment,
    panorama_ovary_compartment_zh,
    panorama_ovary_annotation,
    panorama_ovary_annotation_zh,
    name = "n_cells"
  ) %>%
  mutate(percent_cells = round(100 * n_cells / sum(n_cells), 4)) %>%
  arrange(desc(n_cells))

write_tsv(manual_map, file.path(config_dir, "panorama_ovary_manual_annotation_map.tsv"))
write_tsv(manual_table, file.path(annotation_layer_dir, "manual_annotation_table.tsv"))
write_tsv(manual_table, file.path(manual_dir, "panorama_ovary_manual_annotation_table.tsv"))
write_tsv(cell_assignments, file.path(manual_dir, "panorama_ovary_manual_cell_assignments.tsv"))
write_tsv(compartment_summary, file.path(manual_dir, "panorama_ovary_manual_compartment_summary.tsv"))
write_tsv(annotation_summary, file.path(manual_dir, "panorama_ovary_manual_annotation_summary.tsv"))

output_checkpoint <- file.path(checkpoint_dir, "05_after_manual_panorama_ovary_annotation.rds")
saveRDS(obj, output_checkpoint)

plot_umap <- function(group_col, title, file_name, width = 9, height = 6) {
  if (!"umap" %in% names(obj@reductions)) return(invisible(FALSE))
  p <- DimPlot(obj, reduction = "umap", group.by = group_col, label = TRUE, repel = TRUE) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(figure_dir, file_name), p, width = width, height = height, dpi = 180, limitsize = FALSE)
  invisible(TRUE)
}

plot_umap("panorama_ovary_annotation", annotation_set, "panorama_ovary_manual_annotation_umap.png", width = 11, height = 7)
plot_umap("panorama_ovary_compartment", annotation_set_zh, "panorama_ovary_manual_compartment_umap.png", width = 9, height = 6)

count_status_summary <- manual_table %>%
  count(count_status, name = "n_clusters") %>%
  arrange(count_status)

report_lines <- c(
  "# Panorama ovary annotation",
  "",
  "卵巢全景注释结果",
  "",
  paste0("- project_root: ", project_root),
  paste0("- input_checkpoint: ", input_checkpoint),
  paste0("- output_checkpoint: ", output_checkpoint),
  paste0("- cells: ", total_cells),
  paste0("- cluster_column: ", cluster_col),
  paste0("- manual_clusters: ", nrow(manual_map)),
  paste0("- object_clusters: ", length(unique(obj_clusters))),
  "",
  "## Count Status",
  "",
  paste(capture.output(print(count_status_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Compartment Summary",
  "",
  paste(capture.output(print(compartment_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Annotation Summary",
  "",
  paste(capture.output(print(annotation_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Output Files",
  "",
  paste0("- ", file.path(config_dir, "panorama_ovary_manual_annotation_map.tsv")),
  paste0("- ", file.path(annotation_layer_dir, "manual_annotation_table.tsv")),
  paste0("- ", file.path(manual_dir, "panorama_ovary_manual_annotation_table.tsv")),
  paste0("- ", file.path(manual_dir, "panorama_ovary_manual_cell_assignments.tsv")),
  paste0("- ", file.path(manual_dir, "panorama_ovary_manual_compartment_summary.tsv")),
  paste0("- ", file.path(manual_dir, "panorama_ovary_manual_annotation_summary.tsv")),
  paste0("- ", output_checkpoint),
  paste0("- ", file.path(figure_dir, "panorama_ovary_manual_annotation_umap.png")),
  paste0("- ", file.path(figure_dir, "panorama_ovary_manual_compartment_umap.png"))
)

writeLines(report_lines, file.path(report_dir, "report.md"))

if (any(manual_table$count_status != "match_expected")) {
  warning(
    "Some cluster counts did not match the supplied expected counts. See ",
    file.path(annotation_layer_dir, "manual_annotation_table.tsv"),
    call. = FALSE
  )
}

log_msg("Done. Manual annotation table: ", file.path(annotation_layer_dir, "manual_annotation_table.tsv"))
log_msg("Done. Manual checkpoint: ", output_checkpoint)
