#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
)

source_utf8 <- function(path) {
  source(path, encoding = "UTF-8")
}

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "clustering_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "cluster_marker_audit_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))

load_required_packages(c("dplyr", "tibble", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03d_marker_risk"
prepare_dirs_03(cfg)

panorama_spec <- panorama_layer_spec(cfg)
manifest_03c3 <- read_manifest_local(cfg$module_03c3_manifest_path)
markers_scored <- read_tsv_optional(resolve_output_local(manifest_03c3, "markers_scored_tsv"))
top20_markers <- read_tsv_optional(resolve_output_local(manifest_03c3, "top20_markers_tsv"))
cluster_qc <- read_tsv_optional(resolve_output_local(manifest_03c3, "cluster_marker_qc_tsv"))

gene_risk <- build_gene_marker_risk(markers_scored)
cluster_marker_risk <- build_cluster_marker_risk(top20_markers, gene_risk, cluster_qc)

ambient_marker_risk <- data.frame(
  gene = character(0),
  sample_id = character(0),
  source_cluster = character(0),
  leakage_before = numeric(0),
  leakage_after = numeric(0),
  off_source_reduction = numeric(0),
  ambient_marker_flag = character(0),
  stringsAsFactors = FALSE
)
if (file.exists(cfg$module_02a_manifest_path)) {
  manifest_02a <- tryCatch(read_manifest_local(cfg$module_02a_manifest_path), error = function(e) NULL)
  if (!is.null(manifest_02a) && !is.null(manifest_02a$outputs$ambient_marker_leakage)) {
    leakage_path <- tryCatch(resolve_output_local(manifest_02a, "ambient_marker_leakage"), error = function(e) "")
    leakage <- read_tsv_optional(leakage_path)
    if (nrow(leakage) > 0 && all(c("gene", "source_cluster", "leakage_before", "leakage_after") %in% colnames(leakage))) {
      leakage$leakage_before <- suppressWarnings(as.numeric(leakage$leakage_before))
      leakage$leakage_after <- suppressWarnings(as.numeric(leakage$leakage_after))
      leakage$off_source_reduction <- ifelse(
        is.finite(leakage$leakage_before) & leakage$leakage_before > 0,
        (leakage$leakage_before - leakage$leakage_after) / leakage$leakage_before,
        NA_real_
      )
      ambient_marker_risk <- leakage %>%
        dplyr::mutate(
          ambient_marker_flag = dplyr::case_when(
            is.finite(off_source_reduction) & off_source_reduction >= 0.25 ~ "ambient_leakage_reduced",
            is.finite(leakage_before) & leakage_before >= 0.25 ~ "ambient_high_leakage",
            TRUE ~ ""
          )
        ) %>%
        dplyr::select(dplyr::any_of(c(
          "gene", "sample_id", "source_cluster", "leakage_before",
          "leakage_after", "off_source_reduction", "ambient_marker_flag"
        )))
    }
  }
}

if (nrow(ambient_marker_risk) > 0 && nrow(cluster_marker_risk) > 0) {
  ambient_genes <- unique(ambient_marker_risk$gene[nzchar(ambient_marker_risk$ambient_marker_flag)])
  ambient_top20 <- top20_markers %>%
    dplyr::mutate(ambient_marker = gene %in% ambient_genes) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(ambient_top20_fraction = mean(ambient_marker, na.rm = TRUE), .groups = "drop")
  cluster_marker_risk <- dplyr::left_join(cluster_marker_risk, ambient_top20, by = "cluster")
  cluster_marker_risk$ambient_top20_fraction[is.na(cluster_marker_risk$ambient_top20_fraction)] <- 0
} else if (nrow(cluster_marker_risk) > 0) {
  cluster_marker_risk$ambient_top20_fraction <- NA_real_
}

gene_risk_tsv <- file.path(cfg$marker_risk_table_dir_layer, "gene_risk.tsv")
cluster_marker_risk_tsv <- file.path(cfg$marker_risk_table_dir_layer, "cluster_marker_risk.tsv")
ambient_marker_risk_tsv <- file.path(cfg$marker_risk_table_dir_layer, "ambient_marker_risk.tsv")
report_path <- file.path(cfg$marker_risk_report_dir_layer, "report.md")

write_tsv_local(gene_risk, gene_risk_tsv)
write_tsv_local(cluster_marker_risk, cluster_marker_risk_tsv)
write_tsv_local(ambient_marker_risk, ambient_marker_risk_tsv)

report_lines <- c(
  "# 03d Cluster Marker Risk",
  "",
  sprintf("- layer_id: `%s`", panorama_spec$layer_id),
  sprintf("- gene_risk_rows: `%s`", nrow(gene_risk)),
  sprintf("- cluster_marker_risk_rows: `%s`", nrow(cluster_marker_risk)),
  sprintf("- ambient_marker_risk_rows: `%s`", nrow(ambient_marker_risk)),
  "",
  "## Cluster Marker Risk",
  render_markdown_table_local(cluster_marker_risk[, intersect(c(
    "cluster", "n_cells", "strict_marker_n", "nuisance_top20_fraction",
    "high_background_top20_fraction", "shared_marker_fraction", "ambient_top20_fraction",
    "deg_quality_flag", "marker_quality_flag", "manual_review_reason"
  ), colnames(cluster_marker_risk)), drop = FALSE]),
  "",
  "This stage does not call any cell type and does not mark clusters as invalid. It only prepares manual-review risk evidence."
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_03d_marker_risk_manifest_path,
  new_outputs = list(
    gene_risk_tsv = build_output_entry(gene_risk_tsv, "tsv", module_name, "one row per marker gene with dataset-level marker risk flags", base_dir = cfg$project_root, schema = infer_schema_from_df(gene_risk)),
    cluster_marker_risk_tsv = build_output_entry(cluster_marker_risk_tsv, "tsv", module_name, "one row per final cluster with marker risk fractions and review reasons", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_marker_risk)),
    ambient_marker_risk_tsv = build_output_entry(ambient_marker_risk_tsv, "tsv", module_name, "ambient leakage evidence joined from module 02 when available", base_dir = cfg$project_root, schema = infer_schema_from_df(ambient_marker_risk)),
    report = build_output_entry(report_path, "md", module_name, "human-readable marker risk report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(module_03c3_manifest = cfg$module_03c3_manifest_path),
  version = cfg$module_version,
  depends_on = list(module_03c3 = cfg$module_03c3_manifest_path)
)

message("03d complete. Cluster marker risk: ", cluster_marker_risk_tsv)
