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
source_utf8(file.path(.script_dir, "helpers", "ortholog_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))

load_required_packages(c("dplyr", "Seurat", "jsonlite"))

cfg <- get_single_script_config()
module_name <- "00c_cc_gene_mapping"

read_csv_required <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("缺少输入文件: %s", path), call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

ensure_dir(cfg$output_dir)

manifest <- read_manifest_local(cfg$manifest_path)
human_all_path <- resolve_output_local(manifest, "human_all")
human_all <- read_csv_required(human_all_path)

# 00a writes a chicken-centric best table (one best target per chicken gene).
# 00c needs the inverse human->chicken direction, so it reselects one2one rows by human symbol.
human_one2one <- human_all %>%
  dplyr::mutate(type_bucket = orthology_bucket(orthology_type)) %>%
  dplyr::filter(
    type_bucket == "one2one",
    !is.na(external_gene_name),
    external_gene_name != "",
    !is.na(target_gene_name),
    target_gene_name != ""
  ) %>%
  dplyr::mutate(human_symbol_upper = normalize_key(target_gene_name)) %>%
  dplyr::arrange(
    dplyr::desc(dplyr::coalesce(orthology_confidence, -1)),
    dplyr::desc(dplyr::coalesce(perc_id, -1)),
    dplyr::desc(dplyr::coalesce(perc_id_r1, -1)),
    external_gene_name
  ) %>%
  dplyr::group_by(human_symbol_upper) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

h2c <- setNames(human_one2one$external_gene_name, human_one2one$human_symbol_upper)
cc_genes <- Seurat::cc.genes.updated.2019

map_phase_genes <- function(human_genes, h2c) {
  mapped <- unname(h2c[normalize_key(human_genes)])
  unique(mapped[!is.na(mapped) & nzchar(mapped)])
}

chicken_s <- map_phase_genes(cc_genes$s.genes, h2c)
chicken_g2m <- map_phase_genes(cc_genes$g2m.genes, h2c)

unmapped_s <- cc_genes$s.genes[is.na(h2c[normalize_key(cc_genes$s.genes)])]
unmapped_g2m <- cc_genes$g2m.genes[is.na(h2c[normalize_key(cc_genes$g2m.genes)])]

cc_rds_path <- file.path(cfg$output_dir, "chicken_cc_genes.rds")
summary_path <- file.path(cfg$output_dir, "chicken_cc_genes_summary.csv")
report_path <- file.path(cfg$output_dir, "chicken_cc_genes_report.md")

saveRDS(
  list(
    s.genes = chicken_s,
    g2m.genes = chicken_g2m
  ),
  cc_rds_path
)

summary_df <- data.frame(
  phase = c("S", "G2M"),
  total_human_genes = c(length(cc_genes$s.genes), length(cc_genes$g2m.genes)),
  mapped_chicken_genes = c(length(chicken_s), length(chicken_g2m)),
  mapped_fraction = c(length(chicken_s) / length(cc_genes$s.genes), length(chicken_g2m) / length(cc_genes$g2m.genes)),
  stringsAsFactors = FALSE
)
write.csv(summary_df, summary_path, row.names = FALSE)

report_lines <- c(
  "# Chicken Cell-Cycle Gene Mapping Report",
  "",
  sprintf("Generated: %s", timestamp_now()),
  "",
  "## Summary",
  "",
  sprintf("- S phase: %s / %s mapped (%s)", fmt_int(length(chicken_s)), fmt_int(length(cc_genes$s.genes)), fmt_pct(length(chicken_s) / length(cc_genes$s.genes))),
  sprintf("- G2M phase: %s / %s mapped (%s)", fmt_int(length(chicken_g2m)), fmt_int(length(cc_genes$g2m.genes)), fmt_pct(length(chicken_g2m) / length(cc_genes$g2m.genes))),
  "",
  "## Unmapped Human S Genes",
  "",
  if (length(unmapped_s) > 0) paste(unmapped_s, collapse = ", ") else "None",
  "",
  "## Unmapped Human G2M Genes",
  "",
  if (length(unmapped_g2m) > 0) paste(unmapped_g2m, collapse = ", ") else "None"
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$manifest_path,
  new_outputs = list(
    cc_genes = build_output_entry(cc_rds_path, "rds", module_name, "RDS list with s.genes and g2m.genes", base_dir = cfg$output_dir),
    cc_genes_summary = build_output_entry(summary_path, "csv", module_name, "cell-cycle mapping summary by phase", base_dir = cfg$output_dir, schema = infer_schema_from_df(summary_df)),
    cc_genes_report = build_output_entry(report_path, "md", module_name, "cell-cycle mapping markdown report", base_dir = cfg$output_dir)
  ),
  module_name = cfg$module_contract,
  base_dir = cfg$output_dir,
  inputs = list(),
  version = cfg$module_version,
  depends_on = list()
)

message("00c 完成。输出目录: ", cfg$output_dir)
