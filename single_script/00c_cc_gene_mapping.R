#!/usr/bin/env Rscript

load_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      sprintf("缺少 R 包: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

load_required_packages(c("dplyr", "Seurat", "jsonlite"))

output_dir <- "/home/user_test/syf_f5/05projects/scrna_improve/ortholog_cache"
manifest_path <- file.path(output_dir, "_manifest.json")
module_name <- "00c_cc_gene_mapping"
module_contract <- "00_ortholog_module"
module_version <- "1.0"

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

normalize_key <- function(x) {
  toupper(trimws(as.character(x)))
}

orthology_bucket <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    grepl("one2one$", x) ~ "one2one",
    grepl("one2many$", x) ~ "one2many",
    grepl("many2many$", x) ~ "many2many",
    TRUE ~ "other"
  )
}

read_manifest_local <- function(manifest_path) {
  if (!file.exists(manifest_path)) {
    stop(sprintf("缺少 manifest: %s", manifest_path), call. = FALSE)
  }
  jsonlite::read_json(manifest_path, simplifyVector = FALSE)
}

resolve_output_local <- function(manifest, key) {
  entry <- manifest$outputs[[key]]
  if (is.null(entry) || is.null(entry$path)) {
    stop(sprintf("manifest 缺少输出键: %s", key), call. = FALSE)
  }
  raw_path <- as.character(entry$path)
  if (grepl("^/", raw_path)) {
    normalizePath(raw_path, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(as.character(manifest$base_dir), raw_path), winslash = "/", mustWork = FALSE)
  }
}

write_manifest_local <- function(manifest_path, new_outputs, module_name = NULL, base_dir = NULL, inputs = NULL) {
  if (file.exists(manifest_path)) {
    existing <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
    existing$outputs <- modifyList(existing$outputs %||% list(), new_outputs)
    existing$timestamp <- timestamp_now()
    manifest <- existing
  } else {
    if (is.null(module_name) || is.null(base_dir) || is.null(inputs)) {
      stop("manifest 不存在且缺少初始化参数", call. = FALSE)
    }
    manifest <- list(
      module = module_name,
      version = module_version,
      timestamp = timestamp_now(),
      base_dir = base_dir,
      inputs = inputs,
      outputs = new_outputs
    )
  }
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

relative_path_local <- function(path, base_dir) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  normalized_base <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
  prefix <- paste0(normalized_base, "/")
  if (startsWith(normalized_path, prefix)) {
    substring(normalized_path, nchar(prefix) + 1L)
  } else {
    normalized_path
  }
}

build_output_entry <- function(path, type, produced_by, row_semantics) {
  list(
    path = relative_path_local(path, output_dir),
    type = type,
    produced_by = produced_by,
    row_semantics = row_semantics
  )
}

write_markdown_local <- function(lines, path) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

read_csv_required <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("缺少输入文件: %s", path), call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

fmt_int <- function(x) {
  if (is.na(x)) {
    return("NA")
  }
  prettyNum(round(x), big.mark = ",", scientific = FALSE)
}

fmt_pct <- function(x) {
  if (is.na(x)) {
    return("NA")
  }
  sprintf("%.1f%%", 100 * x)
}

ensure_dir(output_dir)

manifest <- read_manifest_local(manifest_path)
human_all_path <- resolve_output_local(manifest, "human_all")
human_all <- read_csv_required(human_all_path)

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

cc_rds_path <- file.path(output_dir, "chicken_cc_genes.rds")
summary_path <- file.path(output_dir, "chicken_cc_genes_summary.csv")
report_path <- file.path(output_dir, "chicken_cc_genes_report.md")

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
  manifest_path = manifest_path,
  new_outputs = list(
    cc_genes = build_output_entry(cc_rds_path, "rds", module_name, "RDS list with s.genes and g2m.genes"),
    cc_genes_summary = build_output_entry(summary_path, "csv", module_name, "cell-cycle mapping summary by phase"),
    cc_genes_report = build_output_entry(report_path, "md", module_name, "cell-cycle mapping markdown report")
  ),
  module_name = module_contract,
  base_dir = output_dir,
  inputs = list()
)

message("00c 完成。输出目录: ", output_dir)
