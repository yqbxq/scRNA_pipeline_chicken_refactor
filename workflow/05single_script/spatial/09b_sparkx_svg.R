#!/usr/bin/env Rscript
.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})
PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_enrichment_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_svg_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_09b_sparkx_svg"
prepare_dirs_spatial(cfg)
out_dir <- file.path(cfg$spatial_table_dir, "09_svg", "sparkx")
manifest_path <- file.path(cfg$manifest_dir, module_name, "_manifest.json")
ensure_dir(out_dir)

sparkx_allow_proxy <- tolower(Sys.getenv("SPARKX_ALLOW_PROXY", unset = "no")) %in% c("yes", "true", "1", "on")
section_id <- Sys.getenv("SPATIAL_SVG_SECTION_ID", unset = "all")
spatial_rds <- Sys.getenv("SPARKX_SPATIAL_RDS", unset = spatial_first_existing(c(cfg$spatial_panorama_niched_rds, cfg$spatial_panorama_annotated_rds, cfg$spatial_panorama_clustered_rds)))
fallback_h5ad <- Sys.getenv("SPARKX_SPATIAL_H5AD", unset = Sys.getenv("SPATIAL_SVG_H5AD", unset = ""))
svg_tsv <- file.path(out_dir, paste0(ifelse(nzchar(section_id), section_id, "all"), "_svg.tsv"))
manifest_tsv <- file.path(out_dir, "sparkx_manifest.tsv")

write_sparkx_manifest <- function(status, reason, input_mode = "", input_path = "", method = "SPARK-X") {
  st09_write_tsv(data.frame(
    section_id = section_id,
    status = status,
    reason = reason,
    method = method,
    proxy_allowed = ifelse(sparkx_allow_proxy, "yes", "no"),
    input_mode = input_mode,
    input_path = input_path,
    svg_tsv = svg_tsv,
    stringsAsFactors = FALSE
  ), manifest_tsv)
}

write_empty <- function(status, reason, method = "SPARK-X") {
  st09_write_tsv(st09_empty_method_row(section_id, "all", method, status, reason), svg_tsv)
  write_sparkx_manifest(status, reason, "", "", method)
}

extract_seurat_expression <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) return(NULL)
  obj <- readRDS(path)
  counts <- tryCatch(spatial_counts_matrix(obj), error = function(e) NULL)
  if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) return(NULL)
  meta <- tryCatch(obj@meta.data, error = function(e) data.frame(row.names = colnames(counts)))
  condition <- if ("condition" %in% colnames(meta)) st09_scalar(unique(meta$condition), "all") else "all"
  means <- Matrix::rowMeans(counts)
  vars <- vapply(seq_len(nrow(counts)), function(i) stats::var(as.numeric(counts[i, , drop = TRUE])), numeric(1))
  data.frame(
    section_id = section_id,
    condition = condition,
    method = "SPARK-X_proxy",
    gene_id = rownames(counts),
    gene_symbol = rownames(counts),
    score = vars,
    pvalue = 1,
    qvalue = 1,
    mean_expression = as.numeric(means),
    n_spots = ncol(counts),
    n_genes_tested = nrow(counts),
    status = "ok_proxy",
    reason = "SPARKX_ALLOW_PROXY enabled; ranking by Seurat expression variance proxy",
    is_ligand_candidate = "no",
    is_receptor_candidate = "no",
    is_receiver_target_candidate = "no",
    is_region_marker_candidate = "no",
    is_regulator_target_candidate = "no",
    linked_question_ids = "",
    stringsAsFactors = FALSE
  )
}

load_seurat_sparkx_input <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) return(NULL)
  obj <- readRDS(path)
  counts <- tryCatch(spatial_counts_matrix(obj), error = function(e) NULL)
  if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) return(NULL)
  meta <- tryCatch(obj@meta.data, error = function(e) data.frame(row.names = colnames(counts)))
  coord_cols <- c("x", "y")
  if (!all(coord_cols %in% colnames(meta))) coord_cols <- c("imagecol", "imagerow")
  if (!all(coord_cols %in% colnames(meta))) coord_cols <- c("spatial_1", "spatial_2")
  if (!all(coord_cols %in% colnames(meta))) return(NULL)
  common <- intersect(colnames(counts), rownames(meta))
  if (length(common) < 3L) return(NULL)
  counts <- counts[, common, drop = FALSE]
  coords <- as.matrix(meta[common, coord_cols, drop = FALSE])
  storage.mode(coords) <- "numeric"
  condition <- if ("condition" %in% colnames(meta)) st09_scalar(unique(meta[common, "condition"]), "all") else "all"
  list(counts = counts, coords = coords, condition = condition)
}

run_sparkx_formal <- function(path) {
  input <- load_seurat_sparkx_input(path)
  if (is.null(input)) return(NULL)
  sparkx_fun <- get("sparkx", envir = asNamespace("SPARK"))
  result <- tryCatch(
    sparkx_fun(input$counts, input$coords),
    error = function(e) tryCatch(sparkx_fun(count_in = input$counts, locus = input$coords), error = function(e2) e2)
  )
  if (inherits(result, "error")) stop(conditionMessage(result), call. = FALSE)
  df <- as.data.frame(result, stringsAsFactors = FALSE)
  gene_col <- st09_pick_col(df, c("gene_id", "gene", "genes", "Gene", "gene_symbol"))
  p_col <- st09_pick_col(df, c("adjustedPval", "pvalue", "p.value", "pval", "combinedPval"))
  q_col <- st09_pick_col(df, c("adjustedPval", "qvalue", "qval", "padj", "fdr"))
  score_col <- st09_pick_col(df, c("score", "statistic", "combinedPval", p_col))
  if (!nzchar(gene_col)) {
    df$gene_id <- rownames(df)
    gene_col <- "gene_id"
  }
  means <- Matrix::rowMeans(input$counts)
  out <- data.frame(
    section_id = section_id,
    condition = input$condition,
    method = "SPARK-X",
    gene_id = as.character(df[[gene_col]]),
    gene_symbol = as.character(df[[gene_col]]),
    score = if (nzchar(score_col)) suppressWarnings(as.numeric(df[[score_col]])) else NA_real_,
    pvalue = if (nzchar(p_col)) suppressWarnings(as.numeric(df[[p_col]])) else NA_real_,
    qvalue = if (nzchar(q_col)) suppressWarnings(as.numeric(df[[q_col]])) else NA_real_,
    mean_expression = as.numeric(means[as.character(df[[gene_col]])]),
    n_spots = ncol(input$counts),
    n_genes_tested = nrow(input$counts),
    status = "ok",
    reason = "",
    is_ligand_candidate = "no",
    is_receptor_candidate = "no",
    is_receiver_target_candidate = "no",
    is_region_marker_candidate = "no",
    is_regulator_target_candidate = "no",
    linked_question_ids = "",
    stringsAsFactors = FALSE
  )
  sort_col <- if (any(is.finite(out$qvalue))) "qvalue" else if (any(is.finite(out$pvalue))) "pvalue" else "score"
  out <- out[order(out[[sort_col]], out$gene_symbol, decreasing = sort_col == "score"), , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  st09_normalize_method_df(out, "SPARK-X")
}

extract_h5ad_expression <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  py <- cfg$py_spatial_bin
  adapter <- tempfile(fileext = ".py")
  out <- tempfile(fileext = ".tsv")
  writeLines(c(
    "import sys, anndata as ad, numpy as np, pandas as pd",
    "adata = ad.read_h5ad(sys.argv[1])",
    "x = adata.X",
    "means = np.asarray(x.mean(axis=0)).reshape(-1)",
    "vars = np.asarray(x.var(axis=0)).reshape(-1) if hasattr(x, 'var') else np.zeros_like(means)",
    "gene_symbol = adata.var['gene_symbol'].astype(str).to_numpy() if 'gene_symbol' in adata.var else adata.var_names.astype(str)",
    "gene_id = adata.var['gene_id'].astype(str).to_numpy() if 'gene_id' in adata.var else adata.var_names.astype(str)",
    "condition = 'all'",
    "if 'condition' in adata.obs:",
    "    vals = sorted(set(adata.obs['condition'].dropna().astype(str)))",
    "    condition = vals[0] if len(vals) == 1 else ('mixed' if vals else 'unknown')",
    "pd.DataFrame({'condition': condition, 'gene_id': gene_id, 'gene_symbol': gene_symbol, 'score': vars, 'mean_expression': means, 'n_spots': adata.n_obs, 'n_genes_tested': adata.n_vars}).to_csv(sys.argv[2], sep='\\t', index=False)"
  ), adapter)
  run <- tryCatch(system2(py, args = c(adapter, path, out), stdout = TRUE, stderr = TRUE), error = function(e) structure(conditionMessage(e), status = 127))
  if (!file.exists(out)) return(NULL)
  df <- st09_read_tsv(out)
  if (nrow(df) == 0) return(NULL)
  data.frame(
    section_id = section_id,
    condition = df$condition,
    method = "SPARK-X_proxy",
    gene_id = df$gene_id,
    gene_symbol = df$gene_symbol,
    score = suppressWarnings(as.numeric(df$score)),
    pvalue = 1,
    qvalue = 1,
    mean_expression = suppressWarnings(as.numeric(df$mean_expression)),
    n_spots = suppressWarnings(as.integer(df$n_spots)),
    n_genes_tested = suppressWarnings(as.integer(df$n_genes_tested)),
    status = "ok_proxy",
    reason = "SPARKX_ALLOW_PROXY enabled; ranking by H5AD adapter expression variance proxy",
    is_ligand_candidate = "no",
    is_receptor_candidate = "no",
    is_receiver_target_candidate = "no",
    is_region_marker_candidate = "no",
    is_regulator_target_candidate = "no",
    linked_question_ids = "",
    stringsAsFactors = FALSE
  )
}

if (!requireNamespace("SPARK", quietly = TRUE)) {
  if (!sparkx_allow_proxy) {
    write_empty("skipped_no_sparkx", "SPARK package unavailable; set SPARKX_ALLOW_PROXY=yes for exploratory proxy.")
  } else {
    proxy <- extract_seurat_expression(spatial_rds)
    input_mode <- "spatial_rds"
    input_path <- spatial_rds
    if (is.null(proxy)) {
      proxy <- extract_h5ad_expression(fallback_h5ad)
      input_mode <- "spatial_h5ad_adapter"
      input_path <- fallback_h5ad
    }
    if (is.null(proxy)) {
      write_empty("skipped_no_input", "No readable Seurat/RDS or H5AD adapter input for SPARK-X proxy.", "SPARK-X_proxy")
    } else {
      proxy <- proxy[order(-proxy$score, -proxy$mean_expression, proxy$gene_symbol), , drop = FALSE]
      proxy$rank <- seq_len(nrow(proxy))
      proxy <- st09_normalize_method_df(proxy, "SPARK-X_proxy")
      st09_write_tsv(proxy, svg_tsv)
      write_sparkx_manifest("ok_proxy", "SPARK unavailable; proxy enabled.", input_mode, input_path, "SPARK-X_proxy")
    }
  }
} else {
  if (!file.exists(spatial_rds)) {
    write_empty("skipped_no_input", "SPARK-X formal run requires spatial Seurat/RDS input.")
  } else {
    formal <- tryCatch(run_sparkx_formal(spatial_rds), error = function(e) {
      write_empty("failed_sparkx_adapter", conditionMessage(e))
      NULL
    })
    if (!is.null(formal)) {
      st09_write_tsv(formal, svg_tsv)
      write_sparkx_manifest("ok", "", "spatial_rds", spatial_rds, "SPARK-X")
    }
  }
}

df <- st09_read_tsv(manifest_tsv)
st06_write_manifest_local(
  manifest_path,
  list(sparkx_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "SPARK-X SVG manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(df))),
  module_name,
  cfg$project_root,
  inputs = list(spatial_rds = spatial_rds, spatial_h5ad = fallback_h5ad),
  version = cfg$module_07_version
)
message(sprintf("%s complete: %s", module_name, st09_scalar(df$status, "unknown")))
