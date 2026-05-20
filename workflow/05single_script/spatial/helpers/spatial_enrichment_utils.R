st06_empty_df <- function(cols) {
  as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols), stringsAsFactors = FALSE)
}

st06_scalar <- function(x, default = "") {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]])) {
    return(default)
  }
  value <- trimws(as.character(x[[1]]))
  if (!nzchar(value)) default else value
}

st06_status_summary <- function(df) {
  if (nrow(df) == 0 || !"status" %in% colnames(df)) {
    return(st06_empty_df(c("status", "n")))
  }
  out <- as.data.frame(table(status = df$status), stringsAsFactors = FALSE)
  colnames(out) <- c("status", "n")
  out
}

st06_abs_path <- function(path, cfg) {
  path <- st06_scalar(path)
  if (!nzchar(path)) {
    return("")
  }
  if (grepl("^/", path)) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(cfg$project_root, path), winslash = "/", mustWork = FALSE)
  }
}

st06_read_tsv <- function(path) {
  if (!nzchar(path) || !file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
}

st06_write_tsv <- function(df, path) {
  ensure_dir(dirname(path))
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

st06_manifest_output <- function(manifest_path, keys, cfg) {
  if (!file.exists(manifest_path)) {
    return("")
  }
  manifest <- read_manifest_local(manifest_path)
  for (key in keys) {
    if (!is.null(manifest$outputs[[key]]) && !is.null(manifest$outputs[[key]]$path)) {
      return(resolve_output_local(manifest, key))
    }
  }
  ""
}

st06_upstream_manifest_rows <- function(cfg) {
  specs <- list(
    list(path = cfg$module_05a_region_pseudobulk_manifest_path, key = "pseudobulk_de_manifest", source = "spatial_05a_region_pseudobulk"),
    list(path = cfg$module_05_spatial_de_manifest_path, key = "spatial_de_manifest", source = "spatial_05_spatial_de"),
    list(path = cfg$module_05_region_marker_manifest_path, key = "marker_discovery_manifest", source = "spatial_05_region_marker")
  )
  rows <- list()
  for (spec in specs) {
    manifest_tsv <- st06_manifest_output(spec$path, spec$key, cfg)
    df <- st06_read_tsv(manifest_tsv)
    if (nrow(df) == 0) {
      next
    }
    df$upstream_module <- spec$source
    df$upstream_manifest_tsv <- manifest_tsv
    rows[[length(rows) + 1L]] <- df
  }
  if (length(rows) == 0) {
    return(st06_empty_df(c("upstream_module", "upstream_manifest_tsv")))
  }
  all_cols <- unique(unlist(lapply(rows, colnames), use.names = FALSE))
  rows <- lapply(rows, function(df) {
    for (col in setdiff(all_cols, colnames(df))) {
      df[[col]] <- ""
    }
    df[, all_cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

st06_pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) "" else hit[[1]]
}

st06_gene_table <- function(path, cfg, direction = "all") {
  df <- st06_read_tsv(path)
  gene_col <- st06_pick_col(df, c("gene", "Gene", "SYMBOL", "symbol", "gene_symbol", "feature", "feature_id", "gene_id"))
  if (!nzchar(gene_col)) {
    return(list(genes = character(0), status = "skipped_too_few_genes", reason = "input result table has no supported gene column"))
  }
  fc_col <- st06_pick_col(df, c("avg_log2FC", "avg_logFC", "logFC", "log2FC", "logfc", "coef", "estimate"))
  padj_col <- st06_pick_col(df, c("p_val_adj", "p_adj", "padj", "adj.P.Val", "FDR", "fdr", "qval", "qvalue"))
  genes <- trimws(as.character(df[[gene_col]]))
  keep <- nzchar(genes)
  if (nzchar(padj_col)) {
    padj <- suppressWarnings(as.numeric(df[[padj_col]]))
    keep <- keep & (is.na(padj) | padj <= 0.05)
  }
  if (nzchar(fc_col) && direction %in% c("up", "down")) {
    fc <- suppressWarnings(as.numeric(df[[fc_col]]))
    keep <- keep & is.finite(fc) & if (identical(direction, "up")) fc > 0 else fc < 0
  }
  df <- df[keep, , drop = FALSE]
  genes <- genes[keep]
  if (length(genes) == 0) {
    return(list(genes = character(0), status = "skipped_too_few_genes", reason = sprintf("0 genes after direction=%s filter", direction)))
  }
  score <- rep(0, nrow(df))
  if (nzchar(fc_col)) {
    score <- abs(suppressWarnings(as.numeric(df[[fc_col]])))
    score[!is.finite(score)] <- 0
  }
  padj <- rep(1, nrow(df))
  if (nzchar(padj_col)) {
    padj <- suppressWarnings(as.numeric(df[[padj_col]]))
    padj[!is.finite(padj)] <- 1
  }
  ordered <- unique(genes[order(-score, padj, genes)])
  ordered <- ordered[seq_len(min(length(ordered), cfg$spatial_enrich_top_n))]
  list(
    genes = ordered,
    status = if (length(ordered) >= cfg$spatial_enrich_min_genes) "ok" else "skipped_too_few_genes",
    reason = if (length(ordered) >= cfg$spatial_enrich_min_genes) "" else sprintf("%s genes after filtering; minimum is %s", length(ordered), cfg$spatial_enrich_min_genes)
  )
}

st06_load_ortholog_map <- function(cfg) {
  candidates <- c(
    st06_manifest_output(cfg$ortholog_manifest_path, c("chicken_human_best_csv", "human_best", "human_best_csv"), cfg),
    file.path(cfg$ortholog_cache_dir, "chicken_human_orthologs.csv")
  )
  path <- candidates[nzchar(candidates) & file.exists(candidates)][1]
  if (is.na(path) || !nzchar(path)) {
    out <- st06_empty_df(c("chicken_symbol", "human_symbol"))
    attr(out, "path") <- candidates[[length(candidates)]]
    return(out)
  }
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  chicken_col <- st06_pick_col(raw, c("external_gene_name", "chicken_symbol", "gene", "symbol"))
  human_col <- st06_pick_col(raw, c("target_gene_name", "human_symbol", "ortholog_symbol"))
  if (!nzchar(chicken_col) || !nzchar(human_col)) {
    out <- st06_empty_df(c("chicken_symbol", "human_symbol"))
    attr(out, "path") <- path
    return(out)
  }
  out <- unique(data.frame(
    chicken_symbol = trimws(as.character(raw[[chicken_col]])),
    human_symbol = trimws(as.character(raw[[human_col]])),
    stringsAsFactors = FALSE
  ))
  out <- out[nzchar(out$chicken_symbol) & nzchar(out$human_symbol), , drop = FALSE]
  attr(out, "path") <- path
  out
}

st06_map_to_human <- function(genes, ortholog_map) {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0 || nrow(ortholog_map) == 0) {
    return(list(genes = character(0), input_n = length(genes), mapped_n = 0L, mapping_rate = NA_real_))
  }
  idx <- match(toupper(genes), toupper(ortholog_map$chicken_symbol))
  mapped <- unique(ortholog_map$human_symbol[stats::na.omit(idx)])
  mapped <- mapped[nzchar(mapped)]
  list(genes = mapped, input_n = length(genes), mapped_n = length(mapped), mapping_rate = if (length(genes) == 0) NA_real_ else length(mapped) / length(genes))
}

st06_symbols_to_entrez <- function(symbols) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(list(ids = character(0), status = "skipped_no_packages", reason = "AnnotationDbi/org.Hs.eg.db is not installed"))
  }
  mapped <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(symbols), keytype = "SYMBOL", columns = "ENTREZID"),
    error = function(e) data.frame()
  )
  if (nrow(mapped) == 0 || !"ENTREZID" %in% colnames(mapped)) {
    return(list(ids = character(0), status = "skipped_no_entrez_map", reason = "0 human symbols mapped to ENTREZID"))
  }
  ids <- unique(trimws(as.character(mapped$ENTREZID)))
  ids <- ids[nzchar(ids)]
  list(ids = ids, status = "ok", reason = "")
}

st06_enrichment_manifest_cols <- function() {
  c(
    "region_value", "parent_region", "group_by", "layer_id", "comparison_id",
    "gene_direction", "analysis_type", "ontology", "source_species",
    "upstream_module", "upstream_results_tsv", "status", "reason",
    "input_gene_n", "mapped_symbol_n", "entrez_gene_n", "mapping_rate",
    "significant_term_n", "enrichment_tsv", "enrichment_significant_tsv"
  )
}

st06_enrichment_row <- function(row, direction, analysis_type, ontology, upstream_path, status, reason = "",
                                input_gene_n = 0L, mapped_symbol_n = 0L, entrez_gene_n = 0L,
                                mapping_rate = NA_real_, significant_term_n = 0L,
                                enrichment_tsv = "", enrichment_significant_tsv = "") {
  data.frame(
    region_value = st06_scalar(row$region_value, st06_scalar(row$group_value)),
    parent_region = st06_scalar(row$parent_region),
    group_by = st06_scalar(row$group_by),
    layer_id = st06_scalar(row$layer_id),
    comparison_id = st06_scalar(row$comparison_id),
    gene_direction = direction,
    analysis_type = analysis_type,
    ontology = ontology,
    source_species = "human_ortholog",
    upstream_module = st06_scalar(row$upstream_module),
    upstream_results_tsv = upstream_path,
    status = status,
    reason = reason,
    input_gene_n = input_gene_n,
    mapped_symbol_n = mapped_symbol_n,
    entrez_gene_n = entrez_gene_n,
    mapping_rate = mapping_rate,
    significant_term_n = significant_term_n,
    enrichment_tsv = enrichment_tsv,
    enrichment_significant_tsv = enrichment_significant_tsv,
    stringsAsFactors = FALSE
  )
}

st06_region_output_prefix <- function(cfg, row, analysis_type, ontology, direction) {
  parts <- c(st06_scalar(row$group_by, "region"), st06_scalar(row$comparison_id, "comparison"), st06_scalar(row$layer_id, "panorama"), st06_scalar(row$region_value, st06_scalar(row$group_value, "all")), ontology, direction)
  base_dir <- if (identical(analysis_type, "go")) cfg$spatial_region_go_table_dir else cfg$spatial_region_kegg_table_dir
  do.call(file.path, as.list(c(base_dir, vapply(parts, spatial_safe_id, character(1)))))
}

st06_write_enrichment_outputs <- function(result_df, prefix, qvalue_cutoff) {
  ensure_dir(prefix)
  result_path <- file.path(prefix, "enrichment.tsv")
  sig_path <- file.path(prefix, "enrichment_significant.tsv")
  if (is.null(result_df) || nrow(result_df) == 0) {
    result_df <- st06_empty_df(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"))
  }
  st06_write_tsv(result_df, result_path)
  sig <- result_df
  if ("qvalue" %in% colnames(sig)) {
    q <- suppressWarnings(as.numeric(sig$qvalue))
    sig <- sig[is.finite(q) & q <= qvalue_cutoff, , drop = FALSE]
  } else if ("p.adjust" %in% colnames(sig)) {
    p <- suppressWarnings(as.numeric(sig[["p.adjust"]]))
    sig <- sig[is.finite(p) & p <= qvalue_cutoff, , drop = FALSE]
  } else {
    sig <- sig[FALSE, , drop = FALSE]
  }
  st06_write_tsv(sig, sig_path)
  list(result = result_path, significant = sig_path, significant_n = nrow(sig))
}

st06_run_go <- function(entrez_ids, ont, cfg) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(list(result = data.frame(), status = "skipped_no_packages", reason = "clusterProfiler is not installed"))
  }
  res <- tryCatch(
    clusterProfiler::enrichGO(
      gene = unique(entrez_ids),
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(e) structure(list(error = conditionMessage(e)), class = "st06_error")
  )
  if (inherits(res, "st06_error")) {
    return(list(result = data.frame(), status = "failed_enrichment", reason = res$error))
  }
  list(result = as.data.frame(res), status = "ok", reason = "")
}

st06_run_kegg <- function(entrez_ids) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(list(result = data.frame(), status = "skipped_no_packages", reason = "clusterProfiler is not installed"))
  }
  res <- tryCatch(
    clusterProfiler::enrichKEGG(gene = unique(entrez_ids), organism = "hsa", keyType = "kegg", pvalueCutoff = 1, qvalueCutoff = 1),
    error = function(e) structure(list(error = conditionMessage(e)), class = "st06_error")
  )
  if (inherits(res, "st06_error")) {
    return(list(result = data.frame(), status = "failed_enrichment", reason = res$error))
  }
  list(result = as.data.frame(res), status = "ok", reason = "")
}

run_spatial_region_enrichment <- function(cfg, analysis_type = c("go", "kegg")) {
  analysis_type <- match.arg(analysis_type)
  prepare_dirs_spatial(cfg)
  upstream <- st06_upstream_manifest_rows(cfg)
  ortholog_map <- st06_load_ortholog_map(cfg)
  rows <- list()
  if (nrow(upstream) == 0) {
    rows[[1L]] <- data.frame(
      as.list(setNames(rep("", length(st06_enrichment_manifest_cols())), st06_enrichment_manifest_cols())),
      stringsAsFactors = FALSE
    )
    rows[[1L]]$analysis_type <- analysis_type
    rows[[1L]]$status <- "skipped_no_upstream_manifest"
    rows[[1L]]$reason <- "No ST05 marker, pseudobulk, or spot-level DE manifest was available"
  } else if (nrow(ortholog_map) == 0) {
    for (idx in seq_len(nrow(upstream))) {
      rows[[length(rows) + 1L]] <- st06_enrichment_row(upstream[idx, , drop = FALSE], "all", analysis_type, if (identical(analysis_type, "go")) "BP" else "KEGG", "", "skipped_no_ortholog_map", sprintf("ortholog map unavailable: %s", attr(ortholog_map, "path") %||% ""))
    }
  } else {
    for (idx in seq_len(nrow(upstream))) {
      item <- upstream[idx, , drop = FALSE]
      result_path <- st06_abs_path(st06_scalar(item$formal_results_tsv, st06_scalar(item$exploratory_results_tsv)), cfg)
      if (!nzchar(result_path) || !file.exists(result_path)) {
        next
      }
      for (direction in c("all", "up", "down")) {
        gene_set <- st06_gene_table(result_path, cfg, direction)
        onts <- if (identical(analysis_type, "go")) c("BP", "CC", "MF") else "KEGG"
        for (ont in onts) {
          prefix <- st06_region_output_prefix(cfg, item, analysis_type, ont, direction)
          if (!identical(gene_set$status, "ok")) {
            paths <- st06_write_enrichment_outputs(data.frame(), prefix, cfg$spatial_enrich_qvalue)
            rows[[length(rows) + 1L]] <- st06_enrichment_row(item, direction, analysis_type, ont, result_path, gene_set$status, gene_set$reason, input_gene_n = length(gene_set$genes), enrichment_tsv = paths$result, enrichment_significant_tsv = paths$significant)
            next
          }
          mapped <- st06_map_to_human(gene_set$genes, ortholog_map)
          entrez <- st06_symbols_to_entrez(mapped$genes)
          if (!identical(entrez$status, "ok") || length(entrez$ids) < cfg$spatial_enrich_min_genes) {
            status <- if (identical(entrez$status, "ok")) "skipped_no_entrez_map" else entrez$status
            reason <- if (identical(entrez$status, "ok")) sprintf("%s Entrez IDs after ortholog mapping; minimum is %s", length(entrez$ids), cfg$spatial_enrich_min_genes) else entrez$reason
            paths <- st06_write_enrichment_outputs(data.frame(), prefix, cfg$spatial_enrich_qvalue)
            rows[[length(rows) + 1L]] <- st06_enrichment_row(item, direction, analysis_type, ont, result_path, status, reason, mapped$input_n, mapped$mapped_n, length(entrez$ids), mapped$mapping_rate, enrichment_tsv = paths$result, enrichment_significant_tsv = paths$significant)
            next
          }
          run <- if (identical(analysis_type, "go")) st06_run_go(entrez$ids, ont, cfg) else st06_run_kegg(entrez$ids)
          paths <- st06_write_enrichment_outputs(run$result, prefix, cfg$spatial_enrich_qvalue)
          rows[[length(rows) + 1L]] <- st06_enrichment_row(item, direction, analysis_type, ont, result_path, run$status, run$reason, mapped$input_n, mapped$mapped_n, length(entrez$ids), mapped$mapping_rate, paths$significant_n, paths$result, paths$significant)
        }
      }
    }
  }
  if (length(rows) == 0) {
    return(st06_empty_df(st06_enrichment_manifest_cols()))
  }
  out <- do.call(rbind, rows)
  out[, st06_enrichment_manifest_cols(), drop = FALSE]
}
