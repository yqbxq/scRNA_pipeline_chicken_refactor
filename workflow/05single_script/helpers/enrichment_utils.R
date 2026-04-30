empty_df_06 <- function(cols = character()) {
  if (exists("empty_df_05", mode = "function")) {
    return(empty_df_05(cols))
  }
  as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols), stringsAsFactors = FALSE)
}

empty_enrichment_manifest_06 <- function() {
  empty_df_06(c(
    "target_id", "layer_id", "comparison_id", "cluster_id", "gene_direction",
    "analysis_type", "ontology", "source_species", "deg_source",
    "deg_inference_status", "formal_status", "result_level", "biological_replicates",
    "background_tsv", "background_gene_n", "background_mapped_n",
    "status", "reason", "input_gene_n",
    "mapped_gene_n", "mapping_rate", "significant_term_n",
    "enrichment_tsv", "enrichment_significant_tsv", "dotplot_png", "barplot_png"
  ))
}

empty_enrichment_result_06 <- function(extra_cols = character()) {
  standard_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "p_adjust_method", "geneID", "Count")
  empty_df_06(unique(c(extra_cols, standard_cols)))
}

normalize_path_06 <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

build_enrichment_manifest_row_06 <- function(
    target_id = "",
    layer_id,
    comparison_id,
    cluster_id,
    gene_direction,
    analysis_type,
    ontology,
    source_species,
    deg_source,
    deg_inference_status,
    formal_status = "",
    result_level = "",
    biological_replicates = "",
    background_tsv = "",
    background_gene_n = 0L,
    background_mapped_n = 0L,
    status,
    reason = "",
    input_gene_n = 0L,
    mapped_gene_n = 0L,
    mapping_rate = NA_real_,
    significant_term_n = 0L,
    paths = NULL,
    enrichment_tsv = "",
    enrichment_significant_tsv = "",
    dotplot_png = "",
    barplot_png = "") {
  if (!is.null(paths)) {
    enrichment_tsv <- paths$enrichment_tsv %||% enrichment_tsv
    enrichment_significant_tsv <- paths$enrichment_significant_tsv %||% enrichment_significant_tsv
  }
  data.frame(
    target_id = target_id,
    layer_id = layer_id,
    comparison_id = comparison_id,
    cluster_id = cluster_id,
    gene_direction = gene_direction,
    analysis_type = analysis_type,
    ontology = ontology,
    source_species = source_species,
    deg_source = deg_source,
    deg_inference_status = deg_inference_status,
    formal_status = formal_status,
    result_level = result_level,
    biological_replicates = biological_replicates,
    background_tsv = normalize_scalar_value(background_tsv),
    background_gene_n = background_gene_n,
    background_mapped_n = background_mapped_n,
    status = status,
    reason = reason,
    input_gene_n = input_gene_n,
    mapped_gene_n = mapped_gene_n,
    mapping_rate = mapping_rate,
    significant_term_n = significant_term_n,
    enrichment_tsv = if (nzchar(enrichment_tsv)) normalize_path_06(enrichment_tsv) else "",
    enrichment_significant_tsv = if (nzchar(enrichment_significant_tsv)) normalize_path_06(enrichment_significant_tsv) else "",
    dotplot_png = normalize_scalar_value(dotplot_png),
    barplot_png = normalize_scalar_value(barplot_png),
    stringsAsFactors = FALSE
  )
}

human_strategy_enabled_06 <- function(cfg) {
  strategy <- tolower(normalize_scalar_value(cfg$enrichment_species_strategy, "chicken_primary"))
  any(grepl("human|dual|both|mapped", strategy))
}

get_org_db_06 <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Missing OrgDb package: %s", pkg), call. = FALSE)
  }
  get(pkg, envir = asNamespace(pkg))
}

clusterprofiler_available_06 <- function() {
  requireNamespace("clusterProfiler", quietly = TRUE)
}

gprofiler_available_06 <- function() {
  requireNamespace("gprofiler2", quietly = TRUE)
}

optional_org_db_06 <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(NULL)
  }
  get(pkg, envir = asNamespace(pkg))
}

resolve_manifest_output_optional_06 <- function(manifest_path, keys) {
  if (!file.exists(manifest_path)) {
    return("")
  }
  manifest <- read_manifest_local(manifest_path)
  for (key in keys) {
    entry <- manifest$outputs[[key]]
    if (!is.null(entry) && !is.null(entry$path) && nzchar(as.character(entry$path))) {
      return(resolve_output_local(manifest, key))
    }
  }
  ""
}

load_ortholog_map <- function(cfg) {
  manifest_path <- cfg$ortholog_manifest_path %||% file.path(cfg$ortholog_cache_dir, "_manifest.json")
  map_path <- resolve_manifest_output_optional_06(
    manifest_path,
    c("chicken_human_best_csv", "human_best", "human_best_csv")
  )
  if (!nzchar(map_path)) {
    map_path <- file.path(cfg$ortholog_cache_dir, "chicken_human_orthologs.csv")
  }
  if (!file.exists(map_path)) {
    out <- empty_df_06(c("chicken_symbol", "human_symbol"))
    attr(out, "stats") <- list(map_path = map_path, available = FALSE, retained_n = 0L)
    return(out)
  }

  raw <- read.csv(map_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("external_gene_name", "target_gene_name")
  if (!all(required %in% colnames(raw))) {
    out <- empty_df_06(c("chicken_symbol", "human_symbol"))
    attr(out, "stats") <- list(map_path = map_path, available = FALSE, retained_n = 0L)
    return(out)
  }

  raw$external_gene_name <- trimws(as.character(raw$external_gene_name))
  raw$target_gene_name <- trimws(as.character(raw$target_gene_name))
  raw <- raw[nzchar(raw$external_gene_name) & nzchar(raw$target_gene_name), , drop = FALSE]

  preferred <- raw
  if ("target_species" %in% colnames(preferred)) {
    preferred <- preferred[tolower(as.character(preferred$target_species)) == "human", , drop = FALSE]
  }
  if ("orthology_type" %in% colnames(preferred)) {
    preferred <- preferred[tolower(as.character(preferred$orthology_type)) %in% c("ortholog_one2one", "apparent_ortholog_one2one"), , drop = FALSE]
  }
  if ("pair_quality" %in% colnames(preferred)) {
    quality_keep <- tolower(as.character(preferred$pair_quality)) %in% c("gold", "silver")
    if (any(quality_keep)) {
      preferred <- preferred[quality_keep, , drop = FALSE]
    }
  }
  if (nrow(preferred) == 0) {
    preferred <- raw
  }

  out <- unique(data.frame(
    chicken_symbol = preferred$external_gene_name,
    human_symbol = preferred$target_gene_name,
    stringsAsFactors = FALSE
  ))
  attr(out, "stats") <- list(map_path = map_path, available = TRUE, retained_n = nrow(out), raw_n = nrow(raw))
  out
}

map_genes_to_human <- function(gene_list, ortholog_map) {
  genes <- unique(trimws(as.character(gene_list)))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0 || nrow(ortholog_map) == 0) {
    return(list(genes = character(0), input_n = length(genes), mapped_n = 0L, mapping_rate = NA_real_))
  }

  lookup <- ortholog_map
  lookup$key <- toupper(lookup$chicken_symbol)
  hit <- lookup[lookup$key %in% toupper(genes), , drop = FALSE]
  mapped <- unique(trimws(as.character(hit$human_symbol)))
  mapped <- mapped[nzchar(mapped)]
  list(
    genes = mapped,
    input_n = length(genes),
    mapped_n = length(unique(hit$key)),
    mapping_rate = safe_rate(length(unique(hit$key)), length(genes))
  )
}

read_enrichment_input_grid_06 <- function(cfg) {
  targets <- read_tsv_optional(cfg$enrichment_targets_sheet)
  registry <- read_tsv_optional(cfg$gene_program_registry_tsv %||% deg_report_paths_05(cfg)$gene_program_registry_tsv)
  target_cols <- c(
    "target_id", "source_question_id", "comparison_id", "layer_scope",
    "analysis_mode", "gene_program_role", "organism", "database",
    "min_genes", "enabled", "notes"
  )
  registry_cols <- c(
    "comparison_id", "source_question_id", "layer_id", "analysis_mode",
    "gene_program_role", "result_level", "preferred_for_downstream",
    "formal_status", "deg_tsv", "marker_tsv", "top_gene_tsv",
    "n_significant", "warning"
  )
  for (col in target_cols) {
    if (!col %in% colnames(targets)) targets[[col]] <- character(nrow(targets))
  }
  for (col in registry_cols) {
    if (!col %in% colnames(registry)) registry[[col]] <- character(nrow(registry))
  }
  if (nrow(targets) == 0) {
    return(empty_df_06(c(
      "target_id", "layer_id", "comparison_id", "database", "organism",
      "deg_tsv", "deg_source", "deg_inference_status", "formal_status",
      "result_level", "biological_replicates", "background_tsv",
      "background_gene_n", "marker_results_tsv", "pseudobulk_results_tsv",
      "registry_warning"
    )))
  }
  targets <- targets[, target_cols, drop = FALSE]
  registry <- registry[, registry_cols, drop = FALSE]
  targets$enabled <- tolower(vapply(targets$enabled, normalize_scalar_value, character(1), default = "yes"))
  targets <- targets[targets$enabled %in% c("yes", "true", "1", "on") & nzchar(targets$comparison_id), , drop = FALSE]
  if (nrow(targets) == 0) {
    return(empty_df_06(c(
      "target_id", "layer_id", "comparison_id", "database", "organism",
      "deg_tsv", "deg_source", "deg_inference_status", "formal_status",
      "result_level", "biological_replicates", "background_tsv",
      "background_gene_n", "marker_results_tsv", "pseudobulk_results_tsv",
      "registry_warning"
    )))
  }

  rows <- lapply(seq_len(nrow(targets)), function(idx) {
    target <- targets[idx, , drop = FALSE]
    comparison_id <- normalize_scalar_value(target$comparison_id[[1]])
    layer_id <- normalize_scalar_value(target$layer_scope[[1]])
    hit <- registry[registry$comparison_id == comparison_id, , drop = FALSE]
    if (nrow(hit) > 0) {
      exact <- hit[nzchar(hit$layer_id) & hit$layer_id == layer_id, , drop = FALSE]
      if (nrow(exact) > 0) hit <- exact
      hit <- hit[1, , drop = FALSE]
    }

    marker_path <- if (nrow(hit) > 0) normalize_scalar_value(hit$marker_tsv[[1]]) else ""
    deg_path <- if (nrow(hit) > 0) normalize_scalar_value(hit$deg_tsv[[1]]) else ""
    top_path <- if (nrow(hit) > 0) normalize_scalar_value(hit$top_gene_tsv[[1]]) else ""
    deg_tsv <- ""
    deg_source <- "none"
    if (nzchar(top_path) && file.exists(top_path)) {
      deg_tsv <- top_path
      deg_source <- "gene_program_registry:top_gene_tsv"
    } else if (nzchar(deg_path) && file.exists(deg_path)) {
      deg_tsv <- deg_path
      deg_source <- "gene_program_registry:deg_tsv"
    } else if (nzchar(marker_path) && file.exists(marker_path)) {
      deg_tsv <- marker_path
      deg_source <- "gene_program_registry:marker_tsv"
    }

    result_level <- if (nrow(hit) > 0) normalize_scalar_value(hit$result_level[[1]], "unavailable") else "unavailable"
    formal_status_raw <- if (nrow(hit) > 0) normalize_scalar_value(hit$formal_status[[1]], "missing") else "missing"
    warning <- if (nrow(hit) > 0) normalize_scalar_value(hit$warning[[1]]) else "comparison_id missing from gene_program_registry.tsv"
    downstream_status <- if (identical(result_level, "pseudobulk_formal") && identical(formal_status_raw, "formal")) {
      "formal"
    } else if (grepl("formal preferred", warning, ignore.case = TRUE)) {
      "exploratory_forced"
    } else if (identical(result_level, "cell_level_exploratory")) {
      "exploratory_only"
    } else if (identical(result_level, "unavailable")) {
      "unavailable"
    } else {
      formal_status_raw
    }
    bg <- enrichment_background_for_layer_06(cfg, layer_id)

    data.frame(
      target_id = normalize_scalar_value(target$target_id[[1]]),
      source_question_id = normalize_scalar_value(target$source_question_id[[1]]),
      layer_id = layer_id,
      comparison_id = comparison_id,
      database = toupper(normalize_scalar_value(target$database[[1]], "GO")),
      organism = normalize_scalar_value(target$organism[[1]], "chicken_primary"),
      analysis_mode = normalize_scalar_value(target$analysis_mode[[1]], if (nrow(hit) > 0) hit$analysis_mode[[1]] else ""),
      gene_program_role = normalize_scalar_value(target$gene_program_role[[1]], if (nrow(hit) > 0) hit$gene_program_role[[1]] else ""),
      min_genes = normalize_scalar_value(target$min_genes[[1]], as.character(cfg$enrichment_min_input_genes)),
      deg_tsv = deg_tsv,
      deg_source = deg_source,
      deg_inference_status = downstream_status,
      formal_status = downstream_status,
      result_level = result_level,
      biological_replicates = "unknown",
      background_tsv = bg$path,
      background_gene_n = length(bg$genes),
      marker_results_tsv = marker_path,
      pseudobulk_results_tsv = deg_path,
      registry_warning = warning,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

enrichment_layer_status_06 <- function(cfg) {
  layer_df <- read_tsv_optional(cfg$layer_status_file)
  rows <- list()
  if (nrow(layer_df) > 0) {
    expected <- c("layer_id", "layer_role", "annotated_rds", "clustered_rds", "status")
    for (col in expected) {
      if (!col %in% colnames(layer_df)) layer_df[[col]] <- ""
    }
    for (idx in seq_len(nrow(layer_df))) {
      row <- layer_df[idx, , drop = FALSE]
      obj_path <- normalize_scalar_value(row$annotated_rds[[1]])
      if (!nzchar(obj_path) || !file.exists(obj_path)) {
        obj_path <- normalize_scalar_value(row$clustered_rds[[1]])
      }
      if (nzchar(row$layer_id[[1]]) && nzchar(obj_path) && file.exists(obj_path)) {
        rows[[length(rows) + 1L]] <- data.frame(
          layer_id = normalize_scalar_value(row$layer_id[[1]]),
          object_rds = obj_path,
          source = "layer_status",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  panorama_path <- resolve_manifest_output_optional_06(cfg$module_03d_manifest_path, c("annotated_object"))
  if (nzchar(panorama_path) && file.exists(panorama_path)) {
    rows[[length(rows) + 1L]] <- data.frame(layer_id = cfg$panorama_layer_id, object_rds = panorama_path, source = "03d_manifest", stringsAsFactors = FALSE)
  }
  if (file.exists(cfg$module_04b_manifest_path)) {
    manifest <- read_manifest_local(cfg$module_04b_manifest_path)
    keys <- names(manifest$outputs %||% list())
    keys <- keys[startsWith(keys, "annotated_")]
    for (key in keys) {
      obj_path <- resolve_output_local(manifest, key)
      if (file.exists(obj_path)) {
        rows[[length(rows) + 1L]] <- data.frame(layer_id = sub("^annotated_", "", key), object_rds = obj_path, source = "04b_manifest", stringsAsFactors = FALSE)
      }
    }
  }
  if (length(rows) == 0) {
    return(empty_df_06(c("layer_id", "object_rds", "source")))
  }
  out <- dplyr::bind_rows(rows)
  out <- out[nzchar(out$layer_id) & nzchar(out$object_rds), , drop = FALSE]
  out[!duplicated(out$layer_id), , drop = FALSE]
}

enrichment_background_path_06 <- function(cfg, layer_id) {
  file.path(cfg$enrichment_background_table_dir, paste0("background_genes_", safe_id_06(layer_id), ".tsv"))
}

enrichment_background_for_layer_06 <- function(cfg, layer_id) {
  path <- enrichment_background_path_06(cfg, layer_id)
  existing <- read_tsv_optional(path)
  if (nrow(existing) > 0 && "gene" %in% colnames(existing)) {
    genes <- unique(trimws(as.character(existing$gene)))
    genes <- genes[nzchar(genes)]
    return(list(path = path, genes = genes, status = "cached"))
  }

  status_df <- enrichment_layer_status_06(cfg)
  hit <- status_df[status_df$layer_id == layer_id, , drop = FALSE]
  if (nrow(hit) == 0 || !file.exists(hit$object_rds[[1]]) || !requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
    write_tsv_local(empty_df_06(c("layer_id", "gene", "background_source", "status")), path)
    return(list(path = path, genes = character(0), status = "unavailable"))
  }

  genes <- tryCatch(
    {
      obj <- readRDS(hit$object_rds[[1]])
      if (exists("maybe_join_layers", mode = "function")) {
        obj <- maybe_join_layers(obj)
      }
      mat <- if (exists("get_assay_matrix", mode = "function")) {
        get_assay_matrix(obj, assay = "RNA", type = "counts")
      } else {
        Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
      }
      rownames(mat)[Matrix::rowSums(mat > 0) > 0]
    },
    error = function(e) character(0)
  )
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes)]
  out <- if (length(genes) == 0) {
    empty_df_06(c("layer_id", "gene", "background_source", "status"))
  } else {
    data.frame(layer_id = layer_id, gene = genes, background_source = normalize_scalar_value(hit$source[[1]]), status = "ok", stringsAsFactors = FALSE)
  }
  write_tsv_local(out, path)
  list(path = path, genes = genes, status = if (length(genes) == 0) "empty" else "ok")
}

read_background_genes_06 <- function(path) {
  df <- read_tsv_optional(path)
  if (nrow(df) == 0 || !"gene" %in% colnames(df)) {
    return(character(0))
  }
  genes <- unique(trimws(as.character(df$gene)))
  genes[nzchar(genes)]
}

detect_deg_columns_06 <- function(df) {
  pick <- function(candidates) {
    hit <- candidates[candidates %in% colnames(df)]
    if (length(hit) == 0) "" else hit[[1]]
  }
  list(
    gene_col = pick(c("gene", "Gene", "symbol", "SYMBOL", "gene_symbol", "feature", "feature_id", "gene_id")),
    cluster_col = pick(c("cluster_id", "cluster", "cluster_name", "group_id", "label")),
    logfc_col = pick(c("avg_log2FC", "avg_logFC", "log2FC", "logFC", "lfc", "coef", "estimate")),
    padj_col = pick(c("p_val_adj", "p_adj", "padj", "adj.P.Val", "FDR", "fdr", "p_adj.loc", "p_adj.glb", "qval", "qvalue"))
  )
}

standardize_deg_table_06 <- function(df, alpha = 0.05) {
  cols <- detect_deg_columns_06(df)
  if (!nzchar(cols$gene_col)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  out <- as.data.frame(df, stringsAsFactors = FALSE)
  out$.gene <- trimws(as.character(out[[cols$gene_col]]))
  out$.cluster <- if (nzchar(cols$cluster_col)) trimws(as.character(out[[cols$cluster_col]])) else "all"
  out$.cluster[!nzchar(out$.cluster)] <- "all"
  out$.logfc <- if (nzchar(cols$logfc_col)) suppressWarnings(as.numeric(out[[cols$logfc_col]])) else NA_real_
  out$.padj <- if (nzchar(cols$padj_col)) suppressWarnings(as.numeric(out[[cols$padj_col]])) else NA_real_
  out$.has_padj <- nzchar(cols$padj_col)
  out$.has_logfc <- nzchar(cols$logfc_col)
  out$.significant <- if (nzchar(cols$padj_col)) is.finite(out$.padj) & out$.padj <= alpha else TRUE
  out[nzchar(out$.gene), , drop = FALSE]
}

cluster_ids_from_deg_table_06 <- function(deg_df) {
  std <- standardize_deg_table_06(deg_df)
  if (nrow(std) == 0) {
    return(character(0))
  }
  sort(unique(std$.cluster[nzchar(std$.cluster)]))
}

prepare_deg_gene_list_from_df_06 <- function(deg_df, cfg, cluster_id = NULL) {
  std <- standardize_deg_table_06(deg_df, alpha = cfg$deg_alpha)
  if (nrow(std) == 0) {
    empty <- list(genes = character(0), input_n = 0L, status = "missing_gene_column", reason = "DEG table has no supported gene column")
    return(list(all = empty, up = empty, down = empty))
  }
  if (!is.null(cluster_id) && nzchar(cluster_id)) {
    std <- std[std$.cluster == cluster_id, , drop = FALSE]
  }
  total_n <- nrow(std)
  std <- std[std$.significant, , drop = FALSE]
  if (nrow(std) == 0) {
    empty <- list(genes = character(0), input_n = 0L, status = "no_significant_deg", reason = sprintf("0 significant DEG among %s rows", total_n))
    return(list(all = empty, up = empty, down = empty))
  }

  order_genes <- function(x) {
    if (nrow(x) == 0) {
      return(character(0))
    }
    score <- abs(x$.logfc)
    score[!is.finite(score)] <- 0
    padj <- x$.padj
    padj[!is.finite(padj)] <- 1
    x <- x[order(-score, padj, x$.gene), , drop = FALSE]
    unique(x$.gene)[seq_len(min(cfg$enrichment_top_n_genes, length(unique(x$.gene))))]
  }

  split_rows <- list(
    all = std,
    up = if (any(std$.has_logfc)) std[is.finite(std$.logfc) & std$.logfc > 0, , drop = FALSE] else std[FALSE, , drop = FALSE],
    down = if (any(std$.has_logfc)) std[is.finite(std$.logfc) & std$.logfc < 0, , drop = FALSE] else std[FALSE, , drop = FALSE]
  )
  out <- lapply(names(split_rows), function(direction) {
    genes <- order_genes(split_rows[[direction]])
    list(
      genes = genes,
      input_n = length(genes),
      status = if (length(genes) >= cfg$enrichment_min_input_genes) "ok" else "too_few_genes",
      reason = if (length(genes) >= cfg$enrichment_min_input_genes) "" else sprintf("gene_direction=%s has %s genes after filtering", direction, length(genes))
    )
  })
  names(out) <- names(split_rows)
  out
}

prepare_deg_gene_list <- function(deg_tsv_path, cfg, cluster_id = NULL) {
  path <- normalize_scalar_value(deg_tsv_path)
  if (!nzchar(path) || !file.exists(path)) {
    empty <- list(genes = character(0), input_n = 0L, status = "missing_deg_tsv", reason = sprintf("missing DEG TSV: %s", path))
    return(list(all = empty, up = empty, down = empty))
  }
  prepare_deg_gene_list_from_df_06(read_tsv_optional(path), cfg, cluster_id)
}

convert_symbols_to_entrezid <- function(genes, org_db, keytype = "SYMBOL") {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0) {
    return(list(ids = character(0), input_n = 0L, mapped_n = 0L, mapping_rate = NA_real_, keytype = keytype))
  }

  available <- tryCatch(AnnotationDbi::keytypes(org_db), error = function(e) character(0))
  candidate_keytypes <- unique(c(keytype, "SYMBOL", "ALIAS", "ENSEMBL"))
  candidate_keytypes <- candidate_keytypes[candidate_keytypes %in% available]
  if (length(candidate_keytypes) == 0) {
    candidate_keytypes <- keytype
  }

  best <- list(ids = character(0), mapped_n = 0L, keytype = candidate_keytypes[[1]])
  for (kt in candidate_keytypes) {
    query_genes <- if (identical(kt, "ENSEMBL")) sub("\\..*$", "", genes) else genes
    mapped <- tryCatch(
      AnnotationDbi::mapIds(org_db, keys = query_genes, column = "ENTREZID", keytype = kt, multiVals = "first"),
      error = function(e) setNames(rep(NA_character_, length(query_genes)), query_genes)
    )
    ids <- unique(trimws(as.character(mapped)))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    mapped_n <- sum(!is.na(mapped) & nzchar(as.character(mapped)))
    if (mapped_n > best$mapped_n) {
      best <- list(ids = ids, mapped_n = mapped_n, keytype = kt)
    }
  }

  list(
    ids = best$ids,
    input_n = length(genes),
    mapped_n = best$mapped_n,
    mapping_rate = safe_rate(best$mapped_n, length(genes)),
    keytype = best$keytype
  )
}

safe_run_enrichment <- function(fn, ...) {
  tryCatch(
    fn(...),
    error = function(e) list(result = NULL, status = "error", reason = conditionMessage(e))
  )
}

run_go_enrichment_single <- function(gene_list, org_db, ont, pvalue_cutoff, qvalue_cutoff, min_gs_size = 10L, max_gs_size = 500L, universe_ids = character(0)) {
  if (!clusterprofiler_available_06()) {
    return(list(result = NULL, status = "clusterprofiler_unavailable", reason = "clusterProfiler is not loadable in this runtime"))
  }
  if (length(gene_list) == 0) {
    return(list(result = NULL, status = "too_few_mapped_genes", reason = "0 mapped ENTREZ IDs"))
  }
  res <- tryCatch(
    clusterProfiler::enrichGO(
      gene = unique(as.character(gene_list)),
      universe = if (length(universe_ids) > 0) unique(as.character(universe_ids)) else NULL,
      OrgDb = org_db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      readable = TRUE
    ),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(list(result = NULL, status = "error", reason = conditionMessage(res)))
  }
  df <- as.data.frame(res)
  if (nrow(df) == 0) {
    return(list(result = res, status = "no_terms", reason = "clusterProfiler returned 0 terms"))
  }
  list(result = res, status = "ok", reason = "")
}

is_timeout_error_06 <- function(message) {
  grepl("time.?out|timed out|Timeout", message, ignore.case = TRUE)
}

run_kegg_enrichment_single <- function(gene_list, organism, pvalue_cutoff, qvalue_cutoff, min_gs_size = 10L, max_gs_size = 500L, timeout_sec = 60L, universe_ids = character(0)) {
  if (!clusterprofiler_available_06()) {
    return(list(result = NULL, status = "clusterprofiler_unavailable", reason = "clusterProfiler is not loadable in this runtime"))
  }
  if (length(gene_list) == 0) {
    return(list(result = NULL, status = "too_few_mapped_genes", reason = "0 mapped ENTREZ IDs"))
  }

  run_attempt <- function(args) {
    do.call(clusterProfiler::enrichKEGG, args)
  }
  base_args <- list(
    gene = unique(as.character(gene_list)),
    organism = organism,
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size
  )
  if (length(universe_ids) > 0) {
    base_args$universe <- unique(as.character(universe_ids))
  }
  attempts <- list(
    c(base_args, list(keyType = "ncbi-geneid")),
    c(base_args, list(keyType = "kegg")),
    base_args
  )
  timeout_sec <- suppressWarnings(as.numeric(timeout_sec))
  if (!is.finite(timeout_sec) || timeout_sec <= 0) {
    timeout_sec <- 60
  }
  old_timeout <- getOption("timeout")
  options(timeout = timeout_sec)
  on.exit(options(timeout = old_timeout), add = TRUE)

  last_error <- NULL
  for (args in attempts) {
    res <- tryCatch(run_attempt(args), error = function(e) e)
    if (!inherits(res, "error")) {
      df <- as.data.frame(res)
      if (nrow(df) == 0) {
        return(list(result = res, status = "no_terms", reason = "clusterProfiler returned 0 terms"))
      }
      return(list(result = res, status = "ok", reason = ""))
    }
    last_error <- conditionMessage(res)
    if (is_timeout_error_06(last_error)) {
      return(list(result = NULL, status = "timeout", reason = sprintf("enrichKEGG timed out after %s seconds: %s", timeout_sec, last_error)))
    }
  }
  list(result = NULL, status = "error", reason = last_error %||% "enrichKEGG failed")
}

format_gprofiler_result_06 <- function(gost_result, extra_cols = list()) {
  extra_names <- names(extra_cols)
  if (is.null(gost_result) || is.null(gost_result$result) || nrow(gost_result$result) == 0) {
    return(empty_enrichment_result_06(extra_names))
  }
  raw <- as.data.frame(gost_result$result, stringsAsFactors = FALSE)
  value_or <- function(col, default = "") {
    if (col %in% colnames(raw)) raw[[col]] else rep(default, nrow(raw))
  }
  num_or <- function(col) {
    suppressWarnings(as.numeric(value_or(col, NA_real_)))
  }
  intersection <- value_or("intersection", "")
  if (is.list(intersection)) {
    intersection <- vapply(intersection, function(x) paste(as.character(x), collapse = "/"), character(1))
  }
  query_size <- num_or("query_size")
  intersection_size <- num_or("intersection_size")
  term_size <- num_or("term_size")
  effective_domain_size <- num_or("effective_domain_size")
  p <- num_or("p_value")
  df <- data.frame(
    ID = value_or("term_id", ""),
    Description = value_or("term_name", ""),
    GeneRatio = ifelse(is.finite(intersection_size) & is.finite(query_size), paste0(intersection_size, "/", query_size), ""),
    BgRatio = ifelse(is.finite(term_size) & is.finite(effective_domain_size), paste0(term_size, "/", effective_domain_size), ""),
    pvalue = NA_real_,
    p.adjust = p,
    qvalue = NA_real_,
    p_adjust_method = "g_SCS",
    geneID = intersection,
    Count = intersection_size,
    stringsAsFactors = FALSE
  )
  for (nm in extra_names) {
    df[[nm]] <- extra_cols[[nm]]
  }
  df[, unique(c(extra_names, colnames(df))), drop = FALSE]
}

run_gprofiler_enrichment_single <- function(gene_symbols, organism, sources, pvalue_cutoff, background_symbols = character(0)) {
  genes <- unique(trimws(as.character(gene_symbols)))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0) {
    return(list(result = NULL, status = "too_few_genes", reason = "0 submitted gene symbols"))
  }
  if (!gprofiler_available_06()) {
    return(list(result = NULL, status = "enrichment_backend_missing", reason = "Neither clusterProfiler nor gprofiler2 is loadable"))
  }
  res <- tryCatch(
    gprofiler2::gost(
      query = genes,
      organism = organism,
      sources = sources,
      custom_bg = if (length(background_symbols) > 0) unique(as.character(background_symbols)) else NULL,
      significant = FALSE,
      correction_method = "g_SCS"
    ),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(list(result = NULL, status = "gprofiler_error", reason = conditionMessage(res)))
  }
  result_df <- format_gprofiler_result_06(res)
  if (nrow(result_df) == 0) {
    return(list(result = result_df, status = "no_terms", reason = "gprofiler2 returned 0 terms"))
  }
  result_df <- result_df[suppressWarnings(as.numeric(result_df$p.adjust)) <= pvalue_cutoff | is.na(suppressWarnings(as.numeric(result_df$p.adjust))), , drop = FALSE]
  if (nrow(result_df) == 0) {
    return(list(result = result_df, status = "no_terms", reason = "0 gprofiler2 terms passed pvalue_cutoff"))
  }
  list(result = result_df, status = "ok", reason = "gprofiler2 fallback")
}

format_enrichment_result_tsv <- function(enrich_result, extra_cols = list()) {
  extra_names <- names(extra_cols)
  if (is.null(enrich_result)) {
    return(empty_enrichment_result_06(extra_names))
  }
  df <- as.data.frame(enrich_result)
  if (nrow(df) == 0) {
    return(empty_enrichment_result_06(extra_names))
  }
  for (nm in extra_names) {
    df[[nm]] <- extra_cols[[nm]]
  }
  standard_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  for (col in standard_cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  if (!"p_adjust_method" %in% colnames(df)) {
    df$p_adjust_method <- ""
  }
  extra_first <- extra_names[extra_names %in% colnames(df)]
  df[, unique(c(extra_first, standard_cols, "p_adjust_method", setdiff(colnames(df), c(extra_first, standard_cols, "p_adjust_method")))), drop = FALSE]
}

significant_enrichment_df_06 <- function(df, pvalue_cutoff, qvalue_cutoff) {
  if (nrow(df) == 0) {
    return(df)
  }
  padj <- suppressWarnings(as.numeric(df$`p.adjust`))
  qval <- suppressWarnings(as.numeric(df$qvalue))
  keep <- is.finite(padj) & padj <= pvalue_cutoff
  if (any(is.finite(qval))) {
    keep <- keep & is.finite(qval) & qval <= qvalue_cutoff
  }
  df[keep, , drop = FALSE]
}

write_enrichment_result_06 <- function(result_df, paths, cfg) {
  write_tsv_local(result_df, paths$enrichment_tsv)
  sig_df <- significant_enrichment_df_06(result_df, cfg$enrichment_pvalue_cutoff, cfg$enrichment_qvalue_cutoff)
  write_tsv_local(sig_df, paths$enrichment_significant_tsv)
  list(all_n = nrow(result_df), significant_n = nrow(sig_df))
}

top_plot_df_06 <- function(enrich_result, top_n) {
  df <- as.data.frame(enrich_result)
  if (nrow(df) == 0) {
    return(df)
  }
  padj <- suppressWarnings(as.numeric(df$`p.adjust`))
  padj[!is.finite(padj)] <- 1
  df <- df[order(padj, -suppressWarnings(as.numeric(df$Count))), , drop = FALSE]
  df[seq_len(min(top_n, nrow(df))), , drop = FALSE]
}

is_enrichplot_result_06 <- function(enrich_result) {
  if (is.null(enrich_result) || is.data.frame(enrich_result)) {
    return(FALSE)
  }
  any(vapply(c("enrichResult", "gseaResult", "compareClusterResult"), function(cls) {
    inherits(enrich_result, cls) || isTRUE(tryCatch(methods::is(enrich_result, cls), error = function(e) FALSE))
  }, logical(1)))
}

plot_enrichment_fallback_06 <- function(enrich_result, title, top_n, output_png, geom = c("dot", "bar")) {
  geom <- match.arg(geom)
  df <- top_plot_df_06(enrich_result, top_n)
  if (nrow(df) == 0) {
    return("")
  }
  df$Description <- as.character(df$Description)
  df$Description <- ifelse(nchar(df$Description) > 70, paste0(substr(df$Description, 1, 67), "..."), df$Description)
  df$score <- -log10(pmax(suppressWarnings(as.numeric(df$`p.adjust`)), .Machine$double.xmin))
  df$Count <- suppressWarnings(as.numeric(df$Count))
  ensure_dir(dirname(output_png))
  if (identical(geom, "dot")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = score, y = stats::reorder(Description, score), size = Count, color = score)) +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::scale_color_viridis_c(option = "C", end = 0.9) +
      ggplot2::labs(title = title, x = "-log10 adjusted P", y = NULL, color = "-log10 padj") +
      ggplot2::theme_bw(base_size = 10)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = score, y = stats::reorder(Description, score))) +
      ggplot2::geom_col(fill = "#4C78A8") +
      ggplot2::labs(title = title, x = "-log10 adjusted P", y = NULL) +
      ggplot2::theme_bw(base_size = 10)
  }
  ggplot2::ggsave(output_png, p, width = 7.5, height = 5.5, dpi = 300, bg = "white")
  output_png
}

plot_enrichment_dotplot <- function(enrich_result, title, top_n = 20L, output_png) {
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) {
    return("")
  }
  ensure_dir(dirname(output_png))
  plotted <- FALSE
  if (is_enrichplot_result_06(enrich_result) && requireNamespace("enrichplot", quietly = TRUE)) {
    plotted <- tryCatch({
      p <- enrichplot::dotplot(enrich_result, showCategory = min(top_n, nrow(as.data.frame(enrich_result)))) +
        ggplot2::ggtitle(title) +
        ggplot2::theme_bw(base_size = 10)
      ggplot2::ggsave(output_png, p, width = 7.5, height = 5.5, dpi = 300, bg = "white")
      TRUE
    }, error = function(e) FALSE)
  }
  if (!plotted) {
    return(plot_enrichment_fallback_06(enrich_result, title, top_n, output_png, geom = "dot"))
  }
  output_png
}

plot_enrichment_barplot <- function(enrich_result, title, top_n = 15L, output_png) {
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) {
    return("")
  }
  ensure_dir(dirname(output_png))
  plotted <- FALSE
  if (is_enrichplot_result_06(enrich_result) && requireNamespace("enrichplot", quietly = TRUE)) {
    plotted <- tryCatch({
      p <- graphics::barplot(enrich_result, showCategory = min(top_n, nrow(as.data.frame(enrich_result)))) +
        ggplot2::ggtitle(title) +
        ggplot2::theme_bw(base_size = 10)
      ggplot2::ggsave(output_png, p, width = 7.5, height = 5.5, dpi = 300, bg = "white")
      TRUE
    }, error = function(e) FALSE)
  }
  if (!plotted) {
    return(plot_enrichment_fallback_06(enrich_result, title, top_n, output_png, geom = "bar"))
  }
  output_png
}

plot_empty_enrichment_06 <- function(output_png, title = "No significant enrichment") {
  ensure_dir(dirname(output_png))
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = title, size = 4) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void()
  ggplot2::ggsave(output_png, p, width = 7, height = 4.5, dpi = 300, bg = "white")
  output_png
}

read_enrichment_tables_from_manifest_06 <- function(manifest_df) {
  if (is.null(manifest_df) || nrow(manifest_df) == 0 || !"enrichment_tsv" %in% colnames(manifest_df)) {
    return(empty_enrichment_result_06())
  }
  paths <- unique(vapply(manifest_df$enrichment_tsv, normalize_scalar_value, character(1)))
  paths <- paths[nzchar(paths) & file.exists(paths)]
  rows <- lapply(paths, function(path) {
    df <- read_tsv_optional(path)
    if (nrow(df) == 0) {
      return(NULL)
    }
    df
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(empty_enrichment_result_06())
  }
  dplyr::bind_rows(rows)
}

plot_cross_cluster_heatmap <- function(combined_results, value_col, output_png, title = "Cross-cluster enrichment", top_n = 30L) {
  if (is.null(combined_results) || nrow(combined_results) == 0 || !"Description" %in% colnames(combined_results)) {
    return(plot_empty_enrichment_06(output_png, title))
  }
  df <- combined_results
  for (col in c("layer_id", "comparison_id", "cluster_id", "source_species")) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  if (!value_col %in% colnames(df)) {
    value_col <- "p.adjust"
  }
  values <- suppressWarnings(as.numeric(df[[value_col]]))
  score <- if (identical(value_col, "p.adjust")) -log10(pmax(values, .Machine$double.xmin)) else values
  df$score <- score
  df <- df[is.finite(df$score) & nzchar(as.character(df$Description)), , drop = FALSE]
  if (nrow(df) == 0) {
    return(plot_empty_enrichment_06(output_png, title))
  }
  df$cluster_axis <- paste(df$layer_id, df$comparison_id, df$cluster_id, df$source_species, sep = " / ")
  term_order <- aggregate(score ~ Description, df, max)
  term_order <- term_order[order(-term_order$score), , drop = FALSE]
  keep_terms <- head(term_order$Description, top_n)
  df <- df[df$Description %in% keep_terms, , drop = FALSE]
  df$Description <- factor(df$Description, levels = rev(keep_terms))
  ensure_dir(dirname(output_png))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cluster_axis, y = Description, fill = score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_viridis_c(option = "C", end = 0.9) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = if (identical(value_col, "p.adjust")) "-log10 padj" else value_col) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  ggplot2::ggsave(output_png, p, width = 10, height = 7, dpi = 300, bg = "white")
  output_png
}

build_shared_pathways_06 <- function(combined_results, min_cluster_n = 2L) {
  if (is.null(combined_results) || nrow(combined_results) == 0) {
    return(empty_df_06(c("analysis_type", "ontology", "source_species", "ID", "Description", "cluster_n", "clusters", "min_p_adjust")))
  }
  df <- combined_results
  for (col in c("analysis_type", "ontology", "source_species", "ID", "Description", "layer_id", "comparison_id", "cluster_id", "p.adjust")) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  df$cluster_axis <- paste(df$layer_id, df$comparison_id, df$cluster_id, sep = "/")
  split_key <- paste(df$analysis_type, df$ontology, df$source_species, df$ID, df$Description, sep = "\r")
  rows <- lapply(split(df, split_key), function(x) {
    clusters <- sort(unique(x$cluster_axis[nzchar(x$cluster_axis)]))
    data.frame(
      analysis_type = x$analysis_type[[1]],
      ontology = x$ontology[[1]],
      source_species = x$source_species[[1]],
      ID = x$ID[[1]],
      Description = x$Description[[1]],
      cluster_n = length(clusters),
      clusters = paste(clusters, collapse = ";"),
      min_p_adjust = suppressWarnings(min(as.numeric(x$`p.adjust`), na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(rows)
  out <- out[out$cluster_n >= min_cluster_n, , drop = FALSE]
  out[order(-out$cluster_n, out$min_p_adjust, out$Description), , drop = FALSE]
}
