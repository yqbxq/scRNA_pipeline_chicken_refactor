empty_df_07 <- function(cols = character()) {
  if (exists("empty_df_05", mode = "function")) {
    return(empty_df_05(cols))
  }
  as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols), stringsAsFactors = FALSE)
}

normalize_path_07 <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

resolve_manifest_output_optional_07 <- function(manifest_path, keys) {
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

pick_first_col_07 <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) "" else hit[[1]]
}

quality_rank_07 <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(99L, length(x))
  out[x == "gold"] <- 1L
  out[x == "silver"] <- 2L
  out[x == "ambiguous"] <- 3L
  out
}

build_chicken_to_human_lut <- function(ortholog_csv, prefer = c("gold", "silver")) {
  if (!file.exists(ortholog_csv)) {
    stop(sprintf("Ortholog CSV does not exist: %s", ortholog_csv), call. = FALSE)
  }

  raw <- read.csv(ortholog_csv, stringsAsFactors = FALSE, check.names = FALSE)
  chicken_col <- pick_first_col_07(raw, c(
    "chicken_symbol", "external_gene_name", "source_gene_name",
    "source_symbol", "gene_name", "gene_symbol"
  ))
  human_col <- pick_first_col_07(raw, c(
    "human_symbol", "target_gene_name", "target_symbol",
    "hsapiens_homolog_associated_gene_name"
  ))
  if (!nzchar(chicken_col) || !nzchar(human_col)) {
    stop(
      sprintf(
        "Ortholog CSV missing supported chicken/human symbol columns: %s",
        paste(colnames(raw), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  quality_col <- pick_first_col_07(raw, c("pair_quality", "quality", "ortholog_quality"))
  if (!nzchar(quality_col)) {
    raw$pair_quality <- "silver"
    quality_col <- "pair_quality"
  }

  if ("target_species" %in% colnames(raw)) {
    raw <- raw[tolower(as.character(raw$target_species)) %in% c("", "human", "hsapiens", "homo_sapiens"), , drop = FALSE]
  }

  lut <- data.frame(
    chicken_symbol = trimws(as.character(raw[[chicken_col]])),
    human_symbol = trimws(as.character(raw[[human_col]])),
    pair_quality = tolower(trimws(as.character(raw[[quality_col]]))),
    stringsAsFactors = FALSE
  )
  lut <- lut[nzchar(lut$chicken_symbol) & nzchar(lut$human_symbol), , drop = FALSE]
  lut$pair_quality[!nzchar(lut$pair_quality)] <- "silver"
  lut$quality_rank <- quality_rank_07(lut$pair_quality)
  lut$drop_for_cellchat <- !lut$pair_quality %in% prefer

  lut <- lut[order(lut$quality_rank, lut$chicken_symbol, lut$human_symbol), , drop = FALSE]
  lut <- lut[!duplicated(toupper(lut$chicken_symbol)), , drop = FALSE]
  rownames(lut) <- NULL
  lut
}

mapping_coverage_summary <- function(input_genes, lut) {
  genes <- unique(trimws(as.character(input_genes)))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0 || nrow(lut) == 0) {
    return(data.frame(
      n_genes_input = length(genes),
      n_mapped_gold = 0L,
      n_mapped_silver = 0L,
      n_dropped_ambiguous = 0L,
      n_unmapped = length(genes),
      coverage_pct = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  lut$key <- toupper(lut$chicken_symbol)
  hit <- lut[lut$key %in% toupper(genes), , drop = FALSE]
  mapped_preferred <- hit[!hit$drop_for_cellchat, , drop = FALSE]
  mapped_keys <- unique(mapped_preferred$key)
  data.frame(
    n_genes_input = length(genes),
    n_mapped_gold = length(unique(hit$key[hit$pair_quality == "gold"])),
    n_mapped_silver = length(unique(hit$key[hit$pair_quality == "silver"])),
    n_dropped_ambiguous = length(unique(hit$key[hit$pair_quality == "ambiguous"])),
    n_unmapped = length(setdiff(toupper(genes), unique(hit$key))),
    coverage_pct = safe_rate(length(mapped_keys), length(genes)) * 100,
    stringsAsFactors = FALSE
  )
}

collapse_matrix_rows_07 <- function(mat, groups, collapse = c("max", "sum", "mean", "first")) {
  collapse <- match.arg(collapse)
  groups <- as.character(groups)
  unique_groups <- unique(groups[nzchar(groups)])
  rows <- vector("list", length(unique_groups))
  names(rows) <- unique_groups

  for (group in unique_groups) {
    idx <- which(groups == group)
    if (length(idx) == 1L || identical(collapse, "first")) {
      rows[[group]] <- mat[idx[[1]], , drop = FALSE]
    } else {
      sub_mat <- mat[idx, , drop = FALSE]
      vec <- switch(
        collapse,
        sum = Matrix::colSums(sub_mat),
        mean = Matrix::colMeans(sub_mat),
        max = apply(as.matrix(sub_mat), 2, max),
        first = mat[idx[[1]], , drop = FALSE]
      )
      rows[[group]] <- Matrix::Matrix(vec, nrow = 1, sparse = TRUE)
      colnames(rows[[group]]) <- colnames(mat)
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- names(rows)
  out
}

map_expression_matrix_chicken_to_human <- function(mat, lut, collapse = "max") {
  if (is.null(rownames(mat))) {
    stop("Expression matrix has no rownames for ortholog mapping.", call. = FALSE)
  }
  if (nrow(lut) == 0) {
    return(mat[FALSE, , drop = FALSE])
  }

  map <- lut[!lut$drop_for_cellchat, , drop = FALSE]
  map$key <- toupper(map$chicken_symbol)
  map <- map[map$key %in% toupper(rownames(mat)), , drop = FALSE]
  if (nrow(map) == 0) {
    return(mat[FALSE, , drop = FALSE])
  }

  row_key <- toupper(rownames(mat))
  idx <- match(map$key, row_key)
  keep <- !is.na(idx)
  idx <- idx[keep]
  map <- map[keep, , drop = FALSE]
  sub_mat <- mat[idx, , drop = FALSE]
  collapse_matrix_rows_07(sub_mat, map$human_symbol, collapse = collapse)
}

rename_network_human_to_chicken <- function(network_df, lut, gene_cols) {
  if (nrow(network_df) == 0 || nrow(lut) == 0) {
    return(network_df[FALSE, , drop = FALSE])
  }
  missing_cols <- setdiff(gene_cols, colnames(network_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Network table missing gene columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  map <- lut[!lut$drop_for_cellchat, , drop = FALSE]
  map <- map[order(map$quality_rank, map$human_symbol, map$chicken_symbol), , drop = FALSE]
  map <- map[!duplicated(toupper(map$human_symbol)), , drop = FALSE]
  lookup <- setNames(map$chicken_symbol, toupper(map$human_symbol))

  out <- as.data.frame(network_df, stringsAsFactors = FALSE)
  keep <- rep(TRUE, nrow(out))
  for (col in gene_cols) {
    renamed <- unname(lookup[toupper(as.character(out[[col]]))])
    keep <- keep & !is.na(renamed) & nzchar(renamed)
    out[[col]] <- renamed
  }
  unique(out[keep, , drop = FALSE])
}

rename_ligand_target_matrix_human_to_chicken <- function(ligand_target_matrix, lut) {
  map <- lut[!lut$drop_for_cellchat, , drop = FALSE]
  map <- map[order(map$quality_rank, map$human_symbol, map$chicken_symbol), , drop = FALSE]
  map <- map[!duplicated(toupper(map$human_symbol)), , drop = FALSE]
  lookup <- setNames(map$chicken_symbol, toupper(map$human_symbol))

  row_chicken <- unname(lookup[toupper(rownames(ligand_target_matrix))])
  col_chicken <- unname(lookup[toupper(colnames(ligand_target_matrix))])
  row_keep <- !is.na(row_chicken) & nzchar(row_chicken)
  col_keep <- !is.na(col_chicken) & nzchar(col_chicken)
  mat <- ligand_target_matrix[row_keep, col_keep, drop = FALSE]
  rownames(mat) <- row_chicken[row_keep]
  colnames(mat) <- col_chicken[col_keep]
  mat <- mat[!duplicated(rownames(mat)), !duplicated(colnames(mat)), drop = FALSE]
  mat
}

empty_cellchat_index_07 <- function() {
  empty_df_07(c(
    "pair_id", "layer_id", "condition_value", "condition_split_var", "sender", "receiver", "cell_type_col",
    "n_cells", "n_cell_types", "status", "reason", "cellchat_rds_path",
    "lr_table_path", "pathway_table_path", "bubble_png", "network_png", "heatmap_png"
  ))
}

empty_mapping_summary_07 <- function() {
  empty_df_07(c(
    "pair_id", "layer_id", "condition_value", "n_genes_input",
    "n_mapped_gold", "n_mapped_silver", "n_dropped_ambiguous",
    "n_unmapped", "coverage_pct"
  ))
}

empty_nichenet_index_07 <- function() {
  empty_df_07(c(
    "pair_id", "layer_id", "condition_value", "condition_split_var", "sender_set", "receiver_set",
    "sender_cell_types", "receiver_cell_types", "deg_source", "deg_status", "status", "reason",
    "ligand_activity_tsv", "ligand_target_links_tsv", "nichenet_rds",
    "ligand_activity_heatmap_png", "ligand_target_heatmap_png", "circos_png"
  ))
}

empty_roles_resolved_07 <- function() {
  empty_df_07(c(
    "pair_id", "layer_id", "condition_value", "sender_set", "receiver_set",
    "sender_cell_types", "receiver_cell_types", "missing_cell_types", "status"
  ))
}

empty_cross_validation_07 <- function() {
  empty_df_07(c(
    "pair_id", "layer_id", "condition_value", "cellchat_lr_n",
    "nichenet_top_ligand_n", "intersect_ligand_n", "jaccard", "status"
  ))
}

empty_consensus_lr_07 <- function() {
  empty_df_07(c(
    "pair_id", "layer_id", "condition_value", "source", "target",
    "ligand", "receptor", "prob", "pval", "pathway_name",
    "ligand_activity_rank", "ligand_activity_score"
  ))
}

read_communication_comparisons_07 <- function(cfg) {
  df <- read_tsv_optional(cfg$comparison_sheet)
  if (nrow(df) == 0) {
    return(empty_df_07(c(
      "comparison_id", "ident_1", "ident_2", "enabled", "group_var",
      "batch_var", "layer_scope", "min_biological_replicates", "subset_column",
      "subset_value", "force_exploratory", "min_cells_per_group",
      "logfc_threshold", "nichenet_sender", "nichenet_receiver"
    )))
  }
  if (exists("read_deg_comparison_sheet", mode = "function")) {
    base <- read_deg_comparison_sheet(cfg)
  } else {
    base <- read_comparison_sheet_local(cfg)
  }
  for (col in c("nichenet_sender", "nichenet_receiver")) {
    if (col %in% colnames(df)) {
      base[[col]] <- vapply(df[[col]], normalize_scalar_value, character(1))
    } else {
      base[[col]] <- ""
    }
  }
  base
}

communication_layer_status_07 <- function(cfg) {
  layer_df <- read_tsv_optional(cfg$layer_status_file)
  rows <- list()

  if (nrow(layer_df) > 0) {
    expected <- c("layer_id", "layer_role", "cluster_column", "annotated_rds", "clustered_rds", "status")
    for (col in expected) {
      if (!col %in% colnames(layer_df)) {
        layer_df[[col]] <- ""
      }
    }
    for (idx in seq_len(nrow(layer_df))) {
      row <- layer_df[idx, , drop = FALSE]
      obj_path <- normalize_scalar_value(row$annotated_rds[[1]])
      if (!nzchar(obj_path) || !file.exists(obj_path)) {
        obj_path <- normalize_scalar_value(row$clustered_rds[[1]])
      }
      if (!nzchar(obj_path) || !file.exists(obj_path)) {
        next
      }
      layer_id <- normalize_scalar_value(row$layer_id[[1]])
      if (!nzchar(layer_id)) {
        next
      }
      rows[[length(rows) + 1L]] <- data.frame(
        layer_id = layer_id,
        layer_role = normalize_scalar_value(row$layer_role[[1]], ifelse(layer_id == cfg$panorama_layer_id, "panorama", "subcluster")),
        object_rds = obj_path,
        cluster_column = normalize_scalar_value(row$cluster_column[[1]]),
        source = "layer_status",
        stringsAsFactors = FALSE
      )
    }
  }

  panorama_path <- resolve_manifest_output_optional_07(cfg$module_03d_manifest_path, c("annotated_object"))
  if (nzchar(panorama_path) && file.exists(panorama_path)) {
    rows[[length(rows) + 1L]] <- data.frame(
      layer_id = cfg$panorama_layer_id,
      layer_role = "panorama",
      object_rds = panorama_path,
      cluster_column = "",
      source = "03d_manifest",
      stringsAsFactors = FALSE
    )
  }

  if (file.exists(cfg$module_04b_manifest_path)) {
    manifest <- read_manifest_local(cfg$module_04b_manifest_path)
    keys <- names(manifest$outputs %||% list())
    keys <- keys[startsWith(keys, "annotated_")]
    for (key in keys) {
      obj_path <- resolve_output_local(manifest, key)
      if (!file.exists(obj_path)) {
        next
      }
      layer_id <- sub("^annotated_", "", key)
      rows[[length(rows) + 1L]] <- data.frame(
        layer_id = layer_id,
        layer_role = "subcluster",
        object_rds = obj_path,
        cluster_column = "",
        source = "04b_manifest",
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    return(empty_df_07(c("layer_id", "layer_role", "object_rds", "cluster_column", "source")))
  }

  out <- dplyr::bind_rows(rows)
  out <- out[nzchar(out$layer_id) & nzchar(out$object_rds), , drop = FALSE]
  out <- out[!duplicated(out$layer_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

resolve_cell_type_col_07 <- function(seu, layer_id, layer_role = "") {
  configured <- Sys.getenv("COMMUNICATION_CELL_TYPE_COL", unset = "")
  if (nzchar(configured) && configured != "cell_subtype" && configured != "auto") {
    requested <- unique(split_csv_local(configured))
    hit <- requested[requested %in% colnames(seu@meta.data)]
    if (length(hit) == 0) {
      return("")
    }
    return(hit[[1]])
  }

  candidates <- if (identical(configured, "auto")) c() else c("cell_subtype")
  if (!identical(layer_id, "panorama")) {
    candidates <- c(candidates, paste0(layer_id, "_cell_type"))
  }
  candidates <- c(candidates, "cell_type", "annotation_label", "panorama_cell_type", "cluster_id", "seurat_clusters")
  hit <- candidates[candidates %in% colnames(seu@meta.data)]
  if (length(hit) == 0) {
    return("")
  }
  hit[[1]]
}

standardize_comm_metadata_07 <- function(seu) {
  if (exists("standardize_design_metadata_05", mode = "function")) {
    return(standardize_design_metadata_05(seu))
  }
  seu
}

load_comm_layer_object_07 <- function(layer_row) {
  seu <- readRDS(layer_row$object_rds[[1]])
  if (exists("maybe_join_layers", mode = "function")) {
    seu <- maybe_join_layers(seu)
  }
  standardize_comm_metadata_07(seu)
}

comparison_rows_for_layer_07 <- function(comparisons, layer_id) {
  if (nrow(comparisons) == 0) {
    return(comparisons)
  }
  comparisons <- comparisons[comparisons$enabled != "no", , drop = FALSE]
  if (nrow(comparisons) == 0) {
    return(comparisons)
  }
  keep <- vapply(seq_len(nrow(comparisons)), function(i) {
    comparison_applies_to_layer(comparisons[i, , drop = FALSE], layer_id)
  }, logical(1))
  comparisons[keep, , drop = FALSE]
}

resolve_comparison_vars_07 <- function(seu, comparison_row) {
  if (exists("resolve_comparison_vars_05", mode = "function")) {
    return(resolve_comparison_vars_05(seu, comparison_row))
  }
  list(
    comparison_id = normalize_scalar_value(comparison_row$comparison_id[[1]]),
    ident_1 = normalize_scalar_value(comparison_row$ident_1[[1]]),
    ident_2 = normalize_scalar_value(comparison_row$ident_2[[1]]),
    group_var = normalize_scalar_value(comparison_row$group_var[[1]], "group_id"),
    batch_var = normalize_scalar_value(comparison_row$batch_var[[1]], "batch"),
    subset_column = normalize_scalar_value(comparison_row$subset_column[[1]]),
    subset_value = normalize_scalar_value(comparison_row$subset_value[[1]]),
    subset_values = split_csv_local(comparison_row$subset_value[[1]]),
    force_exploratory = FALSE,
    min_biological_replicates = 2L,
    min_cells_per_group = 3L,
    logfc_threshold = 0
  )
}

subset_cells_for_communication_07 <- function(seu, vars) {
  if (exists("subset_cells_for_comparison", mode = "function")) {
    return(subset_cells_for_comparison(seu, vars))
  }
  subset_obj <- apply_comparison_subset(seu, data.frame(
    subset_column = vars$subset_column,
    subset_value = vars$subset_value,
    stringsAsFactors = FALSE
  ))
  if (ncol(subset_obj) == 0) {
    return(list(object = NULL, status = "empty_subset", reason = "comparison matched 0 cells", subset_n = 0L))
  }
  list(object = subset_obj, status = "ok", reason = "", subset_n = ncol(subset_obj))
}

communication_triage_row_07 <- function(layer_id, comparison_id, severity, signal_id, evidence, recommended_action = "") {
  row <- make_triage_row(
    sample_id = sprintf("%s/%s", layer_id, comparison_id),
    severity = severity,
    signal_id = signal_id,
    evidence = evidence,
    recommended_action = recommended_action,
    manual_review_required = "yes"
  )
  row
}

read_deg_for_nichenet_07 <- function(cfg, layer_id, comparison_id) {
  input_grid <- tryCatch(read_enrichment_input_grid_06(cfg), error = function(e) empty_df_07())
  if (nrow(input_grid) > 0 && all(c("layer_id", "comparison_id", "deg_tsv") %in% colnames(input_grid))) {
    hit <- input_grid[input_grid$layer_id == layer_id & input_grid$comparison_id == comparison_id, , drop = FALSE]
    if (nrow(hit) > 0 && nzchar(hit$deg_tsv[[1]]) && file.exists(hit$deg_tsv[[1]])) {
      deg_df <- read_tsv_optional(hit$deg_tsv[[1]])
      prep <- prepare_deg_gene_list_from_df_06(deg_df, cfg, cluster_id = NULL)
      return(list(
        genes = prep$up$genes,
        source = hit$deg_source[[1]],
        status = hit$deg_inference_status[[1]],
        table = deg_df
      ))
    }
  }
  list(genes = character(0), source = "inline_findmarkers", status = "pending_inline", table = data.frame(stringsAsFactors = FALSE))
}

inline_findmarkers_for_nichenet_07 <- function(seu, vars, receiver_cells, alpha = 0.05) {
  if (length(receiver_cells) == 0) {
    return(list(genes = character(0), status = "empty_receiver", table = data.frame(stringsAsFactors = FALSE)))
  }
  obj <- subset(seu, cells = receiver_cells)
  if (exists("maybe_join_layers", mode = "function")) {
    obj <- maybe_join_layers(obj)
  }
  groups <- as.character(obj@meta.data[[vars$group_var]])
  n1 <- sum(groups == vars$ident_1, na.rm = TRUE)
  n2 <- sum(groups == vars$ident_2, na.rm = TRUE)
  if (n1 < vars$min_cells_per_group || n2 < vars$min_cells_per_group) {
    return(list(
      genes = character(0),
      status = "too_few_cells_inline_findmarkers",
      table = data.frame(stringsAsFactors = FALSE)
    ))
  }
  Seurat::Idents(obj) <- vars$group_var
  res <- tryCatch(
    Seurat::FindMarkers(
      obj,
      ident.1 = vars$ident_1,
      ident.2 = vars$ident_2,
      logfc.threshold = vars$logfc_threshold,
      test.use = "wilcox",
      verbose = FALSE
    ),
    error = function(e) e
  )
  if (inherits(res, "error") || is.null(res) || nrow(res) == 0) {
    return(list(genes = character(0), status = "findmarkers_error", table = data.frame(stringsAsFactors = FALSE)))
  }
  out <- tibble::rownames_to_column(as.data.frame(res), "gene")
  std <- standardize_deg_table_06(out, alpha = alpha)
  genes <- unique(std$.gene[std$.significant & (!is.finite(std$.logfc) | std$.logfc > 0)])
  list(genes = genes, status = "inline_findmarkers", table = out)
}

write_empty_png_07 <- function(path, title = "No plot available") {
  ensure_dir(dirname(path))
  grDevices::png(path, width = 1400, height = 900, res = 150)
  graphics::plot.new()
  graphics::text(0.5, 0.5, title)
  grDevices::dev.off()
  path
}
