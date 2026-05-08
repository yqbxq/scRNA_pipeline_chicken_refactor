trajectory_empty_df_09 <- function(cols) {
  as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols), stringsAsFactors = FALSE)
}

trajectory_input_index_cols_09 <- c(
  "pair_id", "source_question_id", "layer_scope", "method", "tools_to_run",
  "split_mode", "split_var", "split_value", "root_group", "terminal_group",
  "coarse_label_var", "fine_label_var", "object_rds", "input_rds",
  "cell_n_before", "outlier_n", "cell_n_after", "retention_fraction",
  "root_group_present", "terminal_group_present", "cc_status", "pca_npcs",
  "umap_reduction", "status", "reason", "produced_at"
)

trajectory_method_index_cols_09 <- c(
  "pair_id", "split_value", "method", "methods_enabled", "input_rds",
  "output_path", "extra_path", "figure_path", "n_cells", "status",
  "reason", "runtime_s"
)

trajectory_read_pairs_09 <- function(cfg) {
  pairs <- read_tsv_optional(cfg$trajectory_pairs_sheet)
  if (nrow(pairs) == 0) {
    return(pairs)
  }
  required <- c(
    "trajectory_id", "source_question_id", "layer_scope", "root_group",
    "terminal_group", "condition_split_var", "condition_split_values",
    "method", "tools_to_run", "methods_extra", "outlier_qc_policy",
    "regress_cell_cycle", "coarse_label_var", "fine_label_var",
    "split_mode", "enabled"
  )
  for (col in required) {
    if (!col %in% colnames(pairs)) {
      pairs[[col]] <- ""
    }
  }
  pairs$enabled <- normalize_flag(pairs$enabled, "yes")
  pairs$method <- tolower(normalize_flag(pairs$method, "trajectory"))
  pairs
}

trajectory_active_pairs_09 <- function(cfg) {
  pairs <- trajectory_read_pairs_09(cfg)
  if (nrow(pairs) == 0) {
    return(pairs)
  }
  pairs[pairs$enabled != "no" & pairs$method == "trajectory", , drop = FALSE]
}

trajectory_manifest_output_optional_09 <- function(manifest_path, key) {
  if (!file.exists(manifest_path)) {
    return("")
  }
  manifest <- read_manifest_local(manifest_path)
  tryCatch(resolve_output_local(manifest, key), error = function(e) "")
}

trajectory_read_method_index_09 <- function(manifest_path, output_key, fallback_path, cols = trajectory_method_index_cols_09) {
  path <- trajectory_manifest_output_optional_09(manifest_path, output_key)
  if (!nzchar(path)) {
    path <- fallback_path
  }
  idx <- read_tsv_optional(path)
  if (nrow(idx) == 0) {
    return(trajectory_empty_df_09(cols))
  }
  for (col in cols) {
    if (!col %in% colnames(idx)) {
      idx[[col]] <- ""
    }
  }
  idx
}

trajectory_read_pseudotime_csv_09 <- function(path, method) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path) || !file.exists(path) || isTRUE(file.info(path)$size == 0)) {
    return(data.frame(cell_id = character(0), pseudotime = numeric(0), method = character(0), stringsAsFactors = FALSE))
  }
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"cell_id" %in% colnames(df) || !"pseudotime" %in% colnames(df)) {
    return(data.frame(cell_id = character(0), pseudotime = numeric(0), method = character(0), stringsAsFactors = FALSE))
  }
  df$cell_id <- as.character(df$cell_id)
  df$pseudotime <- suppressWarnings(as.numeric(df$pseudotime))
  df$method <- method
  df[nzchar(df$cell_id) & is.finite(df$pseudotime), , drop = FALSE]
}

trajectory_read_09a_index <- function(cfg) {
  path <- trajectory_manifest_output_optional_09(cfg$module_09a_manifest_path, "trajectory_input_index")
  if (!nzchar(path)) {
    path <- cfg$trajectory_input_index_tsv
  }
  idx <- read_tsv_optional(path)
  if (nrow(idx) == 0) {
    return(trajectory_empty_df_09(trajectory_input_index_cols_09))
  }
  for (col in trajectory_input_index_cols_09) {
    if (!col %in% colnames(idx)) {
      idx[[col]] <- ""
    }
  }
  idx
}

trajectory_method_extra_tokens_09 <- function(x) {
  x <- normalize_scalar_value(x)
  if (!nzchar(x)) {
    return(character(0))
  }
  parts <- trimws(unlist(strsplit(x, "[,;[:space:]]+", perl = TRUE), use.names = FALSE))
  parts[nzchar(parts)]
}

trajectory_method_enabled_09 <- function(pair_row, method, default = TRUE) {
  method <- tolower(normalize_scalar_value(method))
  tokens <- tolower(trajectory_method_extra_tokens_09(pair_row$methods_extra[[1]]))
  normalized <- gsub("_", "-", tokens, fixed = TRUE)
  method_dash <- gsub("_", "-", method, fixed = TRUE)
  disabled_tokens <- c(paste0("-", method), paste0("-", method_dash), paste0("no-", method_dash), paste0("no_", method))
  enabled_tokens <- c(paste0("+", method), paste0("+", method_dash), method, method_dash)
  if (any(normalized %in% gsub("_", "-", disabled_tokens, fixed = TRUE))) {
    return(FALSE)
  }
  if (any(normalized %in% gsub("_", "-", enabled_tokens, fixed = TRUE))) {
    return(TRUE)
  }
  isTRUE(default)
}

trajectory_execution_units_09 <- function(cfg, include_not_ok = FALSE) {
  pairs <- trajectory_active_pairs_09(cfg)
  idx <- trajectory_read_09a_index(cfg)
  if (nrow(pairs) == 0 || nrow(idx) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  if (!isTRUE(include_not_ok)) {
    idx <- idx[idx$status == "ok", , drop = FALSE]
  }
  if (nrow(idx) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  colnames(pairs)[colnames(pairs) == "trajectory_id"] <- "pair_id"
  merged <- merge(
    idx,
    pairs,
    by = "pair_id",
    all.x = TRUE,
    suffixes = c("", ".pair"),
    sort = FALSE
  )
  for (base_col in c("source_question_id", "layer_scope", "root_group", "terminal_group", "coarse_label_var", "fine_label_var", "split_mode")) {
    pair_col <- paste0(base_col, ".pair")
    if (pair_col %in% colnames(merged)) {
      empty <- !nzchar(normalize_flag(merged[[base_col]], ""))
      merged[[base_col]][empty] <- merged[[pair_col]][empty]
    }
  }
  merged
}

trajectory_load_root_index_09 <- function(cfg) {
  path <- trajectory_manifest_output_optional_09(cfg$module_09c_manifest_path, "trajectory_root_index")
  if (!nzchar(path)) {
    path <- cfg$trajectory_root_index_tsv
  }
  read_tsv_optional(path)
}

trajectory_selected_root_09 <- function(cfg, pair_id, split_value, default_root = "") {
  roots <- trajectory_load_root_index_09(cfg)
  if (nrow(roots) == 0) {
    return(default_root)
  }
  split_value <- display_scalar_value(split_value, "pooled")
  hit <- roots[roots$pair_id == pair_id & roots$split_value == split_value, , drop = FALSE]
  if (nrow(hit) == 0) {
    return(default_root)
  }
  root <- normalize_scalar_value(hit$selected_root[[1]], default_root)
  root
}

trajectory_assay_09 <- function(seu) {
  assays <- tryCatch(names(seu@assays), error = function(e) character(0))
  if ("RNA" %in% assays) {
    return("RNA")
  }
  Seurat::DefaultAssay(seu)
}

trajectory_counts_matrix_09 <- function(seu) {
  assay <- trajectory_assay_09(seu)
  tryCatch(
    get_assay_matrix(seu, assay = assay, type = "counts"),
    error = function(e) get_assay_matrix(seu, assay = assay, type = "data")
  )
}

trajectory_umap_09 <- function(seu) {
  reductions <- tryCatch(names(seu@reductions), error = function(e) character(0))
  hit <- reductions[grepl("umap", reductions, ignore.case = TRUE)]
  if (length(hit) == 0) {
    return(NULL)
  }
  emb <- tryCatch(Seurat::Embeddings(seu, reduction = hit[[1]]), error = function(e) NULL)
  if (is.null(emb) || ncol(emb) < 2) {
    return(NULL)
  }
  emb[, 1:2, drop = FALSE]
}

trajectory_root_cells_09 <- function(seu, label_var, root_label) {
  root_label <- normalize_scalar_value(root_label)
  if (!nzchar(root_label) || root_label %in% c("auto", "*") || !label_var %in% colnames(seu@meta.data)) {
    return(character(0))
  }
  cells <- rownames(seu@meta.data)[as.character(seu@meta.data[[label_var]]) == root_label]
  cells[nzchar(cells)]
}

trajectory_terminal_cells_09 <- function(seu, label_var, terminal_label) {
  terminal_label <- normalize_scalar_value(terminal_label)
  if (!nzchar(terminal_label) || terminal_label %in% c("auto", "*") || !label_var %in% colnames(seu@meta.data)) {
    return(character(0))
  }
  cells <- rownames(seu@meta.data)[as.character(seu@meta.data[[label_var]]) == terminal_label]
  cells[nzchar(cells)]
}

trajectory_scale01_09 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x))
  finite <- is.finite(x)
  if (!any(finite)) {
    return(out)
  }
  rng <- range(x[finite])
  if (diff(rng) == 0) {
    out[finite] <- 0
  } else {
    out[finite] <- (x[finite] - rng[[1]]) / diff(rng)
  }
  out
}

trajectory_write_pseudotime_09 <- function(path, cell_id, pseudotime, label = NULL, extra = NULL) {
  df <- data.frame(
    cell_id = as.character(cell_id),
    pseudotime = suppressWarnings(as.numeric(pseudotime)),
    stringsAsFactors = FALSE
  )
  if (!is.null(label)) {
    df$label <- as.character(label)
  }
  if (!is.null(extra) && nrow(extra) == nrow(df)) {
    df <- cbind(df, extra)
  }
  write.csv(df, path, row.names = FALSE)
  df
}

trajectory_plot_pseudotime_09 <- function(seu, pseudotime, path, title = "") {
  emb <- trajectory_umap_09(seu)
  ensure_dir(dirname(path))
  if (is.null(emb)) {
    plot_obj <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No UMAP") +
      ggplot2::theme_void()
  } else {
    plot_df <- data.frame(
      UMAP_1 = emb[, 1],
      UMAP_2 = emb[, 2],
      pseudotime = suppressWarnings(as.numeric(pseudotime[rownames(emb)])),
      stringsAsFactors = FALSE
    )
    plot_obj <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
      ggplot2::geom_point(size = 0.25, alpha = 0.8) +
      ggplot2::scale_color_viridis_c(na.value = "grey85") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(title = title, x = "UMAP 1", y = "UMAP 2", color = "pseudotime")
  }
  save_plot_local(plot_obj, path, width = 6, height = 4.8)
  invisible(path)
}

trajectory_export_python_input_09 <- function(seu, out_dir, label_var, root_label = "", terminal_label = "") {
  ensure_dir(out_dir)
  counts <- trajectory_counts_matrix_09(seu)
  Matrix::writeMM(counts, file.path(out_dir, "counts.mtx"))
  writeLines(rownames(counts), file.path(out_dir, "features.tsv"), useBytes = TRUE)
  writeLines(colnames(counts), file.path(out_dir, "barcodes.tsv"), useBytes = TRUE)
  meta <- seu@meta.data
  meta$cell_id <- rownames(meta)
  write_tsv_local(meta, file.path(out_dir, "metadata.tsv"))
  emb <- trajectory_umap_09(seu)
  if (!is.null(emb)) {
    umap <- data.frame(cell_id = rownames(emb), UMAP_1 = emb[, 1], UMAP_2 = emb[, 2], stringsAsFactors = FALSE)
    write_tsv_local(umap, file.path(out_dir, "umap.tsv"))
  }
  root_cells <- trajectory_root_cells_09(seu, label_var, root_label)
  terminal_cells <- trajectory_terminal_cells_09(seu, label_var, terminal_label)
  writeLines(root_cells, file.path(out_dir, "root_cells.txt"), useBytes = TRUE)
  writeLines(terminal_cells, file.path(out_dir, "terminal_cells.txt"), useBytes = TRUE)
  list(
    counts_mtx = file.path(out_dir, "counts.mtx"),
    features_tsv = file.path(out_dir, "features.tsv"),
    barcodes_tsv = file.path(out_dir, "barcodes.tsv"),
    metadata_tsv = file.path(out_dir, "metadata.tsv"),
    umap_tsv = file.path(out_dir, "umap.tsv"),
    root_cells_txt = file.path(out_dir, "root_cells.txt"),
    terminal_cells_txt = file.path(out_dir, "terminal_cells.txt")
  )
}

trajectory_python_exe_09 <- function() {
  prefix <- normalize_scalar_value(Sys.getenv("VELOCITY_ENV_PREFIX", ""))
  if (nzchar(prefix)) {
    candidate <- file.path(prefix, "bin", "python")
    if (file.exists(candidate)) {
      return(candidate)
    }
  }
  exe <- normalize_scalar_value(unname(Sys.which("python3")[[1]]))
  if (nzchar(exe) && file.exists(exe)) exe else ""
}

trajectory_run_python_bridge_09 <- function(script_path, jobs_tsv, index_tsv) {
  python <- trajectory_python_exe_09()
  if (!nzchar(python) || !file.exists(python)) {
    return(FALSE)
  }
  message(sprintf(
    "Running Python trajectory bridge: %s %s --jobs %s --index %s",
    shQuote(python),
    shQuote(script_path),
    shQuote(jobs_tsv),
    shQuote(index_tsv)
  ))
  status <- system2(
    python,
    args = c(script_path, "--jobs", jobs_tsv, "--index", index_tsv),
    stdout = "",
    stderr = ""
  )
  exit_status <- suppressWarnings(as.integer(status))
  if (!is.finite(exit_status)) {
    exit_status <- 1L
  }
  if (exit_status != 0L) {
    warning(
      sprintf("Python trajectory bridge failed with exit status %s: %s", exit_status, script_path),
      call. = FALSE
    )
    return(FALSE)
  }
  TRUE
}

trajectory_method_manifest_09 <- function(cfg, manifest_path, module_name, outputs, inputs, depends_on = list()) {
  if (file.exists(manifest_path)) {
    unlink(manifest_path)
  }
  write_manifest_local(
    manifest_path = manifest_path,
    new_outputs = outputs,
    module_name = module_name,
    base_dir = cfg$project_root,
    inputs = inputs,
    version = cfg$module_version,
    depends_on = depends_on
  )
}
