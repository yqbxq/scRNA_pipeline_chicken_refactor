st07_empty_df <- function(cols) {
  as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols), stringsAsFactors = FALSE)
}

st07_scalar <- function(x, default = "") {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]])) {
    return(default)
  }
  value <- trimws(as.character(x[[1]]))
  if (!nzchar(value)) default else value
}

st07_write_tsv <- function(df, path) {
  ensure_dir(dirname(path))
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

st07_read_tsv <- function(path) {
  if (!nzchar(path) || !file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[!grepl("^\\s*#", lines)]
  if (length(lines) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  read.delim(text = paste(lines, collapse = "\n"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
}

st07_abs_path <- function(path, cfg) {
  path <- st07_scalar(path)
  if (!nzchar(path)) {
    return("")
  }
  if (grepl("^/", path)) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(cfg$project_root, path), winslash = "/", mustWork = FALSE)
  }
}

st07_manifest_output <- function(manifest_path, keys, cfg) {
  if (!file.exists(manifest_path) || !requireNamespace("jsonlite", quietly = TRUE)) {
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

st07_json_escape <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\"', x)
  x <- gsub("\n", "\\\\n", x, fixed = TRUE)
  x
}

st07_simple_manifest_value <- function(x, indent = "    ") {
  if (is.null(x)) {
    return("null")
  }
  if (is.list(x) && !is.data.frame(x)) {
    names_x <- names(x)
    if (is.null(names_x)) {
      names_x <- rep("", length(x))
    }
    parts <- character()
    for (idx in seq_along(x)) {
      name <- names_x[[idx]]
      if (!nzchar(name)) {
        next
      }
      parts <- c(parts, sprintf('%s"%s": %s', indent, st07_json_escape(name), st07_simple_manifest_value(x[[idx]], paste0(indent, "  "))))
    }
    return(paste0("{\n", paste(parts, collapse = ",\n"), "\n", sub("  $", "", indent), "}"))
  }
  if (is.numeric(x) || is.integer(x)) {
    return(as.character(x[[1]]))
  }
  sprintf('"%s"', st07_json_escape(x[[1]]))
}

st07_write_manifest_local <- function(manifest_path, new_outputs, module_name, base_dir, inputs = list(), version = "1.0", depends_on = list()) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(write_manifest_local(manifest_path, new_outputs, module_name = module_name, base_dir = base_dir, inputs = inputs, version = version, depends_on = depends_on))
  }
  ensure_dir(dirname(manifest_path))
  manifest <- list(
    module = module_name,
    version = version,
    timestamp = timestamp_now(),
    base_dir = base_dir,
    inputs = inputs,
    outputs = new_outputs,
    depends_on = depends_on
  )
  json <- st07_simple_manifest_value(manifest, "  ")
  writeLines(json, manifest_path, useBytes = TRUE)
  invisible(manifest_path)
}

st07_deconv_cols <- function() {
  c(
    "deconv_id", "section", "tool", "status", "reason", "panorama_input_rds", "n_spots",
    "n_celltypes", "runtime_sec", "proportion_tsv", "proportion_wide_tsv",
    "spot_metadata_tsv", "summary_tsv", "method_object_rds", "method_version"
  )
}

st07_deconv_row <- function(pair, section = "", tool, status, reason = "", n_spots = 0L,
                            n_celltypes = 0L, runtime_sec = NA_real_, outputs = list(),
                            method_version = "1.0", panorama_input_rds = "") {
  data.frame(
    deconv_id = st07_scalar(pair$deconv_id, "default_deconv"),
    section = section,
    tool = tool,
    status = status,
    reason = reason,
    panorama_input_rds = panorama_input_rds,
    n_spots = n_spots,
    n_celltypes = n_celltypes,
    runtime_sec = runtime_sec,
    proportion_tsv = outputs$proportion_tsv %||% "",
    proportion_wide_tsv = outputs$proportion_wide_tsv %||% "",
    spot_metadata_tsv = outputs$spot_metadata_tsv %||% "",
    summary_tsv = outputs$summary_tsv %||% "",
    method_object_rds = outputs$method_object_rds %||% "",
    method_version = method_version,
    stringsAsFactors = FALSE
  )
}

st07_status_summary <- function(df) {
  if (nrow(df) == 0 || !"status" %in% colnames(df)) {
    return(st07_empty_df(c("status", "n")))
  }
  out <- as.data.frame(table(status = df$status), stringsAsFactors = FALSE)
  colnames(out) <- c("status", "n")
  out
}

read_deconv_pairs_st <- function(cfg) {
  pairs <- st07_read_tsv(cfg$deconv_pairs_sheet)
  cols <- c("deconv_id", "source_question_id", "st_scope", "reference_scope", "section_filter", "condition_split_var", "condition_split_values", "tool", "enabled", "notes")
  for (col in cols) {
    if (!col %in% colnames(pairs)) {
      pairs[[col]] <- character(nrow(pairs))
    }
  }
  pairs[, cols, drop = FALSE]
}

deconv_pair_enabled_st <- function(pair) {
  tolower(st07_scalar(pair$enabled, "yes")) %in% c("yes", "true", "1", "on")
}

deconv_pair_allows_tool_st <- function(pair, tool) {
  raw <- tolower(st07_scalar(pair$tool, "all"))
  tools <- trimws(unlist(strsplit(raw, "[,;]+", perl = TRUE), use.names = FALSE))
  "all" %in% tools || tool %in% tools
}

load_spatial_reference_inventory_st <- function(cfg) {
  inv <- st07_read_tsv(cfg$spatial_reference_inventory_file)
  cols <- c("selected_reference", "selected_annotation_col", "reference_sha256")
  for (col in cols) {
    if (!col %in% colnames(inv)) {
      inv[[col]] <- character(nrow(inv))
    }
  }
  if (nrow(inv) == 0) {
    return(list(status = "skipped_no_reference", reason = "spatial_reference_inventory.tsv has no rows", reference_path = "", annotation_col = "", sha256 = ""))
  }
  row <- inv[1, , drop = FALSE]
  reference_path <- st07_abs_path(row$selected_reference, cfg)
  annotation_col <- st07_scalar(row$selected_annotation_col)
  if (!nzchar(reference_path) || !file.exists(reference_path)) {
    return(list(status = "skipped_no_reference", reason = sprintf("selected_reference does not exist: %s", reference_path), reference_path = reference_path, annotation_col = annotation_col, sha256 = st07_scalar(row$reference_sha256)))
  }
  if (!nzchar(annotation_col)) {
    return(list(status = "failed_reference_load", reason = "selected_annotation_col is empty", reference_path = reference_path, annotation_col = annotation_col, sha256 = st07_scalar(row$reference_sha256)))
  }
  list(status = "ok", reason = "", reference_path = reference_path, annotation_col = annotation_col, sha256 = st07_scalar(row$reference_sha256))
}

verify_reference_hash_st <- function(reference_info) {
  expected <- st07_scalar(reference_info$sha256)
  if (!nzchar(expected)) {
    return(list(status = "ok", reason = "reference_sha256 empty; hash verification skipped"))
  }
  actual <- ""
  if (requireNamespace("digest", quietly = TRUE)) {
    actual <- tryCatch(digest::digest(file = reference_info$reference_path, algo = "sha256", serialize = FALSE), error = function(e) "")
  }
  if (!nzchar(actual)) {
    actual <- tryCatch(system2("sha256sum", reference_info$reference_path, stdout = TRUE, stderr = TRUE), error = function(e) "")
    actual <- sub("\\s+.*$", "", actual[[1]] %||% "")
  }
  if (!nzchar(actual)) {
    actual <- tryCatch(system2("shasum", c("-a", "256", reference_info$reference_path), stdout = TRUE, stderr = TRUE), error = function(e) "")
    actual <- sub("\\s+.*$", "", actual[[1]] %||% "")
  }
  if (!nzchar(actual)) {
    return(list(status = "failed_reference_load", reason = "no sha256 implementation returned a digest"))
  }
  if (!identical(tolower(actual), tolower(expected))) {
    return(list(status = "failed_reference_load", reason = sprintf("reference_sha256 mismatch: expected %s got %s", expected, actual)))
  }
  list(status = "ok", reason = "")
}

resolve_panorama_for_deconv_st <- function(cfg) {
  candidates <- c(
    st07_manifest_output(cfg$module_04c_subcluster_eda_manifest_path, c("spatial_panorama_subannotated", "subannotated_object", "panorama_rds"), cfg),
    st07_manifest_output(cfg$module_03a_region_annotation_eda_manifest_path, c("spatial_panorama_annotated", "annotated_object", "panorama_rds"), cfg),
    cfg$spatial_panorama_subannotated_rds,
    cfg$spatial_panorama_annotated_rds,
    cfg$spatial_panorama_clustered_rds
  )
  hit <- candidates[nzchar(candidates) & file.exists(candidates)][1]
  if (is.na(hit)) "" else hit
}

write_empty_deconv_outputs_st <- function(cfg, table_dir, pair, section, tool, status, reason) {
  out_dir <- file.path(table_dir, spatial_safe_id(st07_scalar(pair$deconv_id, "default_deconv")), spatial_safe_id(section))
  ensure_dir(out_dir)
  proportion_tsv <- file.path(out_dir, "spot_celltype_proportions.tsv")
  proportion_wide_tsv <- file.path(out_dir, "spot_celltype_proportions_wide.tsv")
  spot_metadata_tsv <- file.path(out_dir, "spot_metadata.tsv")
  summary_tsv <- file.path(out_dir, "method_summary.tsv")
  st07_write_tsv(st07_empty_df(c("spot_id", "cell_type", "proportion")), proportion_tsv)
  st07_write_tsv(st07_empty_df(c("spot_id")), proportion_wide_tsv)
  st07_write_tsv(st07_empty_df(c("spot_id", "dominant_celltype", "mixing_entropy")), spot_metadata_tsv)
  st07_write_tsv(data.frame(tool = tool, status = status, reason = reason, n_spots = 0L, n_celltypes = 0L, stringsAsFactors = FALSE), summary_tsv)
  list(proportion_tsv = proportion_tsv, proportion_wide_tsv = proportion_wide_tsv, spot_metadata_tsv = spot_metadata_tsv, summary_tsv = summary_tsv)
}

st07_assay_counts <- function(obj, preferred = c("Spatial", "RNA")) {
  assays <- intersect(preferred, names(obj@assays))
  if (length(assays) == 0) {
    assays <- Seurat::DefaultAssay(obj)
  }
  assay <- assays[[1]]
  counts <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = "counts")
  )
  list(assay = assay, counts = counts)
}

st07_ncount <- function(obj, assay, counts) {
  col <- paste0("nCount_", assay)
  if (col %in% colnames(obj@meta.data)) {
    out <- suppressWarnings(as.numeric(obj@meta.data[[col]]))
    names(out) <- rownames(obj@meta.data)
    return(out)
  }
  Matrix::colSums(counts)
}

st07_section_cells <- function(obj, section_filter = "all") {
  meta <- obj@meta.data
  all_cells <- rownames(meta)
  section_col <- if ("section_id" %in% colnames(meta)) "section_id" else if ("spatial_section_id" %in% colnames(meta)) "spatial_section_id" else if ("sample_id" %in% colnames(meta)) "sample_id" else ""
  requested <- st07_scalar(section_filter, "all")
  if (requested %in% c("*", "all", "ALL", "__ALL__") || !nzchar(section_col)) {
    sections <- if (nzchar(section_col)) split(all_cells, as.character(meta[[section_col]])) else list(all = all_cells)
  } else {
    keep <- all_cells[as.character(meta[[section_col]]) == requested]
    sections <- stats::setNames(list(keep), requested)
  }
  sections <- sections[vapply(sections, length, integer(1)) > 0]
  if (length(sections) == 0) {
    sections <- stats::setNames(list(character(0)), requested)
  }
  sections
}

st07_spatial_coords <- function(obj) {
  meta <- obj@meta.data
  x_col <- c("x", "X", "imagecol", "pxl_col_in_fullres", "array_col", "col")[c("x", "X", "imagecol", "pxl_col_in_fullres", "array_col", "col") %in% colnames(meta)][1]
  y_col <- c("y", "Y", "imagerow", "pxl_row_in_fullres", "array_row", "row")[c("y", "Y", "imagerow", "pxl_row_in_fullres", "array_row", "row") %in% colnames(meta)][1]
  if (!is.na(x_col) && !is.na(y_col)) {
    coords <- data.frame(x = suppressWarnings(as.numeric(meta[[x_col]])), y = suppressWarnings(as.numeric(meta[[y_col]])), row.names = rownames(meta))
    return(coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE])
  }
  coords <- tryCatch(Seurat::GetTissueCoordinates(obj), error = function(e) data.frame())
  if (nrow(coords) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  if ("cell" %in% colnames(coords)) {
    rownames(coords) <- as.character(coords$cell)
  } else if ("barcode" %in% colnames(coords)) {
    rownames(coords) <- as.character(coords$barcode)
  }
  x_col <- c("x", "imagecol", "pxl_col_in_fullres", "array_col", "col")[c("x", "imagecol", "pxl_col_in_fullres", "array_col", "col") %in% colnames(coords)][1]
  y_col <- c("y", "imagerow", "pxl_row_in_fullres", "array_row", "row")[c("y", "imagerow", "pxl_row_in_fullres", "array_row", "row") %in% colnames(coords)][1]
  if (is.na(x_col) || is.na(y_col)) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  out <- data.frame(x = suppressWarnings(as.numeric(coords[[x_col]])), y = suppressWarnings(as.numeric(coords[[y_col]])), row.names = rownames(coords))
  out[is.finite(out$x) & is.finite(out$y), , drop = FALSE]
}

st07_normalize_proportions <- function(mat) {
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- 0
  mat[mat < 0] <- 0
  rs <- rowSums(mat)
  keep <- rs > 0
  mat[keep, ] <- sweep(mat[keep, , drop = FALSE], 1, rs[keep], "/")
  mat
}

st07_write_deconv_result_outputs <- function(cfg, table_dir, pair, section, tool, prop_mat, method_object = NULL, runtime_sec = NA_real_) {
  out_dir <- file.path(table_dir, spatial_safe_id(st07_scalar(pair$deconv_id, "default_deconv")), spatial_safe_id(section))
  ensure_dir(out_dir)
  prop_mat <- st07_normalize_proportions(prop_mat)
  prop_mat <- prop_mat[rowSums(prop_mat) > 0, , drop = FALSE]
  if (nrow(prop_mat) == 0 || ncol(prop_mat) == 0) {
    stop("deconvolution returned no nonzero spot-celltype proportions", call. = FALSE)
  }
  proportion_tsv <- file.path(out_dir, "spot_celltype_proportions.tsv")
  proportion_wide_tsv <- file.path(out_dir, "spot_celltype_proportions_wide.tsv")
  spot_metadata_tsv <- file.path(out_dir, "spot_metadata.tsv")
  summary_tsv <- file.path(out_dir, "method_summary.tsv")
  method_object_rds <- file.path(out_dir, "method_object.rds")

  wide <- data.frame(spot_id = rownames(prop_mat), prop_mat, check.names = FALSE, stringsAsFactors = FALSE)
  long <- do.call(rbind, lapply(seq_len(nrow(prop_mat)), function(i) {
    data.frame(spot_id = rownames(prop_mat)[[i]], cell_type = colnames(prop_mat), proportion = as.numeric(prop_mat[i, ]), stringsAsFactors = FALSE)
  }))
  entropy <- -rowSums(ifelse(prop_mat > 0, prop_mat * log(prop_mat), 0))
  dominant <- colnames(prop_mat)[max.col(prop_mat, ties.method = "first")]
  spot_meta <- data.frame(spot_id = rownames(prop_mat), dominant_celltype = dominant, mixing_entropy = entropy, stringsAsFactors = FALSE)
  summary <- data.frame(tool = tool, status = "ok", reason = "", n_spots = nrow(prop_mat), n_celltypes = ncol(prop_mat), runtime_sec = runtime_sec, stringsAsFactors = FALSE)

  st07_write_tsv(long, proportion_tsv)
  st07_write_tsv(wide, proportion_wide_tsv)
  st07_write_tsv(spot_meta, spot_metadata_tsv)
  st07_write_tsv(summary, summary_tsv)
  if (!is.null(method_object)) {
    saveRDS(method_object, method_object_rds)
  } else {
    method_object_rds <- ""
  }
  list(proportion_tsv = proportion_tsv, proportion_wide_tsv = proportion_wide_tsv, spot_metadata_tsv = spot_metadata_tsv, summary_tsv = summary_tsv, method_object_rds = method_object_rds)
}

st07_prepare_deconv_inputs <- function(panorama, reference, annotation_col, cells) {
  st_section <- if (length(cells) > 0) subset(panorama, cells = cells) else panorama[, FALSE]
  if (ncol(st_section) == 0) {
    stop("selected section contains 0 spots", call. = FALSE)
  }
  if (!annotation_col %in% colnames(reference@meta.data)) {
    stop(sprintf("reference metadata missing annotation column: %s", annotation_col), call. = FALSE)
  }
  ref_labels <- as.factor(reference@meta.data[[annotation_col]])
  keep_ref <- !is.na(ref_labels) & nzchar(as.character(ref_labels))
  if (sum(keep_ref) < 2 || length(unique(ref_labels[keep_ref])) < 2) {
    stop("reference needs at least two cells and two cell types", call. = FALSE)
  }
  list(st_section = st_section, reference = subset(reference, cells = rownames(reference@meta.data)[keep_ref]), ref_labels = droplevels(ref_labels[keep_ref]))
}

st07_run_rctd_method <- function(st_section, reference, ref_labels, cfg) {
  st_counts <- st07_assay_counts(st_section, c("Spatial", "SCT", "RNA"))
  ref_counts <- st07_assay_counts(reference, c("RNA", "SCT"))
  coords <- st07_spatial_coords(st_section)
  common_spots <- intersect(colnames(st_counts$counts), rownames(coords))
  if (length(common_spots) < 2) {
    stop("RCTD requires at least two spots with spatial coordinates", call. = FALSE)
  }
  st_counts_mat <- st_counts$counts[, common_spots, drop = FALSE]
  coords <- coords[common_spots, c("x", "y"), drop = FALSE]
  ref_numi <- st07_ncount(reference, ref_counts$assay, ref_counts$counts)
  st_numi <- st07_ncount(st_section, st_counts$assay, st_counts$counts)
  reference_obj <- spacexr::Reference(ref_counts$counts, ref_labels, nUMI = ref_numi[colnames(ref_counts$counts)])
  spatial_rna <- spacexr::SpatialRNA(coords, st_counts_mat, nUMI = st_numi[common_spots])
  rctd <- spacexr::create.RCTD(spatial_rna, reference_obj, max_cores = cfg$spatial_rctd_cores, CELL_MIN_INSTANCE = 25)
  rctd <- spacexr::run.RCTD(rctd, doublet_mode = cfg$spatial_rctd_doublet_mode)
  weights <- rctd@results$weights
  list(prop = as.matrix(weights), object = rctd)
}

st07_run_transfer_method <- function(st_section, reference, ref_labels, cfg) {
  if (!"pca" %in% names(st_section@reductions)) {
    st_section <- Seurat::NormalizeData(st_section, verbose = FALSE)
    st_section <- Seurat::FindVariableFeatures(st_section, verbose = FALSE)
    st_section <- Seurat::ScaleData(st_section, verbose = FALSE)
    st_section <- Seurat::RunPCA(st_section, npcs = 30, verbose = FALSE)
  }
  if (!"pca" %in% names(reference@reductions)) {
    reference <- Seurat::NormalizeData(reference, verbose = FALSE)
    reference <- Seurat::FindVariableFeatures(reference, verbose = FALSE)
    reference <- Seurat::ScaleData(reference, verbose = FALSE)
    reference <- Seurat::RunPCA(reference, npcs = 30, verbose = FALSE)
  }
  reference$.st07_celltype <- ref_labels
  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = st_section,
    normalization.method = "LogNormalize",
    dims = 1:30,
    reduction = "pcaproject",
    verbose = FALSE
  )
  predictions <- Seurat::TransferData(
    anchorset = anchors,
    refdata = reference$.st07_celltype,
    weight.reduction = st_section[["pca"]],
    dims = 1:30,
    verbose = FALSE
  )
  if (is.data.frame(predictions)) {
    score_cols <- grep("^prediction\\.score\\.", colnames(predictions), value = TRUE)
    if (length(score_cols) == 0) {
      stop("TransferData returned no prediction.score.* columns", call. = FALSE)
    }
    pred_mat <- as.matrix(predictions[, score_cols, drop = FALSE])
    colnames(pred_mat) <- sub("^prediction\\.score\\.", "", colnames(pred_mat))
  } else {
    pred_mat <- as.matrix(predictions@data)
  }
  if (ncol(pred_mat) == ncol(st_section) && nrow(pred_mat) != ncol(st_section)) {
    pred_mat <- t(pred_mat)
  }
  if (is.null(rownames(pred_mat)) || !all(rownames(pred_mat) %in% colnames(st_section))) {
    rownames(pred_mat) <- colnames(st_section)[seq_len(nrow(pred_mat))]
  }
  list(prop = pred_mat, object = anchors)
}

st07_run_card_method <- function(st_section, reference, ref_labels, cfg) {
  st_counts <- st07_assay_counts(st_section, c("Spatial", "SCT", "RNA"))
  ref_counts <- st07_assay_counts(reference, c("RNA", "SCT"))
  coords <- st07_spatial_coords(st_section)
  common_spots <- intersect(colnames(st_counts$counts), rownames(coords))
  if (length(common_spots) < 2) {
    stop("CARD requires at least two spots with spatial coordinates", call. = FALSE)
  }
  sc_meta <- data.frame(cellID = colnames(ref_counts$counts), cellType = as.character(ref_labels), stringsAsFactors = FALSE)
  rownames(sc_meta) <- sc_meta$cellID
  card_obj <- CARD::createCARDObject(
    sc_count = ref_counts$counts,
    sc_meta = sc_meta,
    spatial_count = st_counts$counts[, common_spots, drop = FALSE],
    spatial_location = coords[common_spots, c("x", "y"), drop = FALSE],
    ct.varname = "cellType",
    ct.select = unique(sc_meta$cellType)
  )
  card_obj <- CARD::CARD_deconvolution(CARD_object = card_obj)
  list(prop = as.matrix(card_obj@Proportion_CARD), object = card_obj)
}

st07_run_deconv_tool <- function(tool, st_section, reference, ref_labels, cfg) {
  if (identical(tool, "rctd")) {
    return(st07_run_rctd_method(st_section, reference, ref_labels, cfg))
  }
  if (identical(tool, "transfer")) {
    return(st07_run_transfer_method(st_section, reference, ref_labels, cfg))
  }
  if (identical(tool, "card")) {
    return(st07_run_card_method(st_section, reference, ref_labels, cfg))
  }
  stop(sprintf("tool %s is not implemented in the R deconvolution runner", tool), call. = FALSE)
}

run_deconv_entry_st <- function(cfg, tool, package_name, table_dir, fail_status) {
  prepare_dirs_spatial(cfg)
  pairs <- read_deconv_pairs_st(cfg)
  rows <- list()
  if (nrow(pairs) == 0) {
    pair <- data.frame(deconv_id = "", enabled = "yes", tool = tool, stringsAsFactors = FALSE)
    outputs <- write_empty_deconv_outputs_st(cfg, table_dir, pair, "all", tool, "skipped_no_deconv_pairs", "deconv_pairs.tsv has no enabled rows")
    return(st07_deconv_row(pair, "all", tool, "skipped_no_deconv_pairs", "deconv_pairs.tsv has no enabled rows", outputs = outputs))
  }
  reference <- load_spatial_reference_inventory_st(cfg)
  hash <- if (identical(reference$status, "ok")) verify_reference_hash_st(reference) else list(status = reference$status, reason = reference$reason)
  panorama <- resolve_panorama_for_deconv_st(cfg)
  panorama_obj <- NULL
  reference_obj <- NULL
  if (identical(reference$status, "ok") && identical(hash$status, "ok") && nzchar(panorama) && file.exists(panorama) && requireNamespace(package_name, quietly = TRUE) && requireNamespace("Seurat", quietly = TRUE) && requireNamespace("Matrix", quietly = TRUE)) {
    panorama_obj <- readRDS(panorama)
    reference_obj <- readRDS(reference$reference_path)
  }
  for (idx in seq_len(nrow(pairs))) {
    pair <- pairs[idx, , drop = FALSE]
    section <- st07_scalar(pair$section_filter, "all")
    if (!deconv_pair_enabled_st(pair) || !deconv_pair_allows_tool_st(pair, tool)) {
      outputs <- write_empty_deconv_outputs_st(cfg, table_dir, pair, section, tool, "skipped_disabled", "pair disabled or tool filter does not include this method")
      rows[[length(rows) + 1L]] <- st07_deconv_row(pair, section, tool, "skipped_disabled", "pair disabled or tool filter does not include this method", outputs = outputs, panorama_input_rds = panorama)
      next
    }
    if (!identical(reference$status, "ok")) {
      outputs <- write_empty_deconv_outputs_st(cfg, table_dir, pair, section, tool, reference$status, reference$reason)
      rows[[length(rows) + 1L]] <- st07_deconv_row(pair, section, tool, reference$status, reference$reason, outputs = outputs, panorama_input_rds = panorama)
      next
    }
    if (!identical(hash$status, "ok")) {
      outputs <- write_empty_deconv_outputs_st(cfg, table_dir, pair, section, tool, "failed_reference_load", hash$reason)
      rows[[length(rows) + 1L]] <- st07_deconv_row(pair, section, tool, "failed_reference_load", hash$reason, outputs = outputs, panorama_input_rds = panorama)
      next
    }
    if (!requireNamespace(package_name, quietly = TRUE)) {
      outputs <- write_empty_deconv_outputs_st(cfg, table_dir, pair, section, tool, "skipped_no_packages", sprintf("%s is not installed", package_name))
      rows[[length(rows) + 1L]] <- st07_deconv_row(pair, section, tool, "skipped_no_packages", sprintf("%s is not installed", package_name), outputs = outputs, panorama_input_rds = panorama)
      next
    }
    if (!nzchar(panorama) || !file.exists(panorama) || !requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
      reason <- sprintf("panorama/Seurat/Matrix unavailable for %s method execution", tool)
      outputs <- write_empty_deconv_outputs_st(cfg, table_dir, pair, section, tool, fail_status, reason)
      rows[[length(rows) + 1L]] <- st07_deconv_row(pair, section, tool, fail_status, reason, outputs = outputs, panorama_input_rds = panorama)
      next
    }
    sections <- st07_section_cells(panorama_obj, section)
    for (section_name in names(sections)) {
      run <- tryCatch({
        prepared <- st07_prepare_deconv_inputs(panorama_obj, reference_obj, reference$annotation_col, sections[[section_name]])
        timing <- system.time(result <- st07_run_deconv_tool(tool, prepared$st_section, prepared$reference, prepared$ref_labels, cfg))
        outputs <- st07_write_deconv_result_outputs(cfg, table_dir, pair, section_name, tool, result$prop, result$object, unname(timing[["elapsed"]]))
        list(status = "ok", reason = "", outputs = outputs, n_spots = nrow(result$prop), n_celltypes = ncol(result$prop), runtime_sec = unname(timing[["elapsed"]]))
      }, error = function(e) {
        reason <- conditionMessage(e)
        outputs <- write_empty_deconv_outputs_st(cfg, table_dir, pair, section_name, tool, fail_status, reason)
        list(status = fail_status, reason = reason, outputs = outputs, n_spots = 0L, n_celltypes = 0L, runtime_sec = NA_real_)
      })
      rows[[length(rows) + 1L]] <- st07_deconv_row(pair, section_name, tool, run$status, run$reason, run$n_spots, run$n_celltypes, run$runtime_sec, run$outputs, panorama_input_rds = panorama)
    }
  }
  do.call(rbind, rows)
}

run_deconv_method_scaffold_st <- run_deconv_entry_st

write_deconv_module_manifest_st <- function(cfg, manifest_df, manifest_tsv, report_path, module_name, manifest_path, output_key, depends_on = list()) {
  st07_write_tsv(manifest_df, manifest_tsv)
  report_lines <- c(
    sprintf("# %s", module_name),
    "",
    "## Status Summary",
    render_markdown_table_local(st07_status_summary(manifest_df)),
    "",
    "## Preview",
    render_markdown_table_local(utils::head(manifest_df, 50))
  )
  write_markdown_local(report_lines, report_path)
  st07_write_manifest_local(
    manifest_path = manifest_path,
    new_outputs = c(
      stats::setNames(list(build_output_entry(manifest_tsv, "tsv", module_name, "spatial deconvolution method manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df))), output_key),
      report = build_output_entry(report_path, "md", module_name, "spatial deconvolution report", base_dir = cfg$project_root)
    ),
    module_name = module_name,
    base_dir = cfg$project_root,
    inputs = list(spatial_reference_inventory = cfg$spatial_reference_inventory_file, deconv_pairs = cfg$deconv_pairs_sheet),
    version = cfg$module_07_version,
    depends_on = depends_on
  )
}

collect_deconv_method_manifests_st <- function(cfg) {
  specs <- list(
    rctd = list(path = cfg$module_07a_deconvolution_rctd_manifest_path, key = "rctd_manifest", fallback = file.path(cfg$spatial_rctd_table_dir, "rctd_manifest.tsv")),
    transfer = list(path = cfg$module_07b_deconvolution_transfer_manifest_path, key = "transfer_manifest", fallback = file.path(cfg$spatial_transfer_table_dir, "transfer_manifest.tsv")),
    card = list(path = cfg$module_07c_deconvolution_card_manifest_path, key = "card_manifest", fallback = file.path(cfg$spatial_card_table_dir, "card_manifest.tsv")),
    cell2location = list(path = cfg$module_07d_deconvolution_cell2location_manifest_path, key = "cell2location_manifest", fallback = file.path(cfg$spatial_c2l_table_dir, "cell2location_manifest.tsv"))
  )
  rows <- list()
  for (method in names(specs)) {
    manifest_tsv <- st07_manifest_output(specs[[method]]$path, specs[[method]]$key, cfg)
    if (!nzchar(manifest_tsv) && file.exists(specs[[method]]$fallback)) {
      manifest_tsv <- specs[[method]]$fallback
    }
    df <- st07_read_tsv(manifest_tsv)
    if (nrow(df) == 0) {
      next
    }
    df$method <- method
    df$method_manifest_tsv <- manifest_tsv
    rows[[length(rows) + 1L]] <- df
  }
  if (length(rows) == 0) {
    return(st07_empty_df(c(st07_deconv_cols(), "method", "method_manifest_tsv")))
  }
  cols <- unique(unlist(lapply(rows, colnames), use.names = FALSE))
  rows <- lapply(rows, function(df) {
    for (col in setdiff(cols, colnames(df))) {
      df[[col]] <- ""
    }
    df[, cols, drop = FALSE]
  })
  do.call(rbind, rows)
}
