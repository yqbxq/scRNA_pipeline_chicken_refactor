spatial_cell <- function(row, key, default = "") {
  if (is.null(row) || !key %in% colnames(row)) {
    return(default)
  }
  value <- row[[key]][[1]]
  if (is.na(value) || !nzchar(trimws(as.character(value)))) {
    return(default)
  }
  trimws(as.character(value))
}

spatial_bool <- function(value, default = TRUE) {
  value <- tolower(trimws(as.character(value)))
  if (!nzchar(value)) {
    return(default)
  }
  value %in% c("yes", "true", "1", "on", "auto")
}

spatial_first_existing <- function(paths) {
  paths <- paths[nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) "" else hit[[1]]
}

spatial_read_tsv <- function(path) {
  if (exists("read_tsv_optional", mode = "function")) {
    return(read_tsv_optional(path))
  }
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
}

spatial_write_tsv <- function(df, path) {
  ensure_dir(dirname(path))
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

spatial_read_matrix_market <- function(root) {
  matrix_path <- spatial_first_existing(file.path(root, c("counts.mtx", "counts.mtx.gz", "matrix.mtx", "matrix.mtx.gz")))
  feature_path <- spatial_first_existing(file.path(root, c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz")))
  barcode_path <- spatial_first_existing(file.path(root, c("barcodes.tsv", "barcodes.tsv.gz")))
  if (!nzchar(matrix_path) || !nzchar(feature_path) || !nzchar(barcode_path)) {
    stop(sprintf("generic spatial matrix is missing matrix/features/barcodes under %s", root), call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required to read generic spatial matrix inputs", call. = FALSE)
  }
  counts <- Matrix::readMM(matrix_path)
  features <- read.delim(feature_path, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  barcodes <- read.delim(barcode_path, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  feature_names <- if (ncol(features) >= 2) features[[2]] else features[[1]]
  feature_names <- make.unique(as.character(feature_names))
  rownames(counts) <- feature_names
  colnames(counts) <- make.unique(as.character(barcodes[[1]]))
  counts
}

spatial_read_coords <- function(root, cells) {
  coords_path <- spatial_first_existing(file.path(root, c("coords.csv", "coords.tsv", "tissue_positions.csv", "tissue_positions_list.csv")))
  if (!nzchar(coords_path)) {
    return(data.frame(row.names = cells, stringsAsFactors = FALSE))
  }
  sep <- if (grepl("\\.tsv$", coords_path)) "\t" else ","
  coords <- read.delim(coords_path, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"barcode" %in% colnames(coords)) {
    colnames(coords)[[1]] <- "barcode"
  }
  rownames(coords) <- make.unique(as.character(coords$barcode))
  coords <- coords[intersect(cells, rownames(coords)), , drop = FALSE]
  coords[cells, setdiff(colnames(coords), "barcode"), drop = FALSE]
}

contract_seurat_metadata <- function(obj) {
  meta <- tryCatch(obj@meta.data, error = function(e) NULL)
  if (!is.data.frame(meta)) {
    stop("loaded spatial object does not expose Seurat metadata", call. = FALSE)
  }
  required <- c("section_id", "chip_id", "orig.ident", "platform", "stage", "condition", "bundle_layout")
  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0) {
    stop(sprintf("spatial object metadata missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

resolve_spatial_data_dir <- function(samples_row, sections_row) {
  candidates <- c(
    spatial_cell(samples_row, "standardized_outs"),
    spatial_cell(samples_row, "spatial_standardized_outs"),
    spatial_cell(sections_row, "bundle_root"),
    spatial_cell(samples_row, "source_path")
  )
  hit <- spatial_first_existing(candidates)
  if (!nzchar(hit)) {
    stop(sprintf("no spatial bundle directory exists for sample_id=%s", spatial_cell(samples_row, "sample_id")), call. = FALSE)
  }
  normalizePath(hit, winslash = "/", mustWork = FALSE)
}

load_generic_spatial_object <- function(data_dir, samples_row, sections_row) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required to create generic spatial objects", call. = FALSE)
  }
  counts <- spatial_read_matrix_market(data_dir)
  obj <- Seurat::CreateSeuratObject(counts = counts, assay = "Spatial", project = spatial_cell(samples_row, "sample_id", "spatial"))
  coords <- spatial_read_coords(data_dir, colnames(obj))
  if (nrow(coords) > 0) {
    obj <- Seurat::AddMetaData(obj, coords)
  }
  obj
}

load_visium_spatial_object <- function(data_dir, samples_row, sections_row) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required to load Visium spatial objects", call. = FALSE)
  }
  section_id <- spatial_cell(samples_row, "section_id", spatial_cell(sections_row, "section_id", "slice1"))
  obj <- tryCatch(
    Seurat::Load10X_Spatial(data.dir = data_dir, assay = "Spatial", slice = section_id, filter.matrix = TRUE),
    error = function(e) {
      Seurat::Load10X_Spatial(data.dir = data_dir, assay = "Spatial", filter.matrix = TRUE)
    }
  )
  image_names <- tryCatch(names(obj@images), error = function(e) character(0))
  if (length(image_names) > 0 && !section_id %in% image_names) {
    names(obj@images)[[1]] <- section_id
  }
  obj
}

load_spatial_object <- function(samples_row, sections_row, intake_contract_row = NULL) {
  layout <- tolower(spatial_cell(samples_row, "bundle_layout", spatial_cell(intake_contract_row, "bundle_layout", "outs_visium")))
  data_dir <- resolve_spatial_data_dir(samples_row, sections_row)
  sample_id <- spatial_cell(samples_row, "sample_id", spatial_cell(sections_row, "section_id", "spatial_sample"))
  section_id <- spatial_cell(samples_row, "section_id", spatial_cell(sections_row, "section_id", sample_id))
  chip_id <- spatial_cell(samples_row, "chip_id", spatial_cell(sections_row, "chip_id", section_id))
  platform <- tolower(spatial_cell(samples_row, "platform", spatial_cell(sections_row, "platform", "visium")))

  if (layout %in% c("saw_bin50", "saw_bin100", "saw_cellbin", "stomics_native")) {
    stop("intake deferred to M5 for SAW/STOmics native spatial bundles", call. = FALSE)
  }

  obj <- if (layout %in% c("outs_visium", "visium")) {
    load_visium_spatial_object(data_dir, samples_row, sections_row)
  } else if (layout %in% c("generic_spatial_matrix", "spatial_matrix", "generic")) {
    load_generic_spatial_object(data_dir, samples_row, sections_row)
  } else if (layout %in% c("anndata_h5ad", "h5ad")) {
    stop("anndata_h5ad spatial intake requires the future reticulate/SeuratDisk bridge", call. = FALSE)
  } else {
    stop(sprintf("unsupported spatial bundle_layout=%s", layout), call. = FALSE)
  }

  obj$section_id <- section_id
  obj$chip_id <- chip_id
  obj$orig.ident <- sample_id
  obj$sample_id <- sample_id
  obj$platform <- platform
  obj$stage <- spatial_cell(samples_row, "timepoint")
  obj$condition <- spatial_cell(samples_row, "condition", spatial_cell(sections_row, "condition"))
  obj$bundle_layout <- layout
  obj@misc$spatial_source <- list(data_dir = data_dir, image_path = spatial_cell(samples_row, "image_path"))
  contract_seurat_metadata(obj)
  obj
}

manual_chicken_mito_genes <- function() {
  c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
}

read_manual_mito_override <- function(manual_override) {
  manual_override <- trimws(as.character(manual_override %||% ""))
  if (!nzchar(manual_override)) {
    return(character(0))
  }
  if (file.exists(manual_override)) {
    values <- readLines(manual_override, warn = FALSE, encoding = "UTF-8")
  } else {
    values <- unlist(strsplit(manual_override, "[,;[:space:]]+", perl = TRUE), use.names = FALSE)
  }
  values <- trimws(values)
  unique(values[nzchar(values) & !startsWith(values, "#")])
}

mito_genes_from_gtf <- function(features, gtf_path) {
  if (!nzchar(gtf_path) || !file.exists(gtf_path)) {
    return(character(0))
  }
  lines <- readLines(gtf_path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[grepl("^(chrM|MT)\\b|^M\\t", lines, ignore.case = TRUE)]
  if (length(lines) == 0) {
    return(character(0))
  }
  attrs <- sub("^([^\\t]*\\t){8}", "", lines, perl = TRUE)
  gene_names <- sub('.*gene_name "([^"]+)".*', "\\1", attrs)
  gene_ids <- sub('.*gene_id "([^"]+)".*', "\\1", attrs)
  hits <- unique(c(gene_names, gene_ids))
  intersect(features, hits[nzchar(hits)])
}

detect_spatial_mito_features <- function(obj, gtf_path = NULL, manual_override = NULL) {
  features <- rownames(obj)
  order_raw <- Sys.getenv("SPATIAL_MITO_DETECTION_ORDER", unset = "prefix_upper,prefix_lower,gtf_chrM,manual")
  detection_order <- trimws(unlist(strsplit(order_raw, ",", fixed = TRUE), use.names = FALSE))
  manual_values <- read_manual_mito_override(manual_override)

  for (method in detection_order[nzchar(detection_order)]) {
    hits <- switch(
      method,
      prefix_upper = grep("^MT-", features, value = TRUE),
      prefix_lower = grep("^mt-", features, value = TRUE),
      gtf_chrM = mito_genes_from_gtf(features, gtf_path %||% ""),
      manual = {
        if (length(manual_values) > 0) {
          intersect(features, manual_values)
        } else {
          base <- manual_chicken_mito_genes()
          features[toupper(features) %in% toupper(base)]
        }
      },
      character(0)
    )
    hits <- unique(hits[nzchar(hits)])
    if (length(hits) > 0) {
      source <- if (identical(method, "manual") && length(manual_values) > 0) "manual_override" else if (identical(method, "manual")) "manual_chicken_default" else if (identical(method, "gtf_chrM")) "chrM_gtf" else method
      return(list(mito_genes = hits, source = source, n_features = length(hits), valid = TRUE))
    }
  }

  list(mito_genes = character(0), source = "none", n_features = 0L, valid = FALSE)
}

spatial_counts_matrix <- function(obj) {
  assay <- if ("Spatial" %in% names(obj@assays)) "Spatial" else Seurat::DefaultAssay(obj)
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, slot = "counts"),
    error = function(e) Seurat::GetAssayData(obj, assay = assay, layer = "counts")
  )
}

inject_spatial_mito_qc <- function(obj, mito_detection) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required for spatial mito QC injection", call. = FALSE)
  }
  counts <- spatial_counts_matrix(obj)
  mito_genes <- intersect(mito_detection$mito_genes, rownames(counts))
  if (length(mito_genes) == 0) {
    obj$percent.mito <- NA_real_
    obj$mito_qc_valid <- FALSE
  } else {
    total <- Matrix::colSums(counts)
    mito <- Matrix::colSums(counts[mito_genes, , drop = FALSE])
    obj$percent.mito <- ifelse(total > 0, 100 * mito / total, NA_real_)
    obj$mito_qc_valid <- TRUE
  }
  obj@misc$mito_detection <- mito_detection
  obj
}

compute_spot_qc_metrics <- function(obj) {
  meta <- obj@meta.data
  counts <- spatial_counts_matrix(obj)
  nfeature <- if ("nFeature_Spatial" %in% colnames(meta)) meta$nFeature_Spatial else Matrix::colSums(counts > 0)
  ncount <- if ("nCount_Spatial" %in% colnames(meta)) meta$nCount_Spatial else Matrix::colSums(counts)
  percent_mito <- if ("percent.mito" %in% colnames(meta)) meta$percent.mito else rep(NA_real_, nrow(meta))
  in_tissue <- if ("in_tissue" %in% colnames(meta)) meta$in_tissue else rep(TRUE, nrow(meta))
  row_value <- if ("array_row" %in% colnames(meta)) meta$array_row else if ("row" %in% colnames(meta)) meta$row else if ("imagerow" %in% colnames(meta)) meta$imagerow else NA_real_
  col_value <- if ("array_col" %in% colnames(meta)) meta$array_col else if ("col" %in% colnames(meta)) meta$col else if ("imagecol" %in% colnames(meta)) meta$imagecol else NA_real_
  data.frame(
    spot_id = rownames(meta),
    sample_id = if ("sample_id" %in% colnames(meta)) meta$sample_id else meta$orig.ident,
    section_id = meta$section_id,
    nFeature_Spatial = as.numeric(nfeature),
    nCount_Spatial = as.numeric(ncount),
    log10_nCount = log10(pmax(as.numeric(ncount), 1)),
    percent.mito = as.numeric(percent_mito),
    mito_qc_valid = if ("mito_qc_valid" %in% colnames(meta)) as.logical(meta$mito_qc_valid) else FALSE,
    in_tissue = as.logical(in_tissue),
    row = suppressWarnings(as.numeric(row_value)),
    col = suppressWarnings(as.numeric(col_value)),
    stringsAsFactors = FALSE
  )
}

summarize_load_metrics <- function(seurat_list) {
  rows <- lapply(seurat_list, function(obj) {
    metrics <- compute_spot_qc_metrics(obj)
    data.frame(
      sample_id = unique(metrics$sample_id)[[1]],
      section_id = unique(metrics$section_id)[[1]],
      spots = nrow(metrics),
      in_tissue_spots = sum(metrics$in_tissue %in% TRUE, na.rm = TRUE),
      median_nFeature = stats::median(metrics$nFeature_Spatial, na.rm = TRUE),
      median_nCount = stats::median(metrics$nCount_Spatial, na.rm = TRUE),
      median_percent_mito = stats::median(metrics$percent.mito, na.rm = TRUE),
      mito_qc_valid = any(metrics$mito_qc_valid %in% TRUE, na.rm = TRUE),
      mito_source = as.character(obj@misc$mito_detection$source %||% ""),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

spatial_eda_violin <- function(metrics, metric) {
  ggplot2::ggplot(metrics, ggplot2::aes(x = section_id, y = .data[[metric]], fill = section_id)) +
    ggplot2::geom_violin(scale = "width", na.rm = TRUE) +
    ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.2, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = NULL, y = metric)
}

spatial_eda_scatter <- function(metrics) {
  ggplot2::ggplot(metrics, ggplot2::aes(x = nCount_Spatial, y = nFeature_Spatial, color = percent.mito)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.75, na.rm = TRUE) +
    ggplot2::scale_x_log10() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = "nCount_Spatial", y = "nFeature_Spatial", color = "percent.mito")
}

spatial_eda_spatial_plot <- function(metrics, metric) {
  ggplot2::ggplot(metrics, ggplot2::aes(x = col, y = row, color = .data[[metric]])) +
    ggplot2::geom_point(size = 0.9, alpha = 0.85, na.rm = TRUE) +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::labs(color = metric)
}

spatial_numeric_or <- function(value, default) {
  if (length(value) == 0 || is.null(value)) {
    return(default)
  }
  value <- suppressWarnings(as.numeric(value[[1]]))
  if (is.na(value)) default else value
}

spatial_thresholds_table <- function(cfg) {
  thresholds <- spatial_read_tsv(cfg$spatial_qc_threshold_file)
  if (nrow(thresholds) == 0) {
    thresholds <- data.frame(
      section_id = "__DEFAULT__",
      qc_min_nfeature = cfg$qc_min_nfeature_default,
      qc_max_nfeature = "",
      qc_min_ncount = cfg$qc_min_ncount_default,
      qc_max_ncount = "",
      qc_max_mito_pct = cfg$qc_max_mito_pct_default,
      mito_set_override = "",
      spatial_aware_filter = "false",
      excessive_drop_threshold = cfg$excessive_drop_threshold,
      stringsAsFactors = FALSE
    )
  }
  for (col in c(
    "section_id", "qc_min_nfeature", "qc_max_nfeature", "qc_min_ncount",
    "qc_max_ncount", "qc_max_mito_pct", "mito_set_override",
    "spatial_aware_filter", "excessive_drop_threshold"
  )) {
    if (!col %in% colnames(thresholds)) {
      thresholds[[col]] <- ""
    }
  }
  thresholds
}

spatial_threshold_for_section <- function(thresholds, section_id, sample_id = "", cfg = NULL) {
  hit <- thresholds[thresholds$section_id == section_id, , drop = FALSE]
  if (nrow(hit) == 0 && nzchar(sample_id)) {
    hit <- thresholds[thresholds$section_id == sample_id, , drop = FALSE]
  }
  if (nrow(hit) == 0) {
    hit <- thresholds[thresholds$section_id == "__DEFAULT__", , drop = FALSE]
  }
  if (nrow(hit) == 0) {
    hit <- thresholds[1, , drop = FALSE]
  }
  cfg <- cfg %||% list(
    qc_min_nfeature_default = 200,
    qc_max_nfeature_default = Inf,
    qc_min_ncount_default = 500,
    qc_max_ncount_default = Inf,
    qc_max_mito_pct_default = 20,
    excessive_drop_threshold = 0.5
  )
  list(
    section_id = as.character(hit$section_id[[1]]),
    qc_min_nfeature = spatial_numeric_or(hit$qc_min_nfeature, cfg$qc_min_nfeature_default),
    qc_max_nfeature = spatial_numeric_or(hit$qc_max_nfeature, cfg$qc_max_nfeature_default),
    qc_min_ncount = spatial_numeric_or(hit$qc_min_ncount, cfg$qc_min_ncount_default),
    qc_max_ncount = spatial_numeric_or(hit$qc_max_ncount, cfg$qc_max_ncount_default),
    qc_max_mito_pct = spatial_numeric_or(hit$qc_max_mito_pct, cfg$qc_max_mito_pct_default),
    mito_set_override = spatial_cell(hit, "mito_set_override", ""),
    spatial_aware_filter = spatial_bool(spatial_cell(hit, "spatial_aware_filter", "false"), default = FALSE),
    excessive_drop_threshold = spatial_numeric_or(hit$excessive_drop_threshold, cfg$excessive_drop_threshold)
  )
}

spatial_object_ids <- function(obj, fallback = "") {
  meta <- obj@meta.data
  section_id <- if ("section_id" %in% colnames(meta)) unique(as.character(meta$section_id))[1] else fallback
  sample_id <- if ("sample_id" %in% colnames(meta)) unique(as.character(meta$sample_id))[1] else if ("orig.ident" %in% colnames(meta)) unique(as.character(meta$orig.ident))[1] else fallback
  list(sample_id = spatial_safe_id(sample_id), section_id = spatial_safe_id(section_id))
}

apply_qc_filter <- function(obj, thr) {
  metrics <- compute_spot_qc_metrics(obj)
  min_feature <- spatial_numeric_or(thr$qc_min_nfeature, 200)
  max_feature <- spatial_numeric_or(thr$qc_max_nfeature, Inf)
  min_count <- spatial_numeric_or(thr$qc_min_ncount, 500)
  max_count <- spatial_numeric_or(thr$qc_max_ncount, Inf)
  max_mito <- spatial_numeric_or(thr$qc_max_mito_pct, 20)

  low_feature <- metrics$nFeature_Spatial < min_feature
  high_feature <- is.finite(max_feature) & metrics$nFeature_Spatial > max_feature
  low_count <- metrics$nCount_Spatial < min_count
  high_count <- is.finite(max_count) & metrics$nCount_Spatial > max_count
  high_mito <- !is.na(metrics$percent.mito) & is.finite(max_mito) & metrics$percent.mito > max_mito
  keep <- !(low_feature | high_feature | low_count | high_count | high_mito)
  names(keep) <- metrics$spot_id

  reason_df <- data.frame(
    spot_id = metrics$spot_id,
    keep = keep,
    low_feature = low_feature,
    high_feature = high_feature,
    low_count = low_count,
    high_count = high_count,
    high_mito = high_mito,
    stringsAsFactors = FALSE
  )
  reason_summary <- data.frame(
    n_drop_low_feature = sum(low_feature, na.rm = TRUE),
    n_drop_high_feature = sum(high_feature, na.rm = TRUE),
    n_drop_low_count = sum(low_count, na.rm = TRUE),
    n_drop_high_count = sum(high_count, na.rm = TRUE),
    n_drop_high_mito = sum(high_mito, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  list(keep = keep, reasons = reason_df, summary = reason_summary, metrics = metrics)
}

spatial_filter_coords <- function(obj) {
  metrics <- compute_spot_qc_metrics(obj)
  coords <- data.frame(row = metrics$row, col = metrics$col, row.names = metrics$spot_id)
  if (all(!is.finite(coords$row)) || all(!is.finite(coords$col))) {
    coords$row <- seq_len(nrow(coords))
    coords$col <- 1
  }
  coords
}

detect_spatial_drop_cluster <- function(obj, keep, eps_factor = 1.5, min_pts = 5, max_cluster_frac = 0.1) {
  coords <- spatial_filter_coords(obj)
  keep <- as.logical(keep[rownames(coords)])
  keep[is.na(keep)] <- TRUE
  drop_coords <- coords[!keep & is.finite(coords$row) & is.finite(coords$col), c("col", "row"), drop = FALSE]
  total_spots <- nrow(coords)
  if (nrow(drop_coords) < min_pts) {
    return(list(n_clusters = 0L, max_cluster_frac = 0, summary = "no_drop_cluster", status = "ok"))
  }
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    return(list(n_clusters = NA_integer_, max_cluster_frac = NA_real_, summary = "dbscan_unavailable", status = "dbscan_unavailable"))
  }

  sample_idx <- seq_len(nrow(coords))
  if (length(sample_idx) > 200L) {
    set.seed(42L)
    sample_idx <- sample(sample_idx, 200L)
  }
  dist_values <- as.numeric(stats::dist(coords[sample_idx, c("col", "row"), drop = FALSE]))
  dist_values <- dist_values[is.finite(dist_values) & dist_values > 0]
  eps <- if (length(dist_values) > 0) stats::median(dist_values) * eps_factor else eps_factor
  if (!is.finite(eps) || eps <= 0) {
    eps <- eps_factor
  }

  db <- dbscan::dbscan(as.matrix(drop_coords), eps = eps, minPts = min_pts)
  cluster_sizes <- table(db$cluster[db$cluster > 0])
  max_frac <- if (length(cluster_sizes) > 0) max(cluster_sizes) / total_spots else 0
  list(
    n_clusters = length(cluster_sizes),
    max_cluster_frac = as.numeric(max_frac),
    summary = if (is.finite(max_frac) && max_frac > max_cluster_frac) "boundary_drop_detected" else "scattered_drop_ok",
    status = "ok"
  )
}

write_qc_filter_mask_plot <- function(obj, keep, out_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required to write spatial QC mask plots", call. = FALSE)
  }
  coords <- spatial_filter_coords(obj)
  keep <- as.logical(keep[rownames(coords)])
  keep[is.na(keep)] <- TRUE
  plot_df <- data.frame(
    spot_id = rownames(coords),
    row = coords$row,
    col = coords$col,
    qc_status = ifelse(keep, "keep", "drop"),
    stringsAsFactors = FALSE
  )
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = col, y = row, color = qc_status)) +
    ggplot2::geom_point(size = 0.95, alpha = 0.88, na.rm = TRUE) +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed() +
    ggplot2::scale_color_manual(values = c(drop = "#C8102E", keep = "#B8B8B8")) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(color = NULL)
  ensure_dir(dirname(out_path))
  ggplot2::ggsave(out_path, p, width = 5.2, height = 5.2, dpi = 180, bg = "white")
  invisible(out_path)
}

spatial_assay_name <- function(obj) {
  if ("Spatial" %in% names(obj@assays)) "Spatial" else Seurat::DefaultAssay(obj)
}

spatial_row_vars <- function(mat) {
  if (requireNamespace("matrixStats", quietly = TRUE) && !inherits(mat, "sparseMatrix")) {
    return(matrixStats::rowVars(as.matrix(mat)))
  }
  if (inherits(mat, "sparseMatrix")) {
    row_mean <- Matrix::rowMeans(mat)
    row_sq_mean <- Matrix::rowMeans(mat ^ 2)
  } else {
    mat <- as.matrix(mat)
    row_mean <- rowMeans(mat)
    row_sq_mean <- rowMeans(mat ^ 2)
  }
  pmax(as.numeric(row_sq_mean - row_mean ^ 2), 0)
}

spatial_set_assay_data <- function(obj, assay, slot, new.data) {
  tryCatch(
    Seurat::SetAssayData(obj, assay = assay, layer = slot, new.data = new.data),
    error = function(e) Seurat::SetAssayData(obj, assay = assay, slot = slot, new.data = new.data)
  )
}

run_normalize_m0 <- function(obj, hvg_n = 2000L) {
  assay <- spatial_assay_name(obj)
  Seurat::DefaultAssay(obj) <- assay
  counts <- spatial_counts_matrix(obj)
  obj <- spatial_set_assay_data(obj, assay = assay, slot = "data", new.data = counts)
  Seurat::VariableFeatures(obj) <- character(0)
  obj
}

run_normalize_m1 <- function(obj, hvg_n = 2000L) {
  assay <- spatial_assay_name(obj)
  Seurat::DefaultAssay(obj) <- assay
  obj <- Seurat::NormalizeData(obj, assay = assay, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, assay = assay, selection.method = "vst", nfeatures = hvg_n, verbose = FALSE)
  features <- Seurat::VariableFeatures(obj)
  if (length(features) > 0) {
    obj <- Seurat::ScaleData(obj, assay = assay, features = features, verbose = FALSE)
  }
  obj
}

run_normalize_m2 <- function(obj, hvg_n = 2000L) {
  assay <- spatial_assay_name(obj)
  Seurat::SCTransform(obj, assay = assay, vst.flavor = "v1", variable.features.n = hvg_n, verbose = FALSE)
}

run_normalize_m3 <- function(obj, hvg_n = 2000L) {
  assay <- spatial_assay_name(obj)
  Seurat::SCTransform(obj, assay = assay, vst.flavor = "v2", variable.features.n = hvg_n, verbose = FALSE)
}

run_normalize_m4_py_bridge <- function(obj, py_bin, script_path, hvg_n = 2000L) {
  if (!file.exists(script_path)) {
    stop(sprintf("Pearson residual bridge script is missing: %s", script_path), call. = FALSE)
  }
  if (!nzchar(py_bin) || Sys.which(py_bin) == "") {
    stop(sprintf("Python executable is not available for py_spatial: %s", py_bin), call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required for Pearson residual bridge", call. = FALSE)
  }
  assay <- spatial_assay_name(obj)
  counts <- spatial_counts_matrix(obj)
  tmp <- tempfile("spatial_pearson_")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  counts_mtx <- file.path(tmp, "counts.mtx")
  features_tsv <- file.path(tmp, "features.tsv")
  barcodes_tsv <- file.path(tmp, "barcodes.tsv")
  output_mtx <- file.path(tmp, "pearson_residuals.mtx")
  metadata_json <- file.path(tmp, "pearson_metadata.json")
  Matrix::writeMM(counts, counts_mtx)
  writeLines(rownames(counts), features_tsv, useBytes = TRUE)
  writeLines(colnames(counts), barcodes_tsv, useBytes = TRUE)
  status <- system2(
    py_bin,
    c(
      script_path,
      "--counts-mtx", counts_mtx,
      "--features", features_tsv,
      "--barcodes", barcodes_tsv,
      "--output-mtx", output_mtx,
      "--metadata-json", metadata_json
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_status <- attr(status, "status") %||% 0L
  if (!identical(as.integer(exit_status), 0L) || !file.exists(output_mtx)) {
    msg <- paste(status, collapse = "\n")
    stop(sprintf("Pearson residual bridge failed: %s", msg), call. = FALSE)
  }
  residuals <- Matrix::readMM(output_mtx)
  rownames(residuals) <- rownames(counts)
  colnames(residuals) <- colnames(counts)
  obj <- spatial_set_assay_data(obj, assay = assay, slot = "scale.data", new.data = as.matrix(residuals))
  vars <- spatial_row_vars(as(residuals, "dgCMatrix"))
  vars[!is.finite(vars)] <- 0
  hvg <- rownames(residuals)[order(vars, decreasing = TRUE)]
  Seurat::VariableFeatures(obj) <- head(hvg, hvg_n)
  obj@misc$pearson_residuals_bridge <- list(metadata_json = if (file.exists(metadata_json)) jsonlite::read_json(metadata_json, simplifyVector = TRUE) else list())
  obj
}

normalization_result_obj <- function(x) {
  if (inherits(x, "Seurat")) {
    return(x)
  }
  if (is.list(x) && !is.null(x$obj) && inherits(x$obj, "Seurat")) {
    return(x$obj)
  }
  NULL
}

normalization_result_status <- function(x) {
  if (is.list(x) && !is.null(x$status)) {
    return(as.character(x$status))
  }
  if (inherits(x, "Seurat")) "ok" else "failed"
}

normalization_result_data_slot <- function(x) {
  if (is.list(x) && !is.null(x$data_slot)) {
    return(as.character(x$data_slot))
  }
  ""
}

compute_hvg_iou_matrix <- function(results_list) {
  ok_methods <- names(results_list)[vapply(results_list, function(x) !is.null(normalization_result_obj(x)), logical(1))]
  if (length(ok_methods) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  mat <- matrix(NA_real_, nrow = length(ok_methods), ncol = length(ok_methods), dimnames = list(ok_methods, ok_methods))
  hvgs <- lapply(results_list[ok_methods], function(x) Seurat::VariableFeatures(normalization_result_obj(x)))
  for (i in ok_methods) {
    for (j in ok_methods) {
      union_n <- length(union(hvgs[[i]], hvgs[[j]]))
      mat[i, j] <- if (i == j) 1 else if (union_n == 0) NA_real_ else length(intersect(hvgs[[i]], hvgs[[j]])) / union_n
    }
  }
  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  out$method <- rownames(out)
  out[, c("method", setdiff(colnames(out), "method")), drop = FALSE]
}

ensure_spatial_pca <- function(obj, npcs = 10L) {
  if ("pca" %in% names(obj@reductions) && ncol(obj@reductions$pca@cell.embeddings) >= 1) {
    return(obj)
  }
  features <- Seurat::VariableFeatures(obj)
  if (length(features) == 0) {
    counts <- spatial_counts_matrix(obj)
    vars <- spatial_row_vars(counts)
    vars[!is.finite(vars)] <- 0
    features <- head(rownames(counts)[order(vars, decreasing = TRUE)], min(2000L, nrow(counts)))
    Seurat::VariableFeatures(obj) <- features
  }
  assay <- Seurat::DefaultAssay(obj)
  obj <- tryCatch(
    Seurat::ScaleData(obj, assay = assay, features = features, verbose = FALSE),
    error = function(e) obj
  )
  Seurat::RunPCA(obj, assay = assay, features = features, npcs = min(npcs, length(features), max(1L, ncol(obj) - 1L)), verbose = FALSE)
}

compute_pca_confounder_correlation <- function(results_list, confounders = c("nCount_Spatial", "percent.mito")) {
  rows <- list()
  for (method in names(results_list)) {
    obj <- normalization_result_obj(results_list[[method]])
    if (is.null(obj)) {
      next
    }
    if (identical(method, "m0_no_normalization")) {
      for (conf in confounders) {
        rows[[length(rows) + 1L]] <- data.frame(method = method, pc = NA_integer_, confounder = conf, spearman_rho = NA_real_, status = "m0_no_hvg", stringsAsFactors = FALSE)
      }
      next
    }
    pca_obj <- tryCatch(ensure_spatial_pca(obj, npcs = 10L), error = function(e) e)
    if (inherits(pca_obj, "error")) {
      rows[[length(rows) + 1L]] <- data.frame(method = method, pc = NA_integer_, confounder = NA_character_, spearman_rho = NA_real_, status = "failed_pca", stringsAsFactors = FALSE)
      next
    }
    emb <- pca_obj@reductions$pca@cell.embeddings
    meta <- pca_obj@meta.data[rownames(emb), , drop = FALSE]
    for (pc_idx in seq_len(min(10L, ncol(emb)))) {
      for (conf in confounders) {
        if (!conf %in% colnames(meta)) {
          rho <- NA_real_
          status <- "missing_confounder"
        } else {
          rho <- suppressWarnings(stats::cor(emb[, pc_idx], as.numeric(meta[[conf]]), method = "spearman", use = "complete.obs"))
          status <- "ok"
        }
        rows[[length(rows) + 1L]] <- data.frame(method = method, pc = pc_idx, confounder = conf, spearman_rho = rho, status = status, stringsAsFactors = FALSE)
      }
    }
  }
  if (length(rows) == 0) {
    return(data.frame(method = character(), pc = integer(), confounder = character(), spearman_rho = numeric(), status = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

extract_spatial_expression <- function(obj, features) {
  assay <- Seurat::DefaultAssay(obj)
  mat <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, slot = "data"),
    error = function(e) spatial_counts_matrix(obj)
  )
  features <- intersect(features, rownames(mat))
  if (length(features) == 0) {
    return(matrix(nrow = 0, ncol = ncol(mat)))
  }
  as.matrix(mat[features, , drop = FALSE])
}

compute_feature_spatial_coherence <- function(obj, features, max_cells = 500L) {
  coords <- spatial_filter_coords(obj)
  keep <- rownames(coords)
  if (length(keep) > max_cells) {
    set.seed(42L)
    keep <- sample(keep, max_cells)
  }
  coords <- coords[keep, , drop = FALSE]
  finite <- is.finite(coords$row) & is.finite(coords$col)
  coords <- coords[finite, , drop = FALSE]
  if (nrow(coords) < 4) {
    return(NA_real_)
  }
  expr <- extract_spatial_expression(obj, features)
  expr <- expr[, rownames(coords), drop = FALSE]
  if (nrow(expr) == 0) {
    return(NA_real_)
  }
  d <- as.matrix(stats::dist(coords[, c("col", "row"), drop = FALSE]))
  diag(d) <- Inf
  k <- min(4L, nrow(d) - 1L)
  nn <- t(apply(d, 1, function(x) order(x)[seq_len(k)]))
  values <- apply(expr, 1, function(v) {
    if (stats::sd(v, na.rm = TRUE) == 0) {
      return(NA_real_)
    }
    neigh_mean <- vapply(seq_len(nrow(nn)), function(i) mean(v[nn[i, ]], na.rm = TRUE), numeric(1))
    suppressWarnings(stats::cor(v, neigh_mean, method = "spearman", use = "complete.obs"))
  })
  mean(values, na.rm = TRUE)
}

compute_marker_spatial_coherence <- function(results_list, marker_panel = NULL) {
  rows <- list()
  marker_features <- character(0)
  if (!is.null(marker_panel) && nrow(marker_panel) > 0) {
    marker_col <- intersect(c("gene", "gene_symbol", "marker", "feature"), colnames(marker_panel))[1]
    if (!is.na(marker_col)) {
      marker_features <- unique(as.character(marker_panel[[marker_col]]))
    }
  }
  for (method in names(results_list)) {
    obj <- normalization_result_obj(results_list[[method]])
    if (is.null(obj)) {
      next
    }
    features <- marker_features
    if (length(features) == 0) {
      features <- head(Seurat::VariableFeatures(obj), 20L)
    }
    if (length(features) == 0) {
      features <- head(rownames(spatial_counts_matrix(obj)), 20L)
    }
    score <- tryCatch(compute_feature_spatial_coherence(obj, features), error = function(e) NA_real_)
    rows[[length(rows) + 1L]] <- data.frame(method = method, feature_set = if (length(marker_features) > 0) "marker_panel" else "top_hvg_fallback", morans_i_approx = score, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) {
    return(data.frame(method = character(), feature_set = character(), morans_i_approx = numeric(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

summarize_timing <- function(results_list) {
  rows <- lapply(names(results_list), function(method) {
    x <- results_list[[method]]
    data.frame(
      method = method,
      status = normalization_result_status(x),
      data_slot = normalization_result_data_slot(x),
      timing_sec = if (is.list(x) && !is.null(x$timing_sec)) as.numeric(x$timing_sec) else NA_real_,
      message = if (is.list(x) && !is.null(x$message)) as.character(x$message) else "",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

read_normalization_override <- function(override_file, section_id = "") {
  if (!nzchar(override_file) || !file.exists(override_file) || file.info(override_file)$size == 0) {
    return("")
  }
  overrides <- spatial_read_tsv(override_file)
  if (nrow(overrides) == 0) {
    return("")
  }
  if ("section_id" %in% colnames(overrides)) {
    overrides <- overrides[!startsWith(trimws(as.character(overrides$section_id)), "#"), , drop = FALSE]
  }
  if (nrow(overrides) == 0) {
    return("")
  }
  method_col <- intersect(c("override_choice", "final_choice", "normalization_method", "method"), colnames(overrides))[1]
  if (is.na(method_col)) {
    return("")
  }
  if ("section_id" %in% colnames(overrides) && nzchar(section_id)) {
    for (candidate in c(section_id, "__DEFAULT__", "all", "ALL")) {
      hit <- overrides[overrides$section_id == candidate, , drop = FALSE]
      if (nrow(hit) > 0) {
        return(spatial_cell(hit[1, , drop = FALSE], method_col, ""))
      }
    }
    return("")
  }
  spatial_cell(overrides[1, , drop = FALSE], method_col, "")
}

select_normalization_method <- function(default, override = "", available_methods = character()) {
  selected <- if (nzchar(override)) override else default
  if (length(available_methods) > 0 && !selected %in% available_methods) {
    if (default %in% available_methods) {
      return(default)
    }
    return(available_methods[[1]])
  }
  selected
}

spatial_tokenize <- function(value) {
  value <- trimws(as.character(value %||% ""))
  if (!nzchar(value)) {
    return(character(0))
  }
  tokens <- trimws(unlist(strsplit(value, "[,;[:space:]]+", perl = TRUE), use.names = FALSE))
  unique(tokens[nzchar(tokens)])
}

spatial_parse_dims <- function(value, fallback = 1:30) {
  value <- trimws(as.character(value %||% ""))
  if (!nzchar(value)) {
    return(fallback)
  }
  if (grepl("^[0-9]+:[0-9]+$", value)) {
    parts <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
    return(seq(parts[[1]], parts[[2]]))
  }
  dims <- suppressWarnings(as.integer(spatial_tokenize(value)))
  dims <- dims[!is.na(dims) & dims > 0]
  if (length(dims) == 0) fallback else dims
}

spatial_parse_numeric_vector <- function(value, fallback) {
  values <- suppressWarnings(as.numeric(spatial_tokenize(value)))
  values <- values[is.finite(values)]
  if (length(values) == 0) fallback else values
}

spatial_layer_settings <- function(layer_file, layer_id = "panorama_st") {
  layers <- spatial_read_tsv(layer_file)
  if (nrow(layers) == 0) {
    return(list())
  }
  if ("layer_id" %in% colnames(layers)) {
    hit <- layers[layers$layer_id == layer_id, , drop = FALSE]
    if (nrow(hit) == 0) {
      hit <- layers[tolower(layers$layer_role %||% "") == "panorama", , drop = FALSE]
    }
    if (nrow(hit) == 0) {
      hit <- layers[1, , drop = FALSE]
    }
  } else {
    hit <- layers[1, , drop = FALSE]
  }
  as.list(hit[1, , drop = FALSE])
}

spatial_region_subset_layers <- function(layer_file) {
  layers <- spatial_read_tsv(layer_file)
  if (nrow(layers) == 0 || !"layer_role" %in% colnames(layers)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  layers[tolower(trimws(as.character(layers$layer_role))) == "region_subset", , drop = FALSE]
}

spatial_subset_panorama_by_region <- function(panorama, selection_column, selection_values) {
  selection_column <- trimws(as.character(selection_column %||% ""))
  if (!nzchar(selection_column) || !selection_column %in% colnames(panorama@meta.data)) {
    stop(sprintf("selection column is missing from panorama metadata: %s", selection_column), call. = FALSE)
  }
  values <- spatial_tokenize(selection_values)
  if (length(values) == 0) {
    stop("selection_values is empty for region subset layer", call. = FALSE)
  }
  meta_values <- as.character(panorama@meta.data[[selection_column]])
  keep <- rownames(panorama@meta.data)[meta_values %in% values]
  if (length(keep) == 0) {
    stop(sprintf("region subset selection returned no spots: %s in %s", paste(values, collapse = ","), selection_column), call. = FALSE)
  }
  Seurat::subset(panorama, cells = keep)
}

spatial_layer_normalization_method <- function(layer_settings) {
  method <- spatial_tokenize(layer_settings$normalization_methods %||% "")
  method <- tolower(if (length(method) == 0) "log" else method[[1]])
  if (method %in% c("sct", "sctransform", "m2_sct_v1", "m3_sct_v2")) {
    "sct"
  } else {
    "log"
  }
}

rebuild_layer_normalization_st <- function(subset_obj, layer_settings) {
  rebuild <- spatial_bool(layer_settings$rebuild_normalization %||% "yes", default = TRUE)
  method <- spatial_layer_normalization_method(layer_settings)
  hvg_n <- as.integer(spatial_numeric_or(layer_settings$hvg_nfeatures, 1500L))
  hvg_n <- max(50L, min(hvg_n, nrow(subset_obj)))
  vars_to_regress <- intersect(spatial_tokenize(layer_settings$vars_to_regress %||% ""), colnames(subset_obj@meta.data))
  if (!rebuild) {
    subset_obj@misc$subcluster_normalization <- list(method = "inherited", hvg_nfeatures = length(Seurat::VariableFeatures(subset_obj)), vars_to_regress = vars_to_regress)
    return(subset_obj)
  }
  if (identical(method, "sct")) {
    vst_flavor <- Sys.getenv("SPATIAL_SCT_VST_FLAVOR", unset = "v2")
    regress_arg <- if (length(vars_to_regress) > 0) vars_to_regress else NULL
    subset_obj <- Seurat::SCTransform(
      subset_obj,
      variable.features.n = hvg_n,
      vars.to.regress = regress_arg,
      vst.flavor = vst_flavor,
      verbose = FALSE
    )
    Seurat::DefaultAssay(subset_obj) <- "SCT"
  } else {
    assay <- spatial_assay_name(subset_obj)
    Seurat::DefaultAssay(subset_obj) <- assay
    subset_obj <- Seurat::NormalizeData(subset_obj, assay = assay, verbose = FALSE)
    subset_obj <- Seurat::FindVariableFeatures(subset_obj, assay = assay, nfeatures = hvg_n, verbose = FALSE)
    regress_arg <- if (length(vars_to_regress) > 0) vars_to_regress else NULL
    subset_obj <- Seurat::ScaleData(subset_obj, assay = assay, vars.to.regress = regress_arg, verbose = FALSE)
  }
  subset_obj@misc$subcluster_normalization <- list(method = method, hvg_nfeatures = hvg_n, vars_to_regress = vars_to_regress)
  subset_obj
}

pick_subcluster_resolution <- function(metrics_df, target_clusters) {
  if (nrow(metrics_df) == 0 || !"cluster_N" %in% colnames(metrics_df)) {
    return("")
  }
  ok <- metrics_df[metrics_df$status == "ok" & is.finite(metrics_df$cluster_N), , drop = FALSE]
  if (nrow(ok) == 0) {
    return("")
  }
  ok$distance <- abs(as.numeric(ok$cluster_N) - as.numeric(target_clusters))
  ok$silhouette_sort <- if ("silhouette_pca" %in% colnames(ok)) ifelse(is.finite(ok$silhouette_pca), ok$silhouette_pca, -Inf) else -Inf
  ok <- ok[order(ok$distance, -ok$silhouette_sort, as.numeric(ok$resolution)), , drop = FALSE]
  as.character(ok$resolution[[1]])
}

cluster_subset_snn_st <- function(subset_obj, layer_settings) {
  dims <- spatial_parse_dims(layer_settings$pca_dims %||% "1:20", fallback = 1:20)
  target_clusters <- as.integer(spatial_numeric_or(layer_settings$target_clusters, 5L))
  res_range <- spatial_parse_numeric_vector(layer_settings$res_range %||% "", fallback = c(0.2, 0.4, 0.6))
  res_fine_step <- spatial_numeric_or(layer_settings$res_fine_step, 0.05)
  if (length(res_range) >= 2 && is.finite(res_fine_step) && res_fine_step > 0) {
    res_values <- sort(unique(c(res_range, seq(min(res_range), max(res_range), by = res_fine_step))))
  } else {
    res_values <- res_range
  }
  res_values <- res_values[is.finite(res_values) & res_values >= 0]
  if (length(res_values) == 0) {
    res_values <- c(0.4, 0.6)
  }

  assay <- spatial_panorama_assay(subset_obj)
  Seurat::DefaultAssay(subset_obj) <- assay
  features <- intersect(Seurat::VariableFeatures(subset_obj), rownames(subset_obj))
  if (length(features) < 10) {
    features <- compute_panorama_hvg(subset_obj, hvg_n = min(1500L, nrow(subset_obj)))
  }
  features <- intersect(features, rownames(subset_obj))
  npcs <- min(max(dims), length(features), max(1L, ncol(subset_obj) - 1L))
  if (npcs < 1L || length(features) == 0) {
    stop("subset PCA cannot run with zero features or fewer than two spots", call. = FALSE)
  }
  subset_obj <- Seurat::RunPCA(subset_obj, assay = assay, features = features, npcs = npcs, reduction.name = "pca_subcluster", verbose = FALSE)
  dims <- dims[dims <= npcs]
  if (length(dims) == 0) {
    dims <- seq_len(npcs)
  }
  subset_obj <- spatial_run_umap_or_fallback(subset_obj, "pca_subcluster", "umap_subcluster", dims)
  k_param <- max(5L, min(20L, floor(ncol(subset_obj) / 3L)))
  if (k_param >= ncol(subset_obj)) {
    k_param <- max(1L, ncol(subset_obj) - 1L)
  }
  subset_obj <- tryCatch(
    Seurat::FindNeighbors(subset_obj, reduction = "pca_subcluster", dims = dims, k.param = k_param, graph.name = "subcluster_snn", verbose = FALSE),
    error = function(e) Seurat::FindNeighbors(subset_obj, reduction = "pca_subcluster", dims = dims, k.param = k_param, verbose = FALSE)
  )
  use_named_graph <- "subcluster_snn" %in% names(subset_obj@graphs)

  search_rows <- list()
  candidates <- list()
  for (res in res_values) {
    candidate <- tryCatch({
      if (use_named_graph) {
        Seurat::FindClusters(subset_obj, graph.name = "subcluster_snn", resolution = res, verbose = FALSE)
      } else {
        Seurat::FindClusters(subset_obj, resolution = res, verbose = FALSE)
      }
    }, error = function(e) e)
    if (inherits(candidate, "error")) {
      search_rows[[length(search_rows) + 1L]] <- data.frame(
        resolution = res,
        cluster_N = NA_integer_,
        silhouette_pca = NA_real_,
        spatial_coherence_score = NA_real_,
        status = "failed",
        message = conditionMessage(candidate),
        stringsAsFactors = FALSE
      )
      next
    }
    candidate$cluster_default <- factor(candidate@meta.data$seurat_clusters)
    cluster_n <- length(unique(as.character(candidate$cluster_default)))
    candidates[[as.character(res)]] <- candidate
    search_rows[[length(search_rows) + 1L]] <- data.frame(
      resolution = res,
      cluster_N = cluster_n,
      silhouette_pca = compute_silhouette_pca(candidate, "cluster_default", reduction = "pca_subcluster"),
      spatial_coherence_score = compute_spatial_coherence_score(candidate, "cluster_default", k = 6L),
      status = "ok",
      message = "",
      stringsAsFactors = FALSE
    )
  }
  metrics_df <- if (length(search_rows) > 0) do.call(rbind, search_rows) else data.frame()
  picked <- pick_subcluster_resolution(metrics_df, target_clusters)
  if (!nzchar(picked) || is.null(candidates[[picked]])) {
    stop("subcluster SNN failed for all requested resolutions", call. = FALSE)
  }
  out <- candidates[[picked]]
  out@misc$subcluster_snn <- list(
    selected_resolution = picked,
    target_clusters = target_clusters,
    k_param = k_param,
    dims = dims,
    search_table = metrics_df
  )
  list(obj = out, metrics_df = metrics_df, picked_resolution = picked)
}

merge_small_subclusters <- function(subset_obj, cluster_col, threshold_frac) {
  if (!cluster_col %in% colnames(subset_obj@meta.data)) {
    return(subset_obj)
  }
  threshold_frac <- suppressWarnings(as.numeric(threshold_frac))
  if (!is.finite(threshold_frac) || threshold_frac <= 0) {
    return(subset_obj)
  }
  clusters <- as.character(subset_obj@meta.data[[cluster_col]])
  counts <- table(clusters)
  cutoff <- length(clusters) * threshold_frac
  small <- names(counts)[as.numeric(counts) < cutoff]
  if (length(small) == 0) {
    return(subset_obj)
  }
  clusters[clusters %in% small] <- paste0(clusters[clusters %in% small], "_merged_small")
  subset_obj@meta.data[[cluster_col]] <- factor(clusters)
  subset_obj@misc$subcluster_small_clusters <- list(threshold_frac = threshold_frac, small_clusters = small)
  subset_obj
}

spatial_resolve_path <- function(path, base_dir) {
  path <- trimws(as.character(path %||% ""))
  if (!nzchar(path)) {
    return("")
  }
  if (grepl("^/", path)) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(base_dir, path), winslash = "/", mustWork = FALSE)
  }
}

spatial_get_assay_data <- function(obj, assay = NULL, slot = "data") {
  assay <- assay %||% spatial_assay_name(obj)
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = slot),
    error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = slot)
  )
}

spatial_panorama_assay <- function(panorama) {
  norm <- tryCatch(panorama@misc$normalization, error = function(e) NULL)
  if (!is.null(norm$assay) && norm$assay %in% names(panorama@assays)) {
    return(norm$assay)
  }
  if ("SCT" %in% names(panorama@assays)) {
    return("SCT")
  }
  spatial_assay_name(panorama)
}

spatial_panorama_slot <- function(panorama) {
  norm <- tryCatch(panorama@misc$normalization, error = function(e) NULL)
  slot <- as.character(norm$data_slot %||% "")
  if (nzchar(slot)) slot else "data"
}

spatial_valid_python <- function(py_bin) {
  py_bin <- trimws(as.character(py_bin %||% ""))
  if (!nzchar(py_bin)) {
    return(FALSE)
  }
  if (grepl("/", py_bin, fixed = TRUE)) {
    return(file.exists(py_bin) && file.access(py_bin, mode = 1) == 0)
  }
  nzchar(Sys.which(py_bin))
}

load_post_norm_objects <- function(cfg, selected_method_tsv) {
  selected <- spatial_read_tsv(selected_method_tsv)
  if (nrow(selected) == 0) {
    stop(sprintf("selected normalization table is missing or empty: %s", selected_method_tsv), call. = FALSE)
  }
  required <- c("section_id", "final_choice", "data_slot", "canonical_rds")
  missing <- setdiff(required, colnames(selected))
  if (length(missing) > 0) {
    stop(sprintf("selected_method.tsv missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  objects <- list()
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, , drop = FALSE]
    section_id <- spatial_safe_id(spatial_cell(row, "section_id", sprintf("section_%d", i)))
    rds_path <- spatial_resolve_path(spatial_cell(row, "canonical_rds"), cfg$project_root)
    if (!file.exists(rds_path)) {
      stop(sprintf("selected normalized object is missing for section_id=%s: %s", section_id, rds_path), call. = FALSE)
    }
    obj <- readRDS(rds_path)
    obj$section_id <- if ("section_id" %in% colnames(obj@meta.data)) obj$section_id else section_id
    obj$spatial_section_id <- section_id
    obj$normalization_method <- spatial_cell(row, "final_choice")
    obj$normalization_data_slot <- spatial_cell(row, "data_slot", "data")
    obj@misc$normalization <- modifyList(
      obj@misc$normalization %||% list(),
      list(
        method = spatial_cell(row, "final_choice"),
        data_slot = spatial_cell(row, "data_slot", "data"),
        selected_method_tsv = selected_method_tsv,
        canonical_rds = rds_path
      )
    )
    objects[[section_id]] <- obj
  }
  objects
}

merge_spatial_panorama <- function(obj_list) {
  if (length(obj_list) == 0) {
    stop("no spatial objects supplied for panorama merge", call. = FALSE)
  }
  obj_list <- obj_list[vapply(obj_list, inherits, logical(1), what = "Seurat")]
  if (length(obj_list) == 0) {
    stop("no Seurat objects supplied for panorama merge", call. = FALSE)
  }
  section_ids <- names(obj_list)
  if (is.null(section_ids) || any(!nzchar(section_ids))) {
    section_ids <- paste0("section_", seq_along(obj_list))
  }
  for (i in seq_along(obj_list)) {
    obj_list[[i]]$spot_barcode_raw <- colnames(obj_list[[i]])
    obj_list[[i]]$spatial_section_id <- section_ids[[i]]
  }
  if (length(obj_list) == 1) {
    panorama <- obj_list[[1]]
  } else {
    panorama <- Seurat::merge(
      x = obj_list[[1]],
      y = obj_list[-1],
      add.cell.ids = section_ids,
      merge.data = TRUE
    )
  }
  panorama$panorama_id <- colnames(panorama)
  panorama@misc$spatial_panorama <- list(section_ids = section_ids, n_sections = length(section_ids))
  panorama
}

compute_panorama_hvg <- function(panorama, method = "", hvg_n = 2000L) {
  hvg <- tryCatch(Seurat::VariableFeatures(panorama), error = function(e) character(0))
  hvg <- intersect(hvg, rownames(panorama))
  if (length(hvg) >= min(100L, hvg_n)) {
    return(head(hvg, hvg_n))
  }
  assay <- spatial_panorama_assay(panorama)
  mat <- tryCatch(spatial_get_assay_data(panorama, assay = assay, slot = spatial_panorama_slot(panorama)), error = function(e) spatial_counts_matrix(panorama))
  vars <- spatial_row_vars(mat)
  vars[!is.finite(vars)] <- 0
  head(rownames(mat)[order(vars, decreasing = TRUE)], min(hvg_n, nrow(mat)))
}

spatial_create_umap_fallback <- function(obj, source_reduction, target_reduction) {
  emb <- Seurat::Embeddings(obj, reduction = source_reduction)
  if (ncol(emb) == 1) {
    coords <- cbind(emb[, 1], rep(0, nrow(emb)))
  } else {
    coords <- emb[, seq_len(2), drop = FALSE]
  }
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  key <- paste0(toupper(gsub("[^A-Za-z0-9]+", "", target_reduction)), "_")
  obj[[target_reduction]] <- Seurat::CreateDimReducObject(embeddings = coords, key = key, assay = spatial_panorama_assay(obj))
  obj
}

spatial_run_umap_or_fallback <- function(obj, source_reduction, target_reduction, dims) {
  out <- tryCatch(
    Seurat::RunUMAP(obj, reduction = source_reduction, dims = dims, reduction.name = target_reduction, verbose = FALSE),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    msg <- conditionMessage(out)
    out <- spatial_create_umap_fallback(obj, source_reduction, target_reduction)
    out@misc$umap_fallback <- modifyList(out@misc$umap_fallback %||% list(), stats::setNames(list(msg), target_reduction))
  }
  out
}

run_integration_none <- function(panorama, hvg, npcs = 30L) {
  assay <- spatial_panorama_assay(panorama)
  Seurat::DefaultAssay(panorama) <- assay
  hvg <- intersect(hvg, rownames(panorama))
  if (length(hvg) == 0) {
    hvg <- compute_panorama_hvg(panorama, hvg_n = min(2000L, nrow(panorama)))
  }
  Seurat::VariableFeatures(panorama) <- hvg
  panorama <- tryCatch(Seurat::ScaleData(panorama, assay = assay, features = hvg, verbose = FALSE), error = function(e) panorama)
  pca_n <- min(as.integer(npcs), length(hvg), max(1L, ncol(panorama) - 1L))
  panorama <- Seurat::RunPCA(panorama, assay = assay, features = hvg, npcs = pca_n, reduction.name = "pca_none", verbose = FALSE)
  dims <- seq_len(max(1L, min(pca_n, 30L)))
  panorama <- spatial_run_umap_or_fallback(panorama, "pca_none", "umap_none", dims)
  panorama@misc$integration <- modifyList(panorama@misc$integration %||% list(), list(default_mode = "none", selected_mode = "none"))
  panorama
}

run_integration_harmony <- function(panorama, hvg, npcs = 30L, group_by = "section_id") {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("harmony package is not available", call. = FALSE)
  }
  if (!group_by %in% colnames(panorama@meta.data)) {
    stop(sprintf("harmony group column is missing: %s", group_by), call. = FALSE)
  }
  if (!"pca_none" %in% names(panorama@reductions)) {
    panorama <- run_integration_none(panorama, hvg, npcs = npcs)
  }
  panorama <- harmony::RunHarmony(
    object = panorama,
    group.by.vars = group_by,
    reduction = "pca_none",
    reduction.save = "harmony",
    verbose = FALSE
  )
  dims <- seq_len(min(ncol(Seurat::Embeddings(panorama, "harmony")), npcs))
  spatial_run_umap_or_fallback(panorama, "harmony", "umap_harmony", dims)
}

run_integration_cca <- function(panorama, hvg, npcs = 30L) {
  if (!"IntegrateLayers" %in% getNamespaceExports("Seurat")) {
    stop("Seurat::IntegrateLayers is not available", call. = FALSE)
  }
  if (!"pca_none" %in% names(panorama@reductions)) {
    panorama <- run_integration_none(panorama, hvg, npcs = npcs)
  }
  panorama <- Seurat::IntegrateLayers(
    object = panorama,
    method = Seurat::CCAIntegration,
    orig.reduction = "pca_none",
    new.reduction = "cca_integrated",
    verbose = FALSE
  )
  dims <- seq_len(min(ncol(Seurat::Embeddings(panorama, "cca_integrated")), npcs))
  spatial_run_umap_or_fallback(panorama, "cca_integrated", "umap_cca", dims)
}

compute_lisi <- function(panorama, reduction, group_col, perplexity = 30) {
  emb <- Seurat::Embeddings(panorama, reduction = reduction)
  meta <- panorama@meta.data[rownames(emb), , drop = FALSE]
  if (!group_col %in% colnames(meta)) {
    return(list(status = "missing_group", mean_lisi = NA_real_, implementation = "none", message = sprintf("missing group column: %s", group_col)))
  }
  groups <- as.character(meta[[group_col]])
  perplexity <- max(2, min(perplexity, floor(nrow(emb) / 3)))
  if (requireNamespace("lisi", quietly = TRUE) && nrow(emb) >= 4) {
    result <- tryCatch(
      lisi::compute_lisi(emb, meta[, group_col, drop = FALSE], label_colnames = group_col, perplexity = perplexity),
      error = function(e) e
    )
    if (!inherits(result, "error")) {
      return(list(status = "ok", mean_lisi = mean(result[[group_col]], na.rm = TRUE), implementation = "lisi", message = ""))
    }
  }
  tab <- table(groups)
  p <- as.numeric(tab) / sum(tab)
  proxy <- if (length(p) == 0) NA_real_ else 1 / sum(p ^ 2)
  list(status = "ok", mean_lisi = proxy, implementation = "inverse_simpson_proxy", message = "lisi package unavailable or failed; used group diversity proxy")
}

spatial_match_features <- function(genes, features) {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes) & genes != "*"]
  if (length(genes) == 0) {
    return(character(0))
  }
  idx <- match(toupper(genes), toupper(features))
  unique(features[idx[!is.na(idx)]])
}

read_spatial_region_panel <- function(marker_panel_dir, layer_id = "panorama_st") {
  if (!dir.exists(marker_panel_dir)) {
    stop(sprintf("marker panel directory is missing: %s", marker_panel_dir), call. = FALSE)
  }
  files <- list.files(marker_panel_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(files) == 0) {
    stop(sprintf("no active marker panel TSV files found in %s", marker_panel_dir), call. = FALSE)
  }
  panels <- lapply(files, spatial_read_tsv)
  panels <- panels[vapply(panels, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
  if (length(panels) == 0) {
    stop("active marker panel TSV files are empty", call. = FALSE)
  }
  panel <- do.call(rbind, panels)
  required <- c("layer_id", "celltype", "gene", "evidence_source")
  missing <- setdiff(required, colnames(panel))
  if (length(missing) > 0) {
    stop(sprintf("marker panel missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  panel <- panel[panel$layer_id %in% c(layer_id, "*", "all", "ALL"), , drop = FALSE]
  panel$celltype <- trimws(as.character(panel$celltype))
  panel$gene <- trimws(as.character(panel$gene))
  panel <- panel[nzchar(panel$celltype) & nzchar(panel$gene), , drop = FALSE]
  if (nrow(panel) == 0) {
    stop(sprintf("marker panel has no rows for layer_id=%s", layer_id), call. = FALSE)
  }
  panel
}

spatial_region_levels <- function(panel) {
  levels <- unique(as.character(panel$celltype))
  levels <- levels[nzchar(levels) & !tolower(levels) %in% c("uncertain", "mixed_or_uncertain")]
  unique(c(levels, "mixed_or_uncertain"))
}

compute_marker_module_score <- function(panorama, panel) {
  assay <- spatial_panorama_assay(panorama)
  Seurat::DefaultAssay(panorama) <- assay
  regions <- setdiff(spatial_region_levels(panel), "mixed_or_uncertain")
  score_cols <- character()
  for (region in regions) {
    genes <- panel$gene[panel$celltype == region]
    features <- spatial_match_features(genes, rownames(panorama))
    score_col <- paste0("ms_", spatial_safe_id(region))
    if (length(features) == 0) {
      panorama[[score_col]] <- rep(0, ncol(panorama))
    } else {
      panorama <- Seurat::AddModuleScore(panorama, features = list(features), name = score_col, assay = assay, search = FALSE)
      generated <- paste0(score_col, "1")
      if (generated %in% colnames(panorama@meta.data)) {
        panorama@meta.data[[score_col]] <- panorama@meta.data[[generated]]
        panorama@meta.data[[generated]] <- NULL
      }
    }
    score_cols <- c(score_cols, score_col)
  }
  panorama@misc$region_module_score_cols <- stats::setNames(score_cols, regions)
  panorama
}

compute_biological_consistency <- function(panorama, panel, reduction, target_k = 4L) {
  if (is.null(panel) || nrow(panel) == 0 || !reduction %in% names(panorama@reductions)) {
    return(data.frame(mode = reduction, region = character(), n_top_spots = integer(), marker_aligned_cluster_frac = numeric(), stringsAsFactors = FALSE))
  }
  obj <- tryCatch(compute_marker_module_score(panorama, panel), error = function(e) panorama)
  dims <- seq_len(min(30L, ncol(Seurat::Embeddings(obj, reduction))))
  obj <- tryCatch(Seurat::FindNeighbors(obj, reduction = reduction, dims = dims, verbose = FALSE), error = function(e) obj)
  obj <- tryCatch(Seurat::FindClusters(obj, resolution = 0.4, verbose = FALSE), error = function(e) obj)
  cluster_col <- if ("seurat_clusters" %in% colnames(obj@meta.data)) "seurat_clusters" else ""
  rows <- list()
  for (region in setdiff(spatial_region_levels(panel), "mixed_or_uncertain")) {
    ms_col <- paste0("ms_", spatial_safe_id(region))
    if (!ms_col %in% colnames(obj@meta.data) || !nzchar(cluster_col)) {
      next
    }
    scores <- obj@meta.data[[ms_col]]
    top_q <- stats::quantile(scores, 0.75, na.rm = TRUE)
    top_clusters <- obj@meta.data[scores >= top_q, cluster_col, drop = TRUE]
    aligned <- if (length(top_clusters) == 0) 0 else max(table(top_clusters)) / length(top_clusters)
    rows[[length(rows) + 1L]] <- data.frame(
      mode = reduction,
      region = region,
      n_top_spots = length(top_clusters),
      marker_aligned_cluster_frac = as.numeric(aligned),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) {
    return(data.frame(mode = character(), region = character(), n_top_spots = integer(), marker_aligned_cluster_frac = numeric(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

summarize_integration <- function(results) {
  rows <- lapply(names(results), function(mode) {
    x <- results[[mode]]
    data.frame(
      integration_mode = mode,
      status = as.character(x$status %||% "failed"),
      reduction = as.character(x$reduction %||% ""),
      umap = as.character(x$umap %||% ""),
      mean_lisi = as.numeric(x$mean_lisi %||% NA_real_),
      lisi_implementation = as.character(x$lisi_implementation %||% ""),
      biological_consistency = as.numeric(x$biological_consistency %||% NA_real_),
      message = as.character(x$message %||% ""),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

select_integration_mode <- function(layers_tsv, override_tsv = "", layer_id = "panorama_st", default = "none") {
  settings <- spatial_layer_settings(layers_tsv, layer_id = layer_id)
  mode <- trimws(as.character(settings$integration_mode %||% ""))
  env_mode <- Sys.getenv("SPATIAL_INTEGRATION_MODE", unset = "")
  # The env var is a temporary override only when it changes the configured default;
  # otherwise a non-empty spatial_object_layers.tsv value remains the durable decision.
  if (nzchar(env_mode) && (!identical(env_mode, default) || !nzchar(mode))) {
    return(env_mode)
  }
  if (nzchar(mode)) mode else default
}

cluster_b1_seurat <- function(panorama, reduction, dims, target_k, res_range, res_fine_step = 0.05) {
  if (!reduction %in% names(panorama@reductions)) {
    stop(sprintf("missing reduction for b1 clustering: %s", reduction), call. = FALSE)
  }
  max_dim <- ncol(Seurat::Embeddings(panorama, reduction))
  dims <- dims[dims <= max_dim]
  if (length(dims) == 0) {
    dims <- seq_len(max_dim)
  }
  panorama <- tryCatch(
    Seurat::FindNeighbors(panorama, reduction = reduction, dims = dims, graph.name = "b1_snn", verbose = FALSE),
    error = function(e) Seurat::FindNeighbors(panorama, reduction = reduction, dims = dims, verbose = FALSE)
  )
  use_named_graph <- "b1_snn" %in% names(panorama@graphs)
  if (length(res_range) >= 2 && is.finite(res_fine_step) && res_fine_step > 0) {
    res_values <- sort(unique(c(res_range, seq(min(res_range), max(res_range), by = res_fine_step))))
  } else {
    res_values <- res_range
  }
  best <- NULL
  search_rows <- list()
  for (res in res_values) {
    candidate <- tryCatch({
      if (use_named_graph) {
        Seurat::FindClusters(panorama, graph.name = "b1_snn", resolution = res, verbose = FALSE)
      } else {
        Seurat::FindClusters(panorama, resolution = res, verbose = FALSE)
      }
    }, error = function(e) e)
    if (inherits(candidate, "error")) {
      search_rows[[length(search_rows) + 1L]] <- data.frame(resolution = res, cluster_N = NA_integer_, status = "failed", message = conditionMessage(candidate), stringsAsFactors = FALSE)
      next
    }
    clusters <- as.character(candidate@meta.data$seurat_clusters)
    cluster_n <- length(unique(clusters))
    score <- abs(cluster_n - target_k)
    search_rows[[length(search_rows) + 1L]] <- data.frame(resolution = res, cluster_N = cluster_n, status = "ok", message = "", stringsAsFactors = FALSE)
    if (is.null(best) || score < best$score) {
      best <- list(obj = candidate, clusters = clusters, resolution = res, cluster_n = cluster_n, score = score)
    }
  }
  if (is.null(best)) {
    stop("Seurat SNN clustering failed for all requested resolutions", call. = FALSE)
  }
  out <- best$obj
  out$cluster_b1_seurat_snn <- factor(best$clusters)
  out@misc$cluster_b1_seurat_snn <- list(
    selected_resolution = best$resolution,
    target_clusters = target_k,
    search_table = if (length(search_rows) > 0) do.call(rbind, search_rows) else data.frame()
  )
  out
}

cluster_b2_bayesspace <- function(panorama, target_k, dims = 1:30) {
  if (!requireNamespace("BayesSpace", quietly = TRUE)) {
    stop("BayesSpace package is not available", call. = FALSE)
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package is required for BayesSpace backend", call. = FALSE)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required for BayesSpace backend", call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required for BayesSpace backend", call. = FALSE)
  }

  counts <- spatial_counts_matrix(panorama)
  coords <- spatial_filter_coords(panorama)
  cells <- intersect(colnames(panorama), rownames(coords))
  if (length(cells) < 4) {
    stop("BayesSpace backend requires at least 4 spatial spots with coordinates", call. = FALSE)
  }
  counts <- counts[, cells, drop = FALSE]
  coords <- coords[cells, , drop = FALSE]
  meta <- panorama@meta.data[cells, , drop = FALSE]
  section_col <- if ("section_id" %in% colnames(meta)) "section_id" else if ("spatial_section_id" %in% colnames(meta)) "spatial_section_id" else ""
  sections <- if (nzchar(section_col)) as.character(meta[[section_col]]) else rep("panorama", length(cells))
  names(sections) <- cells
  reduction_name <- if ("pca_none" %in% names(panorama@reductions)) "pca_none" else if ("pca" %in% names(panorama@reductions)) "pca" else ""
  reduction_embeddings <- NULL
  if (nzchar(reduction_name)) {
    reduction_embeddings <- Seurat::Embeddings(panorama, reduction = reduction_name)
    reduction_embeddings <- reduction_embeddings[intersect(cells, rownames(reduction_embeddings)), , drop = FALSE]
  }
  nrep <- suppressWarnings(as.integer(Sys.getenv("SPATIAL_BAYESSPACE_NREP", unset = "1000")))
  if (is.na(nrep) || nrep <= 0) {
    nrep <- 1000L
  }
  gamma <- suppressWarnings(as.numeric(Sys.getenv("SPATIAL_BAYESSPACE_GAMMA", unset = "2")))
  if (!is.finite(gamma) || gamma <= 0) {
    gamma <- 2
  }

  labels <- rep(NA_character_, length(cells))
  names(labels) <- cells
  section_results <- list()
  for (section_id in unique(sections)) {
    section_cells <- names(sections)[sections == section_id]
    if (length(section_cells) < 4) {
      next
    }
    section_counts <- counts[, section_cells, drop = FALSE]
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = section_counts))
    lib_size <- Matrix::colSums(section_counts)
    lib_size[!is.finite(lib_size) | lib_size <= 0] <- 1
    norm_counts <- Matrix::t(Matrix::t(section_counts) / lib_size * 10000)
    SingleCellExperiment::logcounts(sce) <- log1p(norm_counts)
    section_coords <- coords[section_cells, , drop = FALSE]
    SummarizedExperiment::colData(sce)$array_row <- as.integer(round(section_coords$row))
    SummarizedExperiment::colData(sce)$array_col <- as.integer(round(section_coords$col))
    SummarizedExperiment::colData(sce)$row <- section_coords$row
    SummarizedExperiment::colData(sce)$col <- section_coords$col
    SummarizedExperiment::colData(sce)$section_id <- section_id

    has_pca <- !is.null(reduction_embeddings) && all(section_cells %in% rownames(reduction_embeddings))
    if (has_pca) {
      section_emb <- reduction_embeddings[section_cells, , drop = FALSE]
      keep_dims <- dims[dims <= ncol(section_emb)]
      if (length(keep_dims) == 0) {
        keep_dims <- seq_len(min(10L, ncol(section_emb)))
      }
      SingleCellExperiment::reducedDim(sce, "PCA") <- as.matrix(section_emb[, keep_dims, drop = FALSE])
      sce <- tryCatch(
        BayesSpace::spatialPreprocess(sce, platform = "Visium", skip.PCA = TRUE),
        error = function(e) sce
      )
    } else {
      n_pcs <- min(max(dims), ncol(sce) - 1L, nrow(sce) - 1L)
      n_hvgs <- min(2000L, nrow(sce))
      sce <- BayesSpace::spatialPreprocess(sce, platform = "Visium", skip.PCA = FALSE, n.PCs = n_pcs, n.HVGs = n_hvgs, log.normalize = TRUE)
    }

    cluster_k <- min(as.integer(target_k), max(2L, floor(length(section_cells) / 2L)))
    d_value <- if ("PCA" %in% SingleCellExperiment::reducedDimNames(sce)) min(length(dims), ncol(SingleCellExperiment::reducedDim(sce, "PCA"))) else min(length(dims), ncol(sce) - 1L)
    d_value <- max(1L, d_value)
    clustered <- BayesSpace::spatialCluster(
      sce,
      q = cluster_k,
      platform = "Visium",
      d = d_value,
      init.method = "mclust",
      model = "t",
      gamma = gamma,
      nrep = nrep,
      save.chain = FALSE
    )
    section_labels <- as.character(SummarizedExperiment::colData(clustered)$spatial.cluster)
    if (length(section_labels) != length(section_cells) || any(is.na(section_labels))) {
      stop(sprintf("BayesSpace did not return valid spatial.cluster labels for section_id=%s", section_id), call. = FALSE)
    }
    labels[section_cells] <- paste(spatial_safe_id(section_id), section_labels, sep = "_")
    section_results[[section_id]] <- list(n_spots = length(section_cells), q = cluster_k, d = d_value)
  }
  if (any(is.na(labels))) {
    missing_n <- sum(is.na(labels))
    stop(sprintf("BayesSpace backend skipped %d spots because one or more sections were too small", missing_n), call. = FALSE)
  }
  final_labels <- labels[colnames(panorama)]
  if (any(is.na(final_labels))) {
    stop("BayesSpace backend did not return labels for all panorama spots", call. = FALSE)
  }
  panorama$cluster_b2_bayesspace <- factor(final_labels)
  panorama@misc$cluster_b2_bayesspace <- list(
    nrep = nrep,
    gamma = gamma,
    target_clusters = target_k,
    section_results = section_results
  )
  panorama
}

export_panorama_anndata <- function(panorama, h5ad_path, reduction = "pca_none", assay = NULL, py_bin = Sys.getenv("PY_SPATIAL_BIN", unset = "python")) {
  if (!spatial_valid_python(py_bin)) {
    stop(sprintf("Python executable is not available for h5ad export: %s", py_bin), call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required for h5ad export", call. = FALSE)
  }
  assay <- assay %||% spatial_panorama_assay(panorama)
  counts <- spatial_counts_matrix(panorama)
  coords <- spatial_filter_coords(panorama)
  coords <- coords[colnames(panorama), , drop = FALSE]
  emb <- if (reduction %in% names(panorama@reductions)) Seurat::Embeddings(panorama, reduction = reduction) else matrix(nrow = ncol(panorama), ncol = 0, dimnames = list(colnames(panorama), NULL))
  emb <- emb[colnames(panorama), , drop = FALSE]
  work_dir <- dirname(h5ad_path)
  ensure_dir(work_dir)
  counts_mtx <- file.path(work_dir, "counts_cells_by_genes.mtx")
  features_tsv <- file.path(work_dir, "features.tsv")
  barcodes_tsv <- file.path(work_dir, "barcodes.tsv")
  coords_tsv <- file.path(work_dir, "coords.tsv")
  pca_tsv <- file.path(work_dir, "pca.tsv")
  writer_py <- file.path(work_dir, "write_h5ad.py")
  Matrix::writeMM(Matrix::t(counts), counts_mtx)
  writeLines(rownames(counts), features_tsv, useBytes = TRUE)
  writeLines(colnames(counts), barcodes_tsv, useBytes = TRUE)
  spatial_write_tsv(data.frame(barcode = rownames(coords), x = coords$col, y = coords$row, stringsAsFactors = FALSE), coords_tsv)
  spatial_write_tsv(data.frame(barcode = rownames(emb), emb, check.names = FALSE), pca_tsv)
  writeLines(c(
    "import argparse",
    "import anndata",
    "import pandas as pd",
    "from scipy.io import mmread",
    "ap = argparse.ArgumentParser()",
    "ap.add_argument('--counts', required=True)",
    "ap.add_argument('--features', required=True)",
    "ap.add_argument('--barcodes', required=True)",
    "ap.add_argument('--coords', required=True)",
    "ap.add_argument('--pca', required=True)",
    "ap.add_argument('--output', required=True)",
    "args = ap.parse_args()",
    "X = mmread(args.counts).tocsr()",
    "features = [line.strip() for line in open(args.features, encoding='utf-8') if line.strip()]",
    "barcodes = [line.strip() for line in open(args.barcodes, encoding='utf-8') if line.strip()]",
    "adata = anndata.AnnData(X=X)",
    "adata.obs_names = barcodes",
    "adata.var_names = features",
    "coords = pd.read_csv(args.coords, sep='\\t').set_index('barcode').loc[barcodes]",
    "adata.obsm['spatial'] = coords[['x', 'y']].to_numpy()",
    "pca = pd.read_csv(args.pca, sep='\\t').set_index('barcode').loc[barcodes]",
    "if pca.shape[1] > 0:",
    "    adata.obsm['X_pca'] = pca.to_numpy()",
    "adata.write_h5ad(args.output)"
  ), writer_py, useBytes = TRUE)
  status <- system2(py_bin, c(writer_py, "--counts", counts_mtx, "--features", features_tsv, "--barcodes", barcodes_tsv, "--coords", coords_tsv, "--pca", pca_tsv, "--output", h5ad_path), stdout = TRUE, stderr = TRUE)
  exit_status <- attr(status, "status") %||% 0L
  if (!identical(as.integer(exit_status), 0L) || !file.exists(h5ad_path)) {
    stop(sprintf("h5ad export failed: %s", paste(status, collapse = "\n")), call. = FALSE)
  }
  invisible(h5ad_path)
}

spatial_read_bridge_clusters <- function(cluster_tsv, cells) {
  clusters <- spatial_read_tsv(cluster_tsv)
  if (nrow(clusters) == 0 || !"barcode" %in% colnames(clusters) || !"cluster" %in% colnames(clusters)) {
    stop(sprintf("bridge output cluster TSV has invalid schema: %s", cluster_tsv), call. = FALSE)
  }
  lookup <- stats::setNames(as.character(clusters$cluster), as.character(clusters$barcode))
  values <- lookup[cells]
  if (any(is.na(values))) {
    stop("bridge output is missing clusters for one or more panorama cells", call. = FALSE)
  }
  factor(values)
}

cluster_python_bridge <- function(panorama, backend, py_bin, script, target_k, work_dir, reduction = "pca_none") {
  if (!spatial_valid_python(py_bin)) {
    stop(sprintf("Python executable is not available: %s", py_bin), call. = FALSE)
  }
  if (!file.exists(script)) {
    stop(sprintf("Python bridge script is missing: %s", script), call. = FALSE)
  }
  ensure_dir(work_dir)
  h5ad_path <- file.path(work_dir, sprintf("%s_input.h5ad", backend))
  cluster_tsv <- file.path(work_dir, sprintf("%s_clusters.tsv", backend))
  metadata_json <- file.path(work_dir, sprintf("%s_metadata.json", backend))
  export_panorama_anndata(panorama, h5ad_path, reduction = reduction, py_bin = py_bin)
  status <- system2(
    py_bin,
    c(script, "--input", h5ad_path, "--output-cluster-tsv", cluster_tsv, "--output-metadata-json", metadata_json, "--target-clusters", as.character(target_k)),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_status <- attr(status, "status") %||% 0L
  if (!identical(as.integer(exit_status), 0L) || !file.exists(cluster_tsv)) {
    stop(sprintf("%s bridge failed: %s", backend, paste(status, collapse = "\n")), call. = FALSE)
  }
  clusters <- spatial_read_bridge_clusters(cluster_tsv, colnames(panorama))
  list(clusters = clusters, cluster_tsv = cluster_tsv, metadata_json = metadata_json)
}

cluster_b3_spagcn_bridge <- function(panorama, py_bin, script, target_k, work_dir, reduction = "pca_none") {
  cluster_python_bridge(panorama, "b3_spagcn", py_bin, script, target_k, work_dir, reduction = reduction)
}

cluster_b4_stagate_bridge <- function(panorama, py_bin, script, target_k, work_dir, reduction = "pca_none") {
  cluster_python_bridge(panorama, "b4_stagate", py_bin, script, target_k, work_dir, reduction = reduction)
}

compute_spatial_coherence_score <- function(panorama, cluster_col, k = 6L) {
  if (!cluster_col %in% colnames(panorama@meta.data)) {
    return(NA_real_)
  }
  coords <- spatial_filter_coords(panorama)
  coords <- coords[rownames(panorama@meta.data), , drop = FALSE]
  finite <- is.finite(coords$row) & is.finite(coords$col)
  coords <- coords[finite, , drop = FALSE]
  if (nrow(coords) < 3) {
    return(NA_real_)
  }
  meta <- panorama@meta.data[rownames(coords), , drop = FALSE]
  section_col <- if ("section_id" %in% colnames(meta)) "section_id" else if ("spatial_section_id" %in% colnames(meta)) "spatial_section_id" else ""
  section_groups <- if (nzchar(section_col)) split(rownames(coords), as.character(meta[[section_col]])) else list(panorama = rownames(coords))
  scores <- vapply(section_groups, function(idx) {
    idx <- intersect(idx, rownames(coords))
    if (length(idx) < 3) {
      return(NA_real_)
    }
    sub_coords <- coords[idx, , drop = FALSE]
    clusters <- as.character(meta[idx, cluster_col, drop = TRUE])
    d <- as.matrix(stats::dist(sub_coords[, c("col", "row"), drop = FALSE]))
    diag(d) <- Inf
    k_local <- min(as.integer(k), nrow(d) - 1L)
    nn <- t(apply(d, 1, function(x) order(x)[seq_len(k_local)]))
    same <- vapply(seq_len(nrow(nn)), function(i) mean(clusters[nn[i, ]] == clusters[[i]], na.rm = TRUE), numeric(1))
    mean(same, na.rm = TRUE)
  }, numeric(1))
  if (all(!is.finite(scores))) NA_real_ else mean(scores, na.rm = TRUE)
}

compute_marker_consistency_score <- function(panorama, cluster_col, panel) {
  if (!cluster_col %in% colnames(panorama@meta.data)) {
    return(NA_real_)
  }
  obj <- tryCatch(compute_marker_module_score(panorama, panel), error = function(e) panorama)
  score_cols <- unlist(obj@misc$region_module_score_cols %||% list(), use.names = FALSE)
  score_cols <- intersect(score_cols, colnames(obj@meta.data))
  if (length(score_cols) == 0) {
    return(NA_real_)
  }
  clusters <- unique(as.character(obj@meta.data[[cluster_col]]))
  winners <- vapply(clusters, function(cl) {
    meta <- obj@meta.data[obj@meta.data[[cluster_col]] == cl, score_cols, drop = FALSE]
    means <- colMeans(meta, na.rm = TRUE)
    if (all(!is.finite(means))) "" else names(means)[which.max(means)]
  }, character(1))
  max(table(winners)) / length(winners)
}

compute_silhouette_pca <- function(panorama, cluster_col, reduction = "pca_none") {
  if (!requireNamespace("cluster", quietly = TRUE) || !cluster_col %in% colnames(panorama@meta.data) || !reduction %in% names(panorama@reductions)) {
    return(NA_real_)
  }
  cells <- rownames(panorama@meta.data)
  max_cells <- suppressWarnings(as.integer(Sys.getenv("SPATIAL_SILHOUETTE_MAX_CELLS", unset = "2000")))
  if (is.na(max_cells) || max_cells <= 0) {
    max_cells <- 2000L
  }
  if (length(cells) > max_cells) {
    set.seed(42L)
    cells <- unlist(lapply(split(cells, panorama@meta.data[cells, cluster_col, drop = TRUE]), function(idx) {
      n <- max(1L, ceiling(length(idx) / length(rownames(panorama@meta.data)) * max_cells))
      sample(idx, min(length(idx), n))
    }), use.names = FALSE)
    cells <- head(unique(cells), max_cells)
  }
  clusters <- as.factor(panorama@meta.data[cells, cluster_col, drop = TRUE])
  if (length(levels(clusters)) < 2 || length(levels(clusters)) >= length(clusters)) {
    return(NA_real_)
  }
  emb <- Seurat::Embeddings(panorama, reduction = reduction)
  emb <- emb[cells, , drop = FALSE]
  d <- stats::dist(emb)
  sil <- cluster::silhouette(as.integer(clusters), d)
  mean(sil[, "sil_width"], na.rm = TRUE)
}

summarize_clustering <- function(results, panel) {
  rows <- lapply(names(results), function(backend) {
    x <- results[[backend]]
    obj <- x$obj
    cluster_col <- as.character(x$cluster_col %||% "")
    data.frame(
      backend = backend,
      status = as.character(x$status %||% "failed"),
      cluster_col = cluster_col,
      cluster_N = if (!is.null(obj) && nzchar(cluster_col) && cluster_col %in% colnames(obj@meta.data)) length(unique(obj@meta.data[[cluster_col]])) else NA_integer_,
      spatial_coherence_score = if (!is.null(obj) && nzchar(cluster_col)) compute_spatial_coherence_score(obj, cluster_col) else NA_real_,
      marker_consistency_score = if (!is.null(obj) && nzchar(cluster_col)) compute_marker_consistency_score(obj, cluster_col, panel) else NA_real_,
      silhouette_pca = if (!is.null(obj) && nzchar(cluster_col)) compute_silhouette_pca(obj, cluster_col, x$reduction %||% "pca_none") else NA_real_,
      timing_sec = as.numeric(x$timing_sec %||% NA_real_),
      message = as.character(x$message %||% ""),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

read_clustering_override <- function(path, section_id = "__DEFAULT__") {
  if (!nzchar(path) || !file.exists(path) || file.info(path)$size == 0) {
    return("")
  }
  overrides <- spatial_read_tsv(path)
  if (nrow(overrides) == 0) {
    return("")
  }
  if ("section_id" %in% colnames(overrides)) {
    overrides <- overrides[!startsWith(trimws(as.character(overrides$section_id)), "#"), , drop = FALSE]
  }
  choice_col <- intersect(c("override_choice", "final_choice", "selected_clustering_backend", "backend", "method"), colnames(overrides))[1]
  if (is.na(choice_col) || nrow(overrides) == 0) {
    return("")
  }
  if ("section_id" %in% colnames(overrides)) {
    for (candidate in c(section_id, "__DEFAULT__", "all", "ALL")) {
      hit <- overrides[overrides$section_id == candidate, , drop = FALSE]
      if (nrow(hit) > 0) {
        return(spatial_cell(hit[1, , drop = FALSE], choice_col, ""))
      }
    }
    return("")
  }
  spatial_cell(overrides[1, , drop = FALSE], choice_col, "")
}

select_clustering_backend <- function(override_tsv, results, default = "b1_seurat_snn") {
  override <- read_clustering_override(override_tsv)
  available <- names(results)[vapply(results, function(x) identical(x$status, "ok") && !is.null(x$obj), logical(1))]
  selected <- if (nzchar(override)) override else default
  if (selected %in% available) {
    return(selected)
  }
  if (default %in% available) {
    return(default)
  }
  if (length(available) > 0) {
    return(available[[1]])
  }
  ""
}

read_region_annotation_override <- function(path) {
  if (!nzchar(path) || !file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(cluster_id = character(), override_region = character(), stringsAsFactors = FALSE))
  }
  overrides <- spatial_read_tsv(path)
  if (nrow(overrides) == 0) {
    return(data.frame(cluster_id = character(), override_region = character(), stringsAsFactors = FALSE))
  }
  if ("cluster_id" %in% colnames(overrides)) {
    overrides <- overrides[!startsWith(trimws(as.character(overrides$cluster_id)), "#"), , drop = FALSE]
  }
  region_col <- intersect(c("override_region", "region", "final_region"), colnames(overrides))[1]
  if (is.na(region_col) || !"cluster_id" %in% colnames(overrides)) {
    return(data.frame(cluster_id = character(), override_region = character(), stringsAsFactors = FALSE))
  }
  data.frame(
    cluster_id = trimws(as.character(overrides$cluster_id)),
    override_region = trimws(as.character(overrides[[region_col]])),
    stringsAsFactors = FALSE
  )
}

region_three_evidence_chain <- function(panorama, cluster_col, panel, thresholds = list(), override_file = "") {
  if (!cluster_col %in% colnames(panorama@meta.data)) {
    stop(sprintf("cluster column is missing: %s", cluster_col), call. = FALSE)
  }
  panorama <- compute_marker_module_score(panorama, panel)
  Seurat::Idents(panorama) <- panorama@meta.data[[cluster_col]]
  marker_table <- tryCatch(
    Seurat::FindAllMarkers(panorama, assay = spatial_panorama_assay(panorama), only.pos = TRUE, logfc.threshold = 0, min.pct = 0.05, verbose = FALSE),
    error = function(e) data.frame(cluster = character(), gene = character(), avg_log2FC = numeric(), p_val_adj = numeric(), pct.1 = numeric(), pct.2 = numeric(), stringsAsFactors = FALSE)
  )
  if (nrow(marker_table) > 0 && !"avg_log2FC" %in% colnames(marker_table) && "avg_logFC" %in% colnames(marker_table)) {
    marker_table$avg_log2FC <- marker_table$avg_logFC
  }
  levels <- spatial_region_levels(panel)
  regions <- setdiff(levels, "mixed_or_uncertain")
  clusters <- sort(unique(as.character(panorama@meta.data[[cluster_col]])))
  score_cols <- panorama@misc$region_module_score_cols %||% list()
  module_matrix <- matrix(NA_real_, nrow = length(clusters), ncol = length(regions), dimnames = list(clusters, regions))
  rows <- list()
  overrides <- read_region_annotation_override(override_file)
  override_lookup <- stats::setNames(overrides$override_region, overrides$cluster_id)
  for (cluster_id in clusters) {
    cells <- rownames(panorama@meta.data)[as.character(panorama@meta.data[[cluster_col]]) == cluster_id]
    marker_genes <- if (nrow(marker_table) > 0 && "cluster" %in% colnames(marker_table)) as.character(marker_table$gene[as.character(marker_table$cluster) == cluster_id]) else character(0)
    panel_hits <- vapply(regions, function(region) {
      genes <- spatial_match_features(panel$gene[panel$celltype == region], rownames(panorama))
      length(intersect(toupper(marker_genes), toupper(genes)))
    }, integer(1))
    module_scores <- vapply(regions, function(region) {
      col <- score_cols[[region]] %||% ""
      if (!nzchar(col) || !col %in% colnames(panorama@meta.data)) {
        return(NA_real_)
      }
      mean(panorama@meta.data[cells, col], na.rm = TRUE)
    }, numeric(1))
    module_matrix[cluster_id, regions] <- module_scores
    panel_winner <- if (length(panel_hits) == 0 || max(panel_hits) == 0) "" else names(panel_hits)[which.max(panel_hits)]
    module_winner <- if (length(module_scores) == 0 || all(!is.finite(module_scores))) "" else names(module_scores)[which.max(module_scores)]
    winner <- if (cluster_id %in% names(override_lookup) && nzchar(override_lookup[[cluster_id]])) {
      override_lookup[[cluster_id]]
    } else if (nzchar(panel_winner) && identical(panel_winner, module_winner)) {
      module_winner
    } else if (nzchar(module_winner) && is.finite(max(module_scores, na.rm = TRUE))) {
      module_winner
    } else {
      "mixed_or_uncertain"
    }
    if (tolower(winner) %in% c("uncertain", "mixed")) {
      winner <- "mixed_or_uncertain"
    }
    if (!winner %in% levels) {
      winner <- "mixed_or_uncertain"
    }
    agreement <- nzchar(panel_winner) && nzchar(module_winner) && identical(panel_winner, module_winner)
    rows[[length(rows) + 1L]] <- data.frame(
      cluster = cluster_id,
      region = winner,
      n_spots = length(cells),
      module_score_max = if (all(!is.finite(module_scores))) NA_real_ else max(module_scores, na.rm = TRUE),
      module_score_winner = module_winner,
      panel_hit_count_winner = if (nzchar(panel_winner)) as.integer(panel_hits[[panel_winner]]) else 0L,
      panel_hit_winner = panel_winner,
      evidence_agreement = agreement,
      confidence = if (agreement && nzchar(panel_winner)) "high" else if (winner != "mixed_or_uncertain") "medium" else "low",
      stringsAsFactors = FALSE
    )
  }
  evidence_df <- do.call(rbind, rows)
  attr(evidence_df, "region_levels") <- levels
  list(
    panorama = panorama,
    evidence_df = evidence_df,
    marker_table = marker_table,
    module_score_matrix = as.data.frame(module_matrix, stringsAsFactors = FALSE)
  )
}

assign_region_labels <- function(panorama, evidence_df) {
  levels <- attr(evidence_df, "region_levels") %||% unique(c(as.character(evidence_df$region), "mixed_or_uncertain"))
  levels <- unique(c(setdiff(levels, "mixed_or_uncertain"), "mixed_or_uncertain"))
  lookup <- stats::setNames(as.character(evidence_df$region), as.character(evidence_df$cluster))
  cluster_col <- if ("cluster_default" %in% colnames(panorama@meta.data)) "cluster_default" else "seurat_clusters"
  values <- lookup[as.character(panorama@meta.data[[cluster_col]])]
  values[is.na(values) | !nzchar(values)] <- "mixed_or_uncertain"
  panorama$region <- factor(values, levels = levels)
  panorama@misc$region_annotation <- list(
    cluster_col = cluster_col,
    evidence = evidence_df,
    levels = levels,
    timestamp = as.character(Sys.time())
  )
  panorama
}

read_spatial_subregion_panel <- function(marker_panel_dir, layer_id) {
  read_spatial_region_panel(marker_panel_dir, layer_id = layer_id)
}

spatial_prefix_subregion <- function(parent_region, sub_label) {
  parent_region <- spatial_safe_id(parent_region)
  sub_label <- spatial_safe_id(sub_label)
  prefix <- paste0(parent_region, "_")
  if (startsWith(sub_label, prefix)) {
    sub_label
  } else {
    paste(parent_region, sub_label, sep = "_")
  }
}

assign_subregion_labels <- function(subset_obj, evidence_df, parent_region) {
  if (nrow(evidence_df) == 0) {
    subset_obj$sub_region <- factor(rep(spatial_prefix_subregion(parent_region, "mixed_or_uncertain"), ncol(subset_obj)))
    return(subset_obj)
  }
  cluster_col <- if ("cluster_default" %in% colnames(subset_obj@meta.data)) "cluster_default" else "seurat_clusters"
  if (!cluster_col %in% colnames(subset_obj@meta.data)) {
    stop("subset object does not contain cluster_default or seurat_clusters", call. = FALSE)
  }
  lookup <- stats::setNames(as.character(evidence_df$region), as.character(evidence_df$cluster))
  raw <- lookup[as.character(subset_obj@meta.data[[cluster_col]])]
  raw[is.na(raw) | !nzchar(raw)] <- "mixed_or_uncertain"
  labels <- vapply(raw, function(x) spatial_prefix_subregion(parent_region, x), character(1))
  subset_obj$sub_region <- factor(labels, levels = unique(labels))
  evidence_out <- evidence_df
  evidence_out$sub_region <- vapply(as.character(evidence_out$region), function(x) spatial_prefix_subregion(parent_region, x), character(1))
  subset_obj@misc$subregion_annotation <- list(
    parent_region = parent_region,
    cluster_col = cluster_col,
    evidence = evidence_out,
    timestamp = as.character(Sys.time())
  )
  subset_obj
}

project_subregion_to_panorama <- function(panorama, subset_objs_list) {
  values <- rep(NA_character_, ncol(panorama))
  names(values) <- colnames(panorama)
  projected <- list()
  mismatches <- list()
  for (name in names(subset_objs_list)) {
    obj <- subset_objs_list[[name]]
    if (is.null(obj) || !"sub_region" %in% colnames(obj@meta.data)) {
      next
    }
    common <- intersect(colnames(panorama), colnames(obj))
    missing <- setdiff(colnames(obj), colnames(panorama))
    if (length(missing) > 0) {
      mismatch_rate <- length(missing) / ncol(obj)
      mismatches[[name]] <- list(n_missing = length(missing), rate = mismatch_rate, head = head(missing, 5))
      if (mismatch_rate > 0.01) {
        warning(sprintf(
          "project_subregion_to_panorama: layer %s has %d/%d (%.1f%%) barcodes missing from panorama (first: %s)",
          name,
          length(missing),
          ncol(obj),
          mismatch_rate * 100,
          paste(head(missing, 3), collapse = ",")
        ), call. = FALSE)
      }
    }
    values[common] <- as.character(obj@meta.data[common, "sub_region", drop = TRUE])
    projected[[name]] <- length(common)
  }
  levels <- unique(values[!is.na(values) & nzchar(values)])
  panorama$sub_region <- factor(values, levels = levels)
  panorama@misc$subregion_annotation <- list(
    projected_layers = projected,
    barcode_mismatches = mismatches,
    n_projected_spots = sum(!is.na(values)),
    timestamp = as.character(Sys.time())
  )
  panorama
}

write_region_subset_assignment_tsv <- function(layer_results, out_path) {
  rows <- lapply(layer_results, function(x) x$assignment %||% data.frame())
  rows <- rows[vapply(rows, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
  if (length(rows) == 0) {
    out <- data.frame(
      layer_id = character(),
      parent_region = character(),
      cluster = character(),
      sub_region = character(),
      n_spots = integer(),
      module_score_max = numeric(),
      panel_hit_count_winner = integer(),
      evidence_agreement = logical(),
      confidence = character(),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, rows)
  }
  spatial_write_tsv(out, out_path)
  out
}

summarize_subcluster_triage <- function(layer_results, cfg) {
  rows <- list()
  for (layer_id in names(layer_results)) {
    result <- layer_results[[layer_id]]
    status <- as.character(result$status %||% "")
    obj <- result$obj
    assignment <- result$assignment %||% data.frame()
    if (!identical(status, "ok") || is.null(obj) || !"sub_region" %in% colnames(obj@meta.data)) {
      rows[[length(rows) + 1L]] <- data.frame(
        layer_id = layer_id,
        signal = if (nzchar(status)) status else "missing_subregion",
        status = if (status %in% c("skipped_no_enabled_layers", "no_subsets_built", "skipped")) "skipped" else "warn",
        value = status,
        threshold = "",
        message = as.character(result$message %||% ""),
        stringsAsFactors = FALSE
      )
      next
    }
    labels <- as.character(obj@meta.data$sub_region)
    clusters <- if ("cluster_default" %in% colnames(obj@meta.data)) as.character(obj@meta.data$cluster_default) else labels
    cluster_counts <- table(clusters)
    small_any <- any(grepl("_merged_small$", names(cluster_counts)))
    overlap_fraction <- if (nrow(assignment) == 0 || !"panel_hit_count_winner" %in% colnames(assignment)) {
      NA_real_
    } else {
      mean(as.numeric(assignment$panel_hit_count_winner) > 0, na.rm = TRUE)
    }
    coherence <- compute_spatial_coherence_score(obj, "sub_region", k = 6L)
    all_undetermined <- all(grepl("mixed_or_uncertain$", labels))
    single_subcluster <- length(unique(labels[!is.na(labels)])) <= 1L
    section_col <- if ("section_id" %in% colnames(obj@meta.data)) "section_id" else if ("spatial_section_id" %in% colnames(obj@meta.data)) "spatial_section_id" else ""
    imbalance <- NA_real_
    if (nzchar(section_col) && length(unique(obj@meta.data[[section_col]])) >= 2L) {
      tab <- prop.table(table(obj@meta.data$sub_region, obj@meta.data[[section_col]]), margin = 2)
      imbalance <- max(apply(tab, 1, function(x) max(x) - min(x)), na.rm = TRUE)
    }
    signal_df <- data.frame(
      layer_id = layer_id,
      signal = c("small_subcluster", "weak_marker_overlap", "low_spatial_coherence", "all_undetermined", "single_subcluster", "cross_section_imbalance"),
      status = c(
        if (small_any) "warn" else "ok",
        if (is.finite(overlap_fraction) && overlap_fraction < cfg$subcluster_triage_overlap_threshold) "warn" else "ok",
        if (is.finite(coherence) && coherence < cfg$subcluster_triage_coherence_threshold) "warn" else "ok",
        if (all_undetermined) "warn" else "ok",
        if (single_subcluster) "warn" else "ok",
        if (is.finite(imbalance) && imbalance > cfg$subcluster_triage_imbalance_threshold) "warn" else "ok"
      ),
      value = c(
        paste(names(cluster_counts), as.integer(cluster_counts), sep = ":", collapse = ","),
        sprintf("%.3f", overlap_fraction),
        sprintf("%.3f", coherence),
        as.character(all_undetermined),
        as.character(single_subcluster),
        sprintf("%.3f", imbalance)
      ),
      threshold = c(
        sprintf("< %.3f of layer spots", cfg$subcluster_small_cluster_frac),
        sprintf(">= %.3f", cfg$subcluster_triage_overlap_threshold),
        sprintf(">= %.3f", cfg$subcluster_triage_coherence_threshold),
        "FALSE",
        "FALSE",
        sprintf("> %.3f", cfg$subcluster_triage_imbalance_threshold)
      ),
      message = c(
        "cluster sizes after small-cluster suffixing",
        "fraction of subclusters with at least one winning panel marker hit",
        "mean same-sub_region spatial kNN fraction within section",
        "all labels are mixed_or_uncertain",
        "layer produced only one sub_region label",
        if (is.finite(imbalance)) "maximum section-wise fraction difference for any sub_region" else "skipped because fewer than two sections are available"
      ),
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- signal_df
  }
  if (length(rows) == 0) {
    return(data.frame(layer_id = character(), signal = character(), status = character(), value = character(), threshold = character(), message = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

compute_region_triage <- function(panorama, evidence_df) {
  region_counts <- table(panorama$region)
  mixed_n <- if ("mixed_or_uncertain" %in% names(region_counts)) as.numeric(region_counts[["mixed_or_uncertain"]]) else 0
  mixed_fraction <- if (length(region_counts) == 0) NA_real_ else mixed_n / sum(region_counts)
  section_col <- if ("section_id" %in% colnames(panorama@meta.data)) "section_id" else if ("spatial_section_id" %in% colnames(panorama@meta.data)) "spatial_section_id" else ""
  balance_ratio <- NA_real_
  if (nzchar(section_col)) {
    cross <- table(panorama$region, panorama@meta.data[[section_col]])
    ratios <- apply(cross, 1, function(x) {
      x <- as.numeric(x)
      positive <- x[x > 0]
      if (length(positive) <= 1) 1 else max(positive) / min(positive)
    })
    balance_ratio <- max(ratios, na.rm = TRUE)
  }
  coherence <- compute_spatial_coherence_score(panorama, "region")
  low_region_n <- if (length(region_counts) == 0) 0L else sum(region_counts < 5)
  marker_strength <- if (nrow(evidence_df) == 0) NA_real_ else mean(evidence_df$panel_hit_count_winner + evidence_df$module_score_max, na.rm = TRUE)
  data.frame(
    signal_name = c("per_region_spot_count", "cluster_region_mapping_table", "mixed_uncertain_fraction", "region_marker_evidence_strength", "cross_section_balance", "region_spatial_coherence"),
    value = c(
      paste(names(region_counts), as.integer(region_counts), sep = ":", collapse = ","),
      sprintf("%d clusters", nrow(evidence_df)),
      sprintf("%.3f", mixed_fraction),
      sprintf("%.3f", marker_strength),
      sprintf("%.3f", balance_ratio),
      sprintf("%.3f", coherence)
    ),
    threshold = c("warn if any region < 5 spots", "must map every cluster", "<=0.20", ">0 preferred", "<=5", ">0.50 preferred"),
    status = c(
      if (low_region_n > 0) "warn" else "ok",
      if (nrow(evidence_df) > 0) "ok" else "fail",
      if (is.finite(mixed_fraction) && mixed_fraction > 0.2) "warn" else "ok",
      if (is.finite(marker_strength) && marker_strength > 0) "ok" else "warn",
      if (is.finite(balance_ratio) && balance_ratio > 5) "warn" else "ok",
      if (is.finite(coherence) && coherence < 0.5) "warn" else "ok"
    ),
    message = c(
      sprintf("%d regions have fewer than 5 spots", low_region_n),
      "cluster to region mapping available",
      "mixed_or_uncertain spot fraction",
      "mean panel-hit plus module-score evidence",
      "maximum non-zero section imbalance ratio across regions",
      "mean kNN same-region fraction"
    ),
    stringsAsFactors = FALSE
  )
}

spatial_split_csv <- function(value) {
  value <- trimws(as.character(value %||% ""))
  if (!nzchar(value)) {
    return(character(0))
  }
  out <- trimws(unlist(strsplit(value, "[,;]+", perl = TRUE), use.names = FALSE))
  out[nzchar(out)]
}

filter_comparisons_by_modality <- function(comparisons_df, modality = "spatial") {
  if (nrow(comparisons_df) == 0) {
    return(comparisons_df)
  }
  if (!"analysis_modality" %in% colnames(comparisons_df)) {
    comparisons_df$analysis_modality <- "scrna"
  }
  if (!"enabled" %in% colnames(comparisons_df)) {
    comparisons_df$enabled <- "yes"
  }
  keep <- tolower(trimws(as.character(comparisons_df$analysis_modality))) == tolower(modality) &
    tolower(trimws(as.character(comparisons_df$enabled))) != "no"
  comparisons_df[keep, , drop = FALSE]
}

read_spatial_comparisons_05 <- function(cfg) {
  expected <- c(
    "comparison_id", "source_question_id", "display_question_id", "output_alias",
    "report_title", "layer_scope", "contrast_axis", "analysis_mode",
    "analysis_modality", "analysis_unit", "stat_level", "group_var", "ident_1",
    "ident_2", "subset_column", "subset_value", "aggregation_group_var",
    "composition_group_var", "batch_var", "enabled", "min_biological_replicates",
    "force_exploratory", "min_cells_per_group", "logfc_threshold",
    "produces_gene_program", "gene_program_role", "notes"
  )
  df <- spatial_read_tsv(cfg$comparison_sheet)
  for (col in expected) {
    if (!col %in% colnames(df)) {
      df[[col]] <- character(nrow(df))
    }
  }
  if (nrow(df) == 0) {
    return(df[, expected, drop = FALSE])
  }
  df <- df[, expected, drop = FALSE]
  char_cols <- setdiff(expected, c("min_biological_replicates", "min_cells_per_group", "logfc_threshold"))
  for (col in char_cols) {
    df[[col]] <- trimws(as.character(df[[col]]))
    df[[col]][is.na(df[[col]])] <- ""
  }
  df$enabled[!nzchar(df$enabled)] <- "yes"
  df$analysis_modality[!nzchar(df$analysis_modality)] <- "scrna"
  df$layer_scope[!nzchar(df$layer_scope)] <- "panorama_st"
  df$group_var[!nzchar(df$group_var)] <- "condition"
  df$analysis_mode[!nzchar(df$analysis_mode)] <- "condition_pairwise"
  df$analysis_unit[!nzchar(df$analysis_unit)] <- "region"
  df$stat_level[!nzchar(df$stat_level)] <- "spot_level_exploratory"
  df$force_exploratory[!nzchar(df$force_exploratory)] <- "no"
  df$min_biological_replicates <- suppressWarnings(as.integer(df$min_biological_replicates))
  df$min_biological_replicates[is.na(df$min_biological_replicates)] <- 3L
  df$min_cells_per_group <- suppressWarnings(as.integer(df$min_cells_per_group))
  df$min_cells_per_group[is.na(df$min_cells_per_group)] <- 3L
  df$logfc_threshold <- suppressWarnings(as.numeric(df$logfc_threshold))
  df$logfc_threshold[is.na(df$logfc_threshold)] <- cfg$spatial_marker_logfc_threshold %||% 0.25
  filter_comparisons_by_modality(df, "spatial")
}

spatial_comparison_vars_05 <- function(row, cfg) {
  comparison_id <- spatial_cell(row, "comparison_id", "")
  if (!nzchar(comparison_id)) {
    comparison_id <- sprintf("%s_vs_%s", spatial_cell(row, "ident_1", "ident1"), spatial_cell(row, "ident_2", "ident2"))
  }
  list(
    comparison_id = comparison_id,
    ident_1 = spatial_cell(row, "ident_1", ""),
    ident_2 = spatial_cell(row, "ident_2", ""),
    group_var = spatial_cell(row, "group_var", "condition"),
    batch_var = spatial_cell(row, "batch_var", "batch"),
    layer_scope = spatial_cell(row, "layer_scope", "panorama_st"),
    analysis_mode = spatial_cell(row, "analysis_mode", "condition_pairwise"),
    analysis_unit = spatial_cell(row, "analysis_unit", "region"),
    stat_level = spatial_cell(row, "stat_level", "spot_level_exploratory"),
    gene_program_role = spatial_cell(row, "gene_program_role", ""),
    subset_column = spatial_cell(row, "subset_column", ""),
    subset_value = spatial_cell(row, "subset_value", ""),
    force_exploratory = tolower(spatial_cell(row, "force_exploratory", "no")) %in% c("yes", "true", "1", "on"),
    min_biological_replicates = as.integer(row$min_biological_replicates[[1]] %||% 3L),
    min_cells_per_group = as.integer(row$min_cells_per_group[[1]] %||% 3L),
    logfc_threshold = as.numeric(row$logfc_threshold[[1]] %||% (cfg$spatial_marker_logfc_threshold %||% 0.25)),
    min_pct = cfg$spatial_marker_min_pct %||% 0.10,
    top_n = cfg$spatial_marker_top_n %||% 20L
  )
}

load_panorama_for_05 <- function(cfg) {
  path <- if (file.exists(cfg$spatial_panorama_subannotated_rds)) {
    cfg$spatial_panorama_subannotated_rds
  } else {
    cfg$spatial_panorama_annotated_rds
  }
  if (!file.exists(path)) {
    stop(sprintf("missing spatial panorama for module 05: %s", path), call. = FALSE)
  }
  obj <- readRDS(path)
  list(object = obj, path = path, has_sub_region = "sub_region" %in% colnames(obj@meta.data))
}

spatial_layer_ids_for_comparison <- function(vars) {
  values <- spatial_split_csv(vars$layer_scope)
  if (length(values) == 0 || any(values %in% c("*", "all"))) {
    return("panorama_st")
  }
  values
}

subset_panorama_for_comparison_05 <- function(obj, vars) {
  if (!nzchar(vars$subset_column) || !nzchar(vars$subset_value)) {
    return(list(object = obj, status = "ok", reason = ""))
  }
  if (!vars$subset_column %in% colnames(obj@meta.data)) {
    return(list(object = obj[, FALSE], status = "failed_design_resolution", reason = sprintf("missing subset_column=%s", vars$subset_column)))
  }
  keep_values <- spatial_split_csv(vars$subset_value)
  keep <- as.character(obj@meta.data[[vars$subset_column]]) %in% keep_values
  if (!any(keep, na.rm = TRUE)) {
    return(list(object = obj[, FALSE], status = "skipped_too_few_spots", reason = "comparison subset matched no spots"))
  }
  list(object = subset(obj, cells = colnames(obj)[keep]), status = "ok", reason = "")
}

spatial_group_bys_05 <- function(obj) {
  out <- "region"
  if ("sub_region" %in% colnames(obj@meta.data)) {
    out <- c(out, "sub_region")
  }
  out
}

spatial_valid_group_values_05 <- function(meta, group_by_col, include_undetermined = FALSE) {
  values <- unique(as.character(meta[[group_by_col]]))
  values <- values[!is.na(values) & nzchar(values)]
  if (identical(group_by_col, "sub_region") && !include_undetermined) {
    values <- setdiff(values, c("mixed_or_uncertain", "undetermined", "unknown", "NA"))
  }
  sort(values)
}

spatial_parent_region_for_value <- function(meta, group_by_col, group_value) {
  if (!identical(group_by_col, "sub_region") || !"region" %in% colnames(meta)) {
    return("")
  }
  hits <- meta[as.character(meta[[group_by_col]]) == group_value, "region", drop = TRUE]
  hits <- unique(as.character(hits[!is.na(hits) & nzchar(hits)]))
  if (length(hits) == 0) "" else paste(hits, collapse = ",")
}

replicate_gate_check_st <- function(panorama, group_by_col, group_var, vars, sample_col = "sample_id") {
  meta <- panorama@meta.data
  if (!group_var %in% colnames(meta)) {
    return(list(status = "failed_design_resolution", pass = FALSE, reason = sprintf("missing group_var=%s", group_var), table = data.frame()))
  }
  if (!sample_col %in% colnames(meta)) {
    sample_col <- if ("orig.ident" %in% colnames(meta)) "orig.ident" else ""
  }
  if (!nzchar(sample_col)) {
    return(list(status = "skipped_replicate_gate", pass = FALSE, reason = "missing sample_id/orig.ident metadata", table = data.frame()))
  }
  groups <- c(vars$ident_1, vars$ident_2)
  gate_rows <- lapply(groups[nzchar(groups)], function(group_id) {
    sub <- meta[as.character(meta[[group_var]]) == group_id, , drop = FALSE]
    samples <- unique(as.character(sub[[sample_col]]))
    samples <- samples[!is.na(samples) & nzchar(samples)]
    data.frame(group_id = group_id, sample_n = length(samples), spot_n = nrow(sub), stringsAsFactors = FALSE)
  })
  gate_table <- if (length(gate_rows) > 0) do.call(rbind, gate_rows) else data.frame(group_id = character(), sample_n = integer(), spot_n = integer())
  min_n <- if (nrow(gate_table) > 0) min(gate_table$sample_n, na.rm = TRUE) else 0L
  pass <- is.finite(min_n) && min_n >= vars$min_biological_replicates
  list(
    status = if (pass) "pass" else "skip",
    pass = pass,
    reason = if (pass) "" else sprintf("min sample_n=%s < min_biological_replicates=%s", min_n, vars$min_biological_replicates),
    table = gate_table
  )
}

compute_region_markers_st <- function(panorama, group_by_col, vars, include_undetermined = FALSE) {
  meta <- panorama@meta.data
  if (!group_by_col %in% colnames(meta)) {
    return(list(markers = data.frame(), summary = data.frame(), status = "failed_design_resolution", reason = sprintf("missing group_by=%s", group_by_col)))
  }
  values <- spatial_valid_group_values_05(meta, group_by_col, include_undetermined)
  count_df <- data.frame(
    group_value = values,
    n_spots = vapply(values, function(value) sum(as.character(meta[[group_by_col]]) == value, na.rm = TRUE), integer(1)),
    stringsAsFactors = FALSE
  )
  count_df$status <- ifelse(count_df$n_spots >= vars$min_cells_per_group, "ok", "skipped_too_few_spots")
  count_df$reason <- ifelse(count_df$status == "ok", "", sprintf("n_spots < min_cells_per_group=%s", vars$min_cells_per_group))
  ok_values <- count_df$group_value[count_df$status == "ok"]
  if (length(ok_values) < 2L) {
    return(list(markers = data.frame(), summary = count_df, status = "skipped_too_few_spots", reason = "fewer than two groups have enough spots"))
  }
  keep <- as.character(meta[[group_by_col]]) %in% ok_values
  obj <- subset(panorama, cells = colnames(panorama)[keep])
  assay <- spatial_assay_name(obj)
  Seurat::DefaultAssay(obj) <- assay
  Seurat::Idents(obj) <- group_by_col
  res <- tryCatch(
    Seurat::FindAllMarkers(
      obj,
      only.pos = TRUE,
      group.by = group_by_col,
      logfc.threshold = vars$logfc_threshold,
      min.pct = vars$min_pct,
      test.use = "wilcox",
      verbose = FALSE
    ),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    count_df$status[count_df$status == "ok"] <- "failed_findmarkers"
    count_df$reason[count_df$status == "failed_findmarkers"] <- conditionMessage(res)
    return(list(markers = data.frame(), summary = count_df, status = "failed_findmarkers", reason = conditionMessage(res)))
  }
  markers <- as.data.frame(res, stringsAsFactors = FALSE)
  if (nrow(markers) > 0 && !"gene" %in% colnames(markers)) {
    markers$gene <- rownames(markers)
  }
  if (nrow(markers) > 0 && "cluster" %in% colnames(markers)) {
    markers$group_value <- as.character(markers$cluster)
  }
  list(markers = markers, summary = count_df, status = if (nrow(markers) > 0) "ok" else "failed_findmarkers", reason = if (nrow(markers) > 0) "" else "FindAllMarkers returned no rows")
}

compute_spatial_composition <- function(panorama, group_by_col, vars) {
  meta <- panorama@meta.data
  section_col <- if ("section_id" %in% colnames(meta)) "section_id" else if ("spatial_section_id" %in% colnames(meta)) "spatial_section_id" else if ("sample_id" %in% colnames(meta)) "sample_id" else "orig.ident"
  sample_col <- if ("sample_id" %in% colnames(meta)) "sample_id" else if ("orig.ident" %in% colnames(meta)) "orig.ident" else section_col
  if (!group_by_col %in% colnames(meta) || !vars$group_var %in% colnames(meta)) {
    return(data.frame())
  }
  tab <- as.data.frame(table(
    section_id = as.character(meta[[section_col]]),
    sample_id = as.character(meta[[sample_col]]),
    condition = as.character(meta[[vars$group_var]]),
    group_value = as.character(meta[[group_by_col]])
  ), stringsAsFactors = FALSE)
  colnames(tab)[ncol(tab)] <- "n_spots"
  tab <- tab[tab$n_spots > 0, , drop = FALSE]
  totals <- aggregate(n_spots ~ section_id + sample_id + condition, tab, sum)
  colnames(totals)[ncol(totals)] <- "section_total_spots"
  out <- merge(tab, totals, by = c("section_id", "sample_id", "condition"), all.x = TRUE)
  out$proportion <- ifelse(out$section_total_spots > 0, out$n_spots / out$section_total_spots, NA_real_)
  out
}

run_formal_propeller_st <- function(prop_df, vars) {
  if (nrow(prop_df) == 0) {
    return(list(result = data.frame(), status = "failed_propeller", reason = "empty proportion table"))
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    return(list(result = data.frame(), status = "skipped_no_packages", reason = "limma is required for formal composition testing"))
  }
  res <- tryCatch({
    groups <- sort(unique(as.character(prop_df$group_value)))
    samples <- sort(unique(as.character(prop_df$sample_id)))
    mat <- matrix(NA_real_, nrow = length(groups), ncol = length(samples), dimnames = list(groups, samples))
    for (i in seq_len(nrow(prop_df))) {
      mat[as.character(prop_df$group_value[[i]]), as.character(prop_df$sample_id[[i]])] <- as.numeric(prop_df$proportion[[i]])
    }
    mat[is.na(mat)] <- 0
    sample_info <- unique(prop_df[, c("sample_id", "condition"), drop = FALSE])
    sample_info <- sample_info[match(samples, sample_info$sample_id), , drop = FALSE]
    group <- factor(sample_info$condition, levels = c(vars$ident_2, vars$ident_1))
    if (length(unique(group[!is.na(group)])) < 2L) {
      stop("composition design has fewer than two groups", call. = FALSE)
    }
    y <- qlogis(pmin(pmax(mat, 1e-5), 1 - 1e-5))
    design <- stats::model.matrix(~group)
    fit <- limma::eBayes(limma::lmFit(y, design))
    out <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
    out$group_value <- rownames(out)
    out
  }, error = function(e) e)
  if (inherits(res, "error")) {
    return(list(result = data.frame(), status = "failed_propeller", reason = conditionMessage(res)))
  }
  list(result = res, status = if (nrow(res) > 0) "ok" else "failed_propeller", reason = if (nrow(res) > 0) "" else "limma returned no rows")
}

run_spatial_pseudobulk_de <- function(panorama, vars, group_by_col, region_value) {
  if (!requireNamespace("edgeR", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
    return(list(result = data.frame(), status = "skipped_no_packages", reason = "edgeR and Matrix are required for formal pseudobulk DE"))
  }
  meta <- panorama@meta.data
  if (!all(c(group_by_col, vars$group_var) %in% colnames(meta))) {
    return(list(result = data.frame(), status = "failed_design_resolution", reason = "missing group_by or group_var column"))
  }
  sample_col <- if ("sample_id" %in% colnames(meta)) "sample_id" else if ("orig.ident" %in% colnames(meta)) "orig.ident" else ""
  if (!nzchar(sample_col)) {
    return(list(result = data.frame(), status = "failed_aggregation", reason = "missing sample_id/orig.ident metadata"))
  }
  keep <- as.character(meta[[group_by_col]]) == region_value & as.character(meta[[vars$group_var]]) %in% c(vars$ident_1, vars$ident_2)
  cells <- rownames(meta)[keep]
  if (length(cells) == 0) {
    return(list(result = data.frame(), status = "skipped_too_few_spots", reason = "region has no spots for this comparison"))
  }
  counts <- spatial_counts_matrix(panorama)[, cells, drop = FALSE]
  meta_sub <- meta[cells, , drop = FALSE]
  sample_ids <- as.character(meta_sub[[sample_col]])
  samples <- unique(sample_ids[nzchar(sample_ids)])
  if (length(samples) < 2L) {
    return(list(result = data.frame(), status = "failed_singular_design", reason = "fewer than two pseudobulk samples"))
  }
  pb_counts <- do.call(cbind, lapply(samples, function(sample_id) {
    Matrix::rowSums(counts[, sample_ids == sample_id, drop = FALSE])
  }))
  rownames(pb_counts) <- rownames(counts)
  colnames(pb_counts) <- samples
  sample_meta <- do.call(rbind, lapply(samples, function(sample_id) {
    rows <- meta_sub[sample_ids == sample_id, , drop = FALSE]
    group_values <- unique(as.character(rows[[vars$group_var]]))
    group_values <- group_values[nzchar(group_values)]
    group_value <- if (length(group_values) > 0) group_values[[1]] else ""
    data.frame(sample_id = sample_id, group = group_value, stringsAsFactors = FALSE)
  }))
  sample_meta <- sample_meta[sample_meta$group %in% c(vars$ident_1, vars$ident_2), , drop = FALSE]
  pb_counts <- pb_counts[, sample_meta$sample_id, drop = FALSE]
  if (length(unique(sample_meta$group)) < 2L) {
    return(list(result = data.frame(), status = "failed_singular_design", reason = "pseudobulk design has fewer than two groups"))
  }
  res <- tryCatch({
    group <- factor(sample_meta$group, levels = c(vars$ident_2, vars$ident_1))
    y <- edgeR::DGEList(counts = pb_counts, group = group)
    keep_genes <- edgeR::filterByExpr(y, group = group)
    y <- y[keep_genes, , keep.lib.sizes = FALSE]
    if (nrow(y) == 0) {
      stop("no genes retained after edgeR filterByExpr", call. = FALSE)
    }
    y <- edgeR::calcNormFactors(y)
    design <- stats::model.matrix(~group)
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)
    test <- edgeR::glmQLFTest(fit, coef = 2)
    out <- as.data.frame(edgeR::topTags(test, n = Inf), stringsAsFactors = FALSE)
    out$gene <- rownames(out)
    out
  }, error = function(e) e)
  if (inherits(res, "error")) {
    msg <- conditionMessage(res)
    status <- if (grepl("design|singular|coef|contrast", msg, ignore.case = TRUE)) "failed_singular_design" else "failed_aggregation"
    return(list(result = data.frame(), status = status, reason = msg))
  }
  list(result = res, status = if (nrow(res) > 0) "ok" else "failed_aggregation", reason = if (nrow(res) > 0) "" else "edgeR returned no rows")
}

run_spotlevel_findmarkers_st <- function(panorama, group_by_col, region_value, vars) {
  meta <- panorama@meta.data
  if (!group_by_col %in% colnames(meta) || !vars$group_var %in% colnames(meta)) {
    return(list(result = data.frame(), status = "failed_design_resolution", reason = "missing group_by or group_var column", n1 = NA_integer_, n2 = NA_integer_))
  }
  keep <- as.character(meta[[group_by_col]]) == region_value
  obj <- subset(panorama, cells = colnames(panorama)[keep])
  groups <- as.character(obj@meta.data[[vars$group_var]])
  n1 <- sum(groups == vars$ident_1, na.rm = TRUE)
  n2 <- sum(groups == vars$ident_2, na.rm = TRUE)
  if (n1 < vars$min_cells_per_group || n2 < vars$min_cells_per_group) {
    return(list(result = data.frame(), status = "skipped_too_few_spots", reason = sprintf("%s=%s, %s=%s, min=%s", vars$ident_1, n1, vars$ident_2, n2, vars$min_cells_per_group), n1 = n1, n2 = n2))
  }
  assay <- spatial_assay_name(obj)
  Seurat::DefaultAssay(obj) <- assay
  Seurat::Idents(obj) <- vars$group_var
  res <- tryCatch(
    Seurat::FindMarkers(
      obj,
      ident.1 = vars$ident_1,
      ident.2 = vars$ident_2,
      logfc.threshold = vars$logfc_threshold,
      min.pct = vars$min_pct,
      test.use = "wilcox",
      verbose = FALSE
    ),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(list(result = data.frame(), status = "failed_findmarkers", reason = conditionMessage(res), n1 = n1, n2 = n2))
  }
  out <- as.data.frame(res, stringsAsFactors = FALSE)
  if (nrow(out) > 0 && !"gene" %in% colnames(out)) {
    out$gene <- rownames(out)
  }
  list(result = out, status = if (nrow(out) > 0) "ok" else "failed_findmarkers", reason = if (nrow(out) > 0) "" else "FindMarkers returned no rows", n1 = n1, n2 = n2)
}
