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
