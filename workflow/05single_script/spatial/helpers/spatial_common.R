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
  lines <- lines[grepl("^(chrM|MT|M)\\b", lines, ignore.case = TRUE)]
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
          intersect(toupper(features), toupper(base))
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
