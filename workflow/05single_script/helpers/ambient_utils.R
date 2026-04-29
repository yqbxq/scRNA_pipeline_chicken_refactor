empty_df_local <- function() {
  data.frame(stringsAsFactors = FALSE)
}

get_assay_matrix <- function(seu, assay = "RNA", type = c("counts", "data")) {
  type <- match.arg(type)
  tryCatch(
    Seurat::GetAssayData(seu, assay = assay, layer = type),
    error = function(e) Seurat::GetAssayData(seu, assay = assay, slot = type)
  )
}

read_matrix_from_path <- function(path) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path) || !file.exists(path)) {
    stop(sprintf("矩阵路径不存在: %s", path), call. = FALSE)
  }
  if (dir.exists(path)) {
    return(pick_expression_matrix(Seurat::Read10X(data.dir = path)))
  }
  if (grepl("\\.h5$", path, ignore.case = TRUE)) {
    return(pick_expression_matrix(Seurat::Read10X_h5(path)))
  }
  stop(sprintf("暂不支持的矩阵路径类型: %s", path), call. = FALSE)
}

resolve_sample_split_var <- function(seu) {
  if ("sample_id" %in% colnames(seu@meta.data)) {
    return("sample_id")
  }
  "orig.ident"
}

split_object_by_sample <- function(seu) {
  split_by <- resolve_sample_split_var(seu)
  split_objs <- Seurat::SplitObject(seu, split.by = split_by)
  split_objs[order(names(split_objs))]
}

merge_named_objects <- function(object_list) {
  nonempty <- object_list[vapply(object_list, function(obj) as.integer(ncol(obj)), integer(1)) > 0]
  if (length(nonempty) == 0) {
    stop("对象列表为空，无法合并。", call. = FALSE)
  }
  merged <- nonempty[[1]]
  if (length(nonempty) > 1) {
    merged <- merge(x = nonempty[[1]], y = nonempty[2:length(nonempty)])
  }
  maybe_join_layers(merged)
}

build_sample_branch_model <- function(sample_obj,
                                      dims,
                                      resolution,
                                      hvg_nfeatures,
                                      min_cells = 50L) {
  sample_obj <- maybe_join_layers(sample_obj)
  if (ncol(sample_obj) == 0) {
    return(list(
      status = "skipped_empty",
      object = sample_obj,
      usable_dims = integer(0),
      cluster_col = "provisional_cluster"
    ))
  }

  sample_obj$provisional_cluster <- "cluster_1"
  sample_obj$provisional_umap_1 <- NA_real_
  sample_obj$provisional_umap_2 <- NA_real_

  if (ncol(sample_obj) < min_cells) {
    return(list(
      status = "skipped_low_cells",
      object = sample_obj,
      usable_dims = integer(0),
      cluster_col = "provisional_cluster"
    ))
  }

  model_obj <- sample_obj
  model_obj <- Seurat::NormalizeData(model_obj, verbose = FALSE)
  model_obj <- Seurat::FindVariableFeatures(
    model_obj,
    selection.method = "vst",
    nfeatures = hvg_nfeatures,
    verbose = FALSE
  )
  model_obj <- Seurat::ScaleData(model_obj, verbose = FALSE)
  model_obj <- Seurat::RunPCA(model_obj, verbose = FALSE)

  usable_dims <- dims[dims <= ncol(Seurat::Embeddings(model_obj, "pca"))]
  usable_dims <- usable_dims[is.finite(usable_dims)]
  if (length(usable_dims) == 0) {
    return(list(
      status = "skipped_no_pcs",
      object = sample_obj,
      usable_dims = integer(0),
      cluster_col = "provisional_cluster"
    ))
  }

  model_obj <- Seurat::FindNeighbors(model_obj, dims = usable_dims, verbose = FALSE)
  model_obj <- Seurat::FindClusters(model_obj, resolution = resolution, verbose = FALSE)
  model_obj$provisional_cluster <- as.character(Seurat::Idents(model_obj))

  if (length(usable_dims) >= 2) {
    model_obj <- Seurat::RunUMAP(model_obj, dims = usable_dims, verbose = FALSE)
    umap_mat <- Seurat::Embeddings(model_obj, reduction = "umap")
    if (ncol(umap_mat) >= 2) {
      model_obj$provisional_umap_1 <- umap_mat[, 1]
      model_obj$provisional_umap_2 <- umap_mat[, 2]
    }
  }

  list(
    status = "built",
    object = model_obj,
    usable_dims = usable_dims,
    cluster_col = "provisional_cluster"
  )
}

safe_find_all_markers <- function(seu, cluster_col = "provisional_cluster", top_n = 3L) {
  if (!cluster_col %in% colnames(seu@meta.data) || length(unique(seu[[cluster_col, drop = TRUE]])) < 2) {
    return(empty_df_local())
  }
  marker_df <- tryCatch(
    Seurat::FindAllMarkers(
      seu,
      only.pos = TRUE,
      group.by = cluster_col,
      assay = Seurat::DefaultAssay(seu),
      verbose = FALSE
    ),
    error = function(e) empty_df_local()
  )
  if (nrow(marker_df) == 0 || !"cluster" %in% colnames(marker_df)) {
    return(marker_df)
  }

  score_col <- if ("avg_log2FC" %in% colnames(marker_df)) "avg_log2FC" else if ("avg_logFC" %in% colnames(marker_df)) "avg_logFC" else ""
  out <- lapply(split(marker_df, marker_df$cluster), function(df) {
    if (nzchar(score_col)) {
      df <- df[order(df[[score_col]], decreasing = TRUE), , drop = FALSE]
    }
    utils::head(df, top_n)
  })
  do.call(rbind, out)
}

calculate_marker_leakage <- function(before_counts, after_counts, cluster_labels, marker_df) {
  if (nrow(marker_df) == 0 || length(cluster_labels) == 0) {
    return(empty_df_local())
  }
  cluster_labels <- as.character(cluster_labels)
  out_rows <- list()
  for (i in seq_len(nrow(marker_df))) {
    row <- marker_df[i, , drop = FALSE]
    gene <- as.character(row$gene[1])
    cluster_id <- as.character(row$cluster[1])
    if (!gene %in% rownames(before_counts) || !gene %in% rownames(after_counts)) {
      next
    }
    in_cluster <- cluster_labels == cluster_id
    out_cluster <- !in_cluster
    before_gene <- as.numeric(before_counts[gene, , drop = TRUE])
    after_gene <- as.numeric(after_counts[gene, , drop = TRUE])
    before_in <- if (any(in_cluster)) mean(before_gene[in_cluster]) else NA_real_
    before_out <- if (any(out_cluster)) mean(before_gene[out_cluster]) else NA_real_
    after_in <- if (any(in_cluster)) mean(after_gene[in_cluster]) else NA_real_
    after_out <- if (any(out_cluster)) mean(after_gene[out_cluster]) else NA_real_
    out_rows[[length(out_rows) + 1]] <- data.frame(
      gene = gene,
      source_cluster = cluster_id,
      mean_in_source_before = before_in,
      mean_outside_before = before_out,
      leakage_before = ifelse(is.finite(before_in) && before_in > 0, before_out / before_in, NA_real_),
      mean_in_source_after = after_in,
      mean_outside_after = after_out,
      leakage_after = ifelse(is.finite(after_in) && after_in > 0, after_out / after_in, NA_real_),
      stringsAsFactors = FALSE
    )
  }
  if (length(out_rows) == 0) {
    return(empty_df_local())
  }
  dplyr::bind_rows(out_rows)
}

choose_ambient_method <- function(readiness_row, cfg) {
  if (is.null(readiness_row) || nrow(readiness_row) == 0) {
    return(list(
      preferred_method = "none",
      fallback_method = "none",
      apply_policy = cfg$ambient_apply_policy,
      soupx_ready = FALSE,
      decontx_ready = FALSE,
      cellbender_ready = FALSE,
      raw_matrix_available = FALSE,
      raw_matrix_kind = "",
      ambient_notes = "missing_readiness_row"
    ))
  }

  as_bool <- function(col_name) {
    if (!col_name %in% colnames(readiness_row)) {
      return(FALSE)
    }
    tolower(normalize_scalar_value(readiness_row[[col_name]][1], "false")) %in% c("true", "yes", "1", "on")
  }

  preferred_method <- if ("ambient_preferred_method" %in% colnames(readiness_row)) {
    normalize_scalar_value(readiness_row$ambient_preferred_method[1], cfg$ambient_primary_method)
  } else {
    cfg$ambient_primary_method
  }
  fallback_method <- if ("ambient_fallback_method" %in% colnames(readiness_row)) {
    normalize_scalar_value(readiness_row$ambient_fallback_method[1], cfg$ambient_fallback_method)
  } else {
    cfg$ambient_fallback_method
  }
  apply_policy <- if ("ambient_apply_default" %in% colnames(readiness_row)) {
    normalize_scalar_value(readiness_row$ambient_apply_default[1], cfg$ambient_apply_policy)
  } else {
    cfg$ambient_apply_policy
  }

  list(
    preferred_method = tolower(preferred_method),
    fallback_method = tolower(fallback_method),
    apply_policy = tolower(apply_policy),
    soupx_ready = as_bool("ambient_soupx_ready"),
    decontx_ready = as_bool("ambient_decontx_ready"),
    cellbender_ready = as_bool("ambient_cellbender_ready"),
    raw_matrix_available = as_bool("raw_matrix_available"),
    raw_matrix_kind = if ("raw_matrix_kind" %in% colnames(readiness_row)) normalize_scalar_value(readiness_row$raw_matrix_kind[1]) else "",
    ambient_notes = if ("ambient_notes" %in% colnames(readiness_row)) normalize_scalar_value(readiness_row$ambient_notes[1]) else ""
  )
}

ambient_method_ready <- function(method_name, method_plan) {
  method_name <- tolower(normalize_scalar_value(method_name, "none"))
  switch(
    method_name,
    soupx = isTRUE(method_plan$soupx_ready),
    decontx = isTRUE(method_plan$decontx_ready),
    FALSE
  )
}

ambient_recommend_apply <- function(contamination_fraction, cfg) {
  contamination_fraction <- suppressWarnings(as.numeric(contamination_fraction))
  if (!is.finite(contamination_fraction)) {
    return(FALSE)
  }
  contamination_fraction >= cfg$ambient_recommend_min_contamination
}

should_apply_ambient <- function(policy, recommended_flag) {
  policy <- tolower(normalize_scalar_value(policy, "manual"))
  if (policy %in% c("auto", "always")) {
    return(TRUE)
  }
  if (policy %in% c("apply_recommended", "recommended")) {
    return(isTRUE(recommended_flag))
  }
  FALSE
}

resolve_sample_barcodes <- function(sample_obj) {
  if ("barcode_raw" %in% colnames(sample_obj@meta.data)) {
    barcodes <- as.character(sample_obj$barcode_raw)
  } else {
    barcodes <- colnames(sample_obj)
  }
  barcodes[!nzchar(barcodes)] <- colnames(sample_obj)[!nzchar(barcodes)]
  barcodes
}

align_count_matrices <- function(filtered_counts, raw_counts) {
  common_features <- intersect(rownames(filtered_counts), rownames(raw_counts))
  if (length(common_features) == 0) {
    stop("filtered/raw 矩阵没有共同 feature。", call. = FALSE)
  }
  list(
    filtered = filtered_counts[common_features, , drop = FALSE],
    raw = raw_counts[common_features, , drop = FALSE]
  )
}

extract_soupx_contamination <- function(sc) {
  candidate_values <- c()
  candidate_values <- c(candidate_values, tryCatch(sc$metaData$rho, error = function(e) numeric(0)))
  candidate_values <- c(candidate_values, tryCatch(sc$fit$rhoEst, error = function(e) numeric(0)))
  candidate_values <- c(candidate_values, tryCatch(sc$fit$rho, error = function(e) numeric(0)))
  candidate_values <- suppressWarnings(as.numeric(candidate_values))
  candidate_values <- candidate_values[is.finite(candidate_values)]
  if (length(candidate_values) == 0) {
    return(NA_real_)
  }
  stats::median(candidate_values)
}

run_soupx_for_sample <- function(sample_obj, raw_counts, branch_model) {
  if (is.null(raw_counts)) {
    stop("SoupX 需要 raw droplets 矩阵。", call. = FALSE)
  }
  filtered_counts <- get_assay_matrix(sample_obj, assay = "RNA", type = "counts")
  sample_barcodes <- resolve_sample_barcodes(sample_obj)
  colnames(filtered_counts) <- sample_barcodes
  aligned <- align_count_matrices(filtered_counts, raw_counts)
  filtered_counts <- aligned$filtered
  raw_counts <- aligned$raw

  cluster_labels <- as.character(branch_model$object$provisional_cluster)
  names(cluster_labels) <- sample_barcodes
  cluster_labels <- cluster_labels[colnames(filtered_counts)]

  sc <- SoupX::SoupChannel(tod = raw_counts, toc = filtered_counts)
  sc <- SoupX::setClusters(sc, clusters = cluster_labels)
  sc <- SoupX::autoEstCont(sc, doPlot = FALSE)
  corrected_counts <- SoupX::adjustCounts(sc, roundToInt = TRUE)
  colnames(corrected_counts) <- colnames(sample_obj)

  contamination <- extract_soupx_contamination(sc)
  list(
    corrected_counts = corrected_counts,
    contamination_fraction = contamination,
    per_cell_contamination = rep(contamination, ncol(sample_obj)),
    method = "soupx",
    cluster_labels = as.character(branch_model$object$provisional_cluster)
  )
}

run_decontx_for_sample <- function(sample_obj, raw_counts, branch_model) {
  filtered_counts <- get_assay_matrix(sample_obj, assay = "RNA", type = "counts")
  sample_barcodes <- resolve_sample_barcodes(sample_obj)
  colnames(filtered_counts) <- sample_barcodes
  cluster_labels <- as.character(branch_model$object$provisional_cluster)
  names(cluster_labels) <- sample_barcodes

  if (!is.null(raw_counts)) {
    aligned <- align_count_matrices(filtered_counts, raw_counts)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = aligned$filtered))
    cluster_labels <- cluster_labels[colnames(aligned$filtered)]
    bg_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = aligned$raw))
    sce <- celda::decontX(sce, z = cluster_labels, background = bg_sce)
  } else {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = filtered_counts))
    sce <- celda::decontX(sce, z = cluster_labels)
  }

  corrected_counts <- tryCatch(
    celda::decontXcounts(sce),
    error = function(e) SummarizedExperiment::assay(sce, "decontXcounts")
  )
  corrected_counts <- as.matrix(corrected_counts)
  mapped_names <- colnames(sample_obj)[match(colnames(corrected_counts), sample_barcodes)]
  mapped_names[is.na(mapped_names)] <- colnames(corrected_counts)[is.na(mapped_names)]
  colnames(corrected_counts) <- mapped_names
  corrected_counts <- corrected_counts[, colnames(sample_obj), drop = FALSE]

  per_cell <- tryCatch(
    as.numeric(SummarizedExperiment::colData(sce)$decontX_contamination),
    error = function(e) rep(NA_real_, ncol(sample_obj))
  )
  if (length(per_cell) != ncol(sample_obj)) {
    per_cell <- rep(NA_real_, ncol(sample_obj))
  }

  list(
    corrected_counts = corrected_counts,
    contamination_fraction = if (all(is.na(per_cell))) NA_real_ else stats::median(per_cell, na.rm = TRUE),
    per_cell_contamination = per_cell,
    method = "decontx",
    cluster_labels = as.character(branch_model$object$provisional_cluster)
  )
}

dispatch_ambient_method <- function(method_name, sample_obj, raw_counts, branch_model) {
  method_name <- tolower(normalize_scalar_value(method_name, "none"))
  runners <- list(
    soupx = run_soupx_for_sample,
    decontx = run_decontx_for_sample
  )
  if (!method_name %in% names(runners)) {
    return(list(result = NULL, status = "skipped_unknown_method", note = method_name))
  }
  if (identical(method_name, "soupx") && !requireNamespace("SoupX", quietly = TRUE)) {
    return(list(result = NULL, status = "skipped_soupx_package_missing", note = "SoupX unavailable"))
  }
  if (identical(method_name, "soupx") && is.null(raw_counts)) {
    return(list(result = NULL, status = "skipped_soupx_raw_unavailable", note = "raw droplets unavailable"))
  }
  if (identical(method_name, "decontx")) {
    have_decontx <- requireNamespace("celda", quietly = TRUE) &&
      requireNamespace("SingleCellExperiment", quietly = TRUE) &&
      requireNamespace("SummarizedExperiment", quietly = TRUE)
    if (!have_decontx) {
      return(list(result = NULL, status = "skipped_decontx_package_missing", note = "celda/SingleCellExperiment unavailable"))
    }
  }

  tryCatch(
    list(
      result = runners[[method_name]](sample_obj, raw_counts, branch_model),
      status = paste0("completed_", method_name),
      note = ""
    ),
    error = function(e) {
      list(result = NULL, status = paste0("failed_", method_name), note = conditionMessage(e))
    }
  )
}

build_cellbender_stub <- function(sample_id, sample_obj, inventory_row, cfg, raw_counts = NULL) {
  raw_path <- if (!is.null(inventory_row) && "raw_matrix_dir" %in% colnames(inventory_row)) {
    normalize_scalar_value(inventory_row$raw_matrix_dir[1])
  } else {
    ""
  }
  output_path <- file.path(cfg$checkpoint_dir, "ambient", sample_id, paste0(sample_id, "_cellbender_filtered.h5"))
  expected_cells <- ncol(sample_obj)
  total_droplets <- if (!is.null(raw_counts)) ncol(raw_counts) else NA_integer_
  args <- c(
    "cellbender", "remove-background",
    "--input", raw_path,
    "--output", output_path,
    "--expected-cells", as.character(expected_cells)
  )
  if (is.finite(total_droplets) && !is.na(total_droplets)) {
    args <- c(args, "--total-droplets-included", as.character(total_droplets))
  }
  if (tolower(cfg$cellbender_cuda) %in% c("yes", "true", "1", "on")) {
    args <- c(args, "--cuda")
  }
  if (is.finite(cfg$cellbender_fpr) && !is.na(cfg$cellbender_fpr)) {
    args <- c(args, "--fpr", as.character(cfg$cellbender_fpr))
  }
  extra_args <- split_csv_02(cfg$cellbender_extra_args)
  if (length(extra_args) > 0) {
    args <- c(args, extra_args)
  }
  data.frame(
    sample_id = sample_id,
    cellbender_mode = cfg$cellbender_mode,
    raw_input = raw_path,
    output_path = output_path,
    expected_cells = expected_cells,
    total_droplets_included = total_droplets,
    command = paste(args, collapse = " "),
    stringsAsFactors = FALSE
  )
}
