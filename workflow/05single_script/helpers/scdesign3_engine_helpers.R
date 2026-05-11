scd_target_metrics_cols <- c(
  "target_id", "target_type", "layer_id", "n_simulations_requested",
  "n_simulations_done", "best_resolution", "median_ARI", "iqr_ARI",
  "median_NMI", "iqr_NMI", "median_max_jaccard", "median_min_Jaccard",
  "median_purity", "primary_metric", "primary_metric_value",
  "pass_threshold", "warn_threshold", "fail_threshold", "gate_status",
  "runtime_s", "status", "reason"
)

scd_per_simulation_cols <- c(
  "target_id", "simulation_id", "sim_i", "resolution", "n_real_cells",
  "n_synth_cells", "n_truth_labels", "n_recovered_clusters", "ARI", "NMI",
  "max_jaccard_mean", "min_Jaccard", "purity", "runtime_s", "status",
  "reason"
)

scd_per_label_cols <- c(
  "target_id", "simulation_id", "sim_i", "resolution", "truth_label",
  "best_matching_cluster", "jaccard", "precision", "recall", "f1",
  "n_truth_cells", "n_recovered_cells", "n_overlap", "status", "reason"
)

scd_engine_status_cols <- c(
  "target_id", "target_type", "layer_id", "input_object",
  "resolved_input_path", "n_simulations_requested", "n_simulations_done",
  "fit_ok", "sim_ok", "score_ok", "status", "reason", "runtime_s",
  "checkpoint_path"
)

scd_empty_df <- function(cols) {
  out <- as.data.frame(setNames(rep(list(character(0)), length(cols)), cols), stringsAsFactors = FALSE)
  out[, cols, drop = FALSE]
}

scd_scalar <- function(x, default = "") {
  if (length(x) == 0 || is.null(x) || is.na(x[[1]]) || !nzchar(trimws(as.character(x[[1]])))) {
    return(default)
  }
  trimws(as.character(x[[1]]))
}

scd_as_integer <- function(x, default) {
  value <- suppressWarnings(as.integer(scd_scalar(x, as.character(default))))
  if (is.na(value)) default else value
}

scd_as_numeric <- function(x, default = NA_real_) {
  value <- suppressWarnings(as.numeric(scd_scalar(x, as.character(default))))
  if (is.na(value)) default else value
}

scd_bool <- function(x, default = FALSE) {
  value <- tolower(scd_scalar(x, ifelse(default, "yes", "no")))
  value %in% c("1", "true", "yes", "y", "on")
}

scd_split <- function(x, sep = ",") {
  value <- scd_scalar(x)
  if (!nzchar(value)) {
    return(character(0))
  }
  out <- trimws(unlist(strsplit(value, sep, fixed = TRUE), use.names = FALSE))
  out[nzchar(out)]
}

scd_parse_resolution_grid <- function(x, default = c(0.6)) {
  raw <- scd_split(x, ",")
  if (length(raw) == 0) {
    return(default)
  }
  value <- suppressWarnings(as.numeric(raw))
  value <- value[is.finite(value)]
  if (length(value) == 0) default else value
}

scd_safe_file_exists <- function(path) {
  path <- scd_scalar(path)
  nzchar(path) && file.exists(path) && isTRUE(file.info(path)$size > 0)
}

scd_bind_rows <- function(items, cols) {
  items <- items[!vapply(items, is.null, logical(1))]
  if (length(items) == 0) {
    return(scd_empty_df(cols))
  }
  out <- do.call(rbind, lapply(items, function(x) {
    for (col in cols) {
      if (!col %in% colnames(x)) x[[col]] <- ""
    }
    x[, cols, drop = FALSE]
  }))
  rownames(out) <- NULL
  out[, cols, drop = FALSE]
}

scd_package_available <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

scd_missing_engine_packages <- function() {
  pkgs <- c("scDesign3", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors")
  pkgs[!vapply(pkgs, scd_package_available, logical(1))]
}

scd_metric_value <- function(metrics, name) {
  if (is.null(metrics) || !name %in% names(metrics)) {
    return(NA_real_)
  }
  value <- suppressWarnings(as.numeric(metrics[[name]][[1]]))
  if (is.na(value)) NA_real_ else value
}

scd_eval_threshold_condition <- function(condition, metrics) {
  condition <- trimws(condition)
  if (!nzchar(condition)) {
    return(TRUE)
  }
  match <- regexec("^([A-Za-z0-9_.]+)\\s*(>=|<=|>|<|==|=)\\s*(-?[0-9.]+)\\s*$", condition, perl = TRUE)
  parts <- regmatches(condition, match)[[1]]
  if (length(parts) != 4) {
    return(NA)
  }
  metric_name <- parts[[2]]
  op <- parts[[3]]
  threshold <- suppressWarnings(as.numeric(parts[[4]]))
  value <- scd_metric_value(metrics, metric_name)
  if (!is.finite(value) || !is.finite(threshold)) {
    return(NA)
  }
  switch(
    op,
    ">=" = value >= threshold,
    "<=" = value <= threshold,
    ">" = value > threshold,
    "<" = value < threshold,
    "==" = value == threshold,
    "=" = value == threshold,
    NA
  )
}

scd_eval_threshold_expr <- function(expr, metrics) {
  parts <- scd_split(expr, ";")
  if (length(parts) == 0) {
    return(NA)
  }
  values <- vapply(parts, scd_eval_threshold_condition, logical(1), metrics = metrics, USE.NAMES = FALSE)
  if (any(is.na(values))) {
    return(NA)
  }
  all(values)
}

scd_gate_decision <- function(metrics, pass_threshold, warn_threshold) {
  pass <- scd_eval_threshold_expr(pass_threshold, metrics)
  if (isTRUE(pass)) {
    return("PASS")
  }
  warn <- scd_eval_threshold_expr(warn_threshold, metrics)
  if (isTRUE(warn)) {
    return("WARN")
  }
  "FAIL"
}

scd_interpretation_for_gate <- function(gate_status) {
  switch(
    scd_scalar(gate_status),
    PASS = "core",
    WARN = "exploratory",
    FAIL = "rejected",
    NOT_APPLICABLE = "context",
    DERIVED = "inherit_parent",
    "no"
  )
}

scd_choose2 <- function(x) {
  x <- as.numeric(x)
  x * (x - 1) / 2
}

scd_adjusted_rand_index <- function(truth, cluster) {
  truth <- as.character(truth)
  cluster <- as.character(cluster)
  keep <- nzchar(truth) & nzchar(cluster) & !is.na(truth) & !is.na(cluster)
  truth <- truth[keep]
  cluster <- cluster[keep]
  if (length(truth) < 2 || length(unique(truth)) < 2 || length(unique(cluster)) < 2) {
    return(NA_real_)
  }
  if (scd_package_available("aricode")) {
    value <- tryCatch(aricode::ARI(truth, cluster), error = function(e) NA_real_)
    if (is.finite(value)) {
      return(as.numeric(value))
    }
  }
  tab <- table(truth, cluster)
  sum_comb <- sum(scd_choose2(tab))
  row_comb <- sum(scd_choose2(rowSums(tab)))
  col_comb <- sum(scd_choose2(colSums(tab)))
  n_comb <- scd_choose2(sum(tab))
  if (!is.finite(n_comb) || n_comb == 0) {
    return(NA_real_)
  }
  expected <- row_comb * col_comb / n_comb
  maximum <- (row_comb + col_comb) / 2
  denominator <- maximum - expected
  if (!is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  (sum_comb - expected) / denominator
}

scd_entropy <- function(x) {
  p <- as.numeric(table(x))
  p <- p / sum(p)
  p <- p[p > 0]
  -sum(p * log(p))
}

scd_normalized_mutual_info <- function(truth, cluster) {
  truth <- as.character(truth)
  cluster <- as.character(cluster)
  keep <- nzchar(truth) & nzchar(cluster) & !is.na(truth) & !is.na(cluster)
  truth <- truth[keep]
  cluster <- cluster[keep]
  if (length(truth) < 2 || length(unique(truth)) < 2 || length(unique(cluster)) < 2) {
    return(NA_real_)
  }
  if (scd_package_available("aricode")) {
    value <- tryCatch(aricode::NMI(truth, cluster), error = function(e) NA_real_)
    if (is.finite(value)) {
      return(as.numeric(value))
    }
  }
  tab <- table(truth, cluster)
  n <- sum(tab)
  truth_marg <- rowSums(tab)
  cluster_marg <- colSums(tab)
  mi <- 0
  nz <- which(tab > 0, arr.ind = TRUE)
  for (idx in seq_len(nrow(nz))) {
    i <- nz[idx, 1]
    j <- nz[idx, 2]
    mi <- mi + (tab[i, j] / n) * log((tab[i, j] * n) / (truth_marg[i] * cluster_marg[j]))
  }
  denom <- sqrt(scd_entropy(truth) * scd_entropy(cluster))
  if (!is.finite(denom) || denom == 0) {
    return(NA_real_)
  }
  mi / denom
}

compute_per_label_metrics <- function(cluster, truth, target_id, sim_i, resolution) {
  truth <- as.character(truth)
  cluster <- as.character(cluster)
  labels <- sort(unique(truth[nzchar(truth) & !is.na(truth)]))
  if (length(labels) == 0) {
    return(scd_empty_df(scd_per_label_cols))
  }
  rows <- lapply(labels, function(label) {
    truth_hit <- truth == label
    clusters <- sort(unique(cluster[nzchar(cluster) & !is.na(cluster)]))
    if (length(clusters) == 0) {
      best <- ""
      best_stats <- c(jaccard = NA_real_, precision = NA_real_, recall = NA_real_, f1 = NA_real_, n_recovered = 0, n_overlap = 0)
    } else {
      stats <- lapply(clusters, function(cl) {
        cluster_hit <- cluster == cl
        overlap <- sum(truth_hit & cluster_hit, na.rm = TRUE)
        union_n <- sum(truth_hit | cluster_hit, na.rm = TRUE)
        precision <- ifelse(sum(cluster_hit, na.rm = TRUE) > 0, overlap / sum(cluster_hit, na.rm = TRUE), NA_real_)
        recall <- ifelse(sum(truth_hit, na.rm = TRUE) > 0, overlap / sum(truth_hit, na.rm = TRUE), NA_real_)
        f1 <- ifelse(is.finite(precision + recall) && (precision + recall) > 0, 2 * precision * recall / (precision + recall), NA_real_)
        c(
          jaccard = ifelse(union_n > 0, overlap / union_n, NA_real_),
          precision = precision,
          recall = recall,
          f1 = f1,
          n_recovered = sum(cluster_hit, na.rm = TRUE),
          n_overlap = overlap
        )
      })
      mat <- do.call(rbind, stats)
      best_idx <- which.max(ifelse(is.finite(mat[, "jaccard"]), mat[, "jaccard"], -Inf))
      best <- clusters[[best_idx]]
      best_stats <- mat[best_idx, ]
    }
    data.frame(
      target_id = target_id,
      simulation_id = sprintf("%s_SIM%02d", target_id, sim_i),
      sim_i = as.character(sim_i),
      resolution = as.character(resolution),
      truth_label = label,
      best_matching_cluster = best,
      jaccard = as.character(unname(best_stats[["jaccard"]])),
      precision = as.character(unname(best_stats[["precision"]])),
      recall = as.character(unname(best_stats[["recall"]])),
      f1 = as.character(unname(best_stats[["f1"]])),
      n_truth_cells = as.character(sum(truth_hit, na.rm = TRUE)),
      n_recovered_cells = as.character(unname(best_stats[["n_recovered"]])),
      n_overlap = as.character(unname(best_stats[["n_overlap"]])),
      status = "ok",
      reason = "",
      stringsAsFactors = FALSE
    )
  })
  scd_bind_rows(rows, scd_per_label_cols)
}

compute_clustering_metrics <- function(cluster, truth) {
  truth <- as.character(truth)
  cluster <- as.character(cluster)
  per_label <- compute_per_label_metrics(cluster, truth, "target", 1L, 0)
  jaccard <- suppressWarnings(as.numeric(per_label$jaccard))
  tab <- table(cluster, truth)
  purity <- if (sum(tab) > 0) sum(apply(tab, 1, max)) / sum(tab) else NA_real_
  c(
    ARI = scd_adjusted_rand_index(truth, cluster),
    NMI = scd_normalized_mutual_info(truth, cluster),
    max_jaccard_mean = ifelse(length(jaccard[is.finite(jaccard)]) > 0, mean(jaccard, na.rm = TRUE), NA_real_),
    min_Jaccard = ifelse(length(jaccard[is.finite(jaccard)]) > 0, min(jaccard, na.rm = TRUE), NA_real_),
    purity = purity
  )
}

scd_extract_counts_meta <- function(obj) {
  if (scd_package_available("Seurat") && inherits(obj, "Seurat")) {
    assay <- tryCatch(Seurat::DefaultAssay(obj), error = function(e) NULL)
    counts <- tryCatch(
      SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
      error = function(e) NULL
    )
    if (is.null(counts)) {
      counts <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
    }
    meta <- obj@meta.data
    return(list(counts = counts, meta = as.data.frame(meta, stringsAsFactors = FALSE), object_class = "Seurat"))
  }
  if (scd_package_available("SingleCellExperiment") && inherits(obj, "SingleCellExperiment")) {
    assay_names <- SummarizedExperiment::assayNames(obj)
    assay_use <- if ("counts" %in% assay_names) "counts" else assay_names[[1]]
    counts <- SummarizedExperiment::assay(obj, assay_use)
    meta <- as.data.frame(SummarizedExperiment::colData(obj), stringsAsFactors = FALSE)
    return(list(counts = counts, meta = meta, object_class = "SingleCellExperiment"))
  }
  if (is.list(obj)) {
    count_key <- intersect(c("counts", "count", "count_mat", "matrix", "expr", "expression"), names(obj))
    meta_key <- intersect(c("meta", "metadata", "colData", "cell_metadata"), names(obj))
    if (length(count_key) > 0 && length(meta_key) > 0) {
      return(list(
        counts = obj[[count_key[[1]]]],
        meta = as.data.frame(obj[[meta_key[[1]]]], stringsAsFactors = FALSE),
        object_class = "list"
      ))
    }
  }
  stop("Unsupported object type for scDesign3 engine; expected Seurat, SingleCellExperiment, or list(counts, meta).", call. = FALSE)
}

scd_align_counts_meta <- function(counts, meta) {
  if (is.null(counts) || is.null(dim(counts)) || length(dim(counts)) != 2) {
    stop("Input object does not expose a two-dimensional counts matrix.", call. = FALSE)
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))
  }
  if (is.null(rownames(counts))) {
    rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
  }
  if (is.null(rownames(meta)) || any(!nzchar(rownames(meta)))) {
    if (nrow(meta) == ncol(counts)) {
      rownames(meta) <- colnames(counts)
    }
  }
  common <- intersect(colnames(counts), rownames(meta))
  if (length(common) == 0 && nrow(meta) == ncol(counts)) {
    rownames(meta) <- colnames(counts)
    common <- colnames(counts)
  }
  if (length(common) < 2) {
    stop("Counts columns and metadata rows could not be aligned.", call. = FALSE)
  }
  list(counts = counts[, common, drop = FALSE], meta = meta[common, , drop = FALSE])
}

scd_row_vars <- function(counts) {
  if (scd_package_available("Matrix") && inherits(counts, "sparseMatrix")) {
    means <- Matrix::rowMeans(counts)
    means_sq <- Matrix::rowMeans(counts ^ 2)
  } else {
    dense <- as.matrix(counts)
    means <- rowMeans(dense)
    means_sq <- rowMeans(dense ^ 2)
  }
  vars <- means_sq - means ^ 2
  vars[!is.finite(vars)] <- 0
  vars
}

scd_prepare_target_data <- function(object_path, truth_col, max_cells_per_label, n_hvg, seed) {
  obj <- readRDS(object_path)
  extracted <- scd_extract_counts_meta(obj)
  aligned <- scd_align_counts_meta(extracted$counts, extracted$meta)
  counts <- aligned$counts
  meta <- aligned$meta
  if (!truth_col %in% colnames(meta)) {
    stop(sprintf("Truth column `%s` is missing from input metadata.", truth_col), call. = FALSE)
  }
  truth <- as.character(meta[[truth_col]])
  keep <- !is.na(truth) & nzchar(truth)
  counts <- counts[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  truth <- truth[keep]
  if (length(unique(truth)) < 2) {
    stop("Fewer than two truth labels after filtering; cannot score cluster recovery.", call. = FALSE)
  }
  set.seed(seed)
  by_label <- split(seq_along(truth), truth)
  selected <- unlist(lapply(by_label, function(idx) {
    if (length(idx) <= max_cells_per_label) idx else sample(idx, max_cells_per_label)
  }), use.names = FALSE)
  selected <- sort(selected)
  counts <- counts[, selected, drop = FALSE]
  meta <- meta[selected, , drop = FALSE]
  truth <- truth[selected]
  if (nrow(counts) > n_hvg) {
    vars <- scd_row_vars(counts)
    hvg_idx <- order(vars, decreasing = TRUE)[seq_len(min(n_hvg, length(vars)))]
    counts <- counts[hvg_idx, , drop = FALSE]
  }
  list(
    counts = counts,
    meta = meta,
    truth = truth,
    truth_col = truth_col,
    n_real_cells = ncol(counts),
    n_truth_labels = length(unique(truth)),
    object_class = extracted$object_class
  )
}

scd_to_sce <- function(prepared) {
  if (!all(vapply(c("SingleCellExperiment", "S4Vectors"), scd_package_available, logical(1)))) {
    stop("SingleCellExperiment/S4Vectors is required to build scDesign3 input.", call. = FALSE)
  }
  col_data <- S4Vectors::DataFrame(prepared$meta)
  rownames(col_data) <- colnames(prepared$counts)
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = prepared$counts),
    colData = col_data
  )
}

scd_call_supported_args <- function(fn, args) {
  formal_names <- tryCatch(names(formals(fn)), error = function(e) character(0))
  if (length(formal_names) > 0 && !"..." %in% formal_names) {
    args <- args[names(args) %in% formal_names]
  }
  do.call(fn, args)
}

scd_get_exported <- function(pkg, name) {
  if (!scd_package_available(pkg) || !exists(name, envir = asNamespace(pkg), inherits = FALSE)) {
    return(NULL)
  }
  get(name, envir = asNamespace(pkg), inherits = FALSE)
}

fit_scdesign3 <- function(prepared, truth_col, family_use = "nb", n_cores = 4L) {
  if (!scd_package_available("scDesign3")) {
    stop("R package scDesign3 is not installed.", call. = FALSE)
  }
  sce <- scd_to_sce(prepared)
  construct_data <- scd_get_exported("scDesign3", "construct_data")
  fit_marginal <- scd_get_exported("scDesign3", "fit_marginal")
  fit_copula <- scd_get_exported("scDesign3", "fit_copula")
  simu_new <- scd_get_exported("scDesign3", "simu_new")
  high_level <- scd_get_exported("scDesign3", "scdesign3")

  if (!is.null(construct_data) && !is.null(fit_marginal) && !is.null(fit_copula) && !is.null(simu_new)) {
    low_level <- tryCatch({
      constructed <- scd_call_supported_args(construct_data, list(
        sce = sce,
        assay_use = "counts",
        celltype = truth_col,
        pseudotime = NULL,
        spatial = NULL,
        other_covariates = NULL,
        corr_by = truth_col
      ))
      marginal <- scd_call_supported_args(fit_marginal, list(
        data = constructed,
        sce = sce,
        assay_use = "counts",
        mu_formula = paste0("~", truth_col),
        sigma_formula = "1",
        family_use = family_use,
        n_cores = n_cores
      ))
      copula <- scd_call_supported_args(fit_copula, list(
        data = constructed,
        sce = sce,
        assay_use = "counts",
        marginal_list = marginal,
        family_use = family_use,
        copula = "gaussian",
        n_cores = n_cores
      ))
      list(
        engine = "scDesign3_low_level",
        sce = sce,
        prepared = prepared,
        truth_col = truth_col,
        family_use = family_use,
        n_cores = n_cores,
        constructed = constructed,
        marginal = marginal,
        copula = copula,
        simu_new = simu_new
      )
    }, error = function(e) {
      structure(list(message = conditionMessage(e)), class = "scd_low_level_error")
    })
    if (!inherits(low_level, "scd_low_level_error")) {
      return(low_level)
    }
    low_level_message <- low_level$message
  } else {
    low_level_message <- "low-level scDesign3 functions are unavailable"
  }

  if (!is.null(high_level)) {
    return(list(
      engine = "scDesign3_high_level",
      sce = sce,
      prepared = prepared,
      truth_col = truth_col,
      family_use = family_use,
      n_cores = n_cores,
      scdesign3 = high_level,
      fit_note = sprintf("Using high-level scdesign3 fallback after low-level fit path failed: %s", low_level_message)
    ))
  }
  stop(sprintf("No usable scDesign3 fit path found: %s", low_level_message), call. = FALSE)
}

scd_extract_synthetic_counts <- function(x, target_gene_names = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.null(dim(x)) && length(dim(x)) == 2) {
    counts <- x
    if (!is.null(target_gene_names) && ncol(counts) == length(target_gene_names) && nrow(counts) != length(target_gene_names)) {
      counts <- t(counts)
    }
    if (is.null(rownames(counts)) && !is.null(target_gene_names) && nrow(counts) == length(target_gene_names)) {
      rownames(counts) <- target_gene_names
    }
    if (is.null(colnames(counts))) {
      colnames(counts) <- paste0("synth_cell_", seq_len(ncol(counts)))
    }
    return(counts)
  }
  if (is.list(x)) {
    preferred <- c("new_count", "new_counts", "simu_count", "simu_counts", "count_mat", "counts", "sim_count", "sim_counts", "newCount")
    for (key in intersect(preferred, names(x))) {
      counts <- scd_extract_synthetic_counts(x[[key]], target_gene_names)
      if (!is.null(counts)) {
        return(counts)
      }
    }
    for (key in names(x)) {
      counts <- scd_extract_synthetic_counts(x[[key]], target_gene_names)
      if (!is.null(counts)) {
        return(counts)
      }
    }
  }
  NULL
}

scd_extract_synthetic_truth <- function(x, truth_col, expected_n, fallback_truth) {
  candidates <- list()
  if (is.list(x)) {
    for (key in intersect(c("newCovariate", "new_covariate", "new_covariates", "covariate", "covariates", "simu_covariate"), names(x))) {
      candidates[[length(candidates) + 1L]] <- x[[key]]
    }
  }
  for (candidate in candidates) {
    if (is.data.frame(candidate) && truth_col %in% colnames(candidate) && nrow(candidate) == expected_n) {
      return(as.character(candidate[[truth_col]]))
    }
  }
  rep_len(as.character(fallback_truth), expected_n)
}

simulate_synthetic_counts <- function(fit, sim_i, seed) {
  set.seed(seed)
  target_gene_names <- rownames(fit$prepared$counts)
  if (identical(fit$engine, "scDesign3_low_level")) {
    raw <- scd_call_supported_args(fit$simu_new, list(
      data = fit$constructed,
      sce = fit$sce,
      assay_use = "counts",
      marginal_list = fit$marginal,
      copula_list = fit$copula,
      family_use = fit$family_use,
      n_cores = fit$n_cores
    ))
  } else {
    raw <- scd_call_supported_args(fit$scdesign3, list(
      sce = fit$sce,
      assay_use = "counts",
      celltype = fit$truth_col,
      pseudotime = NULL,
      spatial = NULL,
      other_covariates = NULL,
      mu_formula = paste0("~", fit$truth_col),
      sigma_formula = "1",
      family_use = fit$family_use,
      n_cores = fit$n_cores,
      usebam = FALSE,
      corr_formula = "1",
      copula = "gaussian",
      DT = TRUE,
      pseudo_obs = FALSE,
      return_model = TRUE,
      nonzerovar = FALSE
    ))
  }
  counts <- scd_extract_synthetic_counts(raw, target_gene_names)
  if (is.null(counts)) {
    stop("scDesign3 simulation did not return a recognizable count matrix.", call. = FALSE)
  }
  counts <- counts[rownames(counts) %in% target_gene_names, , drop = FALSE]
  if (nrow(counts) == 0) {
    stop("Synthetic counts do not share genes with the fitted matrix.", call. = FALSE)
  }
  truth <- scd_extract_synthetic_truth(raw, fit$truth_col, ncol(counts), fit$prepared$truth)
  list(counts = counts, truth = truth, raw = raw)
}

scd_recluster_seurat <- function(counts, resolution, n_pcs = 30L) {
  obj <- Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = min(2000L, nrow(counts)), verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  pcs <- max(2L, min(n_pcs, ncol(counts) - 1L, nrow(counts) - 1L))
  obj <- Seurat::RunPCA(obj, npcs = pcs, verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, dims = seq_len(pcs), verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = resolution, verbose = FALSE)
  as.character(obj$seurat_clusters)
}

scd_recluster_pca_kmeans <- function(counts, expected_clusters, seed) {
  set.seed(seed)
  dense <- as.matrix(counts)
  mat <- t(log1p(dense))
  pcs <- tryCatch(stats::prcomp(mat, center = TRUE, scale. = FALSE)$x, error = function(e) NULL)
  if (is.null(pcs) || ncol(pcs) == 0) {
    pcs <- mat
  }
  k <- max(2L, min(as.integer(expected_clusters), nrow(pcs) - 1L))
  if (k < 2L) {
    return(rep("0", nrow(pcs)))
  }
  stats::kmeans(pcs[, seq_len(min(10L, ncol(pcs))), drop = FALSE], centers = k, nstart = 10)$cluster
}

recluster_synthetic <- function(counts, resolution, expected_clusters, seed, n_pcs = 30L) {
  if (ncol(counts) < 3 || nrow(counts) < 2) {
    stop("Synthetic matrix is too small for reclustering.", call. = FALSE)
  }
  if (scd_package_available("Seurat")) {
    clusters <- tryCatch(
      scd_recluster_seurat(counts, resolution, n_pcs),
      error = function(e) NULL
    )
    if (!is.null(clusters) && length(clusters) == ncol(counts)) {
      return(as.character(clusters))
    }
  }
  as.character(scd_recluster_pca_kmeans(counts, expected_clusters, seed))
}

scd_target_metrics_row <- function(target, n_sim_requested, n_sim_done, status, reason, runtime_s, gate_status = "PLANNED", metrics = NULL, best_resolution = "") {
  metric_value <- if (!is.null(metrics)) scd_metric_value(metrics, scd_scalar(target$primary_metric, "ARI")) else NA_real_
  data.frame(
    target_id = scd_scalar(target$target_id),
    target_type = scd_scalar(target$target_type),
    layer_id = scd_scalar(target$layer_id),
    n_simulations_requested = as.character(n_sim_requested),
    n_simulations_done = as.character(n_sim_done),
    best_resolution = as.character(best_resolution),
    median_ARI = as.character(if (!is.null(metrics)) scd_metric_value(metrics, "ARI") else NA_real_),
    iqr_ARI = as.character(if (!is.null(metrics)) scd_metric_value(metrics, "iqr_ARI") else NA_real_),
    median_NMI = as.character(if (!is.null(metrics)) scd_metric_value(metrics, "NMI") else NA_real_),
    iqr_NMI = as.character(if (!is.null(metrics)) scd_metric_value(metrics, "iqr_NMI") else NA_real_),
    median_max_jaccard = as.character(if (!is.null(metrics)) scd_metric_value(metrics, "max_jaccard_mean") else NA_real_),
    median_min_Jaccard = as.character(if (!is.null(metrics)) scd_metric_value(metrics, "min_Jaccard") else NA_real_),
    median_purity = as.character(if (!is.null(metrics)) scd_metric_value(metrics, "purity") else NA_real_),
    primary_metric = scd_scalar(target$primary_metric, "ARI"),
    primary_metric_value = as.character(metric_value),
    pass_threshold = scd_scalar(target$pass_threshold),
    warn_threshold = scd_scalar(target$warn_threshold),
    fail_threshold = scd_scalar(target$fail_threshold),
    gate_status = gate_status,
    runtime_s = as.character(round(runtime_s, 2)),
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

scd_engine_status_row <- function(target, resolved_input_path, n_sim_requested, n_sim_done, fit_ok, sim_ok, score_ok, status, reason, runtime_s, checkpoint_path) {
  data.frame(
    target_id = scd_scalar(target$target_id),
    target_type = scd_scalar(target$target_type),
    layer_id = scd_scalar(target$layer_id),
    input_object = scd_scalar(target$input_object),
    resolved_input_path = scd_scalar(resolved_input_path),
    n_simulations_requested = as.character(n_sim_requested),
    n_simulations_done = as.character(n_sim_done),
    fit_ok = ifelse(fit_ok, "yes", "no"),
    sim_ok = ifelse(sim_ok, "yes", "no"),
    score_ok = ifelse(score_ok, "yes", "no"),
    status = status,
    reason = reason,
    runtime_s = as.character(round(runtime_s, 2)),
    checkpoint_path = checkpoint_path,
    stringsAsFactors = FALSE
  )
}

scd_save_placeholder_png <- function(path, title, message) {
  ensure_dir(dirname(path))
  grDevices::png(path, width = 1200, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot.new()
  graphics::title(main = title)
  graphics::text(0.5, 0.55, message, cex = 0.9)
}

scd_write_engine_figures <- function(target_id, figure_root, per_sim, per_label, pseudobulk_df = NULL, reason = "") {
  target_dir <- file.path(figure_root, target_id)
  ensure_dir(target_dir)
  ari_path <- file.path(target_dir, "ARI_distribution.png")
  confusion_path <- file.path(target_dir, "recluster_confusion_matrix.png")
  pseudo_path <- file.path(target_dir, "pseudobulk_synth_vs_real.png")

  ok_sim <- per_sim
  if (nrow(ok_sim) > 0) {
    ok_sim$ARI_num <- suppressWarnings(as.numeric(ok_sim$ARI))
    ok_sim <- ok_sim[is.finite(ok_sim$ARI_num), , drop = FALSE]
  }
  if (nrow(ok_sim) > 0) {
    grDevices::png(ari_path, width = 1200, height = 900, res = 150)
    graphics::boxplot(ARI_num ~ resolution, data = ok_sim, xlab = "Resolution", ylab = "ARI", main = paste(target_id, "ARI distribution"))
    grDevices::dev.off()
  } else {
    scd_save_placeholder_png(ari_path, paste(target_id, "ARI distribution"), reason %||% "No scored simulations.")
  }

  ok_label <- per_label
  if (nrow(ok_label) > 0) {
    ok_label$n_overlap_num <- suppressWarnings(as.numeric(ok_label$n_overlap))
    ok_label <- ok_label[is.finite(ok_label$n_overlap_num) & nzchar(ok_label$best_matching_cluster), , drop = FALSE]
  }
  if (nrow(ok_label) > 0) {
    best_sim <- ok_label$simulation_id[[1]]
    best_res <- ok_label$resolution[[1]]
    mat_df <- ok_label[ok_label$simulation_id == best_sim & ok_label$resolution == best_res, , drop = FALSE]
    mat <- xtabs(n_overlap_num ~ truth_label + best_matching_cluster, data = mat_df)
    grDevices::png(confusion_path, width = 1200, height = 900, res = 150)
    graphics::image(t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, main = paste(target_id, "truth vs recluster"))
    graphics::axis(1, at = seq(0, 1, length.out = ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.6)
    graphics::axis(2, at = seq(0, 1, length.out = nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.6)
    grDevices::dev.off()
  } else {
    scd_save_placeholder_png(confusion_path, paste(target_id, "truth vs recluster"), reason %||% "No label recovery rows.")
  }

  if (!is.null(pseudobulk_df) && nrow(pseudobulk_df) > 0) {
    grDevices::png(pseudo_path, width = 1200, height = 900, res = 150)
    graphics::plot(
      pseudobulk_df$real_mean,
      pseudobulk_df$synth_mean,
      pch = 16,
      cex = 0.5,
      xlab = "Real mean expression",
      ylab = "Synthetic mean expression",
      main = paste(target_id, "pseudobulk synthetic vs real")
    )
    graphics::abline(0, 1, col = "red")
    grDevices::dev.off()
  } else {
    scd_save_placeholder_png(pseudo_path, paste(target_id, "pseudobulk synthetic vs real"), reason %||% "No synthetic counts available.")
  }

  c(ari_distribution_png = ari_path, recluster_confusion_matrix_png = confusion_path, pseudobulk_synth_vs_real_png = pseudo_path)
}

scd_skip_result <- function(target, resolved_input_path, n_sim_requested, status, reason, runtime_s, checkpoint_path, figure_root) {
  metrics <- scd_target_metrics_row(target, n_sim_requested, 0L, status, reason, runtime_s, gate_status = "PLANNED")
  engine <- scd_engine_status_row(target, resolved_input_path, n_sim_requested, 0L, FALSE, FALSE, FALSE, status, reason, runtime_s, checkpoint_path)
  figures <- scd_write_engine_figures(scd_scalar(target$target_id), figure_root, scd_empty_df(scd_per_simulation_cols), scd_empty_df(scd_per_label_cols), NULL, reason)
  list(
    target_metrics = metrics,
    per_simulation = scd_empty_df(scd_per_simulation_cols),
    per_label = scd_empty_df(scd_per_label_cols),
    engine_status = engine,
    figures = figures
  )
}

run_cluster_target <- function(target, engine_cfg) {
  start <- proc.time()[["elapsed"]]
  target_id <- scd_scalar(target$target_id)
  target_type <- scd_scalar(target$target_type)
  target_status <- scd_scalar(target$status, "active")
  resolved_input_path <- scd_scalar(target$resolved_input_path)
  n_sim_requested <- scd_as_integer(Sys.getenv("SCDESIGN3_ENGINE_N_SIM", unset = ""), scd_as_integer(target$n_simulations, engine_cfg$n_simulations_default))
  checkpoint_path <- file.path(engine_cfg$checkpoint_dir, paste0(target_id, ".rds"))

  if (!identical(target_type, "cluster_robustness")) {
    return(scd_skip_result(target, resolved_input_path, n_sim_requested, "unsupported_m1", "M1 implements cluster_robustness only; this target type is reserved for a later milestone.", proc.time()[["elapsed"]] - start, checkpoint_path, engine_cfg$figure_root))
  }
  if (!identical(target_status, "active")) {
    return(scd_skip_result(target, resolved_input_path, n_sim_requested, "target_not_active", "Target is planned or disabled; D7 keeps TC layer targets waiting for formal input.", proc.time()[["elapsed"]] - start, checkpoint_path, engine_cfg$figure_root))
  }
  if (!scd_safe_file_exists(resolved_input_path)) {
    return(scd_skip_result(target, resolved_input_path, n_sim_requested, "waiting_input", sprintf("No formal annotated object resolved for layer `%s`.", scd_scalar(target$layer_id)), proc.time()[["elapsed"]] - start, checkpoint_path, engine_cfg$figure_root))
  }

  missing <- scd_missing_engine_packages()
  if (length(missing) > 0) {
    return(scd_skip_result(target, resolved_input_path, n_sim_requested, "missing_engine_dependency", sprintf("Missing R packages: %s", paste(missing, collapse = ", ")), proc.time()[["elapsed"]] - start, checkpoint_path, engine_cfg$figure_root))
  }

  prepared <- tryCatch(
    scd_prepare_target_data(
      resolved_input_path,
      scd_scalar(target$truth_col),
      engine_cfg$max_cells_per_label,
      engine_cfg$n_hvg,
      engine_cfg$seed
    ),
    error = function(e) e
  )
  if (inherits(prepared, "error")) {
    return(scd_skip_result(target, resolved_input_path, n_sim_requested, "input_prepare_failed", conditionMessage(prepared), proc.time()[["elapsed"]] - start, checkpoint_path, engine_cfg$figure_root))
  }

  fit <- tryCatch(
    fit_scdesign3(prepared, scd_scalar(target$truth_col), engine_cfg$family_use, engine_cfg$n_cores),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(scd_skip_result(target, resolved_input_path, n_sim_requested, "fit_failed", conditionMessage(fit), proc.time()[["elapsed"]] - start, checkpoint_path, engine_cfg$figure_root))
  }

  resolutions <- scd_parse_resolution_grid(target$resolution_grid, default = engine_cfg$resolution_default)
  per_sim_rows <- list()
  per_label_rows <- list()
  pseudobulk_df <- NULL

  for (sim_i in seq_len(n_sim_requested)) {
    sim_start <- proc.time()[["elapsed"]]
    synth <- tryCatch(
      simulate_synthetic_counts(fit, sim_i, engine_cfg$seed + sim_i),
      error = function(e) e
    )
    if (inherits(synth, "error")) {
      per_sim_rows[[length(per_sim_rows) + 1L]] <- data.frame(
        target_id = target_id,
        simulation_id = sprintf("%s_SIM%02d", target_id, sim_i),
        sim_i = as.character(sim_i),
        resolution = "",
        n_real_cells = as.character(prepared$n_real_cells),
        n_synth_cells = "0",
        n_truth_labels = as.character(prepared$n_truth_labels),
        n_recovered_clusters = "0",
        ARI = as.character(NA_real_),
        NMI = as.character(NA_real_),
        max_jaccard_mean = as.character(NA_real_),
        min_Jaccard = as.character(NA_real_),
        purity = as.character(NA_real_),
        runtime_s = as.character(round(proc.time()[["elapsed"]] - sim_start, 2)),
        status = "simulation_failed",
        reason = conditionMessage(synth),
        stringsAsFactors = FALSE
      )
      next
    }

    if (is.null(pseudobulk_df)) {
      common_genes <- intersect(rownames(prepared$counts), rownames(synth$counts))
      common_genes <- common_genes[seq_len(min(length(common_genes), 2000L))]
      if (length(common_genes) > 0) {
        real_mean <- rowMeans(as.matrix(prepared$counts[common_genes, , drop = FALSE]))
        synth_mean <- rowMeans(as.matrix(synth$counts[common_genes, , drop = FALSE]))
        pseudobulk_df <- data.frame(gene = common_genes, real_mean = real_mean, synth_mean = synth_mean, stringsAsFactors = FALSE)
      }
    }

    for (resolution in resolutions) {
      res_start <- proc.time()[["elapsed"]]
      clusters <- tryCatch(
        recluster_synthetic(synth$counts, resolution, prepared$n_truth_labels, engine_cfg$seed + sim_i, engine_cfg$n_pcs),
        error = function(e) e
      )
      if (inherits(clusters, "error")) {
        per_sim_rows[[length(per_sim_rows) + 1L]] <- data.frame(
          target_id = target_id,
          simulation_id = sprintf("%s_SIM%02d", target_id, sim_i),
          sim_i = as.character(sim_i),
          resolution = as.character(resolution),
          n_real_cells = as.character(prepared$n_real_cells),
          n_synth_cells = as.character(ncol(synth$counts)),
          n_truth_labels = as.character(prepared$n_truth_labels),
          n_recovered_clusters = "0",
          ARI = as.character(NA_real_),
          NMI = as.character(NA_real_),
          max_jaccard_mean = as.character(NA_real_),
          min_Jaccard = as.character(NA_real_),
          purity = as.character(NA_real_),
          runtime_s = as.character(round(proc.time()[["elapsed"]] - res_start, 2)),
          status = "recluster_failed",
          reason = conditionMessage(clusters),
          stringsAsFactors = FALSE
        )
        next
      }
      m <- compute_clustering_metrics(clusters, synth$truth)
      per_sim_rows[[length(per_sim_rows) + 1L]] <- data.frame(
        target_id = target_id,
        simulation_id = sprintf("%s_SIM%02d", target_id, sim_i),
        sim_i = as.character(sim_i),
        resolution = as.character(resolution),
        n_real_cells = as.character(prepared$n_real_cells),
        n_synth_cells = as.character(ncol(synth$counts)),
        n_truth_labels = as.character(prepared$n_truth_labels),
        n_recovered_clusters = as.character(length(unique(clusters))),
        ARI = as.character(unname(m[["ARI"]])),
        NMI = as.character(unname(m[["NMI"]])),
        max_jaccard_mean = as.character(unname(m[["max_jaccard_mean"]])),
        min_Jaccard = as.character(unname(m[["min_Jaccard"]])),
        purity = as.character(unname(m[["purity"]])),
        runtime_s = as.character(round(proc.time()[["elapsed"]] - res_start, 2)),
        status = "ok",
        reason = "",
        stringsAsFactors = FALSE
      )
      per_label_rows[[length(per_label_rows) + 1L]] <- compute_per_label_metrics(clusters, synth$truth, target_id, sim_i, resolution)
    }
  }

  per_sim <- scd_bind_rows(per_sim_rows, scd_per_simulation_cols)
  per_label <- scd_bind_rows(per_label_rows, scd_per_label_cols)
  ok <- per_sim[per_sim$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0) {
    reason <- "No simulations produced scoreable reclustering metrics."
    return(list(
      target_metrics = scd_target_metrics_row(target, n_sim_requested, 0L, "score_failed", reason, proc.time()[["elapsed"]] - start, "PLANNED"),
      per_simulation = per_sim,
      per_label = per_label,
      engine_status = scd_engine_status_row(target, resolved_input_path, n_sim_requested, 0L, TRUE, FALSE, FALSE, "score_failed", reason, proc.time()[["elapsed"]] - start, checkpoint_path),
      figures = scd_write_engine_figures(target_id, engine_cfg$figure_root, per_sim, per_label, pseudobulk_df, reason)
    ))
  }

  ok$ARI_num <- suppressWarnings(as.numeric(ok$ARI))
  ok <- ok[is.finite(ok$ARI_num), , drop = FALSE]
  if (nrow(ok) == 0) {
    reason <- "Scored simulations did not produce finite ARI values."
    return(list(
      target_metrics = scd_target_metrics_row(target, n_sim_requested, 0L, "score_failed", reason, proc.time()[["elapsed"]] - start, "PLANNED"),
      per_simulation = per_sim,
      per_label = per_label,
      engine_status = scd_engine_status_row(target, resolved_input_path, n_sim_requested, 0L, TRUE, TRUE, FALSE, "score_failed", reason, proc.time()[["elapsed"]] - start, checkpoint_path),
      figures = scd_write_engine_figures(target_id, engine_cfg$figure_root, per_sim, per_label, pseudobulk_df, reason)
    ))
  }
  median_by_res <- tapply(ok$ARI_num, ok$resolution, median, na.rm = TRUE)
  best_resolution <- names(which.max(median_by_res))[[1]]
  best <- ok[ok$resolution == best_resolution, , drop = FALSE]
  n_sim_done <- length(unique(best$simulation_id[is.finite(best$ARI_num)]))
  aggregate_metrics <- c(
    ARI = median(suppressWarnings(as.numeric(best$ARI)), na.rm = TRUE),
    iqr_ARI = stats::IQR(suppressWarnings(as.numeric(best$ARI)), na.rm = TRUE),
    NMI = median(suppressWarnings(as.numeric(best$NMI)), na.rm = TRUE),
    iqr_NMI = stats::IQR(suppressWarnings(as.numeric(best$NMI)), na.rm = TRUE),
    max_jaccard_mean = median(suppressWarnings(as.numeric(best$max_jaccard_mean)), na.rm = TRUE),
    min_Jaccard = median(suppressWarnings(as.numeric(best$min_Jaccard)), na.rm = TRUE),
    purity = median(suppressWarnings(as.numeric(best$purity)), na.rm = TRUE)
  )
  aggregate_metrics[!is.finite(aggregate_metrics)] <- NA_real_
  gate <- scd_gate_decision(aggregate_metrics, scd_scalar(target$pass_threshold), scd_scalar(target$warn_threshold))
  runtime_s <- proc.time()[["elapsed"]] - start
  list(
    target_metrics = scd_target_metrics_row(target, n_sim_requested, n_sim_done, "ok", "", runtime_s, gate, aggregate_metrics, best_resolution),
    per_simulation = per_sim,
    per_label = per_label,
    engine_status = scd_engine_status_row(target, resolved_input_path, n_sim_requested, n_sim_done, TRUE, TRUE, TRUE, "ok", "", runtime_s, checkpoint_path),
    figures = scd_write_engine_figures(target_id, engine_cfg$figure_root, per_sim, per_label, pseudobulk_df, "")
  )
}

compute_composition_metrics <- function(...) {
  stop("composition_robustness is planned for scDesign3 M2; M1 implements cluster_robustness only.", call. = FALSE)
}

compute_trajectory_metrics <- function(...) {
  stop("trajectory_robustness is planned for scDesign3 M3 after module 09 outputs are formalized.", call. = FALSE)
}

aggregate_by_target <- function(per_simulation_df) {
  if (nrow(per_simulation_df) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  per_simulation_df$ARI_num <- suppressWarnings(as.numeric(per_simulation_df$ARI))
  aggregate(ARI_num ~ target_id + resolution, per_simulation_df, median, na.rm = TRUE)
}

gate_decision <- scd_gate_decision
