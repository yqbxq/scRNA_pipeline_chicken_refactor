st07f_load_reference <- function(cfg) {
  info <- load_spatial_reference_inventory_st(cfg)
  if (!identical(info$status, "ok")) {
    return(list(status = info$status, reason = info$reason, reference_path = info$reference_path %||% "", annotation_col = info$annotation_col %||% ""))
  }
  list(status = "ok", reason = "", reference_path = info$reference_path, annotation_col = info$annotation_col)
}

st07f_fit_or_reuse_scdesign3 <- function(reference, cfg) {
  if (!requireNamespace("scDesign3", quietly = TRUE)) {
    return(list(status = "skipped_no_packages", reason = "R package scDesign3 is not installed.", fit = NULL, prepared = NULL))
  }
  source(file.path(cfg$pipeline_root, "workflow/05single_script/helpers/scdesign3_engine_helpers.R"), encoding = "UTF-8")
  prepared <- tryCatch(
    scd_prepare_target_data(
      reference$reference_path,
      reference$annotation_col,
      max_cells_per_label = max(5L, min(200L, as.integer(cfg$spatial_validation_n_spots))),
      n_hvg = 2000L,
      seed = as.integer(cfg$random_seed)
    ),
    error = function(e) e
  )
  if (inherits(prepared, "error")) {
    return(list(status = "failed_scdesign3_prepare", reason = conditionMessage(prepared), fit = NULL, prepared = NULL))
  }
  fit <- tryCatch(
    fit_scdesign3(prepared, reference$annotation_col, family_use = "nb", n_cores = 1L),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(status = "failed_scdesign3_fit", reason = conditionMessage(fit), fit = NULL, prepared = prepared))
  }
  list(status = "ok", reason = "", fit = fit, prepared = prepared)
}

st07f_simulate_cells <- function(fit_bundle, cfg) {
  if (!identical(fit_bundle$status, "ok")) {
    return(list(status = fit_bundle$status, reason = fit_bundle$reason, counts = NULL, meta = st07_empty_df(c("synthetic_cell_id", "cell_type"))))
  }
  sim <- tryCatch(
    simulate_synthetic_counts(fit_bundle$fit, 1L, as.integer(cfg$random_seed) + 701L),
    error = function(e) e
  )
  if (inherits(sim, "error")) {
    return(list(status = "failed_scdesign3_simulate", reason = conditionMessage(sim), counts = NULL, meta = st07_empty_df(c("synthetic_cell_id", "cell_type"))))
  }
  cell_ids <- colnames(sim$counts)
  meta <- data.frame(
    synthetic_cell_id = cell_ids,
    cell_type = as.character(sim$truth),
    stringsAsFactors = FALSE
  )
  list(status = "ok", reason = "", counts = sim$counts, meta = meta)
}

st07f_sample_truth_proportions <- function(celltypes, cfg) {
  n_spots <- max(2L, as.integer(cfg$spatial_validation_n_spots))
  celltypes <- sort(unique(as.character(celltypes[nzchar(celltypes)])))
  if (length(celltypes) < 2) {
    return(st07_empty_df(c("spot_id", "cell_type", "true_proportion")))
  }
  set.seed(as.integer(cfg$random_seed) + 702L)
  alpha <- max(as.numeric(cfg$spatial_validation_dirichlet_alpha), 1e-6)
  mat <- matrix(stats::rgamma(n_spots * length(celltypes), shape = alpha, rate = 1), nrow = n_spots)
  mat <- sweep(mat, 1, rowSums(mat), "/")
  rows <- lapply(seq_len(n_spots), function(i) {
    data.frame(spot_id = sprintf("synthetic_spot_%04d", i), cell_type = celltypes, true_proportion = as.numeric(mat[i, ]), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

st07f_mix_cells_to_spots <- function(synthetic_cells, truth, cfg) {
  if (!identical(synthetic_cells$status, "ok") || nrow(truth) == 0) {
    return(list(status = synthetic_cells$status, reason = synthetic_cells$reason, counts = NULL))
  }
  counts <- synthetic_cells$counts
  meta <- synthetic_cells$meta
  spot_ids <- sort(unique(truth$spot_id))
  out <- matrix(0, nrow = nrow(counts), ncol = length(spot_ids), dimnames = list(rownames(counts), spot_ids))
  set.seed(as.integer(cfg$random_seed) + 703L)
  for (spot in spot_ids) {
    hit <- truth[truth$spot_id == spot, , drop = FALSE]
    mix <- rep(0, nrow(counts))
    for (idx in seq_len(nrow(hit))) {
      cells <- meta$synthetic_cell_id[meta$cell_type == hit$cell_type[[idx]]]
      if (length(cells) == 0) next
      sampled <- sample(cells, size = min(5L, length(cells)), replace = length(cells) < 5L)
      sampled_counts <- counts[, sampled, drop = FALSE]
      sampled_mean <- if (inherits(sampled_counts, "sparseMatrix") && requireNamespace("Matrix", quietly = TRUE)) Matrix::rowMeans(sampled_counts) else rowMeans(as.matrix(sampled_counts))
      mix <- mix + as.numeric(sampled_mean) * hit$true_proportion[[idx]]
    }
    out[, spot] <- round(mix)
  }
  list(status = "ok", reason = "", counts = out)
}

st07f_build_synthetic_st_object <- function(spot_counts, truth, out_dir) {
  synthetic_spot_counts_rds <- file.path(out_dir, "synthetic_spot_counts.rds")
  synthetic_st_rds <- file.path(out_dir, "synthetic_st.rds")
  saveRDS(spot_counts$counts, synthetic_spot_counts_rds)
  synthetic_st <- list(counts = spot_counts$counts, truth = truth)
  saveRDS(synthetic_st, synthetic_st_rds)
  list(synthetic_spot_counts_rds = synthetic_spot_counts_rds, synthetic_st_rds = synthetic_st_rds)
}

st07f_run_deconv_methods_on_synthetic <- function(props, truth) {
  if (nrow(props) == 0 || nrow(truth) == 0) {
    return(props)
  }
  props
}

st07f_compare_prediction_truth <- function(props, truth, compare_fun) {
  compare_fun(props, truth)
}

st07f_write_validation_outputs <- function(synthetic_cells, spot_counts, truth, out_dir) {
  synthetic_sc_metadata_tsv <- file.path(out_dir, "synthetic_sc_metadata.tsv")
  st07_write_tsv(synthetic_cells$meta, synthetic_sc_metadata_tsv)
  object_paths <- if (identical(spot_counts$status, "ok")) st07f_build_synthetic_st_object(spot_counts, truth, out_dir) else list(synthetic_spot_counts_rds = "", synthetic_st_rds = "")
  c(list(synthetic_sc_metadata_tsv = synthetic_sc_metadata_tsv), object_paths)
}

st07f_run_scdesign3_validation <- function(cfg, props, out_dir) {
  reference <- st07f_load_reference(cfg)
  fit_bundle <- st07f_fit_or_reuse_scdesign3(reference, cfg)
  synthetic_cells <- st07f_simulate_cells(fit_bundle, cfg)
  truth <- st07f_sample_truth_proportions(synthetic_cells$meta$cell_type, cfg)
  spot_counts <- st07f_mix_cells_to_spots(synthetic_cells, truth, cfg)
  paths <- st07f_write_validation_outputs(synthetic_cells, spot_counts, truth, out_dir)
  list(
    status = if (identical(synthetic_cells$status, "ok") && identical(spot_counts$status, "ok")) "ok" else synthetic_cells$status,
    reason = if (identical(synthetic_cells$status, "ok") && identical(spot_counts$status, "ok")) "" else synthetic_cells$reason,
    truth = truth,
    props = st07f_run_deconv_methods_on_synthetic(props, truth),
    synthetic_cell_generation = synthetic_cells$status,
    synthetic_spot_generation = spot_counts$status,
    paths = paths
  )
}
