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
  synthetic_sc_counts_rds <- file.path(out_dir, "synthetic_sc_counts.rds")
  synthetic_spot_counts_rds <- file.path(out_dir, "synthetic_spot_counts.rds")
  synthetic_st_rds <- file.path(out_dir, "synthetic_st.rds")
  saveRDS(spot_counts$counts, synthetic_spot_counts_rds)
  synthetic_st <- list(counts = spot_counts$counts, truth = truth)
  saveRDS(synthetic_st, synthetic_st_rds)
  list(synthetic_sc_counts_rds = synthetic_sc_counts_rds, synthetic_spot_counts_rds = synthetic_spot_counts_rds, synthetic_st_rds = synthetic_st_rds)
}

st07f_synthetic_coords <- function(spot_ids) {
  n <- length(spot_ids)
  side <- ceiling(sqrt(n))
  data.frame(
    spot_id = spot_ids,
    x = ((seq_len(n) - 1L) %% side) * 100,
    y = floor((seq_len(n) - 1L) / side) * 100,
    stringsAsFactors = FALSE
  )
}

st07f_write_mtx_bundle <- function(counts, obs, out_dir, prefix, obs_id_col) {
  if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) {
    return(list(mtx = "", genes = "", obs = ""))
  }
  ensure_dir(out_dir)
  counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is required to write synthetic Matrix Market artifacts.", call. = FALSE)
  }
  mtx <- file.path(out_dir, sprintf("%s_counts.mtx", prefix))
  genes <- file.path(out_dir, sprintf("%s_genes.tsv", prefix))
  obs_path <- file.path(out_dir, sprintf("%s_obs.tsv", prefix))
  Matrix::writeMM(Matrix::Matrix(counts, sparse = TRUE), mtx)
  gene_df <- data.frame(gene_id = rownames(counts), gene_name = rownames(counts), stringsAsFactors = FALSE)
  if (!obs_id_col %in% colnames(obs)) {
    obs[[obs_id_col]] <- colnames(counts)
  }
  st07_write_tsv(gene_df, genes)
  st07_write_tsv(obs, obs_path)
  list(mtx = mtx, genes = genes, obs = obs_path)
}

st07f_write_synthetic_h5ad_artifacts <- function(cfg, synthetic_cells, spot_counts, truth, out_dir, validation_id) {
  manifest_tsv <- file.path(out_dir, "synthetic_h5ad_manifest.tsv")
  ref_h5ad <- file.path(out_dir, "synthetic_reference.h5ad")
  st_h5ad <- file.path(out_dir, "synthetic_spatial.h5ad")
  coords_tsv <- file.path(out_dir, "synthetic_spatial_coordinates.tsv")
  spot_obs_tsv <- file.path(out_dir, "synthetic_spot_metadata.tsv")
  ref_obs <- synthetic_cells$meta
  ref_obs$validation_id <- validation_id
  spot_ids <- colnames(spot_counts$counts)
  spot_obs <- data.frame(spot_id = spot_ids, section_id = "synthetic", validation_id = validation_id, stringsAsFactors = FALSE)
  coords <- st07f_synthetic_coords(spot_ids)
  st07_write_tsv(coords, coords_tsv)
  st07_write_tsv(spot_obs, spot_obs_tsv)
  mtx_dir <- file.path(out_dir, "synthetic_h5ad_inputs")
  ref_bundle <- st07f_write_mtx_bundle(synthetic_cells$counts, ref_obs, mtx_dir, "synthetic_reference", "synthetic_cell_id")
  st_bundle <- st07f_write_mtx_bundle(spot_counts$counts, spot_obs, mtx_dir, "synthetic_spatial", "spot_id")
  writer <- file.path(cfg$pipeline_root, "workflow", "04python", "write_synthetic_h5ad.py")
  args <- c(
    writer,
    "--reference-mtx", ref_bundle$mtx,
    "--reference-genes", ref_bundle$genes,
    "--reference-obs", ref_bundle$obs,
    "--spatial-mtx", st_bundle$mtx,
    "--spatial-genes", st_bundle$genes,
    "--spatial-obs", st_bundle$obs,
    "--spatial-coords", coords_tsv,
    "--truth-tsv", file.path(out_dir, "synthetic_truth.tsv"),
    "--out-dir", out_dir,
    "--validation-id", validation_id
  )
  run <- tryCatch(system2(Sys.getenv("PY_SPATIAL_BIN", unset = "python3"), args = args, stdout = TRUE, stderr = TRUE), error = function(e) structure(conditionMessage(e), status = 127))
  code <- attr(run, "status")
  if (is.null(code)) code <- 0L
  if (!file.exists(manifest_tsv)) {
    status <- if (identical(as.integer(code), 127L)) "skipped_no_python_env" else "failed_write_h5ad"
    reason <- paste(as.character(run), collapse = " ")
    fallback <- data.frame(
      artifact_id = paste(validation_id, c("synthetic_reference", "synthetic_spatial"), sep = ":"),
      artifact_role = c("synthetic_reference", "synthetic_spatial"),
      h5ad_path = c(ref_h5ad, st_h5ad),
      source_object = "scDesign3_synthetic",
      n_obs = 0L,
      n_vars = 0L,
      obs_required_cols = c("synthetic_cell_id,cell_type,validation_id", "spot_id,section_id,validation_id"),
      obsm_required_keys = c("", "spatial"),
      var_required_cols = "gene_id,gene_name",
      fingerprint = "",
      status = status,
      reason = reason,
      stringsAsFactors = FALSE
    )
    st07_write_tsv(fallback, manifest_tsv)
  }
  list(
    synthetic_reference_h5ad = ref_h5ad,
    synthetic_spatial_h5ad = st_h5ad,
    synthetic_h5ad_manifest_tsv = manifest_tsv,
    synthetic_spatial_coordinates_tsv = coords_tsv,
    synthetic_spot_metadata_tsv = spot_obs_tsv,
    h5ad_status = if (file.exists(ref_h5ad) && file.exists(st_h5ad)) "ok" else "not_available"
  )
}

st07f_make_synthetic_seurat_objects <- function(synthetic_cells, spot_counts, validation_id) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for synthetic R deconvolution adapters.", call. = FALSE)
  }
  ref <- Seurat::CreateSeuratObject(counts = synthetic_cells$counts, assay = "RNA", min.cells = 0, min.features = 0)
  ref_meta <- synthetic_cells$meta
  rownames(ref_meta) <- ref_meta$synthetic_cell_id
  ref <- Seurat::AddMetaData(ref, ref_meta[colnames(ref), , drop = FALSE])
  st <- Seurat::CreateSeuratObject(counts = spot_counts$counts, assay = "RNA", min.cells = 0, min.features = 0)
  coords <- st07f_synthetic_coords(colnames(spot_counts$counts))
  rownames(coords) <- coords$spot_id
  st_meta <- data.frame(section_id = "synthetic", validation_id = validation_id, x = coords[colnames(st), "x"], y = coords[colnames(st), "y"], row.names = colnames(st), stringsAsFactors = FALSE)
  st <- Seurat::AddMetaData(st, st_meta)
  list(reference = ref, spatial = st, ref_labels = as.factor(ref$cell_type))
}

st07f_read_synthetic_h5ad_for_r <- function(paths, synthetic_cells, spot_counts, validation_id) {
  if (file.exists(paths$synthetic_reference_h5ad %||% "") && file.exists(paths$synthetic_spatial_h5ad %||% "") && requireNamespace("zellkonverter", quietly = TRUE) && requireNamespace("SingleCellExperiment", quietly = TRUE) && requireNamespace("Seurat", quietly = TRUE)) {
    loaded <- tryCatch({
      ref_sce <- zellkonverter::readH5AD(paths$synthetic_reference_h5ad)
      st_sce <- zellkonverter::readH5AD(paths$synthetic_spatial_h5ad)
      ref <- Seurat::as.Seurat(ref_sce, counts = "X", data = NULL)
      st <- Seurat::as.Seurat(st_sce, counts = "X", data = NULL)
      if (!"cell_type" %in% colnames(ref@meta.data)) stop("synthetic reference H5AD obs missing cell_type", call. = FALSE)
      st@meta.data$section_id <- st@meta.data$section_id %||% "synthetic"
      list(reference = ref, spatial = st, ref_labels = as.factor(ref@meta.data$cell_type))
    }, error = function(e) e)
    if (!inherits(loaded, "error")) {
      return(loaded)
    }
  }
  st07f_make_synthetic_seurat_objects(synthetic_cells, spot_counts, validation_id)
}

st07f_empty_synthetic_deconv_manifest <- function() {
  data.frame(
    validation_id = character(),
    method = character(),
    status = character(),
    reason = character(),
    prediction_source = character(),
    n_spots = integer(),
    n_celltypes = integer(),
    proportion_tsv = character(),
    runtime_sec = numeric(),
    stringsAsFactors = FALSE
  )
}

st07f_run_synthetic_r_method <- function(tool, package_name, fail_status, cfg, synthetic_inputs, table_dir, validation_id) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    return(list(row = data.frame(validation_id = validation_id, method = tool, status = "skipped_no_packages", reason = sprintf("%s is not installed", package_name), prediction_source = "synthetic_r_adapter", n_spots = 0L, n_celltypes = 0L, proportion_tsv = "", runtime_sec = NA_real_, stringsAsFactors = FALSE), props = st07_empty_df(c("method", "deconv_id", "section", "spot_id", "cell_type", "proportion"))))
  }
  run <- tryCatch({
    pair <- data.frame(deconv_id = validation_id, enabled = "yes", tool = tool, stringsAsFactors = FALSE)
    timing <- system.time(result <- st07_run_deconv_tool(tool, synthetic_inputs$spatial, synthetic_inputs$reference, synthetic_inputs$ref_labels, cfg))
    outputs <- st07_write_deconv_result_outputs(cfg, table_dir, pair, "synthetic", tool, result$prop, result$object, unname(timing[["elapsed"]]))
    props <- st07_read_tsv(outputs$proportion_tsv)
    props$method <- tool
    props$deconv_id <- validation_id
    props$section <- "synthetic"
    list(
      row = data.frame(validation_id = validation_id, method = tool, status = "ok", reason = "", prediction_source = "synthetic_r_adapter", n_spots = nrow(result$prop), n_celltypes = ncol(result$prop), proportion_tsv = outputs$proportion_tsv, runtime_sec = unname(timing[["elapsed"]]), stringsAsFactors = FALSE),
      props = props[, c("method", "deconv_id", "section", "spot_id", "cell_type", "proportion"), drop = FALSE]
    )
  }, error = function(e) {
    list(row = data.frame(validation_id = validation_id, method = tool, status = fail_status, reason = conditionMessage(e), prediction_source = "synthetic_r_adapter", n_spots = 0L, n_celltypes = 0L, proportion_tsv = "", runtime_sec = NA_real_, stringsAsFactors = FALSE), props = st07_empty_df(c("method", "deconv_id", "section", "spot_id", "cell_type", "proportion")))
  })
  run
}

st07f_run_synthetic_cell2location <- function(cfg, paths, table_dir, validation_id) {
  ref_h5ad <- paths$synthetic_reference_h5ad %||% ""
  st_h5ad <- paths$synthetic_spatial_h5ad %||% ""
  if (!file.exists(ref_h5ad) || !file.exists(st_h5ad)) {
    return(list(row = data.frame(validation_id = validation_id, method = "cell2location", status = "skipped_no_synthetic_h5ad", reason = "synthetic_reference.h5ad or synthetic_spatial.h5ad is unavailable", prediction_source = "synthetic_h5ad_rerun", n_spots = 0L, n_celltypes = 0L, proportion_tsv = "", runtime_sec = NA_real_, stringsAsFactors = FALSE), props = st07_empty_df(c("method", "deconv_id", "section", "spot_id", "cell_type", "proportion"))))
  }
  out_dir <- file.path(table_dir, "cell2location")
  ensure_dir(out_dir)
  sidecar <- file.path(cfg$pipeline_root, "workflow", "04python", "cell2location_pipeline.py")
  args <- c(
    sidecar,
    "--ref-h5ad", ref_h5ad,
    "--st-h5ad", st_h5ad,
    "--annotation-col", "cell_type",
    "--out-dir", out_dir,
    "--ref-epochs", as.character(cfg$spatial_c2l_ref_epochs),
    "--st-epochs", as.character(cfg$spatial_c2l_st_epochs),
    "--detection-alpha", as.character(cfg$spatial_c2l_detection_alpha),
    "--min-shared-genes", "2"
  )
  timing <- system.time(run <- tryCatch(system2(cfg$py_cell2location_bin, args = args, stdout = TRUE, stderr = TRUE), error = function(e) structure(conditionMessage(e), status = 127)))
  code <- attr(run, "status")
  if (is.null(code)) code <- 0L
  text <- paste(as.character(run), collapse = " ")
  status <- if (identical(as.integer(code), 0L) && file.exists(file.path(out_dir, "spot_celltype_proportions.tsv"))) "ok" else if (grepl("cell2location Python stack import failed", text, fixed = TRUE) || identical(as.integer(code), 20L) || identical(as.integer(code), 127L)) "skipped_no_cell2location" else "failed_sidecar"
  props <- if (identical(status, "ok")) st07_read_tsv(file.path(out_dir, "spot_celltype_proportions.tsv")) else st07_empty_df(c("spot_id", "cell_type", "proportion"))
  if (nrow(props) > 0) {
    props$method <- "cell2location"
    props$deconv_id <- validation_id
    props$section <- "synthetic"
    props <- props[, c("method", "deconv_id", "section", "spot_id", "cell_type", "proportion"), drop = FALSE]
  }
  row <- data.frame(validation_id = validation_id, method = "cell2location", status = status, reason = if (identical(status, "ok")) "" else text, prediction_source = "synthetic_h5ad_rerun", n_spots = length(unique(props$spot_id %||% character())), n_celltypes = length(unique(props$cell_type %||% character())), proportion_tsv = if (identical(status, "ok")) file.path(out_dir, "spot_celltype_proportions.tsv") else "", runtime_sec = unname(timing[["elapsed"]]), stringsAsFactors = FALSE)
  list(row = row, props = props)
}

st07f_run_deconv_methods_on_synthetic <- function(cfg, synthetic_cells, spot_counts, paths, validation_id) {
  manifest_tsv <- file.path(dirname(paths$synthetic_h5ad_manifest_tsv), "synthetic_deconv_manifest.tsv")
  table_dir <- file.path(dirname(paths$synthetic_h5ad_manifest_tsv), "synthetic_deconv")
  rows <- list()
  props <- list()
  synthetic_inputs <- tryCatch(st07f_read_synthetic_h5ad_for_r(paths, synthetic_cells, spot_counts, validation_id), error = function(e) e)
  if (inherits(synthetic_inputs, "error")) {
    for (tool in c("rctd", "transfer", "card")) {
      rows[[length(rows) + 1L]] <- data.frame(validation_id = validation_id, method = tool, status = "skipped_no_synthetic_inputs", reason = conditionMessage(synthetic_inputs), prediction_source = "synthetic_r_adapter", n_spots = 0L, n_celltypes = 0L, proportion_tsv = "", runtime_sec = NA_real_, stringsAsFactors = FALSE)
    }
  } else {
    specs <- list(
      rctd = list(package = "spacexr", fail = "failed_rctd_run"),
      transfer = list(package = "Seurat", fail = "failed_transfer_data"),
      card = list(package = "CARD", fail = "failed_card_deconvolution")
    )
    for (tool in names(specs)) {
      result <- st07f_run_synthetic_r_method(tool, specs[[tool]]$package, specs[[tool]]$fail, cfg, synthetic_inputs, file.path(table_dir, tool), validation_id)
      rows[[length(rows) + 1L]] <- result$row
      if (nrow(result$props) > 0) props[[length(props) + 1L]] <- result$props
    }
  }
  c2l <- st07f_run_synthetic_cell2location(cfg, paths, table_dir, validation_id)
  rows[[length(rows) + 1L]] <- c2l$row
  if (nrow(c2l$props) > 0) props[[length(props) + 1L]] <- c2l$props
  manifest <- if (length(rows) == 0) st07f_empty_synthetic_deconv_manifest() else do.call(rbind, rows)
  st07_write_tsv(manifest, manifest_tsv)
  pred <- if (length(props) == 0) st07_empty_df(c("method", "deconv_id", "section", "spot_id", "cell_type", "proportion")) else do.call(rbind, props)
  h5ad_manifest <- st07_read_tsv(paths$synthetic_h5ad_manifest_tsv %||% "")
  h5ad_manifest_ok <- nrow(h5ad_manifest) > 0 && "status" %in% colnames(h5ad_manifest) && all(h5ad_manifest$status == "ok")
  ok_methods <- manifest$method[manifest$status == "ok"]
  ok_has_native <- any(ok_methods == "cell2location")
  ok_has_adapter <- any(ok_methods %in% c("rctd", "transfer", "card"))
  prediction_source <- if (length(ok_methods) == 0) {
    "none"
  } else if (ok_has_native && ok_has_adapter) {
    "mixed_synthetic_rerun"
  } else if (ok_has_native) {
    "synthetic_h5ad_rerun"
  } else {
    "synthetic_r_adapter_fallback"
  }
  list(
    manifest = manifest,
    manifest_tsv = manifest_tsv,
    props = pred,
    synthetic_prediction_source = prediction_source,
    synthetic_rerun_methods = paste(manifest$method, collapse = ","),
    synthetic_rerun_ok_methods = sum(manifest$status == "ok"),
    synthetic_h5ad_manifest_ok = h5ad_manifest_ok
  )
}

st07f_write_validation_outputs <- function(synthetic_cells, spot_counts, truth, out_dir) {
  synthetic_sc_metadata_tsv <- file.path(out_dir, "synthetic_sc_metadata.tsv")
  st07_write_tsv(synthetic_cells$meta, synthetic_sc_metadata_tsv)
  synthetic_sc_counts_rds <- file.path(out_dir, "synthetic_sc_counts.rds")
  if (identical(synthetic_cells$status, "ok")) saveRDS(synthetic_cells$counts, synthetic_sc_counts_rds)
  object_paths <- if (identical(spot_counts$status, "ok")) st07f_build_synthetic_st_object(spot_counts, truth, out_dir) else list(synthetic_sc_counts_rds = synthetic_sc_counts_rds, synthetic_spot_counts_rds = "", synthetic_st_rds = "")
  object_paths$synthetic_sc_counts_rds <- synthetic_sc_counts_rds
  c(list(synthetic_sc_metadata_tsv = synthetic_sc_metadata_tsv), object_paths)
}

st07f_run_scdesign3_validation <- function(cfg, props, out_dir) {
  validation_id <- basename(out_dir)
  reference <- st07f_load_reference(cfg)
  fit_bundle <- st07f_fit_or_reuse_scdesign3(reference, cfg)
  synthetic_cells <- st07f_simulate_cells(fit_bundle, cfg)
  truth <- st07f_sample_truth_proportions(synthetic_cells$meta$cell_type, cfg)
  spot_counts <- st07f_mix_cells_to_spots(synthetic_cells, truth, cfg)
  paths <- st07f_write_validation_outputs(synthetic_cells, spot_counts, truth, out_dir)
  h5ad_paths <- if (identical(synthetic_cells$status, "ok") && identical(spot_counts$status, "ok")) {
    tryCatch(
      st07f_write_synthetic_h5ad_artifacts(cfg, synthetic_cells, spot_counts, truth, out_dir, validation_id),
      error = function(e) {
        manifest_tsv <- file.path(out_dir, "synthetic_h5ad_manifest.tsv")
        st07_write_tsv(data.frame(artifact_id = paste(validation_id, c("synthetic_reference", "synthetic_spatial"), sep = ":"), artifact_role = c("synthetic_reference", "synthetic_spatial"), h5ad_path = file.path(out_dir, c("synthetic_reference.h5ad", "synthetic_spatial.h5ad")), source_object = "scDesign3_synthetic", n_obs = 0L, n_vars = 0L, obs_required_cols = c("synthetic_cell_id,cell_type,validation_id", "spot_id,section_id,validation_id"), obsm_required_keys = c("", "spatial"), var_required_cols = "gene_id,gene_name", fingerprint = "", status = "failed_write_h5ad", reason = conditionMessage(e), stringsAsFactors = FALSE), manifest_tsv)
        list(synthetic_reference_h5ad = file.path(out_dir, "synthetic_reference.h5ad"), synthetic_spatial_h5ad = file.path(out_dir, "synthetic_spatial.h5ad"), synthetic_h5ad_manifest_tsv = manifest_tsv, h5ad_status = "not_available")
      }
    )
  } else {
    manifest_tsv <- file.path(out_dir, "synthetic_h5ad_manifest.tsv")
    st07_write_tsv(data.frame(artifact_id = paste(validation_id, c("synthetic_reference", "synthetic_spatial"), sep = ":"), artifact_role = c("synthetic_reference", "synthetic_spatial"), h5ad_path = file.path(out_dir, c("synthetic_reference.h5ad", "synthetic_spatial.h5ad")), source_object = "scDesign3_synthetic", n_obs = 0L, n_vars = 0L, obs_required_cols = c("synthetic_cell_id,cell_type,validation_id", "spot_id,section_id,validation_id"), obsm_required_keys = c("", "spatial"), var_required_cols = "gene_id,gene_name", fingerprint = "", status = synthetic_cells$status, reason = synthetic_cells$reason, stringsAsFactors = FALSE), manifest_tsv)
    list(synthetic_reference_h5ad = file.path(out_dir, "synthetic_reference.h5ad"), synthetic_spatial_h5ad = file.path(out_dir, "synthetic_spatial.h5ad"), synthetic_h5ad_manifest_tsv = manifest_tsv, h5ad_status = "not_available")
  }
  paths <- c(paths, h5ad_paths)
  rerun <- if (identical(synthetic_cells$status, "ok") && identical(spot_counts$status, "ok")) {
    st07f_run_deconv_methods_on_synthetic(cfg, synthetic_cells, spot_counts, paths, validation_id)
  } else {
    manifest_tsv <- file.path(out_dir, "synthetic_deconv_manifest.tsv")
    empty <- st07f_empty_synthetic_deconv_manifest()
    st07_write_tsv(empty, manifest_tsv)
    list(manifest = empty, manifest_tsv = manifest_tsv, props = st07_empty_df(c("method", "deconv_id", "section", "spot_id", "cell_type", "proportion")), synthetic_prediction_source = "none", synthetic_rerun_methods = "", synthetic_rerun_ok_methods = 0L, synthetic_h5ad_manifest_ok = FALSE)
  }
  list(
    status = if (identical(synthetic_cells$status, "ok") && identical(spot_counts$status, "ok")) "ok" else synthetic_cells$status,
    reason = if (identical(synthetic_cells$status, "ok") && identical(spot_counts$status, "ok")) "" else synthetic_cells$reason,
    truth = truth,
    props = rerun$props,
    synthetic_cell_generation = synthetic_cells$status,
    synthetic_spot_generation = spot_counts$status,
    synthetic_prediction_source = rerun$synthetic_prediction_source,
    synthetic_rerun_methods = rerun$synthetic_rerun_methods,
    synthetic_rerun_ok_methods = rerun$synthetic_rerun_ok_methods,
    synthetic_h5ad_manifest_ok = rerun$synthetic_h5ad_manifest_ok %||% FALSE,
    synthetic_deconv_manifest_tsv = rerun$manifest_tsv,
    paths = paths
  )
}
