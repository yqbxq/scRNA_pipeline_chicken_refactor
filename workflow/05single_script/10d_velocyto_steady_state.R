#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source(file.path(.script_dir, "helpers", "load_helpers_10.R"), encoding = "UTF-8")

load_required_packages(c("dplyr", "jsonlite", "Matrix"))

cfg <- get_single_script_config_10()
module_name <- "10d_velocyto_steady_state"
prepare_dirs_10(cfg)

normalize_velocity_cell_id_10d <- function(cell_id, sample_id = "") {
  cell_id <- as.character(cell_id)
  if (grepl(":", cell_id, fixed = TRUE)) {
    parts <- strsplit(cell_id, ":", fixed = TRUE)[[1]]
    sample_id <- parts[[1]]
    barcode <- paste(parts[-1], collapse = ":")
  } else {
    barcode <- cell_id
  }
  barcode <- trimws(barcode)
  if (grepl("x$", barcode)) {
    barcode <- sub("x$", "-1", barcode)
  } else if (!grepl("-[0-9]+$", barcode)) {
    barcode <- paste0(barcode, "-1")
  }
  if (!nzchar(sample_id)) {
    return(barcode)
  }
  paste0(sample_id, ":", barcode)
}

empty_method_row_10d <- function(unit, output_path, direction_path, status, reason, runtime_s = 0) {
  data.frame(
    pair_id = unit$pair_id[[1]],
    split_value = display_scalar_value(unit$split_value[[1]], "pooled"),
    method = "velocyto_steady_state",
    methods_enabled = ifelse(velocity_method_enabled_10(unit, "velocyto", default = TRUE), "yes", "no"),
    input_path = normalize_scalar_value(unit$metadata_csv[[1]]),
    output_path = if (identical(status, "ok")) output_path else "",
    extra_path = normalize_scalar_value(unit$umap_csv[[1]]),
    figure_path = "",
    n_cells = 0L,
    status = status,
    reason = reason,
    runtime_s = runtime_s,
    direction_tsv = direction_path,
    stringsAsFactors = FALSE
  )
}

read_velocity_cells_10d <- function(metadata_csv) {
  metadata_csv <- normalize_scalar_value(metadata_csv)
  if (!nzchar(metadata_csv) || !file.exists(metadata_csv)) {
    return(character(0))
  }
  meta <- utils::read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE)
  if ("cell_id" %in% colnames(meta)) {
    return(as.character(meta$cell_id))
  }
  as.character(meta[[1]])
}

sample_ids_from_cells_10d <- function(cells) {
  cells <- as.character(cells)
  has_sample <- grepl(":", cells, fixed = TRUE)
  unique(sub(":.*$", "", cells[has_sample]))
}

align_sparse_rows_10d <- function(mat, genes) {
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  missing <- setdiff(genes, rownames(mat))
  if (length(missing) > 0) {
    zero <- Matrix::Matrix(0, nrow = length(missing), ncol = ncol(mat), sparse = TRUE)
    rownames(zero) <- missing
    colnames(zero) <- colnames(mat)
    mat <- rbind(mat, zero)
  }
  mat[genes, , drop = FALSE]
}

combine_sparse_by_rows_10d <- function(mats) {
  mats <- mats[!vapply(mats, is.null, logical(1))]
  if (length(mats) == 0) {
    return(NULL)
  }
  genes <- Reduce(union, lapply(mats, rownames))
  mats <- lapply(mats, align_sparse_rows_10d, genes = genes)
  do.call(cbind, mats)
}

load_velocity_matrices_10d <- function(unit, cells, loom_index) {
  ok <- loom_index[loom_index$status %in% c("ok", "ok_existing") & file.exists(loom_index$loom_path), , drop = FALSE]
  samples <- sample_ids_from_cells_10d(cells)
  if (length(samples) > 0 && "sample_id" %in% colnames(ok)) {
    ok <- ok[ok$sample_id %in% samples, , drop = FALSE]
  }
  if (nrow(ok) == 0) {
    stop("no matching loom files for velocity reference cells", call. = FALSE)
  }

  spliced <- list()
  unspliced <- list()
  for (i in seq_len(nrow(ok))) {
    loom_path <- ok$loom_path[[i]]
    sample_id <- normalize_scalar_value(ok$sample_id[[i]], tools::file_path_sans_ext(basename(loom_path)))
    mats <- velocyto.R::read.loom.matrices(loom_path)
    if (!"spliced" %in% names(mats) || !"unspliced" %in% names(mats)) {
      next
    }
    colnames(mats$spliced) <- vapply(colnames(mats$spliced), normalize_velocity_cell_id_10d, character(1), sample_id = sample_id)
    colnames(mats$unspliced) <- vapply(colnames(mats$unspliced), normalize_velocity_cell_id_10d, character(1), sample_id = sample_id)
    keep <- intersect(cells, intersect(colnames(mats$spliced), colnames(mats$unspliced)))
    if (length(keep) == 0) {
      next
    }
    spliced[[length(spliced) + 1L]] <- mats$spliced[, keep, drop = FALSE]
    unspliced[[length(unspliced) + 1L]] <- mats$unspliced[, keep, drop = FALSE]
  }

  emat <- combine_sparse_by_rows_10d(spliced)
  nmat <- combine_sparse_by_rows_10d(unspliced)
  if (is.null(emat) || is.null(nmat) || ncol(emat) == 0 || ncol(nmat) == 0) {
    stop("loom matrices have no cells overlapping velocity metadata", call. = FALSE)
  }
  common_genes <- intersect(rownames(emat), rownames(nmat))
  common_cells <- intersect(colnames(emat), colnames(nmat))
  list(
    emat = emat[common_genes, common_cells, drop = FALSE],
    nmat = nmat[common_genes, common_cells, drop = FALSE]
  )
}

velocity_matrix_from_fit_10d <- function(vfit) {
  for (name in c("deltaE", "velocity", "vel")) {
    value <- vfit[[name]]
    if (!is.null(value) && length(dim(value)) == 2) {
      return(value)
    }
  }
  if (!is.null(vfit$current) && !is.null(vfit$projected) && all(dim(vfit$current) == dim(vfit$projected))) {
    return(vfit$projected - vfit$current)
  }
  NULL
}

write_direction_10d <- function(vfit, pair_id, split_value, path) {
  ensure_dir(dirname(path))
  vel <- velocity_matrix_from_fit_10d(vfit)
  if (is.null(vel) || is.null(colnames(vel))) {
    out <- data.frame(
      pair_id = pair_id,
      split_value = display_scalar_value(split_value, "pooled"),
      cell_id = character(0),
      velocity_length = numeric(0),
      status = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    len <- sqrt(Matrix::colSums(Matrix::Matrix(vel, sparse = TRUE)^2))
    out <- data.frame(
      pair_id = pair_id,
      split_value = display_scalar_value(split_value, "pooled"),
      cell_id = names(len),
      velocity_length = as.numeric(len),
      status = "ok",
      stringsAsFactors = FALSE
    )
  }
  write_tsv_local(out, path)
  out
}

run_unit_10d <- function(unit, loom_index) {
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  output_path <- velocity_velocyto_steady_rds_path_10(cfg, pair_id, split_value)
  direction_path <- velocity_velocyto_direction_path_10(cfg, pair_id, split_value)
  start <- proc.time()[["elapsed"]]

  if (!velocity_method_enabled_10(unit, "velocyto", default = TRUE)) {
    reason <- velocity_method_disabled_reason_10(unit, "velocyto")
    return(list(
      index = empty_method_row_10d(unit, output_path, direction_path, "skipped_disabled", reason),
      outputs = list()
    ))
  }
  if (!requireNamespace("velocyto.R", quietly = TRUE)) {
    return(list(
      index = empty_method_row_10d(unit, output_path, direction_path, "failed_no_package", "velocyto.R package is unavailable"),
      outputs = list()
    ))
  }

  result <- tryCatch({
    cells <- read_velocity_cells_10d(unit$metadata_csv[[1]])
    if (length(cells) < 3L) {
      stop("fewer than 3 velocity reference cells", call. = FALSE)
    }
    mats <- load_velocity_matrices_10d(unit, cells, loom_index)
    if (ncol(mats$emat) < 3L || nrow(mats$emat) < 3L) {
      stop("fewer than 3 cells or genes after loom/reference intersection", call. = FALSE)
    }
    k_cells <- min(25L, ncol(mats$emat) - 1L)
    vfit <- velocyto.R::gene.relative.velocity.estimates(
      emat = mats$emat,
      nmat = mats$nmat,
      deltaT = 1,
      kCells = k_cells,
      fit.quantile = 0.02,
      n.cores = max(1L, min(4L, cfg$velocyto_threads %||% 1L))
    )
    ensure_dir(dirname(output_path))
    saveRDS(
      list(
        pair_id = pair_id,
        split_value = split_value,
        emat_dim = dim(mats$emat),
        nmat_dim = dim(mats$nmat),
        result = vfit
      ),
      output_path
    )
    direction <- write_direction_10d(vfit, pair_id, split_value, direction_path)
    runtime_s <- proc.time()[["elapsed"]] - start
    index <- empty_method_row_10d(unit, output_path, direction_path, "ok", "", runtime_s)
    index$output_path <- output_path
    index$n_cells <- nrow(direction)
    list(ok = TRUE, index = index)
  }, error = function(e) {
    runtime_s <- proc.time()[["elapsed"]] - start
    list(ok = FALSE, index = empty_method_row_10d(unit, output_path, direction_path, "failed", conditionMessage(e), runtime_s))
  })

  outputs <- list()
  if (isTRUE(result$ok)) {
    outputs[[sprintf("velocyto_steady__%s", velocity_unit_file_id_10(pair_id, split_value))]] <- build_output_entry(
      output_path,
      "rds",
      module_name,
      "velocyto.R steady-state velocity result",
      base_dir = cfg$project_root
    )
    outputs[[sprintf("velocyto_steady_direction__%s", velocity_unit_file_id_10(pair_id, split_value))]] <- build_output_entry(
      direction_path,
      "tsv",
      module_name,
      "per-cell steady-state velocity direction summary",
      base_dir = cfg$project_root
    )
  }
  list(index = result$index, outputs = outputs)
}

units <- velocity_execution_units_10(cfg)
loom_index <- velocity_read_loom_index_10(cfg)

index_rows <- list()
dynamic_outputs <- list()
for (i in seq_len(nrow(units))) {
  unit <- units[i, , drop = FALSE]
  result <- run_unit_10d(unit, loom_index)
  index_rows[[length(index_rows) + 1L]] <- result$index
  if (length(result$outputs) > 0) {
    dynamic_outputs <- c(dynamic_outputs, result$outputs)
  }
}

index_df <- if (length(index_rows) > 0) {
  dplyr::bind_rows(index_rows)
} else {
  velocity_empty_df_10(velocity_velocyto_steady_index_cols_10)
}
write_tsv_local(index_df, cfg$velocity_velocyto_steady_index_tsv)

velocity_manifest_10(
  cfg,
  cfg$module_10d_manifest_path,
  module_name,
  outputs = c(
    list(
      velocyto_steady_index_tsv = build_output_entry(
        cfg$velocity_velocyto_steady_index_tsv,
        "tsv",
        module_name,
        "velocyto.R steady-state status by velocity unit",
        base_dir = cfg$project_root,
        schema = infer_schema_from_df(index_df)
      )
    ),
    dynamic_outputs
  ),
  inputs = list(
    module_10a = cfg$module_10a_manifest_path,
    module_10b = cfg$module_10b_manifest_path,
    velocity_loom_index = cfg$velocity_loom_index_tsv,
    velocity_reference_index = cfg$velocity_reference_index_tsv
  ),
  depends_on = list(
    module_10a = cfg$module_10a_manifest_path,
    module_10b = cfg$module_10b_manifest_path
  )
)

message("10d completed. velocyto steady-state index: ", cfg$velocity_velocyto_steady_index_tsv)
