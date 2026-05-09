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

load_required_packages(c("Seurat", "dplyr", "jsonlite"))

cfg <- get_single_script_config_10()
module_name <- "10b_prepare_velocity_reference"
prepare_dirs_10(cfg)

write_velocity_csv_10b <- function(df, path) {
  ensure_dir(dirname(path))
  utils::write.csv(df, file = path, row.names = FALSE, na = "")
}

failed_reference_row_10b <- function(pair_row, layer_row = NULL, split_var = "", split_value = "", status = "failed", reason = "") {
  data.frame(
    pair_id = velocity_pair_id_10(pair_row),
    source_question_id = normalize_scalar_value(pair_row$source_question_id[[1]]),
    layer_scope = normalize_scalar_value(pair_row$layer_scope[[1]]),
    split_mode = normalize_scalar_value(pair_row$split_mode[[1]], "auto"),
    split_var = split_var,
    split_value = velocity_display_split_10(split_value),
    input_rds = if (is.null(layer_row)) "" else layer_row$object_rds[[1]],
    umap_csv = "",
    metadata_csv = "",
    n_cells = 0L,
    coarse_label_var = normalize_scalar_value(pair_row$coarse_label_var[[1]]),
    fine_label_var = normalize_scalar_value(pair_row$fine_label_var[[1]]),
    status = status,
    reason = reason,
    produced_at = timestamp_now(),
    stringsAsFactors = FALSE
  )
}

velocity_subset_object_10b <- function(source_obj, split_var, split_value) {
  if (!nzchar(split_var)) {
    return(source_obj)
  }
  values <- trimws(as.character(source_obj@meta.data[[split_var]]))
  cells <- rownames(source_obj@meta.data)[values == split_value]
  if (length(cells) == 0) {
    stop(sprintf("split %s=%s matched zero cells", split_var, split_value), call. = FALSE)
  }
  subset(source_obj, cells = cells)
}

velocity_metadata_export_10b <- function(seu, pair_row, split_var) {
  meta <- seu@meta.data
  meta$cell_id <- rownames(meta)

  coarse <- normalize_scalar_value(pair_row$coarse_label_var[[1]])
  fine <- normalize_scalar_value(pair_row$fine_label_var[[1]])
  required <- c(coarse, fine)
  required <- required[nzchar(required)]
  missing_required <- setdiff(required, colnames(meta))
  if (length(missing_required) > 0) {
    stop(
      sprintf("velocity reference metadata is missing label columns: %s", paste(missing_required, collapse = ",")),
      call. = FALSE
    )
  }

  keep_cols <- unique(c(
    "cell_id", "sample_id", "group_id", "condition", "batch",
    "cell_subtype", "cell_type", "seurat_clusters", "cluster",
    "Phase", "S.Score", "G2M.Score", split_var, coarse, fine
  ))
  keep_cols <- keep_cols[nzchar(keep_cols) & keep_cols %in% colnames(meta)]
  out <- meta[, keep_cols, drop = FALSE]

  if (!"cell_type" %in% colnames(out) && nzchar(coarse) && coarse %in% colnames(meta)) {
    out$cell_type <- as.character(meta[[coarse]])
  }
  if (!"cluster" %in% colnames(out) && nzchar(fine) && fine %in% colnames(meta)) {
    out$cluster <- as.character(meta[[fine]])
  }

  out[] <- lapply(out, function(x) {
    if (is.factor(x)) {
      as.character(x)
    } else {
      x
    }
  })
  out
}

velocity_umap_export_10b <- function(seu) {
  emb <- trajectory_umap_09(seu)
  if (is.null(emb) || nrow(emb) == 0) {
    stop("velocity reference object has no usable UMAP reduction", call. = FALSE)
  }
  cells <- intersect(colnames(seu), rownames(emb))
  if (length(cells) == 0) {
    stop("velocity reference UMAP has no overlap with object cells", call. = FALSE)
  }
  data.frame(
    cell_id = cells,
    UMAP_1 = emb[cells, 1],
    UMAP_2 = emb[cells, 2],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

prepare_one_velocity_reference_10b <- function(pair_row, layer_row, source_obj, split_var, split_value) {
  pair_id <- velocity_pair_id_10(pair_row)
  split_label <- velocity_display_split_10(split_value)
  message("10b velocity reference: ", pair_id, " split=", split_label)

  ref_obj <- velocity_subset_object_10b(source_obj, split_var, split_value)
  n_cells <- ncol(ref_obj)
  if (n_cells < cfg$velocity_min_reference_cells) {
    return(list(
      index = failed_reference_row_10b(
        pair_row,
        layer_row,
        split_var,
        split_value,
        "skipped_low_cells",
        sprintf("fewer than %s cells after split", cfg$velocity_min_reference_cells)
      ),
      outputs = list()
    ))
  }

  umap_df <- velocity_umap_export_10b(ref_obj)
  meta_df <- velocity_metadata_export_10b(ref_obj, pair_row, split_var)
  common_cells <- intersect(umap_df$cell_id, meta_df$cell_id)
  if (length(common_cells) == 0) {
    stop("velocity metadata and UMAP exports have no overlapping cell_id values", call. = FALSE)
  }
  umap_df <- umap_df[match(common_cells, umap_df$cell_id), , drop = FALSE]
  meta_df <- meta_df[match(common_cells, meta_df$cell_id), , drop = FALSE]

  umap_csv <- velocity_reference_umap_path_10(cfg, pair_id, split_value)
  metadata_csv <- velocity_reference_metadata_path_10(cfg, pair_id, split_value)
  write_velocity_csv_10b(umap_df, umap_csv)
  write_velocity_csv_10b(meta_df, metadata_csv)

  index <- data.frame(
    pair_id = pair_id,
    source_question_id = normalize_scalar_value(pair_row$source_question_id[[1]]),
    layer_scope = normalize_scalar_value(pair_row$layer_scope[[1]]),
    split_mode = normalize_scalar_value(pair_row$split_mode[[1]], "auto"),
    split_var = split_var,
    split_value = split_label,
    input_rds = layer_row$object_rds[[1]],
    umap_csv = umap_csv,
    metadata_csv = metadata_csv,
    n_cells = length(common_cells),
    coarse_label_var = normalize_scalar_value(pair_row$coarse_label_var[[1]]),
    fine_label_var = normalize_scalar_value(pair_row$fine_label_var[[1]]),
    status = "ok",
    reason = "",
    produced_at = timestamp_now(),
    stringsAsFactors = FALSE
  )

  outputs <- list()
  outputs[[velocity_reference_output_key_10("velocity_umap", pair_id, split_value)]] <- build_output_entry(
    umap_csv,
    "csv",
    module_name,
    "UMAP coordinates for one velocity reference unit",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(umap_df)
  )
  outputs[[velocity_reference_output_key_10("velocity_metadata", pair_id, split_value)]] <- build_output_entry(
    metadata_csv,
    "csv",
    module_name,
    "metadata for one velocity reference unit",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(meta_df)
  )

  list(index = index, outputs = outputs)
}

pairs <- velocity_read_pairs_10(cfg)
if (nrow(pairs) == 0) {
  message("10b: no active method=velocity rows in trajectory_pairs.tsv")
}

index_rows <- list()
dynamic_outputs <- list()

for (i in seq_len(nrow(pairs))) {
  pair_row <- pairs[i, , drop = FALSE]
  layer_row <- NULL
  source_obj <- NULL
  split_plan <- NULL

  setup <- tryCatch({
    layer_row <- velocity_resolve_layer_10(cfg, pair_row$layer_scope[[1]])
    source_obj <- load_comm_layer_object_07(layer_row)
    split_plan <- velocity_split_values_for_pair_10(source_obj, pair_row)
    list(ok = TRUE, error = NULL)
  }, error = function(e) {
    list(ok = FALSE, error = e)
  })

  if (!setup$ok) {
    index_rows[[length(index_rows) + 1L]] <- failed_reference_row_10b(
      pair_row,
      layer_row,
      "",
      "",
      "failed",
      conditionMessage(setup$error)
    )
    next
  }

  for (split_value in split_plan$values) {
    result <- tryCatch(
      prepare_one_velocity_reference_10b(pair_row, layer_row, source_obj, split_plan$split_var, split_value),
      error = function(e) list(
        index = failed_reference_row_10b(pair_row, layer_row, split_plan$split_var, split_value, "failed", conditionMessage(e)),
        outputs = list()
      )
    )
    index_rows[[length(index_rows) + 1L]] <- result$index
    if (length(result$outputs) > 0) {
      dynamic_outputs <- c(dynamic_outputs, result$outputs)
    }
  }
}

reference_index <- if (length(index_rows) > 0) {
  dplyr::bind_rows(index_rows)
} else {
  velocity_empty_df_10(velocity_reference_index_cols_10)
}
write_tsv_local(reference_index, cfg$velocity_reference_index_tsv)

legacy_outputs <- list()
ok_rows <- reference_index[reference_index$status == "ok", , drop = FALSE]
if (nrow(ok_rows) > 0) {
  pooled <- ok_rows[ok_rows$split_value == "pooled", , drop = FALSE]
  legacy_source <- if (nrow(pooled) > 0) pooled[1, , drop = FALSE] else ok_rows[1, , drop = FALSE]
  legacy_umap <- velocity_reference_legacy_umap_path_10(cfg)
  legacy_meta <- velocity_reference_legacy_metadata_path_10(cfg)
  ensure_dir(dirname(legacy_umap))
  file.copy(legacy_source$umap_csv[[1]], legacy_umap, overwrite = TRUE)
  file.copy(legacy_source$metadata_csv[[1]], legacy_meta, overwrite = TRUE)
  legacy_outputs <- list(
    velocity_umap = build_output_entry(
      legacy_umap,
      "csv",
      module_name,
      "compatibility UMAP coordinates for the default velocity reference unit",
      base_dir = cfg$project_root
    ),
    velocity_metadata = build_output_entry(
      legacy_meta,
      "csv",
      module_name,
      "compatibility metadata for the default velocity reference unit",
      base_dir = cfg$project_root
    )
  )
}

fixed_outputs <- list(
  velocity_reference_index = build_output_entry(
    cfg$velocity_reference_index_tsv,
    "tsv",
    module_name,
    "one row per velocity reference object or split",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(reference_index)
  )
)

velocity_manifest_10(
  cfg,
  cfg$module_10b_manifest_path,
  module_name,
  outputs = c(fixed_outputs, dynamic_outputs, legacy_outputs),
  inputs = list(
    trajectory_pairs = cfg$trajectory_pairs_sheet,
    layer_status = cfg$layer_status_file,
    module_03d = cfg$module_03d_manifest_path,
    module_04b = cfg$module_04b_manifest_path
  ),
  depends_on = list(
    module_03d = cfg$module_03d_manifest_path,
    module_04b = cfg$module_04b_manifest_path
  )
)

message("10b completed. velocity reference index: ", cfg$velocity_reference_index_tsv)
