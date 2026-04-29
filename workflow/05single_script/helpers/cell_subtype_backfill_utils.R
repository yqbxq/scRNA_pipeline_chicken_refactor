empty_cell_subtype_backfill_summary_04 <- function() {
  data.frame(
    layer_id = character(0),
    source_rds = character(0),
    subtype_source_col = character(0),
    candidate_cells = integer(0),
    matched_cells = integer(0),
    updated_cells = integer(0),
    missing_subtype_cells = integer(0),
    subtype_values = character(0),
    status = character(0),
    notes = character(0),
    stringsAsFactors = FALSE
  )
}

metadata_table_04 <- function(seu, object_label) {
  meta <- tryCatch(seu@meta.data, error = function(e) NULL)
  if (!is.data.frame(meta)) {
    stop(sprintf("%s does not expose @meta.data", object_label), call. = FALSE)
  }
  meta
}

metadata_cell_ids_04 <- function(meta) {
  ids <- rownames(meta)
  if (is.null(ids)) {
    return(rep("", nrow(meta)))
  }
  as.character(ids)
}

choose_metadata_column_04 <- function(meta, candidates) {
  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    if (candidate %in% colnames(meta)) {
      values <- as.character(meta[[candidate]])
      if (any(nzchar(values) & !is.na(values))) {
        return(candidate)
      }
    }
  }
  ""
}

resolve_panorama_backfill_rds_04 <- function(cfg) {
  manifest <- tryCatch(read_manifest_local(cfg$module_03d_manifest_path), error = function(e) NULL)
  candidates <- character(0)
  if (!is.null(manifest)) {
    manifest_path <- tryCatch(resolve_output_local(manifest, "annotated_object"), error = function(e) "")
    candidates <- c(candidates, manifest_path)
  }
  candidates <- c(candidates, cfg$panorama_annotated_rds, cfg$compat_annotated_rds)
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0) {
    stop(
      sprintf("cannot locate panorama annotated RDS for cell_subtype backfill; tried: %s", paste(candidates, collapse = ", ")),
      call. = FALSE
    )
  }
  hits[[1]]
}

update_panorama_cell_subtype_from_layer_04 <- function(panorama, layer_id, layer_rds) {
  layer <- readRDS(layer_rds)
  layer_meta <- metadata_table_04(layer, sprintf("layer %s", layer_id))
  panorama_meta <- metadata_table_04(panorama, "panorama")

  source_col <- choose_metadata_column_04(
    layer_meta,
    c("cell_subtype", paste0(layer_id, "_cell_type"), "cell_type", "annotation_label")
  )
  layer_cell_ids <- metadata_cell_ids_04(layer_meta)
  panorama_cell_ids <- metadata_cell_ids_04(panorama_meta)
  matched <- layer_cell_ids %in% panorama_cell_ids

  if (!nzchar(source_col)) {
    return(list(
      panorama = panorama,
      row = data.frame(
        layer_id = layer_id,
        source_rds = layer_rds,
        subtype_source_col = "",
        candidate_cells = nrow(layer_meta),
        matched_cells = sum(matched),
        updated_cells = 0L,
        missing_subtype_cells = nrow(layer_meta),
        subtype_values = "",
        status = "skipped",
        notes = "no subtype source column found",
        stringsAsFactors = FALSE
      )
    ))
  }

  subtype_values <- as.character(layer_meta[[source_col]])
  valid_subtype <- matched & nzchar(subtype_values) & !is.na(subtype_values)
  if (any(valid_subtype)) {
    target_index <- match(layer_cell_ids[valid_subtype], panorama_cell_ids)
    panorama@meta.data$cell_subtype[target_index] <- subtype_values[valid_subtype]
  }

  updated_values <- sort(unique(subtype_values[valid_subtype]))
  list(
    panorama = panorama,
    row = data.frame(
      layer_id = layer_id,
      source_rds = layer_rds,
      subtype_source_col = source_col,
      candidate_cells = nrow(layer_meta),
      matched_cells = sum(matched),
      updated_cells = sum(valid_subtype),
      missing_subtype_cells = sum(matched & (is.na(subtype_values) | !nzchar(subtype_values))),
      subtype_values = paste(updated_values, collapse = ","),
      status = ifelse(sum(valid_subtype) > 0, "updated", "skipped"),
      notes = ifelse(sum(matched) == 0, "no overlapping cell ids", ""),
      stringsAsFactors = FALSE
    )
  )
}

backfill_panorama_cell_subtype_04 <- function(cfg, annotation_summary_df) {
  panorama_rds <- resolve_panorama_backfill_rds_04(cfg)
  panorama <- readRDS(panorama_rds)
  panorama_meta <- metadata_table_04(panorama, "panorama")

  default_col <- choose_metadata_column_04(panorama_meta, c("cell_type", "annotation_label", "panorama_cell_type", "cell_subtype"))
  if (!nzchar(default_col)) {
    stop("panorama metadata lacks cell_type/annotation_label/panorama_cell_type/cell_subtype for default cell_subtype backfill", call. = FALSE)
  }

  panorama@meta.data$cell_subtype <- as.character(panorama_meta[[default_col]])
  rows <- list(data.frame(
    layer_id = cfg$panorama_layer_id,
    source_rds = panorama_rds,
    subtype_source_col = default_col,
    candidate_cells = nrow(panorama_meta),
    matched_cells = nrow(panorama_meta),
    updated_cells = sum(nzchar(panorama@meta.data$cell_subtype) & !is.na(panorama@meta.data$cell_subtype)),
    missing_subtype_cells = sum(is.na(panorama@meta.data$cell_subtype) | !nzchar(panorama@meta.data$cell_subtype)),
    subtype_values = paste(sort(unique(panorama@meta.data$cell_subtype[nzchar(panorama@meta.data$cell_subtype) & !is.na(panorama@meta.data$cell_subtype)])), collapse = ","),
    status = "defaulted",
    notes = "initialized from panorama metadata",
    stringsAsFactors = FALSE
  ))

  if (nrow(annotation_summary_df) > 0 && "annotated_rds" %in% colnames(annotation_summary_df)) {
    for (i in seq_len(nrow(annotation_summary_df))) {
      layer_id <- as.character(annotation_summary_df$layer_id[[i]])
      layer_rds <- as.character(annotation_summary_df$annotated_rds[[i]])
      if (!nzchar(layer_id) || !nzchar(layer_rds) || !file.exists(layer_rds)) {
        next
      }
      payload <- update_panorama_cell_subtype_from_layer_04(panorama, layer_id, layer_rds)
      panorama <- payload$panorama
      rows[[length(rows) + 1L]] <- payload$row
    }
  }

  ensure_dir(dirname(panorama_rds))
  saveRDS(panorama, panorama_rds)

  compat_rds <- cfg$compat_annotated_rds
  if (nzchar(compat_rds) && !identical(
    normalizePath(panorama_rds, winslash = "/", mustWork = FALSE),
    normalizePath(compat_rds, winslash = "/", mustWork = FALSE)
  )) {
    ensure_dir(dirname(compat_rds))
    invisible(file.copy(panorama_rds, compat_rds, overwrite = TRUE))
  }

  summary_df <- if (length(rows) > 0) do.call(rbind, rows) else empty_cell_subtype_backfill_summary_04()
  summary_tsv <- file.path(cfg$subcluster_table_dir, "panorama_cell_subtype_backfill.tsv")
  write_tsv_local(summary_df, summary_tsv)

  final_subtype_values <- sort(unique(as.character(panorama@meta.data$cell_subtype)))
  final_subtype_values <- final_subtype_values[nzchar(final_subtype_values) & !is.na(final_subtype_values)]
  status_df <- read_tsv_optional(cfg$layer_status_file)
  status_row <- if (nrow(status_df) > 0 && any(status_df$layer_id == cfg$panorama_layer_id)) {
    status_df[status_df$layer_id == cfg$panorama_layer_id, , drop = FALSE][1, , drop = FALSE]
  } else {
    data.frame(layer_id = cfg$panorama_layer_id, stringsAsFactors = FALSE)
  }
  status_row$cell_subtype_backfilled <- "yes"
  status_row$cell_subtype_column <- "cell_subtype"
  status_row$cell_subtype_backfill_tsv <- normalizePath(summary_tsv, winslash = "/", mustWork = FALSE)
  status_row$cell_subtype_values <- paste(final_subtype_values, collapse = ",")
  invisible(upsert_layer_status(cfg$layer_status_file, status_row))

  list(
    panorama_rds = panorama_rds,
    compat_rds = compat_rds,
    summary_tsv = summary_tsv,
    summary = summary_df
  )
}
