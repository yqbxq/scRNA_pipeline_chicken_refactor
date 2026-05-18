#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")

load_required_packages(c("Seurat", "Matrix", "ggplot2", "jsonlite"))

cfg <- get_spatial_script_config()
module_name <- "spatial_04b_subcluster_annotate"
prepare_dirs_spatial(cfg)
set.seed(cfg$random_seed)

parse_cli <- function(args) {
  out <- list(marker_panel_dir = file.path(cfg$config_dir, "marker_panels"))
  positional <- character()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--marker-panel-dir") && i < length(args)) {
      out$marker_panel_dir <- args[[i + 1L]]
      i <- i + 2L
    } else {
      positional <- c(positional, arg)
      i <- i + 1L
    }
  }
  if (length(positional) > 0 && dir.exists(positional[[1]])) {
    out$marker_panel_dir <- positional[[1]]
  }
  out
}

assignment_schema <- function() {
  data.frame(
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
}

resolve_index_rds <- function(row) {
  out_rds <- trimws(as.character(row$out_rds %||% ""))
  candidates <- character()
  if (nzchar(out_rds)) {
    candidates <- c(candidates, spatial_resolve_path(out_rds, cfg$project_root))
  }
  layer_id <- spatial_cell(row, "layer_id", "")
  if (nzchar(layer_id)) {
    candidates <- c(candidates, file.path(cfg$spatial_region_dir, sprintf("%s.rds", spatial_safe_id(layer_id))))
  }
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) "" else hit[[1]]
}

parent_region_from_index <- function(row) {
  values <- spatial_tokenize(spatial_cell(row, "selection_values", ""))
  if (length(values) > 0) values[[1]] else spatial_cell(row, "layer_id", "")
}

write_outputs_manifest <- function(assignment_df, evidence_summary, inputs) {
  assignment_path <- file.path(cfg$spatial_subcluster_annotate_table_dir, "region_subset_assignment.tsv")
  evidence_path <- file.path(cfg$spatial_subcluster_annotate_table_dir, "evidence_summary.tsv")
  spatial_write_tsv(assignment_df, assignment_path)
  spatial_write_tsv(evidence_summary, evidence_path)
  outputs <- list(
    region_subset_assignment = build_output_entry(assignment_path, "tsv", module_name, "subcluster to sub-region assignment", base_dir = cfg$project_root, schema = infer_schema_from_df(assignment_df)),
    evidence_summary = build_output_entry(evidence_path, "tsv", module_name, "per-layer subregion annotation status and evidence summary", base_dir = cfg$project_root, schema = infer_schema_from_df(evidence_summary)),
    annotated_subsets_dir = build_output_entry(cfg$spatial_region_annotated_dir, "directory", module_name, "per-layer subregion annotated subset RDS files", base_dir = cfg$project_root)
  )
  if (file.exists(cfg$spatial_panorama_subannotated_rds)) {
    outputs$panorama_subannotated <- build_output_entry(cfg$spatial_panorama_subannotated_rds, "rds", module_name, "panorama with projected sub_region metadata", base_dir = cfg$project_root)
  }
  write_manifest_local(
    manifest_path = cfg$module_04b_subcluster_annotate_manifest_path,
    new_outputs = outputs,
    module_name = module_name,
    base_dir = cfg$project_root,
    inputs = inputs,
    version = cfg$module_04_version,
    depends_on = list(spatial_04a_subcluster_build = cfg$module_04a_subcluster_build_manifest_path)
  )
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (file.exists(cfg$spatial_panorama_subannotated_rds)) {
  unlink(cfg$spatial_panorama_subannotated_rds)
}
index_path <- file.path(cfg$spatial_subcluster_build_table_dir, "region_subset_index.tsv")
index_df <- spatial_read_tsv(index_path)
ok_index <- if (nrow(index_df) > 0 && "status" %in% colnames(index_df)) index_df[index_df$status == "ok", , drop = FALSE] else index_df[0, , drop = FALSE]

if (nrow(ok_index) == 0) {
  assignment_df <- assignment_schema()
  evidence_summary <- data.frame(layer_id = "", parent_region = "", status = "no_subsets_built", n_clusters = 0L, n_assigned = 0L, reason = "04a produced no ok region subset layers", stringsAsFactors = FALSE)
  write_outputs_manifest(assignment_df, evidence_summary, inputs = list(region_subset_index = index_path, marker_panel_dir = args$marker_panel_dir))
  message("spatial subcluster annotate skipped: no built subsets")
  quit(status = 0)
}

if (!file.exists(cfg$spatial_panorama_annotated_rds)) {
  stop(sprintf("missing annotated panorama: %s", cfg$spatial_panorama_annotated_rds), call. = FALSE)
}
panorama <- readRDS(cfg$spatial_panorama_annotated_rds)

layer_results <- list()
evidence_rows <- list()
annotated_subset_objs <- list()
for (i in seq_len(nrow(ok_index))) {
  row <- ok_index[i, , drop = FALSE]
  layer_id <- spatial_cell(row, "layer_id", sprintf("layer_%d", i))
  parent_region <- parent_region_from_index(row)
  panel_path <- file.path(args$marker_panel_dir, sprintf("%s.tsv", layer_id))
  if (!file.exists(panel_path)) {
    layer_results[[layer_id]] <- list(status = "failed_panel_missing", assignment = assignment_schema())
    evidence_rows[[length(evidence_rows) + 1L]] <- data.frame(layer_id = layer_id, parent_region = parent_region, status = "failed_panel_missing", n_clusters = 0L, n_assigned = 0L, reason = sprintf("missing required subregion panel: %s", panel_path), stringsAsFactors = FALSE)
    next
  }
  subset_rds <- resolve_index_rds(row)
  if (!nzchar(subset_rds)) {
    layer_results[[layer_id]] <- list(status = "failed", assignment = assignment_schema())
    evidence_rows[[length(evidence_rows) + 1L]] <- data.frame(layer_id = layer_id, parent_region = parent_region, status = "failed", n_clusters = 0L, n_assigned = 0L, reason = "subset RDS listed by 04a is missing", stringsAsFactors = FALSE)
    next
  }
  result <- tryCatch({
    subset_obj <- readRDS(subset_rds)
    panel <- read_spatial_subregion_panel(args$marker_panel_dir, layer_id = layer_id)
    annotation <- region_three_evidence_chain(
      panorama = subset_obj,
      cluster_col = "cluster_default",
      panel = panel,
      override_file = cfg$selected_region_annotation_override_file
    )
    subset_obj <- assign_subregion_labels(annotation$panorama, annotation$evidence_df, parent_region = parent_region)
    out_rds <- file.path(cfg$spatial_region_annotated_dir, sprintf("%s.rds", spatial_safe_id(layer_id)))
    saveRDS(subset_obj, out_rds)

    layer_table_dir <- file.path(cfg$spatial_subcluster_annotate_table_dir, layer_id)
    ensure_dir(layer_table_dir)
    evidence_df <- annotation$evidence_df
    evidence_df$sub_region <- vapply(as.character(evidence_df$region), function(x) spatial_prefix_subregion(parent_region, x), character(1))
    evidence_out <- evidence_df
    evidence_out$layer_id <- layer_id
    evidence_out$parent_region <- parent_region
    evidence_out <- evidence_out[, c("layer_id", "parent_region", setdiff(colnames(evidence_out), c("layer_id", "parent_region"))), drop = FALSE]
    marker_cols <- c("cluster", "gene", "avg_log2FC", "p_val_adj", "pct.1", "pct.2")
    cluster_markers <- annotation$marker_table
    for (col in marker_cols) {
      if (!col %in% colnames(cluster_markers)) {
        cluster_markers[[col]] <- NA
      }
    }
    cluster_markers <- cluster_markers[, marker_cols, drop = FALSE]
    module_score_matrix <- annotation$module_score_matrix
    module_score_matrix$cluster <- rownames(module_score_matrix)
    module_score_matrix <- module_score_matrix[, c("cluster", setdiff(colnames(module_score_matrix), "cluster")), drop = FALSE]
    spatial_write_tsv(evidence_out, file.path(layer_table_dir, "evidence_summary.tsv"))
    spatial_write_tsv(cluster_markers, file.path(layer_table_dir, "cluster_markers.tsv"))
    spatial_write_tsv(module_score_matrix, file.path(layer_table_dir, "module_score_matrix.tsv"))

    assignment <- data.frame(
      layer_id = layer_id,
      parent_region = parent_region,
      cluster = as.character(evidence_df$cluster),
      sub_region = as.character(evidence_df$sub_region),
      n_spots = as.integer(evidence_df$n_spots),
      module_score_max = as.numeric(evidence_df$module_score_max),
      panel_hit_count_winner = as.integer(evidence_df$panel_hit_count_winner),
      evidence_agreement = as.logical(evidence_df$evidence_agreement),
      confidence = as.character(evidence_df$confidence),
      stringsAsFactors = FALSE
    )
    annotated_subset_objs[[layer_id]] <<- subset_obj
    evidence_rows[[length(evidence_rows) + 1L]] <<- data.frame(layer_id = layer_id, parent_region = parent_region, status = "ok", n_clusters = nrow(evidence_df), n_assigned = nrow(assignment), reason = "", stringsAsFactors = FALSE)
    list(status = "ok", assignment = assignment, obj = subset_obj, out_rds = out_rds)
  }, error = function(e) {
    evidence_rows[[length(evidence_rows) + 1L]] <<- data.frame(layer_id = layer_id, parent_region = parent_region, status = "failed", n_clusters = 0L, n_assigned = 0L, reason = conditionMessage(e), stringsAsFactors = FALSE)
    list(status = "failed", assignment = assignment_schema(), obj = NULL, message = conditionMessage(e))
  })
  layer_results[[layer_id]] <- result
}

assignment_df <- write_region_subset_assignment_tsv(layer_results, file.path(cfg$spatial_subcluster_annotate_table_dir, "region_subset_assignment.tsv"))
evidence_summary <- if (length(evidence_rows) > 0) do.call(rbind, evidence_rows) else data.frame(layer_id = character(), parent_region = character(), status = character(), n_clusters = integer(), n_assigned = integer(), reason = character(), stringsAsFactors = FALSE)
spatial_write_tsv(evidence_summary, file.path(cfg$spatial_subcluster_annotate_table_dir, "evidence_summary.tsv"))
if (length(annotated_subset_objs) > 0) {
  panorama <- project_subregion_to_panorama(panorama, annotated_subset_objs)
  saveRDS(panorama, cfg$spatial_panorama_subannotated_rds)
}

write_outputs_manifest(assignment_df, evidence_summary, inputs = list(region_subset_index = index_path, marker_panel_dir = args$marker_panel_dir, panorama_annotated = cfg$spatial_panorama_annotated_rds))
message(sprintf("spatial subcluster annotate complete: %d ok layer(s)", sum(evidence_summary$status == "ok", na.rm = TRUE)))
