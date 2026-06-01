#!/usr/bin/env Rscript

.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/manifest_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_deconv_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_communication_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_08d_communication_neighborhood"
prepare_dirs_spatial(cfg)

candidates <- st07_read_tsv(file.path(cfg$spatial_communication_table_dir, "spatial_comm_candidates.tsv"))
coloc <- st07_read_tsv(file.path(cfg$spatial_communication_table_dir, "lr_colocalization.tsv"))
interaction <- st07_read_tsv(file.path(cfg$spatial_neighborhood_table_dir, "interaction_matrix.tsv"))
nhood <- st07_read_tsv(file.path(cfg$spatial_neighborhood_table_dir, "nhood_enrichment_zscore.tsv"))

lookup_neighbor_08d <- function(candidate) {
  sender <- st08_scalar(candidate$sender_cell_type)
  receiver <- st08_scalar(candidate$receiver_cell_type)
  score <- NA_real_
  co_score <- NA_real_
  if (nrow(interaction) > 0) {
    cols <- colnames(interaction)
    if (all(c("label_a", "label_b", "interaction") %in% cols)) {
      hit <- interaction[interaction$label_a == sender & interaction$label_b == receiver, , drop = FALSE]
      if (nrow(hit) > 0) co_score <- suppressWarnings(as.numeric(hit$interaction[[1]]))
    } else if (sender %in% interaction[[1]] && receiver %in% colnames(interaction)) {
      hit <- interaction[interaction[[1]] == sender, receiver, drop = TRUE]
      co_score <- suppressWarnings(as.numeric(hit[[1]]))
    }
  }
  if (nrow(nhood) > 0) {
    cols <- colnames(nhood)
    if (all(c("label_a", "label_b", "zscore") %in% cols)) {
      hit <- nhood[nhood$label_a == sender & nhood$label_b == receiver, , drop = FALSE]
      if (nrow(hit) > 0) score <- suppressWarnings(as.numeric(hit$zscore[[1]]))
    } else if (sender %in% nhood[[1]] && receiver %in% colnames(nhood)) {
      hit <- nhood[nhood[[1]] == sender, receiver, drop = TRUE]
      score <- suppressWarnings(as.numeric(hit[[1]]))
    }
  }
  support <- (is.finite(score) && score > 0) || (is.finite(co_score) && co_score > 0)
  data.frame(
    comm_candidate_id = candidate$comm_candidate_id,
    lr_axis_id = candidate$lr_axis_id,
    section_id = candidate$section_id,
    sender_cell_type = sender,
    receiver_cell_type = receiver,
    neighborhood_enrichment_score = score,
    co_occurrence_score = co_score,
    neighborhood_support = if (support) "yes" else "no",
    status = if (nrow(interaction) == 0 && nrow(nhood) == 0) "skipped_no_neighborhood_outputs" else "ok",
    reason = if (nrow(interaction) == 0 && nrow(nhood) == 0) "06d neighborhood outputs were not available." else "",
    stringsAsFactors = FALSE
  )
}

consistency <- if (nrow(candidates) == 0) {
  st08_empty_neighborhood()
} else {
  do.call(rbind, lapply(seq_len(nrow(candidates)), function(i) lookup_neighbor_08d(candidates[i, , drop = FALSE])))
}

out_tsv <- file.path(cfg$spatial_communication_table_dir, "communication_neighborhood_consistency.tsv")
spatial_write_tsv(consistency, out_tsv)

st07_write_manifest_local(
  manifest_path = cfg$module_08d_neighborhood_consistency_manifest_path,
  new_outputs = list(
    communication_neighborhood_consistency = build_output_entry(out_tsv, "tsv", module_name, "sender/receiver communication candidates crossed with 06d neighborhood support", base_dir = cfg$project_root, schema = infer_schema_from_df(consistency))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(spatial_comm_candidates = file.path(cfg$spatial_communication_table_dir, "spatial_comm_candidates.tsv"), spatial_neighborhood = cfg$module_06d_neighborhood_manifest_path),
  version = cfg$module_07_version,
  depends_on = list(spatial_08c_lr_colocalization = cfg$module_08c_lr_colocalization_manifest_path, spatial_06d_neighborhood = cfg$module_06d_neighborhood_manifest_path)
)

message(sprintf("08d communication-neighborhood consistency completed. rows=%d", nrow(consistency)))
