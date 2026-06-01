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
module_name <- "spatial_08c_lr_colocalization"
prepare_dirs_spatial(cfg)

candidates_tsv <- file.path(cfg$spatial_communication_table_dir, "spatial_comm_candidates.tsv")
candidates <- st07_read_tsv(candidates_tsv)
props <- st08_collect_deconv_props(cfg)

score_candidate_08c <- function(candidate, props) {
  if (nrow(props) == 0) {
    return(data.frame(
      comm_candidate_id = candidate$comm_candidate_id,
      lr_axis_id = candidate$lr_axis_id,
      section_id = candidate$section_id,
      sender_cell_type = candidate$sender_cell_type,
      receiver_cell_type = candidate$receiver_cell_type,
      sender_spatial_abundance = NA_real_,
      receiver_spatial_abundance = NA_real_,
      sender_receiver_colocalization = NA_real_,
      colocalization_support = "no",
      status = "skipped_no_deconv_proportions",
      reason = "No ok deconvolution proportion outputs were available.",
      stringsAsFactors = FALSE
    ))
  }
  section <- st08_scalar(candidate$section_id, "all")
  sender <- st08_scalar(candidate$sender_cell_type)
  receiver <- st08_scalar(candidate$receiver_cell_type)
  hit <- props
  if (!section %in% c("all", "*", "__ALL__")) {
    hit <- hit[hit$section %in% c(section, "all"), , drop = FALSE]
  }
  if (nrow(hit) > 0) {
    hit <- aggregate(proportion ~ spot_id + cell_type, data = hit, FUN = mean, na.rm = TRUE)
  }
  wide <- reshape(hit[, c("spot_id", "cell_type", "proportion"), drop = FALSE], idvar = "spot_id", timevar = "cell_type", direction = "wide")
  sender_col <- paste0("proportion.", sender)
  receiver_col <- paste0("proportion.", receiver)
  if (!sender_col %in% colnames(wide) || !receiver_col %in% colnames(wide)) {
    return(data.frame(
      comm_candidate_id = candidate$comm_candidate_id,
      lr_axis_id = candidate$lr_axis_id,
      section_id = section,
      sender_cell_type = sender,
      receiver_cell_type = receiver,
      sender_spatial_abundance = if (sender_col %in% colnames(wide)) mean(wide[[sender_col]], na.rm = TRUE) else NA_real_,
      receiver_spatial_abundance = if (receiver_col %in% colnames(wide)) mean(wide[[receiver_col]], na.rm = TRUE) else NA_real_,
      sender_receiver_colocalization = NA_real_,
      colocalization_support = "no",
      status = "blocked_celltype_not_localized",
      reason = "Sender or receiver cell type was not present in ok deconvolution proportions.",
      stringsAsFactors = FALSE
    ))
  }
  sender_ab <- suppressWarnings(as.numeric(wide[[sender_col]]))
  receiver_ab <- suppressWarnings(as.numeric(wide[[receiver_col]]))
  corr <- if (length(sender_ab) > 1 && stats::sd(sender_ab, na.rm = TRUE) > 0 && stats::sd(receiver_ab, na.rm = TRUE) > 0) {
    suppressWarnings(as.numeric(stats::cor(sender_ab, receiver_ab, use = "pairwise.complete.obs")))
  } else {
    NA_real_
  }
  support <- !is.na(corr) && corr > 0 && mean(sender_ab, na.rm = TRUE) > 0 && mean(receiver_ab, na.rm = TRUE) > 0
  data.frame(
    comm_candidate_id = candidate$comm_candidate_id,
    lr_axis_id = candidate$lr_axis_id,
    section_id = section,
    sender_cell_type = sender,
    receiver_cell_type = receiver,
    sender_spatial_abundance = mean(sender_ab, na.rm = TRUE),
    receiver_spatial_abundance = mean(receiver_ab, na.rm = TRUE),
    sender_receiver_colocalization = corr,
    colocalization_support = if (support) "yes" else "no",
    status = "ok",
    reason = "",
    stringsAsFactors = FALSE
  )
}

coloc <- if (nrow(candidates) == 0) {
  st08_empty_colocalization()
} else {
  do.call(rbind, lapply(seq_len(nrow(candidates)), function(i) score_candidate_08c(candidates[i, , drop = FALSE], props)))
}

lr_coloc_tsv <- file.path(cfg$spatial_communication_table_dir, "lr_colocalization.tsv")
sender_receiver_tsv <- file.path(cfg$spatial_communication_table_dir, "sender_receiver_colocalization.tsv")
spatial_write_tsv(coloc, lr_coloc_tsv)
spatial_write_tsv(coloc[, intersect(colnames(st08_empty_colocalization()), colnames(coloc)), drop = FALSE], sender_receiver_tsv)

st07_write_manifest_local(
  manifest_path = cfg$module_08c_lr_colocalization_manifest_path,
  new_outputs = list(
    lr_colocalization = build_output_entry(lr_coloc_tsv, "tsv", module_name, "LR candidate sender/receiver colocalization proxy from deconvolution proportions", base_dir = cfg$project_root, schema = infer_schema_from_df(coloc)),
    sender_receiver_colocalization = build_output_entry(sender_receiver_tsv, "tsv", module_name, "sender receiver colocalization support table", base_dir = cfg$project_root, schema = infer_schema_from_df(coloc))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(spatial_comm_candidates = candidates_tsv),
  version = cfg$module_07_version,
  depends_on = list(spatial_08a_spatial_communication_io = cfg$module_08a_spatial_communication_io_manifest_path, spatial_07e_deconvolution_compare = cfg$module_07e_deconvolution_compare_manifest_path)
)

message(sprintf("08c LR colocalization completed. rows=%d", nrow(coloc)))
