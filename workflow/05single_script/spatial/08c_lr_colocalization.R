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

load_spatial_expression_08c <- function(cfg) {
  candidates <- c(cfg$spatial_panorama_subannotated_rds, cfg$spatial_panorama_annotated_rds, cfg$spatial_panorama_clustered_rds)
  object_path <- candidates[nzchar(candidates) & file.exists(candidates)][1]
  if (is.na(object_path)) {
    return(list(status = "skipped_missing_expression", reason = "No spatial Seurat/RDS expression object was available.", counts = NULL))
  }
  obj <- tryCatch(readRDS(object_path), error = function(e) e)
  if (inherits(obj, "error")) {
    return(list(status = "skipped_missing_expression", reason = conditionMessage(obj), counts = NULL))
  }
  counts <- NULL
  if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
    assay <- tryCatch(Seurat::DefaultAssay(obj), error = function(e) "Spatial")
    counts <- tryCatch(Seurat::GetAssayData(obj, assay = assay, layer = "data"), error = function(e) {
      tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "data"), error = function(e2) NULL)
    })
  } else if (is.list(obj) && !is.null(obj$counts)) {
    counts <- obj$counts
  } else if (!is.null(dim(obj)) && length(dim(obj)) == 2) {
    counts <- obj
  }
  if (is.null(counts) || is.null(rownames(counts)) || is.null(colnames(counts))) {
    return(list(status = "skipped_missing_expression", reason = "Spatial expression matrix could not be extracted with gene and spot names.", counts = NULL))
  }
  list(status = "ok", reason = "", counts = counts)
}

lr_expression_support_08c <- function(candidate, wide, expr) {
  if (!identical(expr$status, "ok") || is.null(expr$counts)) {
    return(list(
      ligand_spot_expression_mean = NA_real_,
      receptor_spot_expression_mean = NA_real_,
      lr_expression_colocalization = NA_real_,
      lr_expression_support = "no",
      expression_status = expr$status,
      expression_reason = expr$reason
    ))
  }
  ligand <- st08_scalar(candidate$ligand)
  receptor <- st08_scalar(candidate$receptor)
  sender_col <- paste0("proportion.", st08_scalar(candidate$sender_cell_type))
  receiver_col <- paste0("proportion.", st08_scalar(candidate$receiver_cell_type))
  if (!ligand %in% rownames(expr$counts) || !receptor %in% rownames(expr$counts) || !sender_col %in% colnames(wide) || !receiver_col %in% colnames(wide)) {
    return(list(
      ligand_spot_expression_mean = NA_real_,
      receptor_spot_expression_mean = NA_real_,
      lr_expression_colocalization = NA_real_,
      lr_expression_support = "no",
      expression_status = "skipped_missing_expression",
      expression_reason = "Ligand/receptor genes or sender/receiver abundance columns were missing."
    ))
  }
  common_spots <- intersect(as.character(wide$spot_id), colnames(expr$counts))
  if (length(common_spots) < 2) {
    return(list(
      ligand_spot_expression_mean = NA_real_,
      receptor_spot_expression_mean = NA_real_,
      lr_expression_colocalization = NA_real_,
      lr_expression_support = "no",
      expression_status = "skipped_missing_expression",
      expression_reason = "Fewer than two shared spots between expression matrix and deconvolution proportions."
    ))
  }
  wide <- wide[match(common_spots, wide$spot_id), , drop = FALSE]
  ligand_expr <- as.numeric(expr$counts[ligand, common_spots])
  receptor_expr <- as.numeric(expr$counts[receptor, common_spots])
  ligand_score <- ligand_expr * suppressWarnings(as.numeric(wide[[sender_col]]))
  receptor_score <- receptor_expr * suppressWarnings(as.numeric(wide[[receiver_col]]))
  corr <- if (stats::sd(ligand_score, na.rm = TRUE) > 0 && stats::sd(receptor_score, na.rm = TRUE) > 0) {
    suppressWarnings(as.numeric(stats::cor(ligand_score, receptor_score, use = "pairwise.complete.obs")))
  } else {
    NA_real_
  }
  list(
    ligand_spot_expression_mean = mean(ligand_expr, na.rm = TRUE),
    receptor_spot_expression_mean = mean(receptor_expr, na.rm = TRUE),
    lr_expression_colocalization = corr,
    lr_expression_support = if (!is.na(corr) && corr > 0) "yes" else "no",
    expression_status = "ok",
    expression_reason = ""
  )
}

expr <- load_spatial_expression_08c(cfg)

score_candidate_08c <- function(candidate, props, expr) {
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
      ligand_spot_expression_mean = NA_real_,
      receptor_spot_expression_mean = NA_real_,
      lr_expression_colocalization = NA_real_,
      lr_expression_support = "no",
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
      ligand_spot_expression_mean = NA_real_,
      receptor_spot_expression_mean = NA_real_,
      lr_expression_colocalization = NA_real_,
      lr_expression_support = "no",
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
  lr_expr <- lr_expression_support_08c(candidate, wide, expr)
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
    ligand_spot_expression_mean = lr_expr$ligand_spot_expression_mean,
    receptor_spot_expression_mean = lr_expr$receptor_spot_expression_mean,
    lr_expression_colocalization = lr_expr$lr_expression_colocalization,
    lr_expression_support = lr_expr$lr_expression_support,
    status = "ok",
    reason = if (identical(lr_expr$expression_status, "ok")) "" else lr_expr$expression_reason,
    stringsAsFactors = FALSE
  )
}

coloc <- if (nrow(candidates) == 0) {
  st08_empty_colocalization()
} else {
  do.call(rbind, lapply(seq_len(nrow(candidates)), function(i) score_candidate_08c(candidates[i, , drop = FALSE], props, expr)))
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
