#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source_utf8 <- function(path) source(path, encoding = "UTF-8")

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_05.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_06.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_07.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_mapping_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_consensus_utils.R"))

load_required_packages(c("jsonlite"))

cfg <- get_single_script_config_07()
module_name <- "07d_communication_consensus"
prepare_dirs_07(cfg)

resolve_consensus_path_07d <- function(path) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path)) {
    return("")
  }
  if (grepl("^/", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(cfg$project_root, path), winslash = "/", mustWork = FALSE)
}

read_tsv_if_exists_07d <- function(path) {
  path <- resolve_consensus_path_07d(path)
  if (!nzchar(path) || !file.exists(path)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  read_tsv_optional(path)
}

ensure_lr_axis_07d <- function(df) {
  if (nrow(df) == 0) {
    return(df)
  }
  if (!"ligand_human" %in% colnames(df) && "ligand" %in% colnames(df)) {
    df$ligand_human <- df$ligand
  }
  if (!"receptor_human" %in% colnames(df) && "receptor" %in% colnames(df)) {
    df$receptor_human <- df$receptor
  }
  if (!"lr_axis_id" %in% colnames(df)) {
    df$lr_axis_id <- ""
  }
  missing_axis <- !nzchar(as.character(df$lr_axis_id))
  required <- c("ligand_human", "receptor_human", "source", "target")
  if (any(missing_axis) && all(required %in% colnames(df))) {
    df$lr_axis_id[missing_axis] <- standardize_lr_axis_id(
      df$ligand_human[missing_axis],
      df$receptor_human[missing_axis],
      df$source[missing_axis],
      df$target[missing_axis]
    )
  }
  df
}

load_cellchat_lr_07d <- function(index_path) {
  index <- read_tsv_if_exists_07d(index_path)
  if (nrow(index) == 0 || !"lr_table_path" %in% colnames(index)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  rows <- list()
  for (i in seq_len(nrow(index))) {
    lr <- read_tsv_if_exists_07d(index$lr_table_path[[i]])
    if (nrow(lr) == 0) {
      next
    }
    for (col in c("pair_id", "layer_id", "condition_value")) {
      if (!col %in% colnames(lr) || !all(nzchar(as.character(lr[[col]])))) {
        lr[[col]] <- index[[col]][[i]] %||% ""
      }
    }
    lr$method <- "cellchat"
    lr$method_evidence_class <- "hypothesis_only"
    lr$can_be_primary <- "no"
    rows[[length(rows) + 1L]] <- ensure_lr_axis_07d(lr)
  }
  if (length(rows) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

load_liana_lr_07d <- function(path) {
  df <- read_tsv_if_exists_07d(path)
  if (nrow(df) == 0) {
    return(df)
  }
  if (!"liana_consensus_hit" %in% colnames(df)) {
    df$liana_consensus_hit <- FALSE
  }
  ensure_lr_axis_07d(df)
}

label_matches_07d <- function(query, labels) {
  query <- normalize_scalar_value(query)
  labels <- unlist(strsplit(normalize_scalar_value(labels), "[,;]", perl = TRUE), use.names = FALSE)
  labels <- trimws(labels[nzchar(trimws(labels))])
  if (length(labels) == 0 || "*" %in% labels || !nzchar(query)) {
    return(TRUE)
  }
  query %in% labels
}

condition_matches_07d <- function(query, candidate) {
  query <- normalize_scalar_value(query, "all")
  candidate <- normalize_scalar_value(candidate, "all")
  candidate %in% c("", "all", query)
}

load_nichenet_tasks_07d <- function(path) {
  df <- read_tsv_if_exists_07d(path)
  if (nrow(df) == 0) {
    return(df)
  }
  df$nichenet_hit <- as_logical_flag_07d(df$success %||% FALSE) & as.character(df$status %||% "") == "ok"
  if (!"method_evidence_class" %in% colnames(df)) {
    df$method_evidence_class <- "downstream_validation"
  }
  df$source <- df$sender_cell_types %||% df$sender_set %||% ""
  df$target <- df$receiver_cell_types %||% df$receiver_set %||% ""
  df
}

apply_nichenet_task_support_07d <- function(consensus_df, nichenet_tasks) {
  if (nrow(consensus_df) == 0 || nrow(nichenet_tasks) == 0) {
    return(consensus_df)
  }
  hit_tasks <- nichenet_tasks[as_logical_flag_07d(nichenet_tasks$nichenet_hit), , drop = FALSE]
  if (nrow(hit_tasks) == 0) {
    return(consensus_df)
  }
  for (i in seq_len(nrow(consensus_df))) {
    if (isTRUE(consensus_df$nichenet_hit[[i]])) {
      next
    }
    matches <- vapply(seq_len(nrow(hit_tasks)), function(j) {
      condition_matches_07d(consensus_df$condition_value[[i]], hit_tasks$condition_value[[j]]) &&
        label_matches_07d(consensus_df$source[[i]], hit_tasks$source[[j]]) &&
        label_matches_07d(consensus_df$target[[i]], hit_tasks$target[[j]])
    }, logical(1))
    if (any(matches)) {
      first <- which(matches)[[1]]
      consensus_df$nichenet_hit[[i]] <- TRUE
      consensus_df$nichenet_method_evidence_class[[i]] <- hit_tasks$method_evidence_class[[first]]
      consensus_df$pair_id[[i]] <- first_nonempty_07d(consensus_df$pair_id[[i]], hit_tasks$pair_id[[first]])
      consensus_df$layer_id[[i]] <- first_nonempty_07d(consensus_df$layer_id[[i]], hit_tasks$layer_id[[first]])
    }
  }
  consensus_df
}

cellchat_lr <- load_cellchat_lr_07d(cfg$cellchat_index_tsv)
liana_lr <- load_liana_lr_07d(cfg$liana_consensus_lr_tsv)
nichenet_tasks <- load_nichenet_tasks_07d(cfg$nichenet_index_tsv)
nichenet_axis <- if ("lr_axis_id" %in% colnames(nichenet_tasks)) ensure_lr_axis_07d(nichenet_tasks) else data.frame(stringsAsFactors = FALSE)

consensus_df <- full_outer_join_lr_tables(cellchat_lr, liana_lr, nichenet_axis)
consensus_df <- apply_nichenet_task_support_07d(consensus_df, nichenet_tasks)
consensus_df <- assign_evidence_tier_07d(consensus_df)

consensus_keep <- unique(c(
  colnames(empty_method_consensus_07d()),
  setdiff(colnames(consensus_df), colnames(empty_method_consensus_07d()))
))
for (col in consensus_keep) {
  if (!col %in% colnames(consensus_df)) {
    consensus_df[[col]] <- ""
  }
}
consensus_df <- consensus_df[, consensus_keep, drop = FALSE]

evidence_cols <- c(
  "lr_axis_id", "condition_value", "pair_id", "layer_id", "source", "target",
  "ligand", "receptor", "cellchat_hit", "liana_consensus_hit", "nichenet_hit",
  "commot_spatial_hit", "evidence_method_n", "evidence_tier", "can_be_primary",
  "evidence_reason"
)
evidence_df <- consensus_df[, intersect(evidence_cols, colnames(consensus_df)), drop = FALSE]

write_tsv_local(consensus_df, cfg$communication_method_consensus_tsv)
write_tsv_local(evidence_df, cfg$communication_evidence_tier_tsv)

if (file.exists(cfg$module_07d_manifest_path)) {
  unlink(cfg$module_07d_manifest_path)
}
fixed_outputs <- list(
  method_consensus_tsv = build_output_entry(
    cfg$communication_method_consensus_tsv,
    "tsv",
    module_name,
    "one row per ligand-receptor axis after CellChat/LIANA/NicheNet evidence merge",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(consensus_df)
  ),
  communication_evidence_tier_tsv = build_output_entry(
    cfg$communication_evidence_tier_tsv,
    "tsv",
    module_name,
    "one row per ligand-receptor axis with four-level evidence tier",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(evidence_df)
  )
)
write_manifest_local(
  manifest_path = cfg$module_07d_manifest_path,
  new_outputs = fixed_outputs,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    cellchat_index_tsv = cfg$cellchat_index_tsv,
    liana_consensus_lr_tsv = cfg$liana_consensus_lr_tsv,
    nichenet_index_tsv = cfg$nichenet_index_tsv,
    consensus_min_methods_for_primary = Sys.getenv("CONSENSUS_MIN_METHODS_FOR_PRIMARY", "2"),
    consensus_require_nichenet_for_primary = Sys.getenv("CONSENSUS_REQUIRE_NICHENET_FOR_PRIMARY", "yes")
  ),
  version = cfg$module_07d_consensus_version,
  depends_on = list(
    module_07a = cfg$module_07a_manifest_path,
    module_07b = cfg$module_07b_manifest_path,
    module_07c = cfg$module_07c_manifest_path
  )
)

message("07d completed. consensus: ", cfg$communication_method_consensus_tsv)
