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
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "deg_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "enrichment_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_mapping_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_pairs_utils.R"))

load_required_packages(c("dplyr", "tibble", "jsonlite", "ggplot2"))

cfg <- get_single_script_config_07()
module_name <- "07c_communication_eda"
prepare_dirs_07(cfg)

cellchat_index <- read_tsv_optional(cfg$cellchat_index_tsv)
nichenet_index <- read_tsv_optional(cfg$nichenet_index_tsv)
cellchat_triage <- read_tsv_optional(cfg$cellchat_triage_tsv)
nichenet_triage <- read_tsv_optional(cfg$nichenet_triage_tsv)

for (col in c("pair_id", "layer_id", "condition_value", "status", "lr_table_path")) {
  if (!col %in% colnames(cellchat_index)) cellchat_index[[col]] <- character(nrow(cellchat_index))
}
for (col in c("pair_id", "layer_id", "condition_value", "status", "ligand_activity_tsv", "sender_cell_types", "receiver_cell_types")) {
  if (!col %in% colnames(nichenet_index)) nichenet_index[[col]] <- character(nrow(nichenet_index))
}

split_roles_07c <- function(value) {
  parts <- split_csv_local(value)
  parts[nzchar(parts)]
}

read_lr_for_row_07c <- function(row) {
  path <- normalize_scalar_value(row$lr_table_path[[1]])
  lr <- read_tsv_optional(path)
  if (nrow(lr) == 0) {
    return(empty_df_07(c("source", "target", "ligand", "receptor", "prob", "pval", "pathway_name")))
  }
  for (col in c("source", "target", "ligand", "receptor", "prob", "pval", "pathway_name")) {
    if (!col %in% colnames(lr)) {
      lr[[col]] <- if (col %in% c("prob", "pval")) NA_real_ else ""
    }
  }
  lr$prob <- suppressWarnings(as.numeric(lr$prob))
  lr$pval <- suppressWarnings(as.numeric(lr$pval))
  lr
}

read_ligand_activity_07c <- function(row) {
  path <- normalize_scalar_value(row$ligand_activity_tsv[[1]])
  df <- read_tsv_optional(path)
  if (nrow(df) == 0) {
    return(empty_df_07(c("test_ligand", "rank", "score")))
  }
  if (!"test_ligand" %in% colnames(df)) {
    df$test_ligand <- rownames(df)
  }
  score_col <- pick_first_col_07(df, c("pearson", "aupr_corrected", "aupr", "auroc"))
  df$score <- if (nzchar(score_col)) suppressWarnings(as.numeric(df[[score_col]])) else NA_real_
  if (!"rank" %in% colnames(df)) {
    df <- df[order(-df$score), , drop = FALSE]
    df$rank <- seq_len(nrow(df))
  }
  df$rank <- suppressWarnings(as.integer(df$rank))
  df
}

filter_lr_direction_07c <- function(lr, sender_cell_types, receiver_cell_types) {
  senders <- split_roles_07c(sender_cell_types)
  receivers <- split_roles_07c(receiver_cell_types)
  if (length(senders) > 0 && "source" %in% colnames(lr)) {
    lr <- lr[as.character(lr$source) %in% senders, , drop = FALSE]
  }
  if (length(receivers) > 0 && "target" %in% colnames(lr)) {
    lr <- lr[as.character(lr$target) %in% receivers, , drop = FALSE]
  }
  lr
}

keys <- unique(dplyr::bind_rows(
  cellchat_index[, intersect(c("pair_id", "layer_id", "condition_value"), colnames(cellchat_index)), drop = FALSE],
  nichenet_index[, intersect(c("pair_id", "layer_id", "condition_value"), colnames(nichenet_index)), drop = FALSE]
))
if (nrow(keys) == 0) {
  keys <- empty_df_07(c("pair_id", "layer_id", "condition_value"))
}
for (col in c("pair_id", "layer_id", "condition_value")) {
  if (!col %in% colnames(keys)) keys[[col]] <- character(nrow(keys))
}
keys <- keys[nzchar(keys$pair_id) & nzchar(keys$layer_id), , drop = FALSE]

cross_rows <- list()
consensus_rows <- list()
triage_rows <- list()

for (idx in seq_len(nrow(keys))) {
  key <- keys[idx, , drop = FALSE]
  cc_hit <- cellchat_index[
    cellchat_index$pair_id == key$pair_id[[1]] &
      cellchat_index$layer_id == key$layer_id[[1]] &
      cellchat_index$condition_value == key$condition_value[[1]],
    ,
    drop = FALSE
  ]
  nn_hit <- nichenet_index[
    nichenet_index$pair_id == key$pair_id[[1]] &
      nichenet_index$layer_id == key$layer_id[[1]] &
      nichenet_index$condition_value == key$condition_value[[1]],
    ,
    drop = FALSE
  ]
  cc_status <- if (nrow(cc_hit) > 0) normalize_scalar_value(cc_hit$status[[1]], "missing") else "missing"
  nn_status <- if (nrow(nn_hit) > 0) normalize_scalar_value(nn_hit$status[[1]], "missing") else "missing"
  status <- "ok"
  if (!identical(cc_status, "ok") && !identical(nn_status, "ok")) {
    status <- "both_missing_or_failed"
  } else if (!identical(cc_status, "ok")) {
    status <- "cellchat_missing_or_failed"
  } else if (!identical(nn_status, "ok")) {
    status <- "nichenet_missing_or_failed"
  }

  lr <- if (nrow(cc_hit) > 0) read_lr_for_row_07c(cc_hit[1, , drop = FALSE]) else empty_df_07()
  activity <- if (nrow(nn_hit) > 0) read_ligand_activity_07c(nn_hit[1, , drop = FALSE]) else empty_df_07()
  if (nrow(nn_hit) > 0 && nrow(lr) > 0) {
    lr <- filter_lr_direction_07c(lr, nn_hit$sender_cell_types[[1]], nn_hit$receiver_cell_types[[1]])
  }
  lr_sig <- lr[is.finite(lr$prob) & lr$prob >= cfg$cellchat_prob_cutoff, , drop = FALSE]
  top_activity <- activity[is.finite(activity$rank) & activity$rank <= cfg$nichenet_top_ligand_n, , drop = FALSE]
  cellchat_ligands <- unique(as.character(lr_sig$ligand))
  nichenet_ligands <- unique(as.character(top_activity$test_ligand))
  cellchat_ligands <- cellchat_ligands[nzchar(cellchat_ligands)]
  nichenet_ligands <- nichenet_ligands[nzchar(nichenet_ligands)]
  intersect_ligands <- intersect(cellchat_ligands, nichenet_ligands)
  union_ligands <- union(cellchat_ligands, nichenet_ligands)
  jaccard <- safe_rate(length(intersect_ligands), length(union_ligands))

  cross_rows[[length(cross_rows) + 1L]] <- data.frame(
    pair_id = key$pair_id[[1]],
    layer_id = key$layer_id[[1]],
    condition_value = key$condition_value[[1]],
    cellchat_lr_n = nrow(lr_sig),
    nichenet_top_ligand_n = length(nichenet_ligands),
    intersect_ligand_n = length(intersect_ligands),
    jaccard = jaccard,
    status = status,
    stringsAsFactors = FALSE
  )

  if (length(intersect_ligands) == 0 && identical(status, "ok")) {
    triage_rows[[length(triage_rows) + 1L]] <- communication_triage_row_07(
      key$layer_id[[1]], key$pair_id[[1]], "warning", "communication_zero_consensus",
      sprintf("No ligand overlap for condition %s", key$condition_value[[1]]),
      "Review CellChat LR direction filters and NicheNet ligand activity ranking."
    )
  }
  if (identical(status, "ok") && is.finite(jaccard) && jaccard < 0.05) {
    triage_rows[[length(triage_rows) + 1L]] <- communication_triage_row_07(
      key$layer_id[[1]], key$pair_id[[1]], "warning", "communication_low_jaccard",
      sprintf("Jaccard %.3f for condition %s", jaccard, key$condition_value[[1]]),
      "Treat this pair as weakly cross-supported until manual review."
    )
  }

  if (nrow(lr_sig) > 0 && nrow(top_activity) > 0) {
    rank_lookup <- top_activity[, c("test_ligand", "rank", "score"), drop = FALSE]
    colnames(rank_lookup) <- c("ligand", "ligand_activity_rank", "ligand_activity_score")
    consensus <- merge(lr_sig[lr_sig$ligand %in% intersect_ligands, , drop = FALSE], rank_lookup, by = "ligand", all.x = TRUE)
    if (nrow(consensus) > 0) {
      consensus$pair_id <- key$pair_id[[1]]
      consensus$layer_id <- key$layer_id[[1]]
      consensus$condition_value <- key$condition_value[[1]]
      consensus_rows[[length(consensus_rows) + 1L]] <- consensus[, c(
        "pair_id", "layer_id", "condition_value", "source", "target",
        "ligand", "receptor", "prob", "pval", "pathway_name",
        "ligand_activity_rank", "ligand_activity_score"
      ), drop = FALSE]
    }
  }
}

cross_df <- if (length(cross_rows) > 0) dplyr::bind_rows(cross_rows) else empty_cross_validation_07()
consensus_df <- if (length(consensus_rows) > 0) dplyr::bind_rows(consensus_rows) else empty_consensus_lr_07()
combined_triage <- dplyr::bind_rows(
  if (nrow(cellchat_triage) > 0) cellchat_triage else empty_triage_df(include_sample = TRUE),
  if (nrow(nichenet_triage) > 0) nichenet_triage else empty_triage_df(include_sample = TRUE),
  if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)
)

write_tsv_local(cross_df, cfg$communication_cross_validation_tsv)
write_tsv_local(consensus_df, cfg$communication_consensus_lr_tsv)
write_tsv_local(combined_triage, cfg$communication_triage_tsv)

if (nrow(cross_df) > 0) {
  plot_df <- cross_df
  plot_df$jaccard <- suppressWarnings(as.numeric(plot_df$jaccard))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = condition_value, y = pair_id, fill = jaccard)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::facet_wrap(~layer_id, scales = "free_y") +
    ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey85") +
    ggplot2::labs(x = "Condition", y = "Communication pair", fill = "Jaccard", title = "CellChat/NicheNet ligand overlap") +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(cfg$communication_summary_plot_png, p, width = 10, height = 7, dpi = 300, bg = "white")
} else {
  write_empty_png_07(cfg$communication_summary_plot_png, "No communication cross-validation rows")
}

layer_report_outputs <- list()
for (layer_id in sort(unique(cross_df$layer_id))) {
  layer_cross <- cross_df[cross_df$layer_id == layer_id, , drop = FALSE]
  layer_consensus <- consensus_df[consensus_df$layer_id == layer_id, , drop = FALSE]
  layer_report <- communication_layer_report_path_07(cfg, layer_id)
  ensure_dir(dirname(layer_report))
  write_markdown_local(
    c(
      sprintf("# 07 Communication Layer Report: %s", layer_id),
      "",
      "## Cross Validation",
      render_markdown_table_local(layer_cross),
      "",
      "## Consensus LR",
      if (nrow(layer_consensus) == 0) "No consensus ligand-receptor rows." else render_markdown_table_local(head(layer_consensus, 100))
    ),
    layer_report
  )
  layer_report_outputs[[paste0("layer_report_md_", safe_id_07(layer_id))]] <-
    build_output_entry(layer_report, "md", module_name, "per-layer communication EDA report", base_dir = cfg$project_root)
}

status_summary <- if (nrow(cross_df) > 0) as.data.frame(table(cross_df$status), stringsAsFactors = FALSE) else data.frame(status = character(0), row_n = integer(0))
if (ncol(status_summary) == 2) {
  colnames(status_summary) <- c("status", "row_n")
}

report_lines <- c(
  "# 07c Communication EDA",
  "",
  sprintf("- CellChat index: `%s`", cfg$cellchat_index_tsv),
  sprintf("- NicheNet index: `%s`", cfg$nichenet_index_tsv),
  sprintf("- cross validation: `%s`", cfg$communication_cross_validation_tsv),
  sprintf("- consensus LR: `%s`", cfg$communication_consensus_lr_tsv),
  sprintf("- summary plot: `%s`", cfg$communication_summary_plot_png),
  "",
  "## Status Summary",
  render_markdown_table_local(status_summary),
  "",
  "## Cross Validation",
  render_markdown_table_local(cross_df),
  "",
  "## Consensus LR",
  if (nrow(consensus_df) == 0) "No consensus ligand-receptor rows." else render_markdown_table_local(head(consensus_df, 100)),
  "",
  "## Triage",
  render_markdown_table_local(combined_triage)
)
write_markdown_local(report_lines, cfg$communication_report_md)

if (file.exists(cfg$module_07c_manifest_path)) {
  unlink(cfg$module_07c_manifest_path)
}
fixed_outputs <- list(
  cross_validation_tsv = build_output_entry(cfg$communication_cross_validation_tsv, "tsv", module_name, "one row per layer/pair/condition CellChat-NicheNet overlap", base_dir = cfg$project_root, schema = infer_schema_from_df(cross_df)),
  consensus_lr_tsv = build_output_entry(cfg$communication_consensus_lr_tsv, "tsv", module_name, "CellChat LR rows supported by top NicheNet ligands", base_dir = cfg$project_root, schema = infer_schema_from_df(consensus_df)),
  report_md = build_output_entry(cfg$communication_report_md, "md", module_name, "communication EDA report", base_dir = cfg$project_root),
  summary_plot_png = build_output_entry(cfg$communication_summary_plot_png, "png", module_name, "communication cross-validation summary plot", base_dir = cfg$project_root),
  triage_tsv = build_output_entry(cfg$communication_triage_tsv, "tsv", module_name, "combined communication triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(combined_triage))
)
write_manifest_local(
  manifest_path = cfg$module_07c_manifest_path,
  new_outputs = c(fixed_outputs, layer_report_outputs),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    cellchat_index_tsv = cfg$cellchat_index_tsv,
    nichenet_index_tsv = cfg$nichenet_index_tsv,
    communication_pairs_sheet = cfg$communication_pairs_sheet
  ),
  version = cfg$module_version,
  depends_on = list(module_07a = cfg$module_07a_manifest_path, module_07b = cfg$module_07b_manifest_path)
)

message("07c completed. report: ", cfg$communication_report_md)
