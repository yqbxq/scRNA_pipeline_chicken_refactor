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
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "comparison_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "deg_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "enrichment_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_mapping_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_pairs_utils.R"))

load_required_packages(c("Seurat", "CellChat", "dplyr", "tibble", "jsonlite", "Matrix", "ggplot2", "patchwork"))

cfg <- get_single_script_config_07()
module_name <- "07a_cellchat"
prepare_dirs_07(cfg)
set.seed(cfg$random_seed)

ortholog_csv <- resolve_manifest_output_optional_07(
  cfg$ortholog_manifest_path,
  c("human_best", "chicken_human_best_csv", "human_best_csv")
)
if (!nzchar(ortholog_csv)) {
  ortholog_csv <- file.path(cfg$ortholog_cache_dir, "chicken_human_orthologs.csv")
}
lut <- build_chicken_to_human_lut(ortholog_csv)
layers <- communication_layer_status_07(cfg)
pairs <- read_communication_pairs(cfg)

empty_lr_table_07a <- function() {
  empty_df_07(c("source", "target", "ligand", "receptor", "prob", "pval", "pathway_name"))
}

empty_pathway_table_07a <- function() {
  empty_df_07(c("pathway_name", "source", "target", "prob", "pval"))
}

write_cellchat_empty_outputs_07a <- function(paths, reason = "No CellChat result") {
  ensure_dir(paths$table_dir)
  ensure_dir(paths$figure_dir)
  write_tsv_local(empty_lr_table_07a(), paths$lr_table_tsv)
  write_tsv_local(empty_pathway_table_07a(), paths$pathway_table_tsv)
  write_empty_png_07(paths$bubble_png, reason)
  write_empty_png_07(paths$network_png, reason)
  write_empty_png_07(paths$heatmap_png, reason)
}

format_lr_table_07a <- function(cc) {
  lr <- tryCatch(CellChat::subsetCommunication(cc), error = function(e) data.frame(stringsAsFactors = FALSE))
  if (nrow(lr) == 0) {
    return(empty_lr_table_07a())
  }
  for (col in c("source", "target", "ligand", "receptor", "prob", "pval", "pathway_name")) {
    if (!col %in% colnames(lr)) {
      lr[[col]] <- if (col %in% c("prob", "pval")) NA_real_ else ""
    }
  }
  lr[, c("source", "target", "ligand", "receptor", "prob", "pval", "pathway_name"), drop = FALSE]
}

format_pathway_table_07a <- function(cc) {
  pw <- tryCatch(CellChat::subsetCommunication(cc, slot.name = "netP"), error = function(e) data.frame(stringsAsFactors = FALSE))
  if (nrow(pw) == 0) {
    return(empty_pathway_table_07a())
  }
  for (col in c("pathway_name", "source", "target", "prob", "pval")) {
    if (!col %in% colnames(pw)) {
      pw[[col]] <- if (col %in% c("prob", "pval")) NA_real_ else ""
    }
  }
  pw[, c("pathway_name", "source", "target", "prob", "pval"), drop = FALSE]
}

cellchat_db_human_07a <- function() {
  db <- tryCatch(get("CellChatDB.human", envir = asNamespace("CellChat")), error = function(e) NULL)
  if (!is.null(db)) {
    return(db)
  }
  data("CellChatDB.human", package = "CellChat", envir = environment())
  get("CellChatDB.human", envir = environment())
}

save_cellchat_plots_07a <- function(cc, paths) {
  ensure_dir(paths$figure_dir)
  tryCatch(
    {
      grDevices::png(paths$bubble_png, width = 1800, height = 1200, res = 180)
      print(CellChat::netVisual_bubble(cc, remove.isolate = TRUE))
      grDevices::dev.off()
    },
    error = function(e) {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
      write_empty_png_07(paths$bubble_png, paste("CellChat bubble plot failed:", conditionMessage(e)))
    }
  )
  tryCatch(
    {
      grDevices::png(paths$network_png, width = 1400, height = 1200, res = 180)
      group_size <- as.numeric(table(cc@idents))
      CellChat::netVisual_circle(cc@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE)
      grDevices::dev.off()
    },
    error = function(e) {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
      write_empty_png_07(paths$network_png, paste("CellChat network plot failed:", conditionMessage(e)))
    }
  )
  tryCatch(
    {
      grDevices::png(paths$heatmap_png, width = 1400, height = 1200, res = 180)
      print(CellChat::netVisual_heatmap(cc))
      grDevices::dev.off()
    },
    error = function(e) {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
      write_empty_png_07(paths$heatmap_png, paste("CellChat heatmap failed:", conditionMessage(e)))
    }
  )
}

run_cellchat_one_07a <- function(seu, mat_human, cell_type_col, paths, sample_label, pair_row = NULL, roles = NULL) {
  meta <- seu@meta.data
  meta$cellchat_group <- as.character(meta[[cell_type_col]])
  meta$samples <- sample_label
  mat_human <- mat_human[, colnames(mat_human) %in% rownames(meta), drop = FALSE]
  meta <- meta[colnames(mat_human), , drop = FALSE]

  cc <- CellChat::createCellChat(object = mat_human, meta = meta, group.by = "cellchat_group")
  cc@DB <- cellchat_db_human_07a()
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)

  ensure_dir(paths$table_dir)
  saveRDS(cc, paths$cellchat_rds)
  lr <- format_lr_table_07a(cc)
  pathway <- format_pathway_table_07a(cc)
  if (!is.null(pair_row) &&
      exists("communication_pair_direction_filter", mode = "function") &&
      communication_pair_direction_filter(pair_row)) {
    lr <- filter_communication_direction_df_07(lr, roles)
    pathway <- filter_communication_direction_df_07(pathway, roles)
  }
  write_tsv_local(lr, paths$lr_table_tsv)
  write_tsv_local(pathway, paths$pathway_table_tsv)
  save_cellchat_plots_07a(cc, paths)
}

condition_values_for_failed_pair_07a <- function(pair_row) {
  values <- split_csv_local(pair_row$condition_split_values[[1]])
  if (length(values) == 0) "all" else values
}

append_cellchat_row_07a <- function(rows, pair_row, layer_id, condition_value, cell_type_col, n_cells, n_cell_types, status, reason, paths) {
  rows[[length(rows) + 1L]] <- data.frame(
    pair_id = pair_row$pair_id[[1]],
    layer_id = layer_id,
    condition_value = condition_value,
    condition_split_var = normalize_scalar_value(pair_row$condition_split_var[[1]]),
    sender = normalize_scalar_value(pair_row$sender[[1]], "*"),
    receiver = normalize_scalar_value(pair_row$receiver[[1]], "*"),
    cell_type_col = cell_type_col,
    direction_filter = normalize_scalar_value(pair_row$direction_filter[[1]], "yes"),
    requires_cell_subtype = normalize_scalar_value(pair_row$requires_cell_subtype[[1]], "no"),
    n_cells = n_cells,
    n_cell_types = n_cell_types,
    status = status,
    reason = reason,
    cellchat_rds_path = if (file.exists(paths$cellchat_rds)) normalize_path_07(paths$cellchat_rds) else "",
    lr_table_path = normalize_path_07(paths$lr_table_tsv),
    pathway_table_path = normalize_path_07(paths$pathway_table_tsv),
    bubble_png = normalize_path_07(paths$bubble_png),
    network_png = normalize_path_07(paths$network_png),
    heatmap_png = normalize_path_07(paths$heatmap_png),
    stringsAsFactors = FALSE
  )
  rows
}

index_rows <- list()
mapping_rows <- list()
triage_rows <- list()
dynamic_outputs <- list()

if (nrow(layers) == 0) {
  triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
    sample_id = "__PROJECT__",
    severity = "error",
    signal_id = "communication_no_layers",
    evidence = sprintf("No readable annotated layers found via %s, %s, or %s", cfg$layer_status_file, cfg$module_03d_manifest_path, cfg$module_04b_manifest_path),
    recommended_action = "Run 03d annotation first, then rerun 07a.",
    manual_review_required = "yes"
  )
}

for (layer_idx in seq_len(nrow(layers))) {
  layer_row <- layers[layer_idx, , drop = FALSE]
  layer_id <- layer_row$layer_id[[1]]
  message("07a CellChat layer: ", layer_id)
  seu <- load_comm_layer_object_07(layer_row)

  layer_pairs <- communication_pairs_for_layer(pairs, layer_id, "cellchat")
  for (pair_idx in seq_len(nrow(layer_pairs))) {
    pair_row <- layer_pairs[pair_idx, , drop = FALSE]
    pair_id <- pair_row$pair_id[[1]]
    col_result <- resolve_pair_cell_type_col_07(seu, layer_id, layer_row$layer_role[[1]], pair_row)
    cell_type_col <- col_result$column
    if (!identical(col_result$status, "ok")) {
      triage_rows[[length(triage_rows) + 1L]] <- communication_triage_row_07(
        layer_id, pair_id, "error", col_result$status,
        col_result$reason,
        "Run M5 subtype backfill or fix communication metadata before rerunning 07a."
      )
      for (condition_value in condition_values_for_failed_pair_07a(pair_row)) {
        paths <- cellchat_paths_07(cfg, layer_id, pair_id, condition_value)
        write_cellchat_empty_outputs_07a(paths, col_result$reason)
        index_rows <- append_cellchat_row_07a(index_rows, pair_row, layer_id, condition_value, cell_type_col, 0L, 0L, col_result$status, col_result$reason, paths)
      }
      next
    }

    base_subset <- subset_by_pair_filter_07(seu, pair_row)
    if (!identical(base_subset$status, "ok")) {
      for (condition_value in condition_values_for_failed_pair_07a(pair_row)) {
        paths <- cellchat_paths_07(cfg, layer_id, pair_id, condition_value)
        write_cellchat_empty_outputs_07a(paths, base_subset$reason)
        index_rows <- append_cellchat_row_07a(index_rows, pair_row, layer_id, condition_value, cell_type_col, 0L, 0L, base_subset$status, base_subset$reason, paths)
      }
      next
    }

    splits <- split_by_condition(base_subset$object, pair_row)
    for (condition_name in names(splits)) {
      split_item <- splits[[condition_name]]
      condition_value <- normalize_scalar_value(split_item$condition_value, "all")
      paths <- cellchat_paths_07(cfg, layer_id, pair_id, condition_value)

      if (!identical(split_item$status, "ok")) {
        write_cellchat_empty_outputs_07a(paths, split_item$reason)
        index_rows <- append_cellchat_row_07a(index_rows, pair_row, layer_id, condition_value, cell_type_col, 0L, 0L, split_item$status, split_item$reason, paths)
        next
      }

      obj <- split_item$object
      roles <- resolve_sender_receiver_sets(pair_row, obj, cell_type_col)
      role_check <- validate_resolved_roles_07(pair_row, roles, strict = isTRUE(col_result$strict))
      if (!identical(role_check$status, "ok")) {
        triage_rows[[length(triage_rows) + 1L]] <- communication_triage_row_07(
          layer_id, pair_id, role_check$severity, role_check$status,
          sprintf("%s/%s: %s", pair_id, condition_value, role_check$reason),
          "Check communication_pairs.tsv sender/receiver labels against annotation metadata."
        )
      }

      n_cells <- ncol(obj)
      n_cell_types <- length(unique(as.character(obj@meta.data[[cell_type_col]])))
      status <- "ok"
      reason <- ""
      if (isTRUE(col_result$strict) && !identical(role_check$status, "ok")) {
        status <- role_check$status
        reason <- role_check$reason
      }
      if (identical(status, "ok") && n_cells < cfg$cellchat_min_cells_per_group) {
        status <- "skipped_low_n"
        reason <- sprintf("%s cells; minimum is %s", n_cells, cfg$cellchat_min_cells_per_group)
      }

      if (identical(status, "ok")) {
        counts_by_type <- table(as.character(obj@meta.data[[cell_type_col]]))
        keep_types <- names(counts_by_type)[counts_by_type >= cfg$communication_min_cells_per_celltype]
        if (length(keep_types) < 2L) {
          status <- "skipped_low_celltype_n"
          reason <- sprintf("fewer than 2 cell types have >= %s cells", cfg$communication_min_cells_per_celltype)
        } else {
          keep_cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[cell_type_col]]) %in% keep_types]
          obj <- subset(obj, cells = keep_cells)
          expr <- get_assay_matrix(obj, assay = "RNA", type = "data")
          if (nrow(expr) == 0) {
            expr <- get_assay_matrix(obj, assay = "RNA", type = "counts")
          }
          coverage <- mapping_coverage_summary(rownames(expr), lut)
          mapping_rows[[length(mapping_rows) + 1L]] <- cbind(
            data.frame(pair_id = pair_id, layer_id = layer_id, condition_value = condition_value, stringsAsFactors = FALSE),
            coverage
          )
          if (is.finite(coverage$coverage_pct[[1]]) && coverage$coverage_pct[[1]] / 100 < cfg$communication_ortholog_min_coverage) {
            triage_rows[[length(triage_rows) + 1L]] <- communication_triage_row_07(
              layer_id, pair_id, "warning", "low_ortholog_coverage",
              sprintf("%s/%s/%s coverage %.1f%%", layer_id, pair_id, condition_value, coverage$coverage_pct[[1]]),
              "Check gene symbols and ortholog cache before interpreting CellChat output."
            )
          }
          mat_human <- map_expression_matrix_chicken_to_human(expr, lut, collapse = "max")
          if (nrow(mat_human) == 0) {
            status <- "skipped_no_orthologs"
            reason <- "0 genes retained after chicken-to-human mapping"
          } else {
            run_res <- tryCatch(
              run_cellchat_one_07a(
                obj,
                mat_human,
                cell_type_col,
                paths,
                sample_label = paste(layer_id, condition_value, sep = "_"),
                pair_row = pair_row,
                roles = roles
              ),
              error = function(e) e
            )
            if (inherits(run_res, "error")) {
              status <- "error"
              reason <- conditionMessage(run_res)
            }
          }
        }
      }

      if (!identical(status, "ok")) {
        write_cellchat_empty_outputs_07a(paths, reason)
      }
      index_rows <- append_cellchat_row_07a(index_rows, pair_row, layer_id, condition_value, cell_type_col, n_cells, n_cell_types, status, reason, paths)

      suffix <- paste(safe_id_07(layer_id), safe_id_07(pair_id), safe_id_07(condition_value), sep = "_")
      if (file.exists(paths$cellchat_rds)) {
        dynamic_outputs[[paste0("cellchat_rds_", suffix)]] <- build_output_entry(paths$cellchat_rds, "rds", module_name, "CellChat object for one layer/pair/condition", base_dir = cfg$project_root)
      }
      dynamic_outputs[[paste0("lr_table_", suffix)]] <- build_output_entry(paths$lr_table_tsv, "tsv", module_name, "CellChat LR pairs for one layer/pair/condition", base_dir = cfg$project_root)
      dynamic_outputs[[paste0("pathway_table_", suffix)]] <- build_output_entry(paths$pathway_table_tsv, "tsv", module_name, "CellChat pathway table for one layer/pair/condition", base_dir = cfg$project_root)
      dynamic_outputs[[paste0("bubble_png_", suffix)]] <- build_output_entry(paths$bubble_png, "png", module_name, "CellChat bubble plot", base_dir = cfg$project_root)
      dynamic_outputs[[paste0("network_png_", suffix)]] <- build_output_entry(paths$network_png, "png", module_name, "CellChat network plot", base_dir = cfg$project_root)
      dynamic_outputs[[paste0("heatmap_png_", suffix)]] <- build_output_entry(paths$heatmap_png, "png", module_name, "CellChat heatmap plot", base_dir = cfg$project_root)
    }
  }
}

index_df <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else empty_cellchat_index_07()
mapping_df <- if (length(mapping_rows) > 0) dplyr::bind_rows(mapping_rows) else empty_mapping_summary_07()
triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)

write_tsv_local(index_df, cfg$cellchat_index_tsv)
write_tsv_local(mapping_df, cfg$cellchat_mapping_summary_tsv)
write_tsv_local(triage_df, cfg$cellchat_triage_tsv)

if (file.exists(cfg$module_07a_manifest_path)) {
  unlink(cfg$module_07a_manifest_path)
}
fixed_outputs <- list(
  cellchat_index_tsv = build_output_entry(cfg$cellchat_index_tsv, "tsv", module_name, "one row per layer/pair/condition CellChat task", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df)),
  mapping_summary_tsv = build_output_entry(cfg$cellchat_mapping_summary_tsv, "tsv", module_name, "one row per layer/pair/condition ortholog mapping summary", base_dir = cfg$project_root, schema = infer_schema_from_df(mapping_df)),
  triage_tsv = build_output_entry(cfg$cellchat_triage_tsv, "tsv", module_name, "CellChat triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df))
)
write_manifest_local(
  manifest_path = cfg$module_07a_manifest_path,
  new_outputs = c(fixed_outputs, dynamic_outputs),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    ortholog_csv = ortholog_csv,
    layer_status_tsv = cfg$layer_status_file,
    communication_pairs_sheet = cfg$communication_pairs_sheet,
    communication_cell_type_col = cfg$communication_cell_type_col,
    module_03d = cfg$module_03d_manifest_path,
    module_04b = cfg$module_04b_manifest_path
  ),
  version = cfg$module_version,
  depends_on = list(
    module_00 = cfg$ortholog_manifest_path,
    module_03d = cfg$module_03d_manifest_path,
    module_04b = cfg$module_04b_manifest_path
  )
)

message("07a completed. index: ", cfg$cellchat_index_tsv)
