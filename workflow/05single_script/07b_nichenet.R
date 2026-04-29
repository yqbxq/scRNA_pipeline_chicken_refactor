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

load_required_packages(c("Seurat", "nichenetr", "dplyr", "tibble", "jsonlite", "Matrix", "ggplot2", "circlize"))

cfg <- get_single_script_config_07()
module_name <- "07b_nichenet"
prepare_dirs_07(cfg)
set.seed(cfg$random_seed)

required_resources <- c(
  lr_network = cfg$nichenet_lr_network_rds,
  ligand_target_matrix = cfg$nichenet_ligand_target_matrix_rds,
  weighted_networks = cfg$nichenet_weighted_networks_rds
)
missing_resources <- required_resources[!file.exists(required_resources)]
if (length(missing_resources) > 0) {
  stop(
    sprintf(
      "NicheNet resources are missing: %s. Prepare them with workflow/06tools/download_nichenet_resources.sh or set NICHENET_*_RDS.",
      paste(sprintf("%s=%s", names(missing_resources), missing_resources), collapse = "; ")
    ),
    call. = FALSE
  )
}

ortholog_csv <- resolve_manifest_output_optional_07(
  cfg$ortholog_manifest_path,
  c("human_best", "chicken_human_best_csv", "human_best_csv")
)
if (!nzchar(ortholog_csv)) {
  ortholog_csv <- file.path(cfg$ortholog_cache_dir, "chicken_human_orthologs.csv")
}
lut <- build_chicken_to_human_lut(ortholog_csv)

standardize_lr_network_07b <- function(lr_network) {
  df <- as.data.frame(lr_network, stringsAsFactors = FALSE)
  ligand_col <- pick_first_col_07(df, c("from", "ligand", "source", "sender"))
  receptor_col <- pick_first_col_07(df, c("to", "receptor", "target"))
  if (!nzchar(ligand_col) || !nzchar(receptor_col)) {
    stop("NicheNet lr_network has no supported ligand/receptor columns.", call. = FALSE)
  }
  if (!identical(ligand_col, "from")) {
    df$from <- df[[ligand_col]]
  }
  if (!identical(receptor_col, "to")) {
    df$to <- df[[receptor_col]]
  }
  df
}

lr_network_chicken <- rename_network_human_to_chicken(
  standardize_lr_network_07b(readRDS(cfg$nichenet_lr_network_rds)),
  lut,
  c("from", "to")
)
ligand_target_matrix_chicken <- rename_ligand_target_matrix_human_to_chicken(readRDS(cfg$nichenet_ligand_target_matrix_rds), lut)

layers <- communication_layer_status_07(cfg)
pairs <- read_communication_pairs(cfg)

write_nichenet_empty_outputs_07b <- function(paths, reason = "No NicheNet result") {
  ensure_dir(paths$table_dir)
  ensure_dir(paths$figure_dir)
  write_tsv_local(empty_df_07(c("test_ligand", "auroc", "aupr", "pearson", "rank")), paths$ligand_activity_tsv)
  write_tsv_local(empty_df_07(c("ligand", "target", "weight")), paths$ligand_target_links_tsv)
  saveRDS(list(status = "empty", reason = reason), paths$nichenet_rds)
  write_empty_png_07(paths$ligand_activity_heatmap_png, reason)
  write_empty_png_07(paths$ligand_target_heatmap_png, reason)
  write_empty_png_07(paths$circos_png, reason)
}

expressed_genes_for_cells_07b <- function(seu, cells, pct = 0.1) {
  cells <- intersect(cells, colnames(seu))
  if (length(cells) == 0) {
    return(character(0))
  }
  obj <- subset(seu, cells = cells)
  mat <- get_assay_matrix(obj, assay = "RNA", type = "data")
  if (nrow(mat) == 0) {
    mat <- get_assay_matrix(obj, assay = "RNA", type = "counts")
  }
  frac <- Matrix::rowMeans(mat > 0)
  names(frac)[is.finite(frac) & frac >= pct]
}

receiver_marker_genes_for_nichenet_07b <- function(seu, receiver_cells, alpha = 0.05) {
  receiver_cells <- intersect(receiver_cells, colnames(seu))
  other_cells <- setdiff(colnames(seu), receiver_cells)
  if (length(receiver_cells) < 3L || length(other_cells) < 3L) {
    return(list(genes = character(0), status = "too_few_receiver_marker_cells", table = data.frame(stringsAsFactors = FALSE)))
  }
  obj <- subset(seu, cells = c(receiver_cells, other_cells))
  obj$nichenet_receiver_status <- ifelse(colnames(obj) %in% receiver_cells, "receiver", "other")
  Seurat::Idents(obj) <- "nichenet_receiver_status"
  res <- tryCatch(
    Seurat::FindMarkers(
      obj,
      ident.1 = "receiver",
      ident.2 = "other",
      logfc.threshold = 0,
      test.use = "wilcox",
      verbose = FALSE
    ),
    error = function(e) e
  )
  if (inherits(res, "error") || is.null(res) || nrow(res) == 0) {
    return(list(genes = character(0), status = "receiver_findmarkers_error", table = data.frame(stringsAsFactors = FALSE)))
  }
  out <- tibble::rownames_to_column(as.data.frame(res), "gene")
  std <- standardize_deg_table_06(out, alpha = alpha)
  genes <- unique(std$.gene[std$.significant & (!is.finite(std$.logfc) | std$.logfc > 0)])
  list(genes = genes, status = "receiver_marker_findmarkers", table = out)
}

condition_values_for_failed_pair_07b <- function(pair_row) {
  values <- split_csv_local(pair_row$condition_split_values[[1]])
  if (length(values) == 0) "all" else values
}

plot_ligand_activity_07b <- function(activity_df, path) {
  ensure_dir(dirname(path))
  if (nrow(activity_df) == 0) {
    return(write_empty_png_07(path, "No ligand activity"))
  }
  score_col <- pick_first_col_07(activity_df, c("pearson", "aupr_corrected", "aupr", "auroc"))
  if (!nzchar(score_col)) {
    return(write_empty_png_07(path, "No supported ligand activity score"))
  }
  df <- activity_df
  df$score <- suppressWarnings(as.numeric(df[[score_col]]))
  df <- df[is.finite(df$score), , drop = FALSE]
  df <- head(df[order(-df$score), , drop = FALSE], 30)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(test_ligand, score), y = score, fill = score)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_viridis_c(option = "C") +
    ggplot2::labs(x = "Ligand", y = score_col, title = "Top NicheNet ligand activity") +
    ggplot2::theme_bw(base_size = 11)
  ggplot2::ggsave(path, p, width = 8, height = 7, dpi = 300, bg = "white")
  path
}

plot_ligand_target_heatmap_07b <- function(links_df, path) {
  ensure_dir(dirname(path))
  if (nrow(links_df) == 0) {
    return(write_empty_png_07(path, "No ligand-target links"))
  }
  df <- links_df
  df$weight <- suppressWarnings(as.numeric(df$weight))
  df <- df[is.finite(df$weight), , drop = FALSE]
  df <- head(df[order(-df$weight), , drop = FALSE], 120)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = ligand, y = target, fill = weight)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "C") +
    ggplot2::labs(x = "Ligand", y = "Target", title = "Top ligand-target links") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(path, p, width = 9, height = 10, dpi = 300, bg = "white")
  path
}

plot_circos_07b <- function(links_df, path) {
  ensure_dir(dirname(path))
  if (nrow(links_df) == 0) {
    return(write_empty_png_07(path, "No ligand-target links"))
  }
  df <- links_df
  df$weight <- suppressWarnings(as.numeric(df$weight))
  df <- df[is.finite(df$weight), , drop = FALSE]
  df <- head(df[order(-df$weight), , drop = FALSE], 80)
  if (nrow(df) == 0) {
    return(write_empty_png_07(path, "No weighted links"))
  }
  grDevices::png(path, width = 1400, height = 1400, res = 180)
  tryCatch(
    {
      circlize::chordDiagram(df[, c("ligand", "target", "weight")], transparency = 0.35)
      circlize::circos.clear()
      grDevices::dev.off()
    },
    error = function(e) {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
      write_empty_png_07(path, paste("Circos failed:", conditionMessage(e)))
    }
  )
  path
}

run_nichenet_one_07b <- function(geneset, background, potential_ligands, paths) {
  geneset <- unique(intersect(geneset, rownames(ligand_target_matrix_chicken)))
  background <- unique(intersect(background, rownames(ligand_target_matrix_chicken)))
  potential_ligands <- unique(intersect(potential_ligands, colnames(ligand_target_matrix_chicken)))

  if (length(geneset) < 5L) {
    return(list(status = "too_few_geneset_genes", reason = sprintf("%s genes in NicheNet target universe", length(geneset))))
  }
  if (length(background) < 20L) {
    return(list(status = "too_few_background_genes", reason = sprintf("%s background genes in NicheNet target universe", length(background))))
  }
  if (length(potential_ligands) == 0) {
    return(list(status = "no_potential_ligands", reason = "0 potential ligands after sender/receptor filtering"))
  }

  activity <- nichenetr::predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix_chicken,
    potential_ligands = potential_ligands
  )
  activity <- as.data.frame(activity, stringsAsFactors = FALSE)
  if (!"test_ligand" %in% colnames(activity)) {
    activity$test_ligand <- rownames(activity)
  }
  score_col <- pick_first_col_07(activity, c("pearson", "aupr_corrected", "aupr", "auroc"))
  if (nzchar(score_col)) {
    activity <- activity[order(-suppressWarnings(as.numeric(activity[[score_col]]))), , drop = FALSE]
  }
  activity$rank <- seq_len(nrow(activity))
  top_ligands <- head(activity$test_ligand, cfg$nichenet_top_ligand_n)

  links <- tryCatch(
    {
      per_ligand <- lapply(top_ligands, function(lg) {
        tryCatch(
          nichenetr::get_weighted_ligand_target_links(
            ligand = lg,
            geneset = geneset,
            ligand_target_matrix = ligand_target_matrix_chicken,
            n = cfg$nichenet_top_target_n
          ),
          error = function(e) NULL
        )
      })
      per_ligand <- per_ligand[!vapply(per_ligand, is.null, logical(1))]
      if (length(per_ligand) == 0) {
        data.frame(ligand = character(0), target = character(0), weight = numeric(0), stringsAsFactors = FALSE)
      } else {
        dplyr::bind_rows(per_ligand)
      }
    },
    error = function(e) data.frame(ligand = character(0), target = character(0), weight = numeric(0), stringsAsFactors = FALSE)
  )
  links <- as.data.frame(links, stringsAsFactors = FALSE)
  if (nrow(links) == 0) {
    links <- data.frame(ligand = character(0), target = character(0), weight = numeric(0), stringsAsFactors = FALSE)
  } else {
    for (col in c("ligand", "target", "weight")) {
      if (!col %in% colnames(links)) {
        links[[col]] <- if (identical(col, "weight")) NA_real_ else ""
      }
    }
  }
  links <- links[, c("ligand", "target", "weight"), drop = FALSE]

  ensure_dir(paths$table_dir)
  write_tsv_local(activity, paths$ligand_activity_tsv)
  write_tsv_local(links, paths$ligand_target_links_tsv)
  saveRDS(
    list(
      ligand_activities = activity,
      ligand_target_links = links,
      geneset = geneset,
      background = background,
      potential_ligands = potential_ligands,
      weighted_networks_resource = cfg$nichenet_weighted_networks_rds
    ),
    paths$nichenet_rds
  )
  plot_ligand_activity_07b(activity, paths$ligand_activity_heatmap_png)
  plot_ligand_target_heatmap_07b(links, paths$ligand_target_heatmap_png)
  plot_circos_07b(links, paths$circos_png)
  list(status = "ok", reason = "")
}

append_nichenet_row_07b <- function(rows, pair_row, layer_id, condition_value, roles, deg_source, deg_status, status, reason, paths) {
  rows[[length(rows) + 1L]] <- data.frame(
    pair_id = pair_row$pair_id[[1]],
    layer_id = layer_id,
    condition_value = condition_value,
    condition_split_var = normalize_scalar_value(pair_row$condition_split_var[[1]]),
    sender_set = normalize_scalar_value(pair_row$sender[[1]], "*"),
    receiver_set = normalize_scalar_value(pair_row$receiver[[1]], "*"),
    sender_cell_types = roles$sender_label %||% "",
    receiver_cell_types = roles$receiver_label %||% "",
    deg_source = deg_source,
    deg_status = deg_status,
    status = status,
    reason = reason,
    ligand_activity_tsv = normalize_path_07(paths$ligand_activity_tsv),
    ligand_target_links_tsv = normalize_path_07(paths$ligand_target_links_tsv),
    nichenet_rds = normalize_path_07(paths$nichenet_rds),
    ligand_activity_heatmap_png = normalize_path_07(paths$ligand_activity_heatmap_png),
    ligand_target_heatmap_png = normalize_path_07(paths$ligand_target_heatmap_png),
    circos_png = normalize_path_07(paths$circos_png),
    stringsAsFactors = FALSE
  )
  rows
}

index_rows <- list()
roles_rows <- list()
triage_rows <- list()
dynamic_outputs <- list()

for (layer_idx in seq_len(nrow(layers))) {
  layer_row <- layers[layer_idx, , drop = FALSE]
  layer_id <- layer_row$layer_id[[1]]
  message("07b NicheNet layer: ", layer_id)
  seu <- load_comm_layer_object_07(layer_row)
  cell_type_col <- resolve_cell_type_col_07(seu, layer_id, layer_row$layer_role[[1]])
  if (!nzchar(cell_type_col)) {
    triage_rows[[length(triage_rows) + 1L]] <- communication_triage_row_07(
      layer_id, "", "error", "communication_missing_cell_type",
      sprintf("Layer %s has no supported cell type metadata column.", layer_id),
      "Check 03/04 annotation outputs and metadata column names."
    )
    next
  }

  layer_pairs <- communication_pairs_for_layer(pairs, layer_id, "nichenet")
  for (pair_idx in seq_len(nrow(layer_pairs))) {
    pair_row <- layer_pairs[pair_idx, , drop = FALSE]
    pair_id <- pair_row$pair_id[[1]]
    base_subset <- subset_by_pair_filter_07(seu, pair_row)
    if (!identical(base_subset$status, "ok")) {
      for (condition_value in condition_values_for_failed_pair_07b(pair_row)) {
        paths <- nichenet_paths_07(cfg, layer_id, pair_id, condition_value)
        write_nichenet_empty_outputs_07b(paths, base_subset$reason)
        index_rows <- append_nichenet_row_07b(index_rows, pair_row, layer_id, condition_value, list(sender_label = "", receiver_label = ""), "none", "missing", base_subset$status, base_subset$reason, paths)
      }
      next
    }

    splits <- split_by_condition(base_subset$object, pair_row)
    for (condition_name in names(splits)) {
      split_item <- splits[[condition_name]]
      condition_value <- normalize_scalar_value(split_item$condition_value, "all")
      paths <- nichenet_paths_07(cfg, layer_id, pair_id, condition_value)
      roles <- list(sender = character(0), receiver = character(0), sender_label = "", receiver_label = "", missing = character(0))
      status <- split_item$status
      reason <- split_item$reason
      deg_source <- "none"
      deg_status <- "missing"

      if (identical(status, "ok")) {
        obj <- split_item$object
        roles <- resolve_sender_receiver_sets(pair_row, obj, cell_type_col)
        roles_rows[[length(roles_rows) + 1L]] <- data.frame(
          pair_id = pair_id,
          layer_id = layer_id,
          condition_value = condition_value,
          sender_set = normalize_scalar_value(pair_row$sender[[1]], "*"),
          receiver_set = normalize_scalar_value(pair_row$receiver[[1]], "*"),
          sender_cell_types = roles$sender_label,
          receiver_cell_types = roles$receiver_label,
          missing_cell_types = paste(roles$missing, collapse = ","),
          status = ifelse(length(roles$sender) > 0 && length(roles$receiver) > 0, "ok", "missing_roles"),
          stringsAsFactors = FALSE
        )
        if (length(roles$sender) == 0 || length(roles$receiver) == 0) {
          status <- "missing_roles"
          reason <- sprintf("Missing sender/receiver roles: %s", paste(roles$missing, collapse = ","))
        }
      }

      if (identical(status, "ok")) {
        obj <- split_item$object
        meta <- obj@meta.data
        sender_cells <- rownames(meta)[as.character(meta[[cell_type_col]]) %in% roles$sender]
        receiver_cells <- rownames(meta)[as.character(meta[[cell_type_col]]) %in% roles$receiver]
        deg <- read_deg_for_nichenet_07(cfg, layer_id, pair_id)
        geneset <- deg$genes
        deg_source <- deg$source
        deg_status <- deg$status
        if (length(geneset) == 0) {
          receiver_markers <- receiver_marker_genes_for_nichenet_07b(obj, receiver_cells, alpha = cfg$deg_alpha)
          geneset <- receiver_markers$genes
          deg_source <- "receiver_marker_findmarkers"
          deg_status <- receiver_markers$status
        }

        sender_expressed <- expressed_genes_for_cells_07b(obj, sender_cells, cfg$nichenet_expression_pct)
        receiver_expressed <- expressed_genes_for_cells_07b(obj, receiver_cells, cfg$nichenet_expression_pct)
        lr_hit <- lr_network_chicken[
          lr_network_chicken$from %in% sender_expressed & lr_network_chicken$to %in% receiver_expressed,
          ,
          drop = FALSE
        ]
        potential_ligands <- unique(lr_hit$from)
        run_res <- tryCatch(
          run_nichenet_one_07b(
            geneset = geneset,
            background = receiver_expressed,
            potential_ligands = potential_ligands,
            paths = paths
          ),
          error = function(e) e
        )
        if (inherits(run_res, "error")) {
          status <- "error"
          reason <- conditionMessage(run_res)
        } else {
          status <- run_res$status
          reason <- run_res$reason
        }
      }

      if (!identical(status, "ok")) {
        write_nichenet_empty_outputs_07b(paths, reason)
        triage_rows[[length(triage_rows) + 1L]] <- communication_triage_row_07(
          layer_id, pair_id, ifelse(identical(status, "error"), "error", "warning"),
          paste0("nichenet_", status),
          sprintf("%s/%s: %s", condition_value, status, reason),
          "Review communication_pairs.tsv roles, condition split, DEG availability, and NicheNet resources."
        )
      }

      index_rows <- append_nichenet_row_07b(index_rows, pair_row, layer_id, condition_value, roles, deg_source, deg_status, status, reason, paths)
      suffix <- paste(safe_id_07(layer_id), safe_id_07(pair_id), safe_id_07(condition_value), sep = "_")
      dynamic_outputs[[paste0("ligand_activity_", suffix)]] <- build_output_entry(paths$ligand_activity_tsv, "tsv", module_name, "NicheNet ligand activity table", base_dir = cfg$project_root)
      dynamic_outputs[[paste0("ligand_target_links_", suffix)]] <- build_output_entry(paths$ligand_target_links_tsv, "tsv", module_name, "NicheNet ligand-target links", base_dir = cfg$project_root)
      dynamic_outputs[[paste0("nichenet_rds_", suffix)]] <- build_output_entry(paths$nichenet_rds, "rds", module_name, "NicheNet task object", base_dir = cfg$project_root)
      dynamic_outputs[[paste0("ligand_activity_heatmap_png_", suffix)]] <- build_output_entry(paths$ligand_activity_heatmap_png, "png", module_name, "NicheNet ligand activity plot", base_dir = cfg$project_root)
      dynamic_outputs[[paste0("ligand_target_heatmap_png_", suffix)]] <- build_output_entry(paths$ligand_target_heatmap_png, "png", module_name, "NicheNet ligand-target heatmap", base_dir = cfg$project_root)
      dynamic_outputs[[paste0("circos_png_", suffix)]] <- build_output_entry(paths$circos_png, "png", module_name, "NicheNet ligand-target circos plot", base_dir = cfg$project_root)
    }
  }
}

index_df <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else empty_nichenet_index_07()
roles_df <- if (length(roles_rows) > 0) dplyr::bind_rows(roles_rows) else empty_roles_resolved_07()
triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)

write_tsv_local(index_df, cfg$nichenet_index_tsv)
write_tsv_local(roles_df, cfg$nichenet_roles_resolved_tsv)
write_tsv_local(triage_df, cfg$nichenet_triage_tsv)

if (file.exists(cfg$module_07b_manifest_path)) {
  unlink(cfg$module_07b_manifest_path)
}
fixed_outputs <- list(
  nichenet_index_tsv = build_output_entry(cfg$nichenet_index_tsv, "tsv", module_name, "one row per layer/pair/condition NicheNet task", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df)),
  roles_resolved_tsv = build_output_entry(cfg$nichenet_roles_resolved_tsv, "tsv", module_name, "resolved sender/receiver roles", base_dir = cfg$project_root, schema = infer_schema_from_df(roles_df)),
  triage_tsv = build_output_entry(cfg$nichenet_triage_tsv, "tsv", module_name, "NicheNet triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df))
)
write_manifest_local(
  manifest_path = cfg$module_07b_manifest_path,
  new_outputs = c(fixed_outputs, dynamic_outputs),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    ortholog_csv = ortholog_csv,
    lr_network_rds = cfg$nichenet_lr_network_rds,
    ligand_target_matrix_rds = cfg$nichenet_ligand_target_matrix_rds,
    weighted_networks_rds = cfg$nichenet_weighted_networks_rds,
    communication_pairs_sheet = cfg$communication_pairs_sheet,
    communication_cell_type_col = cfg$communication_cell_type_col,
    module_05d = cfg$module_05d_manifest_path
  ),
  version = cfg$module_version,
  depends_on = list(
    module_00 = cfg$ortholog_manifest_path,
    module_05d = cfg$module_05d_manifest_path,
    module_07a = cfg$module_07a_manifest_path
  )
)

message("07b completed. index: ", cfg$nichenet_index_tsv)
