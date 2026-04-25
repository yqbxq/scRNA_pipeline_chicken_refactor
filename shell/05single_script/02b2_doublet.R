#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
)

source_utf8 <- function(path) {
  source(path, encoding = "UTF-8")
}

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "gtf_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "doublet_utils.R"))

load_required_packages(c("Seurat", "dplyr", "tibble", "jsonlite", "scDblFinder", "SingleCellExperiment", "SummarizedExperiment"))

cfg <- get_single_script_config_02()
module_name <- "02b2_doublet"
prepare_dirs_02(cfg)
set.seed(cfg$random_seed)

manifest_02b1 <- read_manifest_local(cfg$module_02b1_manifest_path)
qc_filter_rds <- resolve_output_local(manifest_02b1, "qc_filter_object")
if (!file.exists(qc_filter_rds)) {
  stop(sprintf("缺少 QC filter 输入对象: %s", qc_filter_rds), call. = FALSE)
}

qc_sample_list <- readRDS(qc_filter_rds)
if (!is.list(qc_sample_list) || length(qc_sample_list) == 0) {
  stop("qc_filter_object 必须是非空 per-sample Seurat list。", call. = FALSE)
}

clean_list <- list()
summary_rows <- list()
diagnostic_rows <- list()
concordance_rows <- list()
cluster_risk_rows <- list()

for (sample_id in names(qc_sample_list)) {
  message("运行 doublet detection: ", sample_id)
  sample_obj <- maybe_join_layers(qc_sample_list[[sample_id]])
  if (!"qc_pass" %in% colnames(sample_obj@meta.data)) {
    stop(sprintf("样本 %s 缺少 qc_pass metadata，请先运行 02b1。", sample_id), call. = FALSE)
  }

  sample_meta <- sample_obj@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::mutate(
      sample_id = sample_id,
      barcode = if ("barcode_raw" %in% colnames(.)) barcode_raw else cell_id,
      qc_pass = as.logical(qc_pass)
    )
  qc_cells <- sample_meta$cell_id[!is.na(sample_meta$qc_pass) & sample_meta$qc_pass]
  qc_cells <- qc_cells[nzchar(qc_cells)]
  qc_obj <- if (length(qc_cells) > 0) subset(sample_obj, cells = qc_cells) else sample_obj[, 0]

  platform_resolved <- if ("platform" %in% colnames(sample_obj@meta.data)) {
    normalize_scalar_value(sample_obj$platform[1], "generic_mex")
  } else {
    "generic_mex"
  }

  diag_df <- data.frame(
    sample_id = sample_id,
    barcode = sample_meta$barcode,
    cell_id = sample_meta$cell_id,
    primary_call = NA_character_,
    secondary_call = NA_character_,
    primary_score = NA_real_,
    secondary_score = NA_real_,
    final_call = ifelse(sample_meta$qc_pass, "singlet", "qc_fail"),
    stringsAsFactors = FALSE
  )
  rownames(diag_df) <- diag_df$cell_id

  expected_rate <- estimate_expected_doublet_rate(ncol(qc_obj), cfg$doublet_rate_per_1k)
  expected_doublets <- estimate_expected_doublet_count(ncol(qc_obj), cfg$doublet_rate_per_1k)
  primary_status <- "skipped_low_cells"
  secondary_status <- if (cfg$doublet_secondary_enabled) "not_run" else "disabled"
  branch_status <- "skipped_empty"
  selected_pK <- NA_real_
  cluster_labels <- character(0)
  model_obj <- NULL

  if (ncol(qc_obj) > 0) {
    branch_model <- build_sample_branch_model(
      sample_obj = qc_obj,
      dims = cfg$doublet_dims,
      resolution = cfg$ambient_cluster_resolution,
      hvg_nfeatures = cfg$hvg_nfeatures,
      min_cells = cfg$doublet_min_cells
    )
    branch_status <- branch_model$status
    model_obj <- branch_model$object

    if (identical(branch_model$status, "built")) {
      cluster_labels <- as.character(model_obj$provisional_cluster)
      names(cluster_labels) <- colnames(model_obj)

      primary_result <- run_primary_scdblfinder(model_obj, platform_resolved, expected_rate)
      primary_status <- primary_result$status
      primary_cells <- intersect(names(primary_result$call), rownames(diag_df))
      diag_df[primary_cells, "primary_call"] <- unname(primary_result$call[primary_cells])
      diag_df[primary_cells, "primary_score"] <- unname(primary_result$score[primary_cells])

      if (length(branch_model$usable_dims) >= 5) {
        secondary_result <- run_secondary_doubletfinder(
          model_obj = model_obj,
          usable_dims = branch_model$usable_dims,
          expected_doublets = expected_doublets,
          enabled = cfg$doublet_secondary_enabled
        )
        secondary_status <- secondary_result$status
        selected_pK <- secondary_result$selected_pK
        secondary_cells <- intersect(names(secondary_result$call), rownames(diag_df))
        diag_df[secondary_cells, "secondary_call"] <- unname(secondary_result$call[secondary_cells])
        diag_df[secondary_cells, "secondary_score"] <- unname(secondary_result$score[secondary_cells])
      } else {
        secondary_status <- if (cfg$doublet_secondary_enabled) "skipped_low_dims" else "disabled"
      }
    } else {
      primary_status <- branch_model$status
      secondary_status <- if (cfg$doublet_secondary_enabled) branch_model$status else "disabled"
    }
  }

  qc_diag <- diag_df[diag_df$final_call != "qc_fail", , drop = FALSE]
  primary_available <- any(!is.na(qc_diag$primary_call))
  secondary_available <- any(!is.na(qc_diag$secondary_call))
  if (primary_available) {
    diag_df$final_call <- ifelse(
      diag_df$final_call == "qc_fail",
      "qc_fail",
      ifelse(vapply(diag_df$primary_call, is_primary_doublet, logical(1)), "doublet", "singlet")
    )
  } else if (secondary_available) {
    diag_df$final_call <- ifelse(
      diag_df$final_call == "qc_fail",
      "qc_fail",
      ifelse(vapply(diag_df$secondary_call, is_secondary_doublet, logical(1)), "doublet", "singlet")
    )
  }

  for (col_name in c("primary_call", "secondary_call", "primary_score", "secondary_score", "final_call")) {
    sample_obj@meta.data[rownames(diag_df), col_name] <- diag_df[rownames(diag_df), col_name]
  }
  retained_cells <- diag_df$cell_id[diag_df$final_call == "singlet"]
  retained_obj <- if (length(retained_cells) > 0) subset(sample_obj, cells = retained_cells) else sample_obj[, 0]
  if (ncol(retained_obj) > 0) {
    clean_list[[sample_id]] <- maybe_join_layers(retained_obj)
  }

  comparable <- diag_df[
    !is.na(diag_df$primary_call) & !is.na(diag_df$secondary_call) & diag_df$final_call != "qc_fail",
    ,
    drop = FALSE
  ]
  primary_doublet <- vapply(comparable$primary_call, is_primary_doublet, logical(1))
  secondary_doublet <- vapply(comparable$secondary_call, is_secondary_doublet, logical(1))
  concordance_rate <- if (nrow(comparable) > 0) mean(primary_doublet == secondary_doublet) else NA_real_
  discordant_pct <- if (nrow(comparable) > 0) mean(primary_doublet != secondary_doublet) else NA_real_
  primary_only <- if (nrow(comparable) > 0) sum(primary_doublet & !secondary_doublet) else 0L
  secondary_only <- if (nrow(comparable) > 0) sum(!primary_doublet & secondary_doublet) else 0L
  concordant <- if (nrow(comparable) > 0) sum(primary_doublet == secondary_doublet) else 0L

  qc_only <- diag_df[diag_df$final_call != "qc_fail", , drop = FALSE]
  primary_doublets <- sum(vapply(qc_only$primary_call, is_primary_doublet, logical(1)), na.rm = TRUE)
  secondary_doublets <- sum(vapply(qc_only$secondary_call, is_secondary_doublet, logical(1)), na.rm = TRUE)
  observed_primary_rate <- if (nrow(qc_only) > 0 && primary_available) primary_doublets / nrow(qc_only) else NA_real_

  summary_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    cells_in = nrow(qc_only),
    primary_caller = "scDblFinder",
    primary_doublets = primary_doublets,
    secondary_caller = if (cfg$doublet_secondary_enabled) cfg$doublet_secondary_caller else "",
    secondary_doublets = secondary_doublets,
    concordant = concordant,
    primary_only = primary_only,
    secondary_only = secondary_only,
    expected_rate = expected_rate,
    observed_primary_rate = observed_primary_rate,
    status = sprintf("branch=%s;primary=%s;secondary=%s", branch_status, primary_status, secondary_status),
    stringsAsFactors = FALSE
  )

  concordance_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    concordance_rate = concordance_rate,
    kappa = cohen_kappa_binary(primary_doublet, secondary_doublet),
    primary_secondary_discordant_pct = discordant_pct,
    stringsAsFactors = FALSE
  )

  if (!is.null(model_obj) && length(cluster_labels) > 0) {
    marker_df <- safe_find_all_markers(model_obj, top_n = 1L)
    top_marker <- setNames(rep("", length(unique(cluster_labels))), unique(cluster_labels))
    if (nrow(marker_df) > 0 && all(c("cluster", "gene") %in% colnames(marker_df))) {
      marker_df <- marker_df[!duplicated(marker_df$cluster), , drop = FALSE]
      top_marker[as.character(marker_df$cluster)] <- as.character(marker_df$gene)
    }
    primary_map <- setNames(vapply(diag_df$primary_call, is_primary_doublet, logical(1)), diag_df$cell_id)
    sample_baseline <- if (length(cluster_labels) > 0 && any(!is.na(primary_map[names(cluster_labels)]))) {
      mean(primary_map[names(cluster_labels)], na.rm = TRUE)
    } else {
      NA_real_
    }
    cluster_risk_rows[[sample_id]] <- dplyr::bind_rows(lapply(sort(unique(cluster_labels)), function(cluster_id) {
      cells <- names(cluster_labels)[cluster_labels == cluster_id]
      flags <- primary_map[cells]
      n_doublet <- sum(flags, na.rm = TRUE)
      doublet_fraction <- if (length(flags) > 0) mean(flags, na.rm = TRUE) else NA_real_
      data.frame(
        sample_id = sample_id,
        cluster = cluster_id,
        doublet_fraction = doublet_fraction,
        top_marker = normalize_scalar_value(top_marker[[cluster_id]]),
        risk_level = classify_cluster_doublet_risk(doublet_fraction, sample_baseline, n_doublet),
        stringsAsFactors = FALSE
      )
    }))
  }

  diagnostic_rows[[sample_id]] <- diag_df[, c("sample_id", "barcode", "primary_call", "secondary_call", "primary_score", "secondary_score", "final_call"), drop = FALSE]
}

if (length(clean_list) == 0) {
  stop("所有样本在 QC/doublet 后都没有保留细胞。", call. = FALSE)
}

post_qc_obj <- merge_named_objects(clean_list)
summary_df <- dplyr::bind_rows(summary_rows)
diagnostic_df <- dplyr::bind_rows(diagnostic_rows)
concordance_df <- dplyr::bind_rows(concordance_rows)
cluster_risk_df <- if (length(cluster_risk_rows) > 0) {
  dplyr::bind_rows(cluster_risk_rows)
} else {
  data.frame(
    sample_id = character(0),
    cluster = character(0),
    doublet_fraction = numeric(0),
    top_marker = character(0),
    risk_level = character(0),
    stringsAsFactors = FALSE
  )
}

post_qc_object_path <- file.path(cfg$checkpoint_dir, "02_after_qc_doublet.rds")
doublet_summary_path <- file.path(cfg$table_dir, "doublet_sample_summary.tsv")
doublet_cell_path <- file.path(cfg$table_dir, "doublet_cell_diagnostics.csv")
concordance_path <- file.path(cfg$table_dir, "doublet_concordance_summary.tsv")
cluster_risk_path <- file.path(cfg$table_dir, "doublet_cluster_risk.tsv")

saveRDS(post_qc_obj, post_qc_object_path)
write_tsv_local(summary_df, doublet_summary_path)
write_csv_local(diagnostic_df, doublet_cell_path)
write_tsv_local(concordance_df, concordance_path)
write_tsv_local(cluster_risk_df, cluster_risk_path)

write_manifest_local(
  manifest_path = cfg$module_02b2_manifest_path,
  new_outputs = list(
    post_qc_object = build_output_entry(post_qc_object_path, "rds", module_name, "merged Seurat object after QC and doublet removal", base_dir = cfg$project_root),
    doublet_sample_summary = build_output_entry(doublet_summary_path, "tsv", module_name, "one row per sample doublet summary", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    doublet_cell_diagnostics = build_output_entry(doublet_cell_path, "csv", module_name, "one row per cell doublet diagnostics", base_dir = cfg$project_root, schema = infer_schema_from_df(diagnostic_df)),
    doublet_concordance_summary = build_output_entry(concordance_path, "tsv", module_name, "one row per sample primary/secondary concordance", base_dir = cfg$project_root, schema = infer_schema_from_df(concordance_df)),
    doublet_cluster_risk = build_output_entry(cluster_risk_path, "tsv", module_name, "one row per provisional cluster doublet risk", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_risk_df))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    manifest_02b1 = cfg$module_02b1_manifest_path,
    qc_filter_object = qc_filter_rds
  ),
  version = cfg$module_version,
  depends_on = list(
    module_02b1 = cfg$module_02b1_manifest_path
  )
)

message("02b2 完成。post-QC 对象: ", post_qc_object_path)
