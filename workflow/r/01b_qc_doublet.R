source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

input_rds <- file.path(cfg$checkpoint_dir, "01_after_ambient_branch.rds")
if (!file.exists(input_rds)) {
  stop("缺少输入检查点: 01_after_ambient_branch.rds", call. = FALSE)
}

if (!requireNamespace("scDblFinder", quietly = TRUE) ||
    !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  stop("缺少 scDblFinder / SingleCellExperiment 依赖。请先运行 workflow/01_install_envs.sh 更新主流程环境。", call. = FALSE)
}

secondary_pkg_available <- requireNamespace("DoubletFinder", quietly = TRUE)

input_obj <- readRDS(input_rds)
input_obj <- maybe_join_layers(input_obj)
threshold_overrides <- read_qc_threshold_overrides(cfg)
inventory_df <- read_input_inventory(cfg)

run_secondary_doubletfinder <- function(model_obj, usable_dims, expected_doublets) {
  if (!secondary_pkg_available || !cfg$doublet_secondary_enabled) {
    return(list(
      status = if (!cfg$doublet_secondary_enabled) "disabled" else "package_missing",
      class = rep(NA_character_, ncol(model_obj)),
      pANN = rep(NA_real_, ncol(model_obj)),
      selected_pK = NA_real_
    ))
  }

  doublet_pkg <- asNamespace("DoubletFinder")
  sweep_fun <- if (exists("paramSweep_v3", envir = doublet_pkg, mode = "function")) {
    get("paramSweep_v3", envir = doublet_pkg)
  } else {
    get("paramSweep", envir = doublet_pkg)
  }
  summarize_fun <- get("summarizeSweep", envir = doublet_pkg)
  find_pk_fun <- get("find.pK", envir = doublet_pkg)
  df_fun <- if (exists("doubletFinder_v3", envir = doublet_pkg, mode = "function")) {
    get("doubletFinder_v3", envir = doublet_pkg)
  } else {
    get("doubletFinder", envir = doublet_pkg)
  }

  best_pk <- NA_real_
  status <- "completed"
  sweep_res <- tryCatch(
    sweep_fun(model_obj, PCs = usable_dims, sct = FALSE),
    error = function(e) NULL
  )
  if (!is.null(sweep_res)) {
    sweep_stats <- tryCatch(summarize_fun(sweep_res, GT = FALSE), error = function(e) NULL)
    bcmvn <- tryCatch(find_pk_fun(sweep_stats), error = function(e) NULL)
    if (!is.null(bcmvn) && nrow(bcmvn) > 0) {
      best_pk <- suppressWarnings(as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])))
    }
  }
  if (!is.finite(best_pk)) {
    best_pk <- 0.09
    status <- "completed_fallback_pk"
  }

  df_obj <- tryCatch(
    df_fun(
      model_obj,
      PCs = usable_dims,
      pN = 0.25,
      pK = best_pk,
      nExp = expected_doublets,
      sct = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(df_obj)) {
    return(list(
      status = "failed",
      class = rep(NA_character_, ncol(model_obj)),
      pANN = rep(NA_real_, ncol(model_obj)),
      selected_pK = best_pk
    ))
  }

  class_cols <- grep("DF.classifications", colnames(df_obj@meta.data), value = TRUE)
  pann_cols <- grep("^pANN", colnames(df_obj@meta.data), value = TRUE)
  class_vec <- if (length(class_cols) > 0) as.character(df_obj@meta.data[[tail(class_cols, 1)]]) else rep(NA_character_, ncol(model_obj))
  pann_vec <- if (length(pann_cols) > 0) suppressWarnings(as.numeric(df_obj@meta.data[[tail(pann_cols, 1)]])) else rep(NA_real_, ncol(model_obj))
  names(class_vec) <- colnames(model_obj)
  names(pann_vec) <- colnames(model_obj)

  list(
    status = status,
    class = class_vec,
    pANN = pann_vec,
    selected_pK = best_pk
  )
}

split_objs <- split_object_by_sample(input_obj)
clean_list <- list()
summary_rows <- list()
diagnostic_rows <- list()
concordance_rows <- list()
cluster_risk_rows <- list()

for (sample_id in names(split_objs)) {
  message("运行 sample-wise QC + doublet branch: ", sample_id)
  sample_obj <- maybe_join_layers(split_objs[[sample_id]])
  inventory_row <- resolve_inventory_row(cfg, sample_id, inventory_df = inventory_df)
  platform_resolved <- if (!is.null(inventory_row) && "platform_resolved" %in% colnames(inventory_row)) {
    normalize_scalar_value(inventory_row$platform_resolved[1], "generic_mex")
  } else {
    "generic_mex"
  }
  thresholds <- get_sample_qc_thresholds(cfg, sample_id, threshold_overrides)

  meta_df <- sample_obj@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    mutate(
      sample_id = sample_id,
      qc_pass = nFeature_RNA >= thresholds$qc_min_nfeature &
        nCount_RNA >= thresholds$qc_min_ncount &
        log10GenesPerUMI >= thresholds$qc_min_log10umi &
        percent.mito <= thresholds$qc_max_mito_pct,
      qc_fail_reason = ifelse(qc_pass, "", "basic_qc"),
      provisional_cluster = NA_character_,
      provisional_umap_1 = NA_real_,
      provisional_umap_2 = NA_real_,
      scDblFinder.score = NA_real_,
      scDblFinder.class = NA_character_,
      DoubletFinder.class = NA_character_,
      DoubletFinder.pANN = NA_real_,
      final_doublet_keep = FALSE,
      final_filter_caller = "none"
    )

  qc_cells <- meta_df$cell_id[meta_df$qc_pass]
  qc_obj <- subset(sample_obj, cells = qc_cells)

  primary_status <- "skipped_low_cells"
  secondary_status <- if (cfg$doublet_secondary_enabled) "not_run" else "disabled"
  expected_rate <- estimate_expected_doublet_rate(ncol(qc_obj), cfg$doublet_rate_per_1k)
  expected_doublets <- estimate_expected_doublet_count(ncol(qc_obj), cfg$doublet_rate_per_1k)
  primary_detected_doublets <- 0L
  secondary_detected_doublets <- 0L
  concordant_doublets <- 0L
  discordant_calls <- 0L
  final_filter_caller <- "none"
  selected_pK <- NA_real_

  if (ncol(qc_obj) > 0) {
    branch_model <- build_sample_branch_model(
      sample_obj = qc_obj,
      dims = cfg$doublet_dims,
      resolution = cfg$ambient_cluster_resolution,
      hvg_nfeatures = cfg$hvg_nfeatures,
      min_cells = cfg$doublet_min_cells
    )
    model_obj <- branch_model$object

    if (identical(branch_model$status, "built")) {
      qc_model_meta <- model_obj@meta.data %>%
        tibble::rownames_to_column("cell_id") %>%
        transmute(
          cell_id,
          provisional_cluster = as.character(provisional_cluster),
          provisional_umap_1 = provisional_umap_1,
          provisional_umap_2 = provisional_umap_2
        )
      meta_df <- meta_df %>%
        left_join(qc_model_meta, by = "cell_id", suffix = c("", ".new")) %>%
        mutate(
          provisional_cluster = ifelse(!is.na(provisional_cluster.new), provisional_cluster.new, provisional_cluster),
          provisional_umap_1 = ifelse(!is.na(provisional_umap_1.new), provisional_umap_1.new, provisional_umap_1),
          provisional_umap_2 = ifelse(!is.na(provisional_umap_2.new), provisional_umap_2.new, provisional_umap_2)
        ) %>%
        select(-ends_with(".new"))

      counts_mat <- get_assay_matrix(model_obj, assay = "RNA", type = "counts")
      sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_mat))
      SummarizedExperiment::colData(sce)$cluster <- as.character(model_obj$provisional_cluster)

      primary_args <- list(sce = sce, clusters = SummarizedExperiment::colData(sce)$cluster, verbose = FALSE)
      if (!identical(platform_resolved, "10x_cellranger") && is.finite(expected_rate)) {
        primary_args$dbr <- expected_rate
      }

      sce <- tryCatch(
        do.call(scDblFinder::scDblFinder, primary_args),
        error = function(e) {
          message(sprintf("scDblFinder failed for %s: %s", sample_id, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(sce)) {
        primary_status <- if (!identical(platform_resolved, "10x_cellranger") && is.finite(expected_rate)) "completed_manual_rate" else "completed_auto_rate"
        primary_df <- data.frame(
          cell_id = colnames(model_obj),
          scDblFinder.score = suppressWarnings(as.numeric(SummarizedExperiment::colData(sce)$scDblFinder.score)),
          scDblFinder.class = as.character(SummarizedExperiment::colData(sce)$scDblFinder.class),
          stringsAsFactors = FALSE
        )
        meta_df <- meta_df %>%
          left_join(primary_df, by = "cell_id", suffix = c("", ".new")) %>%
          mutate(
            scDblFinder.score = ifelse(!is.na(scDblFinder.score.new), scDblFinder.score.new, scDblFinder.score),
            scDblFinder.class = ifelse(!is.na(scDblFinder.class.new), scDblFinder.class.new, scDblFinder.class)
          ) %>%
          select(-ends_with(".new"))
        primary_detected_doublets <- sum(tolower(meta_df$scDblFinder.class) == "doublet", na.rm = TRUE)
      } else {
        primary_status <- "failed"
      }

      if (length(branch_model$usable_dims) >= 5) {
        secondary_result <- run_secondary_doubletfinder(
          model_obj = model_obj,
          usable_dims = branch_model$usable_dims,
          expected_doublets = expected_doublets
        )
        secondary_status <- secondary_result$status
        selected_pK <- secondary_result$selected_pK
        secondary_df <- data.frame(
          cell_id = names(secondary_result$class),
          DoubletFinder.class = unname(secondary_result$class),
          DoubletFinder.pANN = unname(secondary_result$pANN),
          stringsAsFactors = FALSE
        )
        meta_df <- meta_df %>%
          left_join(secondary_df, by = "cell_id", suffix = c("", ".new")) %>%
          mutate(
            DoubletFinder.class = ifelse(!is.na(DoubletFinder.class.new), DoubletFinder.class.new, DoubletFinder.class),
            DoubletFinder.pANN = ifelse(!is.na(DoubletFinder.pANN.new), DoubletFinder.pANN.new, DoubletFinder.pANN)
          ) %>%
          select(-ends_with(".new"))
        secondary_detected_doublets <- sum(meta_df$DoubletFinder.class == "Doublet", na.rm = TRUE)
      } else {
        secondary_status <- if (cfg$doublet_secondary_enabled) "skipped_low_dims" else "disabled"
      }
    } else {
      primary_status <- branch_model$status
      secondary_status <- if (cfg$doublet_secondary_enabled) branch_model$status else "disabled"
    }
  }

  primary_available <- any(!is.na(meta_df$scDblFinder.class[meta_df$qc_pass]))
  secondary_available <- any(!is.na(meta_df$DoubletFinder.class[meta_df$qc_pass]))

  comparable_df <- meta_df %>%
    filter(qc_pass, !is.na(scDblFinder.class), !is.na(DoubletFinder.class)) %>%
    mutate(
      primary_doublet = tolower(scDblFinder.class) == "doublet",
      secondary_doublet = DoubletFinder.class == "Doublet"
    )
  if (nrow(comparable_df) > 0) {
    concordant_doublets <- sum(comparable_df$primary_doublet & comparable_df$secondary_doublet)
    discordant_calls <- sum(xor(comparable_df$primary_doublet, comparable_df$secondary_doublet))
  }

  if (primary_available) {
    meta_df$final_doublet_keep <- ifelse(meta_df$qc_pass, tolower(meta_df$scDblFinder.class) != "doublet", FALSE)
    meta_df$final_filter_caller <- ifelse(meta_df$qc_pass, "scDblFinder", "qc")
    final_filter_caller <- "scDblFinder"
  } else if (secondary_available) {
    meta_df$final_doublet_keep <- ifelse(meta_df$qc_pass, meta_df$DoubletFinder.class != "Doublet", FALSE)
    meta_df$final_filter_caller <- ifelse(meta_df$qc_pass, "DoubletFinder", "qc")
    final_filter_caller <- "DoubletFinder"
  } else {
    meta_df$final_doublet_keep <- meta_df$qc_pass
    meta_df$final_filter_caller <- ifelse(meta_df$qc_pass, "none", "qc")
    final_filter_caller <- "none"
  }

  singlet_cells <- meta_df$cell_id[meta_df$final_doublet_keep]
  singlet_obj <- if (length(singlet_cells) > 0) subset(sample_obj, cells = singlet_cells) else sample_obj[, 0]
  if (ncol(singlet_obj) > 0) {
    updated_meta <- meta_df[match(colnames(singlet_obj), meta_df$cell_id), , drop = FALSE]
    rownames(updated_meta) <- updated_meta$cell_id
    for (col_name in setdiff(colnames(updated_meta), "cell_id")) {
      singlet_obj@meta.data[[col_name]] <- updated_meta[[col_name]]
    }
  }
  clean_list[[sample_id]] <- maybe_join_layers(singlet_obj)

  sample_baseline_rate <- if (sum(meta_df$qc_pass) > 0) {
    mean(tolower(meta_df$scDblFinder.class[meta_df$qc_pass]) == "doublet", na.rm = TRUE)
  } else {
    NA_real_
  }
  cluster_risk_df <- meta_df %>%
    filter(qc_pass, !is.na(provisional_cluster)) %>%
    group_by(sample_id, provisional_cluster) %>%
    summarise(
      n_cells = dplyr::n(),
      n_primary_doublet = sum(tolower(scDblFinder.class) == "doublet", na.rm = TRUE),
      primary_doublet_rate = mean(tolower(scDblFinder.class) == "doublet", na.rm = TRUE),
      sample_baseline_rate = sample_baseline_rate,
      enrichment_vs_sample_baseline = ifelse(is.finite(sample_baseline_rate) && sample_baseline_rate > 0, primary_doublet_rate / sample_baseline_rate, NA_real_),
      high_risk_cluster = ifelse(is.finite(enrichment_vs_sample_baseline) && enrichment_vs_sample_baseline >= 1.5 && n_primary_doublet >= 2, "yes", "no"),
      .groups = "drop"
    )
  cluster_risk_rows[[sample_id]] <- cluster_risk_df

  summary_rows[[sample_id]] <- data.frame(
    sample = sample_id,
    cells_raw = ncol(sample_obj),
    cells_after_qc = sum(meta_df$qc_pass),
    singlets_after_doublet = sum(meta_df$final_doublet_keep),
    doublets_removed = sum(meta_df$qc_pass) - sum(meta_df$final_doublet_keep),
    doublet_expected_rate = expected_rate,
    expected_doublets = expected_doublets,
    primary_caller = cfg$doublet_primary_caller,
    primary_status = primary_status,
    primary_detected_doublets = primary_detected_doublets,
    primary_detected_rate = ifelse(sum(meta_df$qc_pass) > 0, primary_detected_doublets / sum(meta_df$qc_pass), NA_real_),
    secondary_caller = cfg$doublet_secondary_caller,
    secondary_enabled = cfg$doublet_secondary_enabled,
    secondary_status = secondary_status,
    secondary_detected_doublets = secondary_detected_doublets,
    secondary_detected_rate = ifelse(sum(meta_df$qc_pass) > 0, secondary_detected_doublets / sum(meta_df$qc_pass), NA_real_),
    concordant_doublets = concordant_doublets,
    discordant_calls = discordant_calls,
    selected_pK = selected_pK,
    final_filter_caller = final_filter_caller,
    stringsAsFactors = FALSE
  )

  concordance_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    comparable_cells = nrow(comparable_df),
    both_doublet = sum(tolower(comparable_df$scDblFinder.class) == "doublet" & comparable_df$DoubletFinder.class == "Doublet"),
    primary_only_doublet = sum(tolower(comparable_df$scDblFinder.class) == "doublet" & comparable_df$DoubletFinder.class != "Doublet"),
    secondary_only_doublet = sum(tolower(comparable_df$scDblFinder.class) != "doublet" & comparable_df$DoubletFinder.class == "Doublet"),
    both_singlet = sum(tolower(comparable_df$scDblFinder.class) != "doublet" & comparable_df$DoubletFinder.class != "Doublet"),
    stringsAsFactors = FALSE
  )

  diagnostic_rows[[sample_id]] <- meta_df %>%
    mutate(
      doublet_expected_rate = expected_rate,
      final_filter_keep = final_doublet_keep
    )
}

nonempty_list <- clean_list[vapply(clean_list, function(obj) as.integer(ncol(obj)), integer(1)) > 0]
if (length(nonempty_list) == 0) {
  stop("所有样本在 QC/doublet 后都没有保留细胞。", call. = FALSE)
}

final_obj <- merge_named_objects(nonempty_list)
summary_df <- bind_rows(summary_rows)
diagnostic_df <- bind_rows(diagnostic_rows)
concordance_df <- bind_rows(concordance_rows)
cluster_risk_df <- bind_rows(cluster_risk_rows)

saveRDS(final_obj, file.path(cfg$checkpoint_dir, "01_after_qc_doublet.rds"))
write.csv(summary_df, file.path(cfg$table_dir, "qc_doublet_summary.csv"), row.names = FALSE)

diag_path <- file.path(cfg$table_dir, "qc_doublet_cell_level.csv.gz")
diag_conn <- gzfile(diag_path, open = "wt")
write.csv(diagnostic_df, diag_conn, row.names = FALSE)
close(diag_conn)

write_tsv(summary_df, file.path(cfg$table_dir, "doublet_sample_summary.tsv"))
write_tsv(concordance_df, file.path(cfg$table_dir, "doublet_concordance_summary.tsv"))
write_tsv(cluster_risk_df, file.path(cfg$table_dir, "doublet_cluster_risk.tsv"))

message("已保存 01_after_qc_doublet.rds、qc_doublet_summary.csv、qc_doublet_cell_level.csv.gz 以及 doublet 汇总表。")
