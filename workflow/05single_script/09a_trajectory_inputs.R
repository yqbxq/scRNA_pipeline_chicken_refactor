#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source(file.path(.script_dir, "helpers", "load_helpers_09.R"), encoding = "UTF-8")

load_required_packages(c("Seurat", "dplyr", "jsonlite", "Matrix"))

cfg <- get_single_script_config_09()
module_name <- "09a_trajectory_inputs"
prepare_dirs_09(cfg)
set.seed(cfg$random_seed)

empty_df_09 <- function(cols) {
  as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols), stringsAsFactors = FALSE)
}

trajectory_input_index_cols_09 <- c(
  "pair_id", "source_question_id", "layer_scope", "method", "tools_to_run",
  "split_mode", "split_var", "split_value", "root_group", "terminal_group",
  "coarse_label_var", "fine_label_var", "object_rds", "input_rds",
  "cell_n_before", "outlier_n", "cell_n_after", "retention_fraction",
  "root_group_present", "terminal_group_present", "cc_status", "pca_npcs",
  "umap_reduction", "status", "reason", "produced_at"
)

outlier_summary_cols_09 <- c(
  "pair_id", "split_value", "policy", "cell_n_before", "outlier_n",
  "cell_n_after", "retention_fraction", "nfeature_low_cutoff",
  "nfeature_high_cutoff", "percent_mt_high_cutoff",
  "neighbor_distance_z_cutoff", "filter_status", "reason"
)

outlier_cell_cols_09 <- c(
  "pair_id", "split_value", "cell_id", "nFeature_RNA", "nCount_RNA",
  "percent.mt", "neighbor_distance_z", "retained", "outlier_reason"
)

cc_summary_cols_09 <- c(
  "pair_id", "split_value", "policy", "status", "reason",
  "s_features_matched", "g2m_features_matched", "mean_s_score",
  "mean_g2m_score", "phase_counts"
)

pca_summary_cols_09 <- c(
  "pair_id", "split_value", "pc", "stdev", "variance_fraction",
  "cumulative_variance_fraction", "n_cells", "n_features", "dims_used"
)

split_index_cols_09 <- c(
  "pair_id", "split_var", "split_value", "label_var", "label_value",
  "n_cells", "fraction", "root_group", "terminal_group", "is_root_group",
  "is_terminal_group"
)

trajectory_default_assay_09 <- function(seu) {
  assays <- tryCatch(names(seu@assays), error = function(e) character(0))
  if ("RNA" %in% assays) {
    return("RNA")
  }
  Seurat::DefaultAssay(seu)
}

trajectory_reset_reductions_09 <- function(seu) {
  if (exists("maybe_join_layers", mode = "function")) {
    seu <- maybe_join_layers(seu)
  }
  if ("reductions" %in% slotNames(seu)) {
    seu@reductions <- list()
  }
  if ("graphs" %in% slotNames(seu)) {
    seu@graphs <- list()
  }
  if ("neighbors" %in% slotNames(seu)) {
    seu@neighbors <- list()
  }
  if ("commands" %in% slotNames(seu)) {
    seu@commands <- list()
  }
  VariableFeatures(seu) <- character(0)
  seu
}

read_trajectory_pairs_09 <- function(cfg) {
  pairs <- read_tsv_optional(cfg$trajectory_pairs_sheet)
  if (nrow(pairs) == 0) {
    stop(sprintf("trajectory_pairs.tsv is empty or missing: %s", cfg$trajectory_pairs_sheet), call. = FALSE)
  }
  required <- c(
    "trajectory_id", "source_question_id", "layer_scope", "root_group",
    "terminal_group", "condition_split_var", "condition_split_values",
    "method", "tools_to_run", "outlier_qc_policy", "regress_cell_cycle",
    "coarse_label_var", "fine_label_var", "split_mode", "enabled"
  )
  require_columns_local(pairs, required, "trajectory_pairs.tsv", cfg$trajectory_pairs_sheet)
  pairs$enabled <- normalize_flag(pairs$enabled, "yes")
  pairs$method <- tolower(normalize_flag(pairs$method, "trajectory"))
  pairs[pairs$enabled != "no" & pairs$method == "trajectory", , drop = FALSE]
}

resolve_trajectory_layer_09 <- function(cfg, layer_scope) {
  layer_scope <- normalize_scalar_value(layer_scope, cfg$panorama_layer_id)
  layers <- communication_layer_status_07(cfg)
  if (nrow(layers) == 0) {
    stop("no usable 03d/04b annotation objects found; complete annotation/subcluster outputs first", call. = FALSE)
  }

  aliases <- unique(c(layer_scope, safe_id_09(layer_scope)))
  if (tolower(layer_scope) %in% c("panorama", tolower(cfg$panorama_layer_id))) {
    aliases <- unique(c(aliases, cfg$panorama_layer_id, "panorama"))
  }

  hit <- layers[layers$layer_id %in% aliases, , drop = FALSE]
  if (nrow(hit) == 0) {
    hit <- layers[tolower(layers$layer_id) %in% tolower(aliases), , drop = FALSE]
  }
  if (nrow(hit) == 0) {
    stop(
      sprintf(
        "layer_scope=%s did not match any usable annotation object; available layers: %s",
        layer_scope,
        paste(layers$layer_id, collapse = ",")
      ),
      call. = FALSE
    )
  }
  hit[1, , drop = FALSE]
}

trajectory_existing_umap_09 <- function(seu) {
  reductions <- tryCatch(names(seu@reductions), error = function(e) character(0))
  if (length(reductions) == 0) {
    return(NULL)
  }
  hit <- reductions[grepl("umap", reductions, ignore.case = TRUE)]
  if (length(hit) == 0) {
    return(NULL)
  }
  emb <- tryCatch(Seurat::Embeddings(seu, reduction = hit[[1]]), error = function(e) NULL)
  if (is.null(emb) || ncol(emb) < 2) {
    return(NULL)
  }
  emb[, 1:2, drop = FALSE]
}

nearest_neighbor_z_09 <- function(seu) {
  coords <- trajectory_existing_umap_09(seu)
  cells <- colnames(seu)
  z <- rep(NA_real_, length(cells))
  names(z) <- cells
  if (is.null(coords) || nrow(coords) < 3) {
    return(z)
  }
  coords <- coords[intersect(rownames(coords), cells), , drop = FALSE]
  if (nrow(coords) < 3) {
    return(z)
  }

  k <- min(15L, nrow(coords) - 1L)
  nn_dist <- NULL
  if (requireNamespace("FNN", quietly = TRUE)) {
    nn <- FNN::get.knn(coords, k = k)
    nn_dist <- rowMeans(nn$nn.dist)
  } else if (nrow(coords) <= 5000L) {
    d <- as.matrix(stats::dist(coords))
    diag(d) <- NA_real_
    nn_dist <- apply(d, 1, function(x) mean(sort(x, na.last = NA)[seq_len(k)]))
  }

  if (is.null(nn_dist)) {
    return(z)
  }
  scaled <- as.numeric(scale(nn_dist))
  z[rownames(coords)] <- scaled
  z
}

metric_from_meta_or_counts_09 <- function(seu, assay) {
  cells <- colnames(seu)
  meta <- seu@meta.data
  counts <- tryCatch(
    get_assay_matrix(seu, assay = assay, type = "counts"),
    error = function(e) get_assay_matrix(seu, assay = assay, type = "data")
  )

  nfeature <- if ("nFeature_RNA" %in% colnames(meta)) {
    suppressWarnings(as.numeric(meta[["nFeature_RNA"]]))
  } else {
    as.numeric(Matrix::colSums(counts > 0))
  }
  ncount <- if ("nCount_RNA" %in% colnames(meta)) {
    suppressWarnings(as.numeric(meta[["nCount_RNA"]]))
  } else {
    as.numeric(Matrix::colSums(counts))
  }

  mt_col <- intersect(c("percent.mt", "percent_mito", "percent.mito", "mito_pct"), colnames(meta))
  percent_mt <- if (length(mt_col) > 0) {
    suppressWarnings(as.numeric(meta[[mt_col[[1]]]]))
  } else {
    gene_names <- rownames(counts)
    mt_genes <- grepl("^MT-|^mt-|^Mt-|^MTRNR|^ND[0-9]|^COX[0-9]|^CYTB$|^ATP[68]$", gene_names)
    if (any(mt_genes)) {
      denom <- Matrix::colSums(counts)
      100 * as.numeric(Matrix::colSums(counts[mt_genes, , drop = FALSE])) / pmax(as.numeric(denom), 1)
    } else {
      rep(NA_real_, length(cells))
    }
  }

  data.frame(
    cell_id = cells,
    nFeature_RNA = as.numeric(nfeature),
    nCount_RNA = as.numeric(ncount),
    percent.mt = as.numeric(percent_mt),
    neighbor_distance_z = nearest_neighbor_z_09(seu)[cells],
    stringsAsFactors = FALSE
  )
}

outlier_thresholds_09 <- function(cfg, policy) {
  if (identical(policy, "strict")) {
    return(list(
      low_prob = cfg$trajectory_outlier_strict_low_prob,
      high_prob = cfg$trajectory_outlier_strict_high_prob,
      neighbor_cutoff = cfg$trajectory_outlier_strict_neighbor_cutoff
    ))
  }
  list(
    low_prob = cfg$trajectory_outlier_low_prob,
    high_prob = cfg$trajectory_outlier_high_prob,
    neighbor_cutoff = cfg$trajectory_outlier_neighbor_cutoff
  )
}

add_outlier_flags_09 <- function(qc_df, policy, cfg) {
  policy <- tolower(normalize_scalar_value(policy, "standard"))
  if (!policy %in% c("off", "none", "standard", "strict")) {
    policy <- "standard"
  }
  if (policy %in% c("off", "none")) {
    qc_df$retained <- TRUE
    qc_df$outlier_reason <- ""
    return(list(
      cells = qc_df,
      summary = list(
        nfeature_low_cutoff = NA_real_,
        nfeature_high_cutoff = NA_real_,
        percent_mt_high_cutoff = NA_real_,
        neighbor_distance_z_cutoff = NA_real_,
        filter_status = "off",
        reason = ""
      )
    ))
  }

  thresholds <- outlier_thresholds_09(cfg, policy)
  low_prob <- thresholds$low_prob
  high_prob <- thresholds$high_prob
  neighbor_cutoff <- thresholds$neighbor_cutoff
  nfeature_low <- safe_quantile(qc_df$nFeature_RNA, low_prob)
  nfeature_high <- safe_quantile(qc_df$nFeature_RNA, high_prob)
  mt_high <- safe_quantile(qc_df$percent.mt, high_prob)

  reasons <- vector("list", nrow(qc_df))
  flag <- rep(FALSE, nrow(qc_df))
  add_reason <- function(flag, reasons, mask, label) {
    mask[is.na(mask)] <- FALSE
    flag <- flag | mask
    idx <- which(mask)
    if (length(idx) > 0) {
      for (i in idx) {
        reasons[[i]] <- c(reasons[[i]], label)
      }
    }
    list(flag = flag, reasons = reasons)
  }

  if (is.finite(nfeature_low)) {
    state <- add_reason(flag, reasons, qc_df$nFeature_RNA < nfeature_low, "nFeature_low")
    flag <- state$flag
    reasons <- state$reasons
  }
  if (is.finite(nfeature_high)) {
    state <- add_reason(flag, reasons, qc_df$nFeature_RNA > nfeature_high, "nFeature_high")
    flag <- state$flag
    reasons <- state$reasons
  }
  if (is.finite(mt_high)) {
    state <- add_reason(flag, reasons, qc_df$percent.mt > mt_high, "percent_mt_high")
    flag <- state$flag
    reasons <- state$reasons
  }
  state <- add_reason(flag, reasons, qc_df$neighbor_distance_z > neighbor_cutoff, "neighbor_distance_high")
  flag <- state$flag
  reasons <- state$reasons

  retained <- !flag
  filter_status <- "applied"
  reason <- ""
  if (sum(retained) < 3L) {
    retained <- rep(TRUE, nrow(qc_df))
    filter_status <- "disabled_min_cells"
    reason <- "outlier filtering would leave fewer than 3 cells"
  }

  qc_df$retained <- retained
  qc_df$outlier_reason <- vapply(reasons, paste, character(1), collapse = ";")
  list(
    cells = qc_df,
    summary = list(
      nfeature_low_cutoff = nfeature_low,
      nfeature_high_cutoff = nfeature_high,
      percent_mt_high_cutoff = mt_high,
      neighbor_distance_z_cutoff = neighbor_cutoff,
      filter_status = filter_status,
      reason = reason
    )
  )
}

resolve_cell_cycle_genes_09 <- function(cfg) {
  if (!file.exists(cfg$module_00c_manifest_path)) {
    return(list(status = "missing_manifest", reason = cfg$module_00c_manifest_path, genes = NULL))
  }
  manifest <- read_manifest_local(cfg$module_00c_manifest_path)
  cc_path <- tryCatch(
    resolve_output_local(manifest, "cc_genes"),
    error = function(e) tryCatch(resolve_output_local(manifest, "chicken_cc_genes"), error = function(e2) "")
  )
  if (!nzchar(cc_path) || !file.exists(cc_path)) {
    return(list(status = "missing_cc_genes", reason = cc_path, genes = NULL))
  }
  list(status = "ok", reason = "", genes = readRDS(cc_path))
}

match_features_case_insensitive_09 <- function(features, genes) {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes)]
  exact <- intersect(genes, features)
  if (length(exact) > 0) {
    return(exact)
  }
  feature_map <- setNames(features, toupper(features))
  matched <- unname(feature_map[toupper(genes)])
  unique(matched[!is.na(matched) & nzchar(matched)])
}

apply_cell_cycle_09 <- function(seu, cfg, policy) {
  policy <- tolower(normalize_scalar_value(policy, "auto"))
  if (policy %in% c("no", "none", "off", "false")) {
    return(list(object = seu, vars_to_regress = character(0), summary = data.frame(
      policy = policy,
      status = "skipped_by_policy",
      reason = "",
      s_features_matched = 0L,
      g2m_features_matched = 0L,
      mean_s_score = NA_real_,
      mean_g2m_score = NA_real_,
      phase_counts = "",
      stringsAsFactors = FALSE
    )))
  }

  cc <- resolve_cell_cycle_genes_09(cfg)
  if (!identical(cc$status, "ok")) {
    if (policy %in% c("yes", "required", "true")) {
      stop(sprintf("cell-cycle regression required but %s: %s", cc$status, cc$reason), call. = FALSE)
    }
    return(list(object = seu, vars_to_regress = character(0), summary = data.frame(
      policy = policy,
      status = paste0("skipped_", cc$status),
      reason = cc$reason,
      s_features_matched = 0L,
      g2m_features_matched = 0L,
      mean_s_score = NA_real_,
      mean_g2m_score = NA_real_,
      phase_counts = "",
      stringsAsFactors = FALSE
    )))
  }

  s_features <- match_features_case_insensitive_09(rownames(seu), cc$genes$s.genes %||% character(0))
  g2m_features <- match_features_case_insensitive_09(rownames(seu), cc$genes$g2m.genes %||% character(0))
  if (length(s_features) == 0 || length(g2m_features) == 0) {
    if (policy %in% c("yes", "required", "true")) {
      stop("cell-cycle regression required but S/G2M genes did not match object rownames.", call. = FALSE)
    }
    return(list(object = seu, vars_to_regress = character(0), summary = data.frame(
      policy = policy,
      status = "skipped_no_feature_match",
      reason = "S/G2M genes did not match object rownames",
      s_features_matched = length(s_features),
      g2m_features_matched = length(g2m_features),
      mean_s_score = NA_real_,
      mean_g2m_score = NA_real_,
      phase_counts = "",
      stringsAsFactors = FALSE
    )))
  }

  seu <- Seurat::CellCycleScoring(
    seu,
    s.features = s_features,
    g2m.features = g2m_features,
    set.ident = FALSE
  )
  phase_counts <- if ("Phase" %in% colnames(seu@meta.data)) {
    paste(sprintf("%s=%s", names(table(seu$Phase)), as.integer(table(seu$Phase))), collapse = ";")
  } else {
    ""
  }
  list(object = seu, vars_to_regress = c("S.Score", "G2M.Score"), summary = data.frame(
    policy = policy,
    status = "applied",
    reason = "",
    s_features_matched = length(s_features),
    g2m_features_matched = length(g2m_features),
    mean_s_score = mean(seu$S.Score, na.rm = TRUE),
    mean_g2m_score = mean(seu$G2M.Score, na.rm = TRUE),
    phase_counts = phase_counts,
    stringsAsFactors = FALSE
  ))
}

split_values_for_pair_09 <- function(seu, pair_row) {
  split_mode <- tolower(normalize_scalar_value(pair_row$split_mode[[1]], "auto"))
  split_var <- normalize_scalar_value(pair_row$condition_split_var[[1]])
  if (split_mode %in% c("pooled", "force_pooled") || !nzchar(split_var)) {
    return(list(split_var = "", values = ""))
  }
  if (!split_var %in% colnames(seu@meta.data)) {
    stop(sprintf("condition_split_var=%s is absent from object metadata", split_var), call. = FALSE)
  }
  values <- split_csv_local(pair_row$condition_split_values[[1]])
  if (length(values) == 0) {
    values <- sort(unique(trimws(as.character(seu@meta.data[[split_var]]))))
    values <- values[nzchar(values)]
  }
  if (length(values) == 0) {
    return(list(split_var = "", values = ""))
  }
  list(split_var = split_var, values = values)
}

label_presence_09 <- function(labels, token) {
  token <- normalize_scalar_value(token)
  if (!nzchar(token) || token %in% c("*", "auto")) {
    return("not_applicable")
  }
  if (token %in% labels) "yes" else "no"
}

build_split_index_09 <- function(seu, pair_row, split_var, split_value) {
  pair_id <- normalize_scalar_value(pair_row$trajectory_id[[1]])
  coarse <- normalize_scalar_value(pair_row$coarse_label_var[[1]])
  fine <- normalize_scalar_value(pair_row$fine_label_var[[1]])
  root <- normalize_scalar_value(pair_row$root_group[[1]])
  terminal <- normalize_scalar_value(pair_row$terminal_group[[1]])
  label_vars <- unique(c(coarse, fine))
  label_vars <- label_vars[nzchar(label_vars) & label_vars %in% colnames(seu@meta.data)]
  if (length(label_vars) == 0) {
    return(empty_df_09(split_index_cols_09))
  }

  rows <- list()
  for (label_var in label_vars) {
    values <- as.character(seu@meta.data[[label_var]])
    tab <- sort(table(values), decreasing = TRUE)
    if (length(tab) == 0) {
      next
    }
    total <- sum(tab)
    rows[[length(rows) + 1L]] <- data.frame(
      pair_id = pair_id,
      split_var = split_var,
      split_value = display_scalar_value(split_value, "pooled"),
      label_var = label_var,
      label_value = names(tab),
      n_cells = as.integer(tab),
      fraction = as.numeric(tab) / total,
      root_group = root,
      terminal_group = terminal,
      is_root_group = names(tab) == root,
      is_terminal_group = names(tab) == terminal,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) {
    return(empty_df_09(split_index_cols_09))
  }
  dplyr::bind_rows(rows)
}

build_pca_summary_09 <- function(seu, pair_id, split_value, dims_used) {
  stdev <- tryCatch(seu[["pca"]]@stdev, error = function(e) numeric(0))
  if (length(stdev) == 0) {
    return(empty_df_09(pca_summary_cols_09))
  }
  variance <- stdev^2
  frac <- variance / sum(variance)
  data.frame(
    pair_id = pair_id,
    split_value = display_scalar_value(split_value, "pooled"),
    pc = seq_along(stdev),
    stdev = stdev,
    variance_fraction = frac,
    cumulative_variance_fraction = cumsum(frac),
    n_cells = ncol(seu),
    n_features = nrow(seu),
    dims_used = paste(dims_used, collapse = ","),
    stringsAsFactors = FALSE
  )
}

prepare_one_split_09 <- function(pair_row, layer_row, source_obj, split_var, split_value) {
  pair_id <- normalize_scalar_value(pair_row$trajectory_id[[1]])
  layer_scope <- normalize_scalar_value(pair_row$layer_scope[[1]])
  split_label <- display_scalar_value(split_value, "pooled")
  message("09a trajectory input: ", pair_id, " split=", split_label)

  split_obj <- source_obj
  if (nzchar(split_var)) {
    values <- trimws(as.character(split_obj@meta.data[[split_var]]))
    cells <- rownames(split_obj@meta.data)[values == split_value]
    if (length(cells) == 0) {
      stop(sprintf("split %s=%s matched zero cells", split_var, split_value), call. = FALSE)
    }
    split_obj <- subset(split_obj, cells = cells)
  }

  coarse <- normalize_scalar_value(pair_row$coarse_label_var[[1]])
  fine <- normalize_scalar_value(pair_row$fine_label_var[[1]])
  missing_labels <- setdiff(c(coarse, fine), colnames(split_obj@meta.data))
  missing_labels <- missing_labels[nzchar(missing_labels)]
  if (length(missing_labels) > 0) {
    stop(sprintf("object metadata is missing trajectory label columns: %s", paste(missing_labels, collapse = ",")), call. = FALSE)
  }

  assay <- trajectory_default_assay_09(split_obj)
  qc_df <- metric_from_meta_or_counts_09(split_obj, assay)
  outlier <- add_outlier_flags_09(qc_df, pair_row$outlier_qc_policy[[1]], cfg)
  retained_cells <- outlier$cells$cell_id[outlier$cells$retained]
  outlier_n <- sum(!outlier$cells$retained)

  outlier_summary <- data.frame(
    pair_id = pair_id,
    split_value = split_label,
    policy = normalize_scalar_value(pair_row$outlier_qc_policy[[1]], "standard"),
    cell_n_before = nrow(outlier$cells),
    outlier_n = outlier_n,
    cell_n_after = length(retained_cells),
    retention_fraction = safe_rate(length(retained_cells), nrow(outlier$cells)),
    nfeature_low_cutoff = outlier$summary$nfeature_low_cutoff,
    nfeature_high_cutoff = outlier$summary$nfeature_high_cutoff,
    percent_mt_high_cutoff = outlier$summary$percent_mt_high_cutoff,
    neighbor_distance_z_cutoff = outlier$summary$neighbor_distance_z_cutoff,
    filter_status = outlier$summary$filter_status,
    reason = outlier$summary$reason,
    stringsAsFactors = FALSE
  )
  outlier_cells <- cbind(
    data.frame(pair_id = pair_id, split_value = split_label, stringsAsFactors = FALSE),
    outlier$cells
  )

  if (length(retained_cells) < 3L) {
    return(list(
      index = data.frame(
        pair_id = pair_id,
        source_question_id = normalize_scalar_value(pair_row$source_question_id[[1]]),
        layer_scope = layer_scope,
        method = "trajectory",
        tools_to_run = normalize_scalar_value(pair_row$tools_to_run[[1]]),
        split_mode = normalize_scalar_value(pair_row$split_mode[[1]], "auto"),
        split_var = split_var,
        split_value = split_label,
        root_group = normalize_scalar_value(pair_row$root_group[[1]]),
        terminal_group = normalize_scalar_value(pair_row$terminal_group[[1]]),
        coarse_label_var = coarse,
        fine_label_var = fine,
        object_rds = layer_row$object_rds[[1]],
        input_rds = "",
        cell_n_before = nrow(outlier$cells),
        outlier_n = outlier_n,
        cell_n_after = length(retained_cells),
        retention_fraction = safe_rate(length(retained_cells), nrow(outlier$cells)),
        root_group_present = "not_evaluated",
        terminal_group_present = "not_evaluated",
        cc_status = "not_run",
        pca_npcs = 0L,
        umap_reduction = "",
        status = "skipped_low_cells",
        reason = "fewer than 3 retained cells",
        produced_at = timestamp_now(),
        stringsAsFactors = FALSE
      ),
      outlier_summary = outlier_summary,
      outlier_cells = outlier_cells,
      cc_summary = empty_df_09(cc_summary_cols_09),
      pca_summary = empty_df_09(pca_summary_cols_09),
      split_index = empty_df_09(split_index_cols_09),
      output_key = "",
      output_path = ""
    ))
  }

  work <- subset(split_obj, cells = retained_cells)
  work <- trajectory_reset_reductions_09(work)
  DefaultAssay(work) <- assay
  work$trajectory_pair_id <- pair_id
  work$trajectory_source_question_id <- normalize_scalar_value(pair_row$source_question_id[[1]])
  work$trajectory_layer_scope <- layer_scope
  work$trajectory_split_var <- split_var
  work$trajectory_split_value <- split_label
  work$trajectory_root_group <- normalize_scalar_value(pair_row$root_group[[1]])
  work$trajectory_terminal_group <- normalize_scalar_value(pair_row$terminal_group[[1]])
  work$trajectory_coarse_label_var <- coarse
  work$trajectory_fine_label_var <- fine

  work <- Seurat::NormalizeData(work, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  cc <- apply_cell_cycle_09(work, cfg, pair_row$regress_cell_cycle[[1]])
  work <- cc$object
  work <- Seurat::FindVariableFeatures(
    work,
    selection.method = "vst",
    nfeatures = min(cfg$trajectory_hvg_nfeatures, nrow(work)),
    verbose = FALSE
  )
  vars_arg <- if (length(cc$vars_to_regress) > 0) cc$vars_to_regress else NULL
  work <- Seurat::ScaleData(work, features = rownames(work), vars.to.regress = vars_arg, verbose = FALSE)

  features <- Seurat::VariableFeatures(work)
  if (length(features) < 2L) {
    stop("trajectory input has fewer than 2 variable features for PCA.", call. = FALSE)
  }
  npcs <- min(max(cfg$trajectory_pca_dims), ncol(work) - 1L, length(features) - 1L)
  if (!is.finite(npcs) || npcs < 2L) {
    stop("trajectory input has fewer than 2 usable PCA dimensions.", call. = FALSE)
  }
  work <- Seurat::RunPCA(work, features = features, npcs = npcs, verbose = FALSE)
  pca_mat <- Seurat::Embeddings(work, reduction = "pca")
  dims_used <- cfg$trajectory_pca_dims[cfg$trajectory_pca_dims <= ncol(pca_mat)]
  dims_used <- dims_used[is.finite(dims_used)]
  if (length(dims_used) < 2L) {
    stop("trajectory input has fewer than 2 requested PCA dims after RunPCA.", call. = FALSE)
  }
  n_neighbors <- min(cfg$trajectory_umap_n_neighbors, ncol(work) - 1L)
  work <- Seurat::RunUMAP(
    work,
    reduction = "pca",
    dims = dims_used,
    n.neighbors = n_neighbors,
    reduction.name = "umap",
    reduction.key = "UMAP_",
    verbose = FALSE
  )

  input_rds <- trajectory_input_rds_path_09(cfg, pair_id, split_value)
  ensure_dir(dirname(input_rds))
  saveRDS(work, input_rds)

  labels <- unique(as.character(work@meta.data[[coarse]]))
  root_present <- label_presence_09(labels, pair_row$root_group[[1]])
  terminal_present <- label_presence_09(labels, pair_row$terminal_group[[1]])
  cc_summary <- cbind(
    data.frame(pair_id = pair_id, split_value = split_label, stringsAsFactors = FALSE),
    cc$summary
  )

  output_key <- if (nzchar(normalize_scalar_value(split_value))) {
    sprintf("trajectory_input_%s__%s_seurat", safe_id_09(pair_id), safe_id_09(split_value))
  } else {
    sprintf("trajectory_input_%s_seurat", safe_id_09(pair_id))
  }

  list(
    index = data.frame(
      pair_id = pair_id,
      source_question_id = normalize_scalar_value(pair_row$source_question_id[[1]]),
      layer_scope = layer_scope,
      method = "trajectory",
      tools_to_run = normalize_scalar_value(pair_row$tools_to_run[[1]]),
      split_mode = normalize_scalar_value(pair_row$split_mode[[1]], "auto"),
      split_var = split_var,
      split_value = split_label,
      root_group = normalize_scalar_value(pair_row$root_group[[1]]),
      terminal_group = normalize_scalar_value(pair_row$terminal_group[[1]]),
      coarse_label_var = coarse,
      fine_label_var = fine,
      object_rds = layer_row$object_rds[[1]],
      input_rds = input_rds,
      cell_n_before = nrow(outlier$cells),
      outlier_n = outlier_n,
      cell_n_after = ncol(work),
      retention_fraction = safe_rate(ncol(work), nrow(outlier$cells)),
      root_group_present = root_present,
      terminal_group_present = terminal_present,
      cc_status = cc$summary$status[[1]],
      pca_npcs = ncol(pca_mat),
      umap_reduction = "umap",
      status = "ok",
      reason = "",
      produced_at = timestamp_now(),
      stringsAsFactors = FALSE
    ),
    outlier_summary = outlier_summary,
    outlier_cells = outlier_cells,
    cc_summary = cc_summary,
    pca_summary = build_pca_summary_09(work, pair_id, split_value, dims_used),
    split_index = build_split_index_09(work, pair_row, split_var, split_value),
    output_key = output_key,
    output_path = input_rds
  )
}

failed_split_09 <- function(pair_row, layer_row, split_var, split_value, err) {
  pair_id <- normalize_scalar_value(pair_row$trajectory_id[[1]])
  split_label <- display_scalar_value(split_value, "pooled")
  data.frame(
    pair_id = pair_id,
    source_question_id = normalize_scalar_value(pair_row$source_question_id[[1]]),
    layer_scope = normalize_scalar_value(pair_row$layer_scope[[1]]),
    method = "trajectory",
    tools_to_run = normalize_scalar_value(pair_row$tools_to_run[[1]]),
    split_mode = normalize_scalar_value(pair_row$split_mode[[1]], "auto"),
    split_var = split_var,
    split_value = split_label,
    root_group = normalize_scalar_value(pair_row$root_group[[1]]),
    terminal_group = normalize_scalar_value(pair_row$terminal_group[[1]]),
    coarse_label_var = normalize_scalar_value(pair_row$coarse_label_var[[1]]),
    fine_label_var = normalize_scalar_value(pair_row$fine_label_var[[1]]),
    object_rds = if (is.null(layer_row)) "" else layer_row$object_rds[[1]],
    input_rds = "",
    cell_n_before = 0L,
    outlier_n = 0L,
    cell_n_after = 0L,
    retention_fraction = NA_real_,
    root_group_present = "not_evaluated",
    terminal_group_present = "not_evaluated",
    cc_status = "not_run",
    pca_npcs = 0L,
    umap_reduction = "",
    status = "failed",
    reason = conditionMessage(err),
    produced_at = timestamp_now(),
    stringsAsFactors = FALSE
  )
}

pairs <- read_trajectory_pairs_09(cfg)
if (nrow(pairs) == 0) {
  message("09a: no active method=trajectory rows in trajectory_pairs.tsv")
}

index_rows <- list()
outlier_summary_rows <- list()
outlier_cell_rows <- list()
cc_rows <- list()
pca_rows <- list()
split_rows <- list()
dynamic_outputs <- list()

for (i in seq_len(nrow(pairs))) {
  pair_row <- pairs[i, , drop = FALSE]
  layer_row <- NULL
  source_obj <- NULL
  split_plan <- NULL

  setup <- tryCatch({
    layer_row <- resolve_trajectory_layer_09(cfg, pair_row$layer_scope[[1]])
    source_obj <- load_comm_layer_object_07(layer_row)
    split_plan <- split_values_for_pair_09(source_obj, pair_row)
    list(ok = TRUE, error = NULL)
  }, error = function(e) {
    list(ok = FALSE, error = e)
  })

  if (!setup$ok) {
    index_rows[[length(index_rows) + 1L]] <- failed_split_09(pair_row, layer_row, "", "", setup$error)
    next
  }

  for (split_value in split_plan$values) {
    result <- tryCatch(
      prepare_one_split_09(pair_row, layer_row, source_obj, split_plan$split_var, split_value),
      error = function(e) list(index = failed_split_09(pair_row, layer_row, split_plan$split_var, split_value, e))
    )
    index_rows[[length(index_rows) + 1L]] <- result$index
    if (!is.null(result$outlier_summary)) outlier_summary_rows[[length(outlier_summary_rows) + 1L]] <- result$outlier_summary
    if (!is.null(result$outlier_cells)) outlier_cell_rows[[length(outlier_cell_rows) + 1L]] <- result$outlier_cells
    if (!is.null(result$cc_summary) && nrow(result$cc_summary) > 0) cc_rows[[length(cc_rows) + 1L]] <- result$cc_summary
    if (!is.null(result$pca_summary) && nrow(result$pca_summary) > 0) pca_rows[[length(pca_rows) + 1L]] <- result$pca_summary
    if (!is.null(result$split_index) && nrow(result$split_index) > 0) split_rows[[length(split_rows) + 1L]] <- result$split_index
    if (!is.null(result$output_key) && nzchar(result$output_key) && file.exists(result$output_path)) {
      dynamic_outputs[[result$output_key]] <- build_output_entry(
        result$output_path,
        "rds",
        module_name,
        "trajectory-ready Seurat object",
        base_dir = cfg$project_root
      )
    }
  }
}

input_index <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else empty_df_09(trajectory_input_index_cols_09)
outlier_summary <- if (length(outlier_summary_rows) > 0) dplyr::bind_rows(outlier_summary_rows) else empty_df_09(outlier_summary_cols_09)
outlier_cells <- if (length(outlier_cell_rows) > 0) dplyr::bind_rows(outlier_cell_rows) else empty_df_09(outlier_cell_cols_09)
cc_summary <- if (length(cc_rows) > 0) dplyr::bind_rows(cc_rows) else empty_df_09(cc_summary_cols_09)
pca_summary <- if (length(pca_rows) > 0) dplyr::bind_rows(pca_rows) else empty_df_09(pca_summary_cols_09)
split_index <- if (length(split_rows) > 0) dplyr::bind_rows(split_rows) else empty_df_09(split_index_cols_09)

write_tsv_local(input_index, cfg$trajectory_input_index_tsv)
write_tsv_local(outlier_summary, cfg$trajectory_outlier_qc_summary_tsv)
write_tsv_local(outlier_cells, cfg$trajectory_outlier_qc_cells_tsv)
write_tsv_local(cc_summary, cfg$trajectory_cc_score_summary_tsv)
write_tsv_local(pca_summary, cfg$trajectory_pca_summary_tsv)
write_tsv_local(split_index, cfg$trajectory_split_index_tsv)

if (file.exists(cfg$module_09a_manifest_path)) {
  unlink(cfg$module_09a_manifest_path)
}
fixed_outputs <- list(
  trajectory_input_index = build_output_entry(cfg$trajectory_input_index_tsv, "tsv", module_name, "one row per trajectory input object or skipped/failed split", base_dir = cfg$project_root, schema = infer_schema_from_df(input_index)),
  outlier_qc_summary = build_output_entry(cfg$trajectory_outlier_qc_summary_tsv, "tsv", module_name, "outlier QC summary by pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(outlier_summary)),
  outlier_qc_cells = build_output_entry(cfg$trajectory_outlier_qc_cells_tsv, "tsv", module_name, "cell-level outlier QC calls before trajectory input filtering", base_dir = cfg$project_root, schema = infer_schema_from_df(outlier_cells)),
  cc_score_summary = build_output_entry(cfg$trajectory_cc_score_summary_tsv, "tsv", module_name, "cell-cycle scoring and regression status by pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(cc_summary)),
  pca_summary = build_output_entry(cfg$trajectory_pca_summary_tsv, "tsv", module_name, "PCA variance summary by pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(pca_summary)),
  split_index = build_output_entry(cfg$trajectory_split_index_tsv, "tsv", module_name, "label composition by pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(split_index))
)
write_manifest_local(
  manifest_path = cfg$module_09a_manifest_path,
  new_outputs = c(fixed_outputs, dynamic_outputs),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    trajectory_pairs = cfg$trajectory_pairs_sheet,
    module_00c = cfg$module_00c_manifest_path,
    module_03d = cfg$module_03d_manifest_path,
    module_04b = cfg$module_04b_manifest_path,
    layer_status = cfg$layer_status_file
  ),
  version = module_version_09(cfg, module_name),
  depends_on = list(
    module_00c = cfg$module_00c_manifest_path,
    module_03d = cfg$module_03d_manifest_path,
    module_04b = cfg$module_04b_manifest_path
  )
)

message("09a completed. trajectory input index: ", cfg$trajectory_input_index_tsv)
