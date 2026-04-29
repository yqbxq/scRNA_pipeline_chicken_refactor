estimate_expected_doublet_rate <- function(n_cells, rate_per_1k) {
  n_cells <- suppressWarnings(as.numeric(n_cells))
  rate_per_1k <- suppressWarnings(as.numeric(rate_per_1k))
  if (!is.finite(n_cells) || n_cells <= 0 || !is.finite(rate_per_1k) || rate_per_1k <= 0) {
    return(NA_real_)
  }
  rate <- rate_per_1k * (n_cells / 1000)
  min(max(rate, 0.001), 0.30)
}

estimate_expected_doublet_count <- function(n_cells, rate_per_1k) {
  rate <- estimate_expected_doublet_rate(n_cells, rate_per_1k)
  if (!is.finite(rate)) {
    return(NA_integer_)
  }
  as.integer(max(1, round(rate * n_cells)))
}

resolve_doubletfinder_fns <- function() {
  pkg <- asNamespace("DoubletFinder")
  list(
    sweep = if (exists("paramSweep_v3", envir = pkg, mode = "function")) {
      get("paramSweep_v3", envir = pkg)
    } else {
      get("paramSweep", envir = pkg)
    },
    summarize = get("summarizeSweep", envir = pkg),
    find_pk = get("find.pK", envir = pkg),
    run = if (exists("doubletFinder_v3", envir = pkg, mode = "function")) {
      get("doubletFinder_v3", envir = pkg)
    } else {
      get("doubletFinder", envir = pkg)
    }
  )
}

run_primary_scdblfinder <- function(model_obj, platform_resolved, expected_rate) {
  counts_mat <- get_assay_matrix(model_obj, assay = "RNA", type = "counts")
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_mat))
  SummarizedExperiment::colData(sce)$cluster <- as.character(model_obj$provisional_cluster)

  primary_args <- list(
    sce = sce,
    clusters = SummarizedExperiment::colData(sce)$cluster,
    verbose = FALSE
  )
  if (!identical(platform_resolved, "10x_cellranger") && is.finite(expected_rate)) {
    primary_args$dbr <- expected_rate
  }

  sce <- tryCatch(
    do.call(scDblFinder::scDblFinder, primary_args),
    error = function(e) {
      message(sprintf("scDblFinder failed: %s", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(sce)) {
    return(list(
      status = "failed",
      call = rep(NA_character_, ncol(model_obj)),
      score = rep(NA_real_, ncol(model_obj))
    ))
  }

  call <- as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)
  score <- suppressWarnings(as.numeric(SummarizedExperiment::colData(sce)$scDblFinder.score))
  names(call) <- colnames(model_obj)
  names(score) <- colnames(model_obj)
  list(
    status = if (!identical(platform_resolved, "10x_cellranger") && is.finite(expected_rate)) "completed_manual_rate" else "completed_auto_rate",
    call = call,
    score = score
  )
}

run_secondary_doubletfinder <- function(model_obj, usable_dims, expected_doublets, enabled = TRUE) {
  if (!isTRUE(enabled) || !requireNamespace("DoubletFinder", quietly = TRUE)) {
    return(list(
      status = if (!isTRUE(enabled)) "disabled" else "package_missing",
      call = rep(NA_character_, ncol(model_obj)),
      score = rep(NA_real_, ncol(model_obj)),
      selected_pK = NA_real_
    ))
  }

  fns <- resolve_doubletfinder_fns()
  expected_doublets <- suppressWarnings(as.integer(expected_doublets))
  if (is.na(expected_doublets) || expected_doublets < 1L) {
    expected_doublets <- max(1L, round(0.008 * ncol(model_obj)))
  }

  best_pk <- NA_real_
  status <- "completed"
  sweep_res <- tryCatch(
    fns$sweep(model_obj, PCs = usable_dims, sct = FALSE),
    error = function(e) NULL
  )
  if (!is.null(sweep_res)) {
    sweep_stats <- tryCatch(fns$summarize(sweep_res, GT = FALSE), error = function(e) NULL)
    bcmvn <- tryCatch(fns$find_pk(sweep_stats), error = function(e) NULL)
    if (!is.null(bcmvn) && nrow(bcmvn) > 0) {
      best_pk <- suppressWarnings(as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])))
    }
  }
  if (!is.finite(best_pk)) {
    best_pk <- 0.09
    status <- "completed_fallback_pk"
  }

  df_obj <- tryCatch(
    fns$run(
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
      call = rep(NA_character_, ncol(model_obj)),
      score = rep(NA_real_, ncol(model_obj)),
      selected_pK = best_pk
    ))
  }

  class_cols <- grep("DF.classifications", colnames(df_obj@meta.data), value = TRUE)
  pann_cols <- grep("^pANN", colnames(df_obj@meta.data), value = TRUE)
  call <- if (length(class_cols) > 0) as.character(df_obj@meta.data[[tail(class_cols, 1)]]) else rep(NA_character_, ncol(model_obj))
  score <- if (length(pann_cols) > 0) suppressWarnings(as.numeric(df_obj@meta.data[[tail(pann_cols, 1)]])) else rep(NA_real_, ncol(model_obj))
  names(call) <- colnames(model_obj)
  names(score) <- colnames(model_obj)

  list(
    status = status,
    call = call,
    score = score,
    selected_pK = best_pk
  )
}

is_primary_doublet <- function(x) {
  tolower(normalize_scalar_value(x)) == "doublet"
}

is_secondary_doublet <- function(x) {
  normalize_scalar_value(x) == "Doublet" || tolower(normalize_scalar_value(x)) == "doublet"
}

cohen_kappa_binary <- function(primary_doublet, secondary_doublet) {
  keep <- !is.na(primary_doublet) & !is.na(secondary_doublet)
  primary_doublet <- primary_doublet[keep]
  secondary_doublet <- secondary_doublet[keep]
  if (length(primary_doublet) == 0) {
    return(NA_real_)
  }
  po <- mean(primary_doublet == secondary_doublet)
  p_primary_true <- mean(primary_doublet)
  p_secondary_true <- mean(secondary_doublet)
  pe <- p_primary_true * p_secondary_true + (1 - p_primary_true) * (1 - p_secondary_true)
  if (!is.finite(pe) || pe >= 1) {
    return(NA_real_)
  }
  (po - pe) / (1 - pe)
}

classify_cluster_doublet_risk <- function(doublet_fraction, sample_baseline, n_doublet) {
  doublet_fraction <- suppressWarnings(as.numeric(doublet_fraction))
  sample_baseline <- suppressWarnings(as.numeric(sample_baseline))
  n_doublet <- suppressWarnings(as.numeric(n_doublet))
  if (!is.finite(doublet_fraction)) {
    return("unknown")
  }
  high_cutoff <- if (is.finite(sample_baseline) && sample_baseline > 0) max(0.15, sample_baseline * 1.5) else 0.15
  if (doublet_fraction >= high_cutoff && is.finite(n_doublet) && n_doublet >= 2) {
    return("high")
  }
  if (doublet_fraction >= 0.10) {
    return("medium")
  }
  "low"
}
