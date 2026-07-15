empty_cluster_marker_df <- function() {
  data.frame(
    p_val = numeric(0),
    avg_log2FC = numeric(0),
    pct.1 = numeric(0),
    pct.2 = numeric(0),
    p_val_adj = numeric(0),
    cluster = character(0),
    gene = character(0),
    pct_diff = numeric(0),
    pct_ratio = numeric(0),
    specificity_score = numeric(0),
    rank_specificity = integer(0),
    stringsAsFactors = FALSE
  )
}

standardize_marker_columns <- function(markers) {
  if (is.null(markers) || nrow(markers) == 0) {
    return(empty_cluster_marker_df())
  }
  markers <- as.data.frame(markers, stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(markers)) {
    markers$gene <- rownames(markers)
  }
  if (!"avg_log2FC" %in% colnames(markers) && "avg_logFC" %in% colnames(markers)) {
    markers$avg_log2FC <- markers$avg_logFC
  }
  for (col in c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) {
    if (!col %in% colnames(markers)) {
      markers[[col]] <- NA_real_
    }
    markers[[col]] <- suppressWarnings(as.numeric(markers[[col]]))
  }
  if (!"cluster" %in% colnames(markers)) {
    markers$cluster <- ""
  }
  markers$cluster <- as.character(markers$cluster)
  markers$gene <- as.character(markers$gene)
  markers$pct_diff <- markers$pct.1 - markers$pct.2
  markers$pct_ratio <- (markers$pct.1 + 1e-6) / (markers$pct.2 + 1e-6)
  markers$specificity_score <- markers$avg_log2FC * pmax(markers$pct_diff, 0)
  markers <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(
      dplyr::desc(specificity_score),
      dplyr::desc(avg_log2FC),
      p_val_adj,
      gene,
      .by_group = TRUE
    ) %>%
    dplyr::mutate(rank_specificity = dplyr::row_number()) %>%
    dplyr::ungroup()
  markers
}

filter_annotation_markers <- function(markers, cfg) {
  if (nrow(markers) == 0) {
    return(markers)
  }
  markers %>%
    dplyr::filter(
      is.finite(p_val_adj),
      p_val_adj < cfg$cluster_marker_anno_padj,
      avg_log2FC >= cfg$cluster_marker_anno_logfc,
      pct.1 >= cfg$cluster_marker_anno_pct1,
      pct_diff >= cfg$cluster_marker_anno_pct_diff
    ) %>%
    dplyr::arrange(cluster, rank_specificity)
}

filter_strict_markers <- function(markers, cfg) {
  if (nrow(markers) == 0) {
    return(markers)
  }
  markers %>%
    dplyr::filter(
      is.finite(p_val_adj),
      p_val_adj < cfg$cluster_marker_strict_padj,
      avg_log2FC >= cfg$cluster_marker_strict_logfc,
      pct.1 >= cfg$cluster_marker_strict_pct1,
      pct_diff >= cfg$cluster_marker_strict_pct_diff
    ) %>%
    dplyr::arrange(cluster, rank_specificity)
}

top_ranked_markers <- function(markers, top_n = 20L) {
  if (nrow(markers) == 0) {
    return(markers)
  }
  markers %>%
    dplyr::arrange(cluster, rank_specificity) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
}

build_cluster_marker_qc <- function(seu, cluster_col, markers_raw, markers_anno, markers_strict, cfg) {
  cluster_values <- as.character(seu@meta.data[[cluster_col]])
  cluster_ids <- sort(unique(cluster_values))
  cluster_sizes <- table(cluster_values)
  rows <- lapply(cluster_ids, function(cluster_id) {
    raw_i <- markers_raw[markers_raw$cluster == cluster_id, , drop = FALSE]
    anno_i <- markers_anno[markers_anno$cluster == cluster_id, , drop = FALSE]
    strict_i <- markers_strict[markers_strict$cluster == cluster_id, , drop = FALSE]
    top20 <- top_ranked_markers(raw_i, 20L)
    top50 <- top_ranked_markers(raw_i, 50L)
    strict_n <- nrow(strict_i)
    top20_mean_pct_diff <- if (nrow(top20) > 0) mean(top20$pct_diff, na.rm = TRUE) else NA_real_
    top20_median_specificity <- if (nrow(top20) > 0) stats::median(top20$specificity_score, na.rm = TRUE) else NA_real_
    top20_high_pct2_n <- if (nrow(top20) > 0) sum(top20$pct.2 >= 0.25, na.rm = TRUE) else 0L
    top20_neg_pct_diff_n <- if (nrow(top20) > 0) sum(top20$pct_diff <= 0, na.rm = TRUE) else 0L
    flag <- "ok"
    if (nrow(raw_i) == 0) {
      flag <- "no_positive_markers"
    } else if (strict_n < cfg$cluster_marker_min_strict_markers) {
      flag <- "few_strict_markers"
    } else if (is.finite(top20_mean_pct_diff) && top20_mean_pct_diff < cfg$cluster_marker_anno_pct_diff) {
      flag <- "weak_top_markers"
    } else if (top20_high_pct2_n >= ceiling(nrow(top20) * 0.5)) {
      flag <- "high_background_top_markers"
    }
    data.frame(
      cluster = cluster_id,
      cluster_col = cluster_col,
      n_cells = unname(cluster_sizes[[cluster_id]]),
      raw_marker_n = nrow(raw_i),
      anno_marker_n = nrow(anno_i),
      strict_marker_n = strict_n,
      top20_mean_pct_diff = top20_mean_pct_diff,
      top20_median_specificity_score = top20_median_specificity,
      top20_neg_pct_diff_n = top20_neg_pct_diff_n,
      top20_high_pct2_n = top20_high_pct2_n,
      top50_mean_pct_diff = if (nrow(top50) > 0) mean(top50$pct_diff, na.rm = TRUE) else NA_real_,
      deg_quality_flag = flag,
      risk_level = ifelse(flag == "ok", "low", ifelse(flag %in% c("few_strict_markers", "weak_top_markers"), "medium", "high")),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

run_cluster_marker_audit <- function(seu, cluster_col, cfg, assay = cfg$cluster_marker_assay) {
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop(sprintf("cluster column not found: %s", cluster_col), call. = FALSE)
  }
  if (!assay %in% Assays(seu)) {
    stop(sprintf("assay `%s` not found for marker audit", assay), call. = FALSE)
  }
  DefaultAssay(seu) <- assay
  Idents(seu) <- cluster_col
  if (length(unique(as.character(Idents(seu)))) < 2L) {
    markers_raw <- empty_cluster_marker_df()
  } else {
    markers_raw <- tryCatch(
      FindAllMarkers(
        seu,
        group.by = cluster_col,
        only.pos = TRUE,
        test.use = "wilcox",
        min.pct = cfg$cluster_marker_min_pct,
        logfc.threshold = cfg$cluster_marker_logfc_threshold,
        verbose = FALSE
      ),
      error = function(e) {
        warning(conditionMessage(e), call. = FALSE)
        empty_cluster_marker_df()
      }
    )
  }
  markers_scored <- standardize_marker_columns(markers_raw)
  markers_anno <- filter_annotation_markers(markers_scored, cfg)
  markers_strict <- filter_strict_markers(markers_scored, cfg)
  cluster_qc <- build_cluster_marker_qc(seu, cluster_col, markers_scored, markers_anno, markers_strict, cfg)
  list(
    markers_raw = markers_raw,
    markers_scored = markers_scored,
    markers_anno = markers_anno,
    markers_strict = markers_strict,
    cluster_qc = cluster_qc,
    top20_markers = top_ranked_markers(markers_scored, 20L),
    top50_markers = top_ranked_markers(markers_scored, 50L)
  )
}

write_cluster_marker_audit_outputs <- function(audit, out_dir) {
  ensure_dir(out_dir)
  paths <- list(
    markers_raw = file.path(out_dir, "markers_raw.tsv"),
    markers_scored = file.path(out_dir, "markers_scored.tsv"),
    markers_anno = file.path(out_dir, "markers_anno.tsv"),
    markers_strict = file.path(out_dir, "markers_strict.tsv"),
    cluster_marker_qc = file.path(out_dir, "cluster_marker_qc.tsv"),
    top20_markers = file.path(out_dir, "top20_markers.tsv"),
    top50_markers = file.path(out_dir, "top50_markers.tsv")
  )
  write_tsv_local(audit$markers_raw, paths$markers_raw)
  write_tsv_local(audit$markers_scored, paths$markers_scored)
  write_tsv_local(audit$markers_anno, paths$markers_anno)
  write_tsv_local(audit$markers_strict, paths$markers_strict)
  write_tsv_local(audit$cluster_qc, paths$cluster_marker_qc)
  write_tsv_local(audit$top20_markers, paths$top20_markers)
  write_tsv_local(audit$top50_markers, paths$top50_markers)
  paths
}

candidate_marker_quality_row <- function(candidate_row, audit, cfg, out_dir) {
  qc <- audit$cluster_qc
  cluster_count <- nrow(qc)
  supported <- if (cluster_count > 0) mean(qc$strict_marker_n >= cfg$cluster_marker_min_strict_markers, na.rm = TRUE) else NA_real_
  data.frame(
    resolution = candidate_row$resolution[[1]],
    candidate_id = candidate_row$candidate_id[[1]],
    candidate_cluster_col = candidate_row$candidate_cluster_col[[1]],
    n_clusters = candidate_row$n_clusters[[1]],
    cluster_marker_dir = out_dir,
    marker_supported_cluster_fraction = supported,
    median_strict_marker_n = if (cluster_count > 0) stats::median(qc$strict_marker_n, na.rm = TRUE) else NA_real_,
    minimum_strict_marker_n = if (cluster_count > 0) min(qc$strict_marker_n, na.rm = TRUE) else NA_real_,
    median_top20_pct_diff = if (cluster_count > 0) stats::median(qc$top20_mean_pct_diff, na.rm = TRUE) else NA_real_,
    fraction_clusters_with_few_strict_markers = if (cluster_count > 0) mean(qc$strict_marker_n < cfg$cluster_marker_min_strict_markers, na.rm = TRUE) else NA_real_,
    fraction_clusters_with_high_background_markers = if (cluster_count > 0) mean(qc$deg_quality_flag == "high_background_top_markers", na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}

marker_quality_score <- function(df) {
  if (nrow(df) == 0) {
    return(numeric(0))
  }
  score <- 0.35 * safe_minmax_score_local(df$marker_supported_cluster_fraction, TRUE) +
    0.25 * safe_minmax_score_local(log1p(df$median_strict_marker_n), TRUE) +
    0.20 * safe_minmax_score_local(df$median_top20_pct_diff, TRUE) +
    0.10 * safe_minmax_score_local(df$fraction_clusters_with_few_strict_markers, FALSE) +
    0.10 * safe_minmax_score_local(df$fraction_clusters_with_high_background_markers, FALSE)
  score[!is.finite(score)] <- 0.5
  pmax(0, pmin(1, score))
}

copy_marker_audit_outputs <- function(from_dir, to_dir) {
  ensure_dir(to_dir)
  files <- c(
    "markers_raw.tsv", "markers_scored.tsv", "markers_anno.tsv",
    "markers_strict.tsv", "cluster_marker_qc.tsv",
    "top20_markers.tsv", "top50_markers.tsv"
  )
  for (file in files) {
    source_path <- file.path(from_dir, file)
    target_path <- file.path(to_dir, file)
    if (file.exists(source_path)) {
      file.copy(source_path, target_path, overwrite = TRUE)
    }
  }
  invisible(TRUE)
}

nuisance_gene_class <- function(genes) {
  gene <- toupper(as.character(genes))
  out <- rep("", length(gene))
  out[grepl("^MT-|^MT\\.", gene)] <- "mitochondrial"
  out[out == "" & grepl("^(RPL|RPS|MRPL|MRPS)", gene)] <- "ribosomal"
  out[out == "" & gene %in% c("MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF", "CENPE", "UBE2C", "BIRC5", "CDC20", "CDK1", "CCNB1", "CCNB2")] <- "cell_cycle"
  out[out == "" & gene %in% c("FOS", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "ATF3", "HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1")] <- "stress_ieg"
  out[out == "" & grepl("^(HB|HBA|HBB|HBE|HBM|HBZ)", gene)] <- "hemoglobin"
  out
}

build_gene_marker_risk <- function(markers_scored) {
  if (nrow(markers_scored) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  markers_scored %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      marker_cluster_n = dplyr::n_distinct(cluster),
      max_pct1 = max(pct.1, na.rm = TRUE),
      max_pct2 = max(pct.2, na.rm = TRUE),
      max_pct_diff = max(pct_diff, na.rm = TRUE),
      max_specificity_score = max(specificity_score, na.rm = TRUE),
      best_rank_specificity = min(rank_specificity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      nuisance_class = nuisance_gene_class(gene),
      high_background_pct2 = max_pct2 >= 0.25,
      low_specificity = max_pct_diff < 0.10,
      multi_cluster_marker = marker_cluster_n >= 3L,
      broad_expression = max_pct1 >= 0.80 & max_pct2 >= 0.50,
      risk_flag = paste(
        ifelse(nzchar(nuisance_class), nuisance_class, ""),
        ifelse(high_background_pct2, "high_background_pct2", ""),
        ifelse(low_specificity, "low_specificity", ""),
        ifelse(multi_cluster_marker, "multi_cluster_marker", ""),
        ifelse(broad_expression, "broad_expression", ""),
        sep = ";"
      ),
      risk_flag = gsub("(^;+|;+$)", "", gsub(";+", ";", risk_flag))
    ) %>%
    dplyr::arrange(dplyr::desc(max_specificity_score), gene)
}

build_cluster_marker_risk <- function(top20_markers, gene_risk, cluster_qc) {
  if (nrow(cluster_qc) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  risk_lookup <- gene_risk[, c("gene", "nuisance_class", "high_background_pct2", "low_specificity", "multi_cluster_marker", "broad_expression"), drop = FALSE]
  merged <- dplyr::left_join(top20_markers, risk_lookup, by = "gene")
  merged$nuisance_class[is.na(merged$nuisance_class)] <- ""
  bool_cols <- c("high_background_pct2", "low_specificity", "multi_cluster_marker", "broad_expression")
  for (col in bool_cols) {
    merged[[col]][is.na(merged[[col]])] <- FALSE
  }
  rows <- lapply(cluster_qc$cluster, function(cluster_id) {
    m <- merged[merged$cluster == cluster_id, , drop = FALSE]
    n <- nrow(m)
    nuisance_fraction <- if (n > 0) mean(nzchar(m$nuisance_class), na.rm = TRUE) else NA_real_
    shared_fraction <- if (n > 0) mean(m$multi_cluster_marker, na.rm = TRUE) else NA_real_
    background_fraction <- if (n > 0) mean(m$high_background_pct2 | m$broad_expression, na.rm = TRUE) else NA_real_
    flag <- "ok"
    reasons <- character(0)
    if (is.finite(nuisance_fraction) && nuisance_fraction >= 0.40) {
      flag <- "nuisance_marker_dominated"
      reasons <- c(reasons, "top20 nuisance marker fraction >= 0.40")
    }
    if (is.finite(background_fraction) && background_fraction >= 0.40) {
      flag <- ifelse(flag == "ok", "high_background_marker_dominated", paste(flag, "high_background_marker_dominated", sep = ";"))
      reasons <- c(reasons, "top20 high-background marker fraction >= 0.40")
    }
    if (is.finite(shared_fraction) && shared_fraction >= 0.50) {
      flag <- ifelse(flag == "ok", "shared_marker_dominated", paste(flag, "shared_marker_dominated", sep = ";"))
      reasons <- c(reasons, "top20 multi-cluster marker fraction >= 0.50")
    }
    data.frame(
      cluster = cluster_id,
      top20_marker_n = n,
      nuisance_top20_fraction = nuisance_fraction,
      high_background_top20_fraction = background_fraction,
      shared_marker_fraction = shared_fraction,
      marker_quality_flag = flag,
      manual_review_reason = paste(reasons, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(rows)
  dplyr::left_join(cluster_qc, out, by = "cluster")
}
