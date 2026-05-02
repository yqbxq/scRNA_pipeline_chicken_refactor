empty_decoupler_activity_08 <- function() {
  empty_df_07(c(
    "layer_id", "resource", "method", "source_id", "group",
    "score", "statistic", "p_value", "target_n", "status"
  ))
}

empty_decoupler_mapping_summary_08 <- function() {
  empty_df_07(c(
    "resource", "source_id", "human_target_n", "mapped_chicken_target_n",
    "mapping_rate", "kept", "drop_reason"
  ))
}

empty_decoupler_resource_fingerprint_08 <- function() {
  empty_df_07(c(
    "resource", "source", "species_origin", "mapped_species",
    "row_n_raw", "row_n_mapped", "package_version", "resource_md5"
  ))
}

resolve_decoupler_ortholog_map_path_08 <- function(cfg) {
  env_path <- normalize_scalar_value(Sys.getenv("SCENIC_ORTHOLOG_MAP_FILE", ""))
  if (nzchar(env_path) && file.exists(env_path)) {
    return(env_path)
  }

  legacy_pipeline_map <- file.path(cfg$ortholog_cache_dir, "chicken_human_orthologs_for_pipeline.csv")
  if (file.exists(legacy_pipeline_map)) {
    return(legacy_pipeline_map)
  }

  map_path <- resolve_manifest_output_optional_07(
    cfg$ortholog_manifest_path,
    c("human_best", "chicken_human_best_csv", "human_best_csv")
  )
  if (nzchar(map_path) && file.exists(map_path)) {
    return(map_path)
  }

  file.path(cfg$ortholog_cache_dir, "chicken_human_orthologs.csv")
}

load_human_to_chicken_map_08 <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("缺少 decoupleR ortholog map: %s", path), call. = FALSE)
  }

  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  chicken_col <- pick_first_col_07(raw, c(
    "external_gene_name", "chicken_symbol", "source_gene_name",
    "source_symbol", "gene_name", "gene_symbol"
  ))
  human_col <- pick_first_col_07(raw, c(
    "target_gene_name", "hsapiens_homolog_associated_gene_name",
    "human_symbol", "target_symbol"
  ))
  if (!nzchar(chicken_col) || !nzchar(human_col)) {
    stop(sprintf("Ortholog map 缺少 chicken/human symbol 列: %s", path), call. = FALSE)
  }

  if ("target_species" %in% colnames(raw)) {
    human_rows <- tolower(trimws(as.character(raw$target_species))) %in% c("", "human", "hsapiens", "homo_sapiens")
    if (any(human_rows)) {
      raw <- raw[human_rows, , drop = FALSE]
    }
  }

  out <- data.frame(
    human_symbol = trimws(as.character(raw[[human_col]])),
    chicken_symbol = trimws(as.character(raw[[chicken_col]])),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$human_symbol) & !is.na(out$human_symbol) &
    nzchar(out$chicken_symbol) & !is.na(out$chicken_symbol), , drop = FALSE]
  out$human_key <- toupper(out$human_symbol)
  out <- unique(out[, c("human_key", "human_symbol", "chicken_symbol"), drop = FALSE])
  rownames(out) <- NULL
  out
}

selected_regulation_layers_08 <- function(cfg) {
  layers <- communication_layer_status_07(cfg)
  requested <- unique(trimws(as.character(cfg$regulation_layers)))
  requested <- requested[nzchar(requested)]
  if (length(requested) == 0) {
    requested <- cfg$panorama_layer_id
  }
  if (!any(tolower(requested) %in% c("all", "*"))) {
    layers <- layers[layers$layer_id %in% requested, , drop = FALSE]
  }
  layers
}

decoupler_cache_paths_08 <- function(cfg, resource) {
  list(
    raw_rds = file.path(cfg$decoupler_resource_dir, paste0(resource, "_human_raw.rds")),
    mapped_tsv = file.path(cfg$decoupler_resource_dir, paste0(resource, "_human_to_chicken.tsv"))
  )
}

load_package_data_08 <- function(package, data_names) {
  data_env <- new.env(parent = emptyenv())
  for (data_name in data_names) {
    suppressWarnings(utils::data(list = data_name, package = package, envir = data_env))
    if (exists(data_name, envir = data_env, inherits = FALSE)) {
      return(get(data_name, envir = data_env, inherits = FALSE))
    }
  }
  NULL
}

load_dorothea_human_08 <- function(levels) {
  raw <- load_package_data_08("dorothea", c("dorothea_hs", "dorothea_human"))
  if (is.null(raw) && exists("get_dorothea", envir = asNamespace("decoupleR"), inherits = FALSE)) {
    raw <- get("get_dorothea", envir = asNamespace("decoupleR"))(organism = "human", levels = levels)
  }
  if (is.null(raw)) {
    stop("无法从 dorothea/decoupleR 载入 human DoRothEA resource", call. = FALSE)
  }

  raw <- as.data.frame(raw, stringsAsFactors = FALSE)
  tf_col <- pick_first_col_07(raw, c("tf", "source", "TF"))
  target_col <- pick_first_col_07(raw, c("target", "gene", "target_gene"))
  mor_col <- pick_first_col_07(raw, c("mor", "weight", "mode_of_regulation"))
  likelihood_col <- pick_first_col_07(raw, c("likelihood", "confidence_score"))
  confidence_col <- pick_first_col_07(raw, c("confidence", "level"))
  if (!nzchar(tf_col) || !nzchar(target_col)) {
    stop("DoRothEA resource 缺少 TF/source 或 target 列", call. = FALSE)
  }
  if (!nzchar(mor_col)) {
    raw$mor <- 1
    mor_col <- "mor"
  }
  if (!nzchar(likelihood_col)) {
    raw$likelihood <- 1
    likelihood_col <- "likelihood"
  }
  if (!nzchar(confidence_col)) {
    raw$confidence <- ""
    confidence_col <- "confidence"
  }

  out <- data.frame(
    resource = "dorothea",
    source = trimws(as.character(raw[[tf_col]])),
    target = trimws(as.character(raw[[target_col]])),
    mor = suppressWarnings(as.numeric(raw[[mor_col]])),
    likelihood = suppressWarnings(as.numeric(raw[[likelihood_col]])),
    confidence = trimws(as.character(raw[[confidence_col]])),
    stringsAsFactors = FALSE
  )
  out$mor[is.na(out$mor)] <- 1
  out$likelihood[is.na(out$likelihood)] <- 1
  out <- out[nzchar(out$source) & nzchar(out$target), , drop = FALSE]
  if (length(levels) > 0 && any(nzchar(out$confidence))) {
    out <- out[out$confidence %in% levels, , drop = FALSE]
  }
  unique(out)
}

load_progeny_human_08 <- function() {
  raw <- load_package_data_08("progeny", c("model_human", "progenyModel"))
  if (is.null(raw) && exists("getModel", envir = asNamespace("progeny"), inherits = FALSE)) {
    raw <- get("getModel", envir = asNamespace("progeny"))(organism = "Human", top = 1000)
  }
  if (is.null(raw) && exists("get_progeny", envir = asNamespace("decoupleR"), inherits = FALSE)) {
    raw <- get("get_progeny", envir = asNamespace("decoupleR"))(organism = "human", top = 1000)
  }
  if (is.null(raw)) {
    stop("无法从 progeny/decoupleR 载入 human PROGENy resource", call. = FALSE)
  }

  raw <- as.data.frame(raw, stringsAsFactors = FALSE)
  pathway_col <- pick_first_col_07(raw, c("pathway", "source", "Pathway"))
  gene_col <- pick_first_col_07(raw, c("gene", "target", "Gene"))
  weight_col <- pick_first_col_07(raw, c("weight", "mor", "coefficient"))
  p_col <- pick_first_col_07(raw, c("p.value", "p_value", "pval"))
  if (!nzchar(pathway_col) || !nzchar(gene_col)) {
    stop("PROGENy resource 缺少 pathway/source 或 gene/target 列", call. = FALSE)
  }
  if (!nzchar(weight_col)) {
    raw$weight <- 1
    weight_col <- "weight"
  }
  if (!nzchar(p_col)) {
    raw$p.value <- NA_real_
    p_col <- "p.value"
  }

  out <- data.frame(
    resource = "progeny",
    source = trimws(as.character(raw[[pathway_col]])),
    target = trimws(as.character(raw[[gene_col]])),
    mor = suppressWarnings(as.numeric(raw[[weight_col]])),
    likelihood = abs(suppressWarnings(as.numeric(raw[[weight_col]]))),
    confidence = "",
    p_value = suppressWarnings(as.numeric(raw[[p_col]])),
    stringsAsFactors = FALSE
  )
  out$mor[is.na(out$mor)] <- 1
  out$likelihood[is.na(out$likelihood)] <- 1
  out <- out[nzchar(out$source) & nzchar(out$target), , drop = FALSE]
  unique(out)
}

load_or_build_decoupler_raw_resource_08 <- function(cfg, resource) {
  paths <- decoupler_cache_paths_08(cfg, resource)
  use_cache <- identical(tolower(cfg$decoupler_use_cache), "yes")
  if (use_cache && file.exists(paths$raw_rds) && file.info(paths$raw_rds)$size > 0) {
    return(list(data = readRDS(paths$raw_rds), source = "cache", path = paths$raw_rds))
  }

  raw <- switch(
    resource,
    dorothea = load_dorothea_human_08(cfg$decoupler_confidence_levels),
    progeny = load_progeny_human_08(),
    stop(sprintf("Unsupported decoupleR resource: %s", resource), call. = FALSE)
  )
  ensure_dir(dirname(paths$raw_rds))
  saveRDS(raw, paths$raw_rds)
  list(data = raw, source = "package", path = paths$raw_rds)
}

map_decoupler_network_to_chicken_08 <- function(raw_network, ortholog_map, min_targets) {
  raw <- as.data.frame(raw_network, stringsAsFactors = FALSE)
  raw$human_target_key <- toupper(raw$target)
  map <- ortholog_map[, c("human_key", "chicken_symbol"), drop = FALSE]
  merged <- merge(raw, map, by.x = "human_target_key", by.y = "human_key", all.x = FALSE, all.y = FALSE)
  if (nrow(merged) == 0) {
    mapped <- raw[FALSE, , drop = FALSE]
    mapped$target <- character(0)
  } else {
    merged$target <- merged$chicken_symbol
    mapped <- unique(merged[, intersect(
      c("resource", "source", "target", "mor", "likelihood", "confidence", "p_value"),
      colnames(merged)
    ), drop = FALSE])
  }

  source_ids <- sort(unique(raw$source))
  summary_rows <- lapply(source_ids, function(source_id) {
    raw_targets <- unique(raw$target[raw$source == source_id])
    mapped_targets <- unique(mapped$target[mapped$source == source_id])
    mapped_n <- length(mapped_targets[nzchar(mapped_targets)])
    human_n <- length(raw_targets[nzchar(raw_targets)])
    kept <- mapped_n >= min_targets
    data.frame(
      resource = raw$resource[[1]],
      source_id = source_id,
      human_target_n = human_n,
      mapped_chicken_target_n = mapped_n,
      mapping_rate = if (human_n > 0) mapped_n / human_n else NA_real_,
      kept = ifelse(kept, "yes", "no"),
      drop_reason = ifelse(kept, "", ifelse(mapped_n == 0, "no_mapped_targets", "too_few_mapped_targets")),
      stringsAsFactors = FALSE
    )
  })
  mapping_summary <- if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else empty_decoupler_mapping_summary_08()
  kept_sources <- mapping_summary$source_id[mapping_summary$kept == "yes"]
  mapped <- mapped[mapped$source %in% kept_sources, , drop = FALSE]
  rownames(mapped) <- NULL
  list(network = mapped, mapping_summary = mapping_summary)
}

write_resource_cache_tsv_08 <- function(network, path) {
  write_tsv_local(as.data.frame(network, stringsAsFactors = FALSE), path)
  normalize_path_07(path)
}

package_version_string_08 <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return("")
  }
  as.character(utils::packageVersion(pkg))
}

resource_fingerprint_row_08 <- function(resource, source, raw_network, mapped_network, mapped_tsv) {
  package_name <- if (identical(resource, "dorothea")) "dorothea" else "progeny"
  data.frame(
    resource = resource,
    source = source,
    species_origin = "human",
    mapped_species = "chicken",
    row_n_raw = nrow(raw_network),
    row_n_mapped = nrow(mapped_network),
    package_version = paste0(
      package_name, "=", package_version_string_08(package_name),
      ";decoupleR=", package_version_string_08("decoupleR")
    ),
    resource_md5 = unname(tools::md5sum(mapped_tsv)),
    stringsAsFactors = FALSE
  )
}

get_decoupler_expression_matrix_08 <- function(obj) {
  data_mat <- get_assay_matrix(obj, assay = "RNA", type = "data")
  if (nrow(data_mat) > 0 && ncol(data_mat) > 0 && Matrix::nnzero(data_mat) > 0) {
    return(list(matrix = data_mat, source = "RNA:data"))
  }

  counts <- get_assay_matrix(obj, assay = "RNA", type = "counts")
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("decoupleR 输入对象缺少 RNA counts/data matrix", call. = FALSE)
  }
  lib_size <- Matrix::colSums(counts)
  scale <- ifelse(lib_size > 0, 10000 / lib_size, 0)
  norm <- t(t(counts) * scale)
  norm <- log1p(norm)
  list(matrix = norm, source = "RNA:counts_log1p_cpm")
}

resolve_decoupler_group_col_08 <- function(obj, cfg, layer_id) {
  configured <- normalize_scalar_value(cfg$decoupler_group_col, "auto")
  meta_cols <- colnames(obj@meta.data)
  if (!identical(configured, "auto")) {
    if (!configured %in% meta_cols) {
      stop(sprintf("DECOUPLER_GROUP_COL=%s 但对象 metadata 中没有该列", configured), call. = FALSE)
    }
    return(configured)
  }

  candidates <- if (identical(layer_id, cfg$panorama_layer_id)) {
    c("cell_type", "annotation_label", "panorama_cell_type", "cell_subtype")
  } else if (identical(layer_id, "GC_subcluster")) {
    c("cell_subtype", paste0(layer_id, "_cell_type"), "cell_type", "annotation_label")
  } else {
    c(paste0(layer_id, "_cell_type"), "cell_subtype", "cell_type", "annotation_label", "panorama_cell_type")
  }
  hit <- candidates[candidates %in% meta_cols]
  if (length(hit) == 0) {
    stop(sprintf("layer %s 缺少 decoupleR group metadata column", layer_id), call. = FALSE)
  }
  hit[[1]]
}

group_average_expression_08 <- function(mat, groups) {
  groups <- trimws(as.character(groups))
  names(groups) <- colnames(mat)
  keep_cells <- names(groups)[nzchar(groups) & !is.na(groups) & names(groups) %in% colnames(mat)]
  if (length(keep_cells) == 0) {
    stop("decoupleR group_average 没有可用细胞", call. = FALSE)
  }
  mat <- mat[, keep_cells, drop = FALSE]
  groups <- groups[keep_cells]
  group_levels <- sort(unique(groups))
  avg <- vapply(group_levels, function(group) {
    cells <- names(groups)[groups == group]
    Matrix::rowMeans(mat[, cells, drop = FALSE])
  }, numeric(nrow(mat)))
  rownames(avg) <- rownames(mat)
  colnames(avg) <- group_levels
  avg
}

run_decoupler_method_08 <- function(mat, network, method, min_targets) {
  if (nrow(network) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  method <- tolower(method)
  fn_name <- switch(
    method,
    wmean = "run_wmean",
    ulm = "run_ulm",
    wsum = "run_wsum",
    mlm = "run_mlm",
    stop(sprintf("Unsupported decoupleR method: %s", method), call. = FALSE)
  )
  fn <- getExportedValue("decoupleR", fn_name)
  formals_names <- names(formals(fn))
  args <- list(
    mat = mat,
    network = network,
    .source = "source",
    .target = "target",
    .mor = "mor"
  )
  if ("minsize" %in% formals_names) {
    args$minsize <- min_targets
  }
  if (".likelihood" %in% formals_names && "likelihood" %in% colnames(network)) {
    args$.likelihood <- "likelihood"
  }
  if ("times" %in% formals_names && identical(method, "wmean")) {
    args$times <- 100L
  }
  as.data.frame(do.call(fn, args), stringsAsFactors = FALSE)
}

standardize_decoupler_result_08 <- function(result, layer_id, resource, method, network) {
  if (nrow(result) == 0) {
    return(empty_decoupler_activity_08())
  }
  source_col <- pick_first_col_07(result, c("source", "tf", "pathway", "source_id"))
  condition_col <- pick_first_col_07(result, c("condition", "sample", "group", "sample_id"))
  score_col <- pick_first_col_07(result, c("score", "statistic", "norm_wmean", "ulm", "estimate"))
  statistic_col <- pick_first_col_07(result, c("statistic", "method"))
  p_col <- pick_first_col_07(result, c("p_value", "p.value", "pval"))
  if (!nzchar(source_col) || !nzchar(condition_col) || !nzchar(score_col)) {
    stop(
      sprintf("无法识别 decoupleR 输出列: %s", paste(colnames(result), collapse = ", ")),
      call. = FALSE
    )
  }

  target_n <- stats::aggregate(
    target ~ source,
    data = unique(network[, c("source", "target"), drop = FALSE]),
    FUN = length
  )
  names(target_n) <- c("source_id", "target_n")
  out <- data.frame(
    layer_id = layer_id,
    resource = resource,
    method = method,
    source_id = trimws(as.character(result[[source_col]])),
    group = trimws(as.character(result[[condition_col]])),
    score = suppressWarnings(as.numeric(result[[score_col]])),
    statistic = if (nzchar(statistic_col)) as.character(result[[statistic_col]]) else method,
    p_value = if (nzchar(p_col)) suppressWarnings(as.numeric(result[[p_col]])) else NA_real_,
    stringsAsFactors = FALSE
  )
  out <- merge(out, target_n, by = "source_id", all.x = TRUE, sort = FALSE)
  out$status <- ifelse(is.na(out$score), "missing_score", "ok")
  out[, colnames(empty_decoupler_activity_08()), drop = FALSE]
}

write_decoupler_heatmap_08 <- function(activity, path, title, top_n = 25L) {
  ensure_dir(dirname(path))
  if (nrow(activity) == 0 || !any(is.finite(activity$score))) {
    write_empty_png_07(path, title = title)
    return(path)
  }

  score_summary <- stats::aggregate(
    abs(score) ~ source_id,
    data = activity[is.finite(activity$score), , drop = FALSE],
    FUN = max
  )
  names(score_summary) <- c("source_id", "max_abs_score")
  score_summary <- score_summary[order(-score_summary$max_abs_score, score_summary$source_id), , drop = FALSE]
  keep_sources <- head(score_summary$source_id, top_n)
  plot_df <- activity[activity$source_id %in% keep_sources, c("source_id", "group", "score"), drop = FALSE]
  wide <- tidyr::pivot_wider(plot_df, names_from = "group", values_from = "score")
  wide <- as.data.frame(wide, stringsAsFactors = FALSE)
  rownames(wide) <- wide$source_id
  wide$source_id <- NULL
  mat <- as.matrix(wide)
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    write_empty_png_07(path, title = title)
    return(path)
  }

  grDevices::png(path, width = 1200, height = 900, res = 150)
  pheatmap::pheatmap(
    mat,
    cluster_rows = nrow(mat) > 1,
    cluster_cols = ncol(mat) > 1,
    scale = "row",
    main = title,
    border_color = NA
  )
  grDevices::dev.off()
  path
}
