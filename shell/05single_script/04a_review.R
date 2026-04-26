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
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "reduction_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble", "jsonlite"))

cfg <- get_single_script_config_04()
module_name <- "04a_review"
prepare_dirs_04(cfg)
set.seed(cfg$random_seed)

empty_compare_04 <- function() {
  data.frame(
    layer_id = character(0),
    candidate_id = character(0),
    normalization = character(0),
    integration = character(0),
    reduction_name = character(0),
    umap_name = character(0),
    max_sample_r2 = numeric(0),
    max_group_r2 = numeric(0),
    silhouette_by_orig_ident = numeric(0),
    same_sample_knn_fraction = numeric(0),
    runtime_sec = numeric(0),
    downgrade_reason = character(0),
    rank_score = numeric(0),
    stringsAsFactors = FALSE
  )
}

r2_by_factor <- function(values, group) {
  group <- as.factor(group)
  if (length(unique(stats::na.omit(group))) < 2) {
    return(NA_real_)
  }
  fit <- stats::lm(values ~ group)
  summary(fit)$r.squared
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  max(x)
}

sample_rows <- function(n, max_n = 1500L) {
  if (n <= max_n) {
    return(seq_len(n))
  }
  sort(sample(seq_len(n), max_n))
}

mean_silhouette_by_group <- function(emb, group) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    return(NA_real_)
  }
  group <- as.character(group)
  keep <- !is.na(group) & nzchar(group)
  if (sum(keep) < 3 || length(unique(group[keep])) < 2) {
    return(NA_real_)
  }
  emb <- emb[keep, , drop = FALSE]
  group <- group[keep]
  idx <- sample_rows(nrow(emb))
  emb <- emb[idx, , drop = FALSE]
  group <- group[idx]
  if (length(unique(group)) < 2) {
    return(NA_real_)
  }
  sil <- cluster::silhouette(as.integer(factor(group)), stats::dist(emb))
  mean(sil[, "sil_width"], na.rm = TRUE)
}

same_sample_knn_fraction <- function(emb, group, k = 5L) {
  group <- as.character(group)
  keep <- !is.na(group) & nzchar(group)
  if (sum(keep) <= k || length(unique(group[keep])) < 2) {
    return(NA_real_)
  }
  emb <- emb[keep, , drop = FALSE]
  group <- group[keep]
  idx <- sample_rows(nrow(emb))
  emb <- emb[idx, , drop = FALSE]
  group <- group[idx]
  d <- as.matrix(stats::dist(emb))
  diag(d) <- Inf
  k <- min(k, nrow(d) - 1L)
  nn <- t(apply(d, 1, function(x) order(x)[seq_len(k)]))
  mean(vapply(seq_len(nrow(nn)), function(i) mean(group[nn[i, ]] == group[i]), numeric(1)), na.rm = TRUE)
}

rank_metric <- function(x, decreasing = FALSE) {
  if (length(x) == 0) {
    return(numeric(0))
  }
  if (all(is.na(x))) {
    return(rep(1, length(x)))
  }
  value <- x
  fill <- if (decreasing) min(value, na.rm = TRUE) - 1 else max(value, na.rm = TRUE) + 1
  value[is.na(value)] <- fill
  rank(if (decreasing) -value else value, ties.method = "min")
}

candidate_metrics_for_layer <- function(candidate_df) {
  rows <- list()
  for (idx in seq_len(nrow(candidate_df))) {
    row <- candidate_df[idx, , drop = FALSE]
    seu <- readRDS(row$out_rds[[1]])
    seu <- maybe_join_layers(seu)
    reduction_name <- row$reduction_name[[1]]
    umap_name <- row$umap_name[[1]]
    emb <- Embeddings(seu, reduction = reduction_name)
    pc_ids <- seq_len(min(20L, ncol(emb)))
    sample_group <- if ("orig.ident" %in% colnames(seu@meta.data)) seu$orig.ident else seu$sample_id
    sample_r2 <- vapply(pc_ids, function(i) r2_by_factor(emb[, i], sample_group), numeric(1))
    group_var <- if ("analysis_group" %in% colnames(seu@meta.data)) "analysis_group" else if ("condition" %in% colnames(seu@meta.data)) "condition" else ""
    group_r2 <- if (nzchar(group_var)) vapply(pc_ids, function(i) r2_by_factor(emb[, i], seu@meta.data[[group_var]]), numeric(1)) else rep(NA_real_, length(pc_ids))
    metric_emb <- if (umap_name %in% Reductions(seu)) Embeddings(seu, reduction = umap_name) else emb[, seq_len(min(2L, ncol(emb))), drop = FALSE]
    rows[[length(rows) + 1]] <- data.frame(
      layer_id = row$layer_id[[1]],
      candidate_id = row$candidate_id[[1]],
      normalization = row$normalization[[1]],
      integration = row$integration[[1]],
      reduction_name = reduction_name,
      umap_name = umap_name,
      max_sample_r2 = safe_max(sample_r2),
      max_group_r2 = safe_max(group_r2),
      silhouette_by_orig_ident = mean_silhouette_by_group(metric_emb, sample_group),
      same_sample_knn_fraction = same_sample_knn_fraction(metric_emb, sample_group),
      runtime_sec = suppressWarnings(as.numeric(row$runtime_sec[[1]])),
      downgrade_reason = normalize_scalar_value(row$downgrade_reason[[1]]),
      stringsAsFactors = FALSE
    )
  }
  out <- if (length(rows) > 0) dplyr::bind_rows(rows) else empty_compare_04()
  out$rank_score <- rank_metric(out$max_sample_r2) +
    rank_metric(out$silhouette_by_orig_ident) +
    rank_metric(out$same_sample_knn_fraction)
  dplyr::arrange(out, rank_score, max_sample_r2, silhouette_by_orig_ident)
}

write_selected_integration_if_needed <- function(path, selected_value, layer_id, mode) {
  ensure_dir(dirname(path))
  existing <- read_selected_integration_value_local(path)
  if (nzchar(existing)) {
    return(existing)
  }
  writeLines(
    c(
      sprintf("# automatic fallback recommendation from 04a_review for layer %s", layer_id),
      "# review subcluster_review_summary.tsv before approving the subcluster gate",
      sprintf("# mode=%s", mode),
      selected_value
    ),
    path,
    useBytes = TRUE
  )
  selected_value
}

manifest_04a <- read_manifest_local(cfg$module_04a_manifest_path)
candidate_index <- read_tsv_optional(resolve_output_local(manifest_04a, "candidate_index_tsv"))
subcluster_index <- read_tsv_optional(resolve_output_local(manifest_04a, "subcluster_index_tsv"))
if (nrow(candidate_index) == 0) {
  candidate_index <- data.frame(
    layer_id = character(0),
    normalization = character(0),
    integration = character(0),
    candidate_id = character(0),
    candidate_key = character(0),
    mode = character(0),
    reduction_name = character(0),
    umap_name = character(0),
    runtime_sec = numeric(0),
    downgrade_reason = character(0),
    out_rds = character(0),
    selected_integration_file = character(0),
    stringsAsFactors = FALSE
  )
}

summary_rows <- list()
compare_outputs <- list()

if (nrow(subcluster_index) > 0) {
  for (layer_id in subcluster_index$layer_id) {
    idx_row <- subcluster_index[subcluster_index$layer_id == layer_id, , drop = FALSE][1, , drop = FALSE]
    layer_candidates <- candidate_index[candidate_index$layer_id == layer_id, , drop = FALSE]
    selected_file <- normalize_scalar_value(idx_row$selected_integration_file[[1]], selected_integration_file_04(cfg, layer_id))
    mode <- normalize_scalar_value(idx_row$mode[[1]])
    recommended <- normalize_scalar_value(idx_row$selected_integration[[1]])
    compare_path <- ""

    if (nrow(layer_candidates) > 1 || identical(mode, "candidate")) {
      compare_df <- candidate_metrics_for_layer(layer_candidates)
      recommended <- compare_df$candidate_id[[1]]
      selected_value <- write_selected_integration_if_needed(selected_file, recommended, layer_id, "candidate")
      compare_path <- file.path(layer_integration_report_dir_04(cfg, layer_id), "subcluster_integration_compare.tsv")
      write_tsv_local(compare_df, compare_path)
      compare_outputs[[paste0("integration_compare_tsv_", layer_id)]] <- build_output_entry(
        compare_path,
        "tsv",
        module_name,
        sprintf("subcluster integration diagnostics for %s", layer_id),
        base_dir = cfg$project_root,
        schema = infer_schema_from_df(compare_df)
      )
    } else if (nrow(layer_candidates) == 1) {
      recommended <- layer_candidates$candidate_id[[1]]
      selected_value <- write_selected_integration_if_needed(selected_file, recommended, layer_id, mode)
    } else {
      selected_value <- ""
    }

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      layer_id = layer_id,
      parent_layer = normalize_scalar_value(idx_row$parent_layer[[1]]),
      mode = mode,
      status = normalize_scalar_value(idx_row$status[[1]]),
      candidate_count = if (nrow(layer_candidates) > 0) nrow(layer_candidates) else suppressWarnings(as.integer(idx_row$candidate_count[[1]])),
      selected_integration = selected_value,
      recommended_integration = recommended,
      selected_integration_file = selected_file,
      integration_compare_tsv = compare_path,
      clustered_key = normalize_scalar_value(idx_row$clustered_key[[1]]),
      clustered_rds = normalize_scalar_value(idx_row$clustered_rds[[1]]),
      reason = normalize_scalar_value(idx_row$reason[[1]]),
      stringsAsFactors = FALSE
    )
  }
}

review_summary <- if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else data.frame(
  layer_id = character(0),
  parent_layer = character(0),
  mode = character(0),
  status = character(0),
  candidate_count = integer(0),
  selected_integration = character(0),
  recommended_integration = character(0),
  selected_integration_file = character(0),
  integration_compare_tsv = character(0),
  clustered_key = character(0),
  clustered_rds = character(0),
  reason = character(0),
  stringsAsFactors = FALSE
)

candidate_layers_count <- length(unique(review_summary$layer_id[
  review_summary$mode == "candidate" & !nzchar(review_summary$clustered_rds)
]))
write_tsv_local(review_summary, cfg$subcluster_review_summary_file)
writeLines(as.character(candidate_layers_count), cfg$subcluster_candidate_layers_count_file, useBytes = TRUE)

review_section_lines <- function(title, df) {
  c(
    sprintf("## %s", title),
    render_markdown_table_local(df),
    ""
  )
}

candidate_review <- review_summary[
  review_summary$mode == "candidate" & !nzchar(review_summary$clustered_rds),
  ,
  drop = FALSE
]
inherited_review <- review_summary[
  review_summary$mode %in% c("inherited", "explicit_single") |
    (review_summary$mode == "candidate" & nzchar(review_summary$clustered_rds)),
  ,
  drop = FALSE
]
skipped_review <- review_summary[
  startsWith(review_summary$mode, "skipped") | startsWith(review_summary$status, "skipped"),
  ,
  drop = FALSE
]

report_path <- file.path(cfg$subcluster_report_dir, "04a_review.md")
ensure_dir(dirname(report_path))
report_lines <- c(
  "# 04a Subcluster Integration Review",
  "",
  sprintf("- review_summary: `%s`", cfg$subcluster_review_summary_file),
  sprintf("- candidate_layers_count: `%s`", candidate_layers_count),
  "",
  review_section_lines("Candidate Layers (Review Needed)", candidate_review),
  review_section_lines("Inherited Or Explicit Single Layers", inherited_review),
  review_section_lines("Skipped Layers", skipped_review),
  "## Review Action",
  "For rows with `mode=candidate`, edit the corresponding `selected_integration_file` if needed, then approve the `subcluster` gate.",
  "",
  "After approval and finalize, `subcluster_candidate_index.tsv` may contain only the selected candidate row for candidate-mode layers. Full multi-candidate diagnostics remain in the per-layer `subcluster_integration_compare.tsv` files listed above."
)
write_markdown_local(report_lines, report_path)

output_entries <- c(
  list(
    subcluster_review_summary_tsv = build_output_entry(cfg$subcluster_review_summary_file, "tsv", module_name, "one row per subcluster review item", base_dir = cfg$project_root, schema = infer_schema_from_df(review_summary)),
    candidate_layers_count_txt = build_output_entry(cfg$subcluster_candidate_layers_count_file, "txt", module_name, "number of layers requiring subcluster integration gate review", base_dir = cfg$project_root),
    report = build_output_entry(report_path, "md", module_name, "subcluster integration review report", base_dir = cfg$project_root)
  ),
  compare_outputs
)

if (file.exists(cfg$module_04a_review_manifest_path)) {
  unlink(cfg$module_04a_review_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_04a_review_manifest_path,
  new_outputs = output_entries,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_04a_manifest = cfg$module_04a_manifest_path,
    candidate_index_tsv = resolve_output_local(manifest_04a, "candidate_index_tsv"),
    subcluster_index_tsv = resolve_output_local(manifest_04a, "subcluster_index_tsv")
  ),
  version = cfg$module_version,
  depends_on = list(module_04a = cfg$module_04a_manifest_path)
)

message("04a_review 完成。candidate_layers_count: ", candidate_layers_count)
