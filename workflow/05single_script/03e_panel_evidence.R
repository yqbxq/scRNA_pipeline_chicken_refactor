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
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "annotation_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))

load_required_packages(c("dplyr", "tidyr", "tibble", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03e_panel_evidence"
prepare_dirs_03(cfg)

panorama_spec <- panorama_layer_spec(cfg)
mode <- cfg$panel_evidence_mode
if (!mode %in% c("off", "audit_only", "validate_manual")) {
  stop(sprintf("Unsupported PANEL_EVIDENCE_MODE: %s", mode), call. = FALSE)
}

policy <- data.frame(
  marker_role = c("core", "supporting", "review_only", "shared_risk", "exclude"),
  scoring_window = c(200L, 500L, 0L, 0L, 0L),
  audit_window = c(500L, 1000L, 2000L, 1000L, 1000L),
  positive_scoring = c(TRUE, TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

empty_hit <- data.frame(
  cluster = character(0),
  celltype = character(0),
  gene = character(0),
  marker_role = character(0),
  rank_specificity = integer(0),
  scoring_hit = logical(0),
  audit_hit = logical(0),
  panel_file = character(0),
  stringsAsFactors = FALSE
)
empty_evidence <- data.frame(
  cluster = character(0),
  celltype = character(0),
  annotation_score = numeric(0),
  core_scoring_n = integer(0),
  supporting_scoring_n = integer(0),
  shared_risk_n = integer(0),
  exclude_n = integer(0),
  review_only_n = integer(0),
  scoring_genes = character(0),
  audit_genes = character(0),
  evidence_call = character(0),
  stringsAsFactors = FALSE
)
empty_conflict <- data.frame(
  cluster = character(0),
  conflict_type = character(0),
  evidence = character(0),
  stringsAsFactors = FALSE
)
empty_sensitivity <- data.frame(
  cluster = character(0),
  window = integer(0),
  best_celltype = character(0),
  best_score = numeric(0),
  support_call = character(0),
  stringsAsFactors = FALSE
)

panel_hit_long <- empty_hit
panel_candidate_evidence <- empty_evidence
panel_conflict <- empty_conflict
panel_window_sensitivity <- empty_sensitivity
manual_validation <- data.frame(stringsAsFactors = FALSE)

manifest_03c3 <- read_manifest_local(cfg$module_03c3_manifest_path)
markers_scored <- read_tsv_optional(resolve_output_local(manifest_03c3, "markers_scored_tsv"))

panel_df <- data.frame(stringsAsFactors = FALSE)
if (!identical(mode, "off")) {
  panel_df <- read_marker_panel_rows(cfg, layer_id = panorama_spec$layer_id)
}

if (!identical(mode, "off") && nrow(panel_df) > 0 && nrow(markers_scored) > 0) {
  for (col in c("cluster", "gene", "rank_specificity")) {
    if (!col %in% colnames(markers_scored)) {
      stop(sprintf("markers_scored.tsv missing required column: %s", col), call. = FALSE)
    }
  }
  markers_scored$rank_specificity <- suppressWarnings(as.integer(markers_scored$rank_specificity))
  panel_df$gene <- as.character(panel_df$gene)
  panel_df$celltype <- as.character(panel_df$celltype)
  panel_df$marker_role <- if ("marker_role" %in% colnames(panel_df)) panel_df$marker_role else "core"
  panel_df <- dplyr::left_join(panel_df, policy, by = "marker_role")
  panel_df$audit_window[is.na(panel_df$audit_window)] <- 1000L
  panel_df$scoring_window[is.na(panel_df$scoring_window)] <- 0L
  panel_df$positive_scoring[is.na(panel_df$positive_scoring)] <- FALSE

  merged <- dplyr::inner_join(
    markers_scored[, c("cluster", "gene", "rank_specificity"), drop = FALSE],
    panel_df,
    by = "gene"
  ) %>%
    dplyr::filter(rank_specificity <= audit_window) %>%
    dplyr::mutate(
      scoring_hit = positive_scoring & rank_specificity <= scoring_window,
      audit_hit = TRUE
    )
  panel_hit_long <- merged %>%
    dplyr::select(dplyr::any_of(c(
      "cluster", "celltype", "gene", "marker_role", "rank_specificity",
      "scoring_hit", "audit_hit", "panel_file", "evidence_source", "evidence_note"
    ))) %>%
    dplyr::arrange(cluster, celltype, rank_specificity, gene)

  panel_candidate_evidence <- panel_hit_long %>%
    dplyr::group_by(cluster, celltype) %>%
    dplyr::summarise(
      annotation_score = 4L * sum(marker_role == "core" & scoring_hit, na.rm = TRUE) +
        sum(marker_role == "supporting" & scoring_hit, na.rm = TRUE),
      core_scoring_n = sum(marker_role == "core" & scoring_hit, na.rm = TRUE),
      supporting_scoring_n = sum(marker_role == "supporting" & scoring_hit, na.rm = TRUE),
      shared_risk_n = sum(marker_role == "shared_risk", na.rm = TRUE),
      exclude_n = sum(marker_role == "exclude", na.rm = TRUE),
      review_only_n = sum(marker_role == "review_only", na.rm = TRUE),
      scoring_genes = paste(unique(gene[scoring_hit]), collapse = ","),
      audit_genes = paste(unique(gene), collapse = ","),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      evidence_call = dplyr::case_when(
        annotation_score > 0 & exclude_n > 0 ~ "positive_with_exclude_conflict",
        annotation_score > 0 ~ "positive_panel_evidence",
        shared_risk_n + exclude_n + review_only_n > 0 ~ "audit_or_conflict_only",
        TRUE ~ "no_panel_support"
      )
    ) %>%
    dplyr::arrange(cluster, dplyr::desc(annotation_score), dplyr::desc(core_scoring_n), celltype)

  conflict_rows <- list()
  for (cluster_id in unique(panel_candidate_evidence$cluster)) {
    e <- panel_candidate_evidence[panel_candidate_evidence$cluster == cluster_id, , drop = FALSE]
    positives <- e[e$annotation_score > 0, , drop = FALSE]
    if (nrow(positives) >= 2) {
      conflict_rows[[length(conflict_rows) + 1]] <- data.frame(
        cluster = cluster_id,
        conflict_type = "multiple_positive_panel_candidates",
        evidence = paste(sprintf("%s(score=%s)", positives$celltype, positives$annotation_score), collapse = "; "),
        stringsAsFactors = FALSE
      )
    }
    exclude_hit <- e[e$exclude_n > 0 & e$annotation_score > 0, , drop = FALSE]
    if (nrow(exclude_hit) > 0) {
      conflict_rows[[length(conflict_rows) + 1]] <- data.frame(
        cluster = cluster_id,
        conflict_type = "positive_candidate_has_exclude_marker",
        evidence = paste(exclude_hit$celltype, collapse = ","),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(conflict_rows) > 0) {
    panel_conflict <- dplyr::bind_rows(conflict_rows)
  }

  windows <- c(20L, 50L, 100L, 200L, 500L, 1000L, 2000L)
  sensitivity_rows <- list()
  for (cluster_id in unique(markers_scored$cluster)) {
    for (window in windows) {
      genes <- markers_scored %>%
        dplyr::filter(cluster == cluster_id, rank_specificity <= window) %>%
        dplyr::pull(gene) %>%
        unique()
      hit <- panel_df[panel_df$gene %in% genes & panel_df$positive_scoring, , drop = FALSE]
      if (nrow(hit) == 0) {
        best_celltype <- ""
        best_score <- 0
      } else {
        score_df <- hit %>%
          dplyr::group_by(celltype) %>%
          dplyr::summarise(best_score = 4L * sum(marker_role == "core") + sum(marker_role == "supporting"), .groups = "drop") %>%
          dplyr::arrange(dplyr::desc(best_score), celltype)
        best_celltype <- score_df$celltype[[1]]
        best_score <- score_df$best_score[[1]]
      }
      sensitivity_rows[[length(sensitivity_rows) + 1]] <- data.frame(
        cluster = cluster_id,
        window = window,
        best_celltype = best_celltype,
        best_score = best_score,
        support_call = ifelse(best_score > 0, best_celltype, "Uncertain"),
        stringsAsFactors = FALSE
      )
    }
  }
  panel_window_sensitivity <- dplyr::bind_rows(sensitivity_rows)

  if (identical(mode, "validate_manual")) {
    manual <- read_tsv_optional(cfg$manual_annotation_file)
    if (nrow(manual) > 0 && all(c("layer_id", "cluster", "manual_annotation") %in% colnames(manual))) {
      manual <- manual[manual$layer_id %in% c(panorama_spec$layer_id, "*", ""), , drop = FALSE]
      best_panel <- panel_candidate_evidence %>%
        dplyr::group_by(cluster) %>%
        dplyr::slice_max(order_by = annotation_score, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()
      manual_validation <- dplyr::left_join(manual, best_panel, by = "cluster") %>%
        dplyr::mutate(
          panel_supports_manual = manual_annotation == celltype & annotation_score > 0,
          panel_validation_flag = dplyr::case_when(
            panel_supports_manual ~ "manual_supported_by_panel",
            annotation_score > 0 & manual_annotation != celltype ~ "manual_panel_mismatch",
            TRUE ~ "no_positive_panel_support"
          )
        )
    }
  }
}

panel_window_policy_tsv <- file.path(cfg$panel_evidence_table_dir_layer, "panel_window_policy.tsv")
panel_hit_long_tsv <- file.path(cfg$panel_evidence_table_dir_layer, "panel_hit_long.tsv")
panel_conflict_tsv <- file.path(cfg$panel_evidence_table_dir_layer, "panel_conflict.tsv")
panel_window_sensitivity_tsv <- file.path(cfg$panel_evidence_table_dir_layer, "panel_window_sensitivity.tsv")
panel_candidate_evidence_tsv <- file.path(cfg$panel_evidence_table_dir_layer, "panel_candidate_evidence.tsv")
manual_validation_tsv <- file.path(cfg$panel_evidence_table_dir_layer, "manual_annotation_panel_validation.tsv")
report_path <- file.path(cfg$panel_evidence_report_dir_layer, "report.md")

write_tsv_local(policy, panel_window_policy_tsv)
write_tsv_local(panel_hit_long, panel_hit_long_tsv)
write_tsv_local(panel_conflict, panel_conflict_tsv)
write_tsv_local(panel_window_sensitivity, panel_window_sensitivity_tsv)
write_tsv_local(panel_candidate_evidence, panel_candidate_evidence_tsv)
write_tsv_local(manual_validation, manual_validation_tsv)

report_lines <- c(
  "# 03e Optional Panel Evidence",
  "",
  sprintf("- layer_id: `%s`", panorama_spec$layer_id),
  sprintf("- PANEL_EVIDENCE_MODE: `%s`", mode),
  sprintf("- panel_rows: `%s`", nrow(panel_df)),
  sprintf("- panel_hit_rows: `%s`", nrow(panel_hit_long)),
  sprintf("- panel_conflict_rows: `%s`", nrow(panel_conflict)),
  "",
  "## Candidate Evidence",
  render_markdown_table_local(head(panel_candidate_evidence[, intersect(c(
    "cluster", "celltype", "annotation_score", "core_scoring_n",
    "supporting_scoring_n", "shared_risk_n", "exclude_n", "evidence_call"
  ), colnames(panel_candidate_evidence)), drop = FALSE], 80)),
  "",
  "This stage is audit-only. It never writes final_annotation and never changes the Seurat object."
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_03e_panel_evidence_manifest_path,
  new_outputs = list(
    panel_window_policy_tsv = build_output_entry(panel_window_policy_tsv, "tsv", module_name, "role-aware panel rank window policy", base_dir = cfg$project_root, schema = infer_schema_from_df(policy)),
    panel_hit_long_tsv = build_output_entry(panel_hit_long_tsv, "tsv", module_name, "one row per cluster/panel marker hit", base_dir = cfg$project_root, schema = infer_schema_from_df(panel_hit_long)),
    panel_conflict_tsv = build_output_entry(panel_conflict_tsv, "tsv", module_name, "panel conflict evidence for manual review", base_dir = cfg$project_root, schema = infer_schema_from_df(panel_conflict)),
    panel_window_sensitivity_tsv = build_output_entry(panel_window_sensitivity_tsv, "tsv", module_name, "panel support across rank windows", base_dir = cfg$project_root, schema = infer_schema_from_df(panel_window_sensitivity)),
    panel_candidate_evidence_tsv = build_output_entry(panel_candidate_evidence_tsv, "tsv", module_name, "cluster by candidate panel evidence", base_dir = cfg$project_root, schema = infer_schema_from_df(panel_candidate_evidence)),
    manual_annotation_panel_validation_tsv = build_output_entry(manual_validation_tsv, "tsv", module_name, "manual annotation validation against panel evidence when requested", base_dir = cfg$project_root, schema = infer_schema_from_df(manual_validation)),
    report = build_output_entry(report_path, "md", module_name, "human-readable optional panel evidence report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(module_03c3_manifest = cfg$module_03c3_manifest_path, marker_panel_dir = cfg$marker_panel_dir),
  version = cfg$module_version,
  depends_on = list(module_03c3 = cfg$module_03c3_manifest_path, module_03d_marker_risk = cfg$module_03d_marker_risk_manifest_path)
)

message("03e complete. Panel evidence mode: ", mode)
