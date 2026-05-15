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

load_required_packages(c("dplyr", "jsonlite"))

cfg <- get_single_script_config_04()
module_name <- "04d_cluster_robustness"
prepare_dirs_04(cfg)
set.seed(cfg$random_seed)

scd_table_dir <- file.path(cfg$table_dir, "04d_cluster_robustness")
scd_report_dir <- file.path(cfg$results_dir, "reports")
ensure_dir(scd_table_dir)
ensure_dir(scd_report_dir)

question_map_path <- env_or_default_03("SCDESIGN3_QUESTION_MAP", file.path(cfg$metadata_dir, "scdesign3_question_map.tsv"))
targets_path <- env_or_default_03("SCDESIGN3_TARGETS_SHEET", file.path(cfg$metadata_dir, "scdesign3_targets.tsv"))
communication_pairs_path <- env_or_default_03("COMMUNICATION_PAIRS_SHEET", file.path(cfg$metadata_dir, "communication_pairs.tsv"))

question_map <- read_tsv_optional(question_map_path)
targets <- read_tsv_optional(targets_path)

required_question_cols <- c(
  "question_id", "section", "question_family", "status",
  "scdesign3_role", "validation_layer", "validation_object",
  "ground_truth_col", "gate_target", "affects_modules",
  "required_before_core_interpretation", "derived_parent_question_id",
  "enabled", "reason"
)
missing_question_cols <- setdiff(required_question_cols, colnames(question_map))
if (length(missing_question_cols) > 0) {
  stop(sprintf("scdesign3_question_map.tsv missing columns: %s", paste(missing_question_cols, collapse = ", ")), call. = FALSE)
}

required_target_cols <- c(
  "target_id", "target_type", "layer_id", "input_object", "truth_col",
  "questions_covered", "n_simulations", "resolution_grid",
  "mixture_design", "primary_metric", "pass_threshold", "warn_threshold",
  "fail_threshold", "max_cells_per_label", "n_hvg", "n_pcs",
  "output_dir", "status"
)
missing_target_cols <- setdiff(required_target_cols, colnames(targets))
if (length(missing_target_cols) > 0) {
  stop(sprintf("scdesign3_targets.tsv missing columns: %s", paste(missing_target_cols, collapse = ", ")), call. = FALSE)
}

empty_df_04d <- function(cols) {
  out <- as.data.frame(setNames(rep(list(character(0)), length(cols)), cols), stringsAsFactors = FALSE)
  out[, cols, drop = FALSE]
}

split_semicolon_04d <- function(x) {
  x <- trimws(as.character(x %||% ""))
  if (!nzchar(x)) {
    return(character(0))
  }
  out <- trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
  out[nzchar(out)]
}

safe_path_exists_04d <- function(path) {
  nzchar(path) && file.exists(path) && isTRUE(file.info(path)$size > 0)
}

manifest_output_or_empty_04d <- function(manifest_path, key) {
  if (!file.exists(manifest_path)) {
    return("")
  }
  out <- tryCatch(resolve_output_local(read_manifest_local(manifest_path), key), error = function(e) "")
  normalize_scalar_value(out)
}

read_layer_status_04d <- function(path) {
  df <- read_tsv_optional(path)
  if (nrow(df) == 0) {
    return(df)
  }
  for (col in c("layer_id", "annotated_rds", "status")) {
    if (!col %in% colnames(df)) df[[col]] <- ""
  }
  df
}

resolve_layer_object_04d <- function(layer_id) {
  layer_id <- normalize_scalar_value(layer_id)
  candidates <- character(0)
  if (layer_id %in% c("panorama", cfg$panorama_layer_id)) {
    candidates <- c(
      manifest_output_or_empty_04d(cfg$module_03d_manifest_path, "annotated_object"),
      cfg$compat_annotated_rds,
      cfg$panorama_annotated_rds
    )
  } else {
    layer_status <- read_layer_status_04d(cfg$layer_status_file)
    if (nrow(layer_status) > 0) {
      hit <- layer_status[layer_status$layer_id == layer_id, , drop = FALSE]
      if (nrow(hit) > 0) {
        candidates <- c(candidates, hit$annotated_rds)
      }
    }
    candidates <- c(
      candidates,
      manifest_output_or_empty_04d(cfg$module_04b_manifest_path, paste0("annotated_", layer_id)),
      file.path(cfg$checkpoint_dir, "layers", layer_id, sprintf("%s_after_annotation.rds", layer_id))
    )
  }
  candidates <- trimws(as.character(candidates))
  candidates[is.na(candidates)] <- ""
  candidates <- unique(candidates[nzchar(candidates)])
  hit <- candidates[vapply(candidates, safe_path_exists_04d, logical(1))]
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  ""
}

target_questions_04d <- function(target) {
  split_semicolon_04d(target$questions_covered[[1]])
}

target_gate_row_04d <- function(target) {
  target_id <- target$target_id[[1]]
  target_type <- target$target_type[[1]]
  layer_id <- target$layer_id[[1]]
  target_status <- target$status[[1]]
  target_questions <- target_questions_04d(target)
  object_path <- resolve_layer_object_04d(layer_id)

  gate_status <- "PLANNED"
  gate_level <- "registered_target"
  interpretation_allowed <- "no"
  recommended_resolution <- ""
  reason <- "Registered in scdesign3_targets.tsv; waiting for runtime validation."

  if (target_status == "planned" || layer_id == "ST" || target_type %in% c("deconv_validation", "spatial_validation")) {
    gate_status <- "PLANNED"
    gate_level <- ifelse(layer_id == "ST", "ST_E2_pending", "planned_target")
    reason <- ifelse(
      layer_id == "ST",
      "ST/deconvolution/spatial support validation is registered here and executed by the future ST-E2 scDesign3 module.",
      "Target is planned until prerequisite cell subtype, split, or velocity inputs are available."
    )
  } else if (!nzchar(object_path)) {
    gate_status <- "PLANNED"
    gate_level <- "waiting_input"
    reason <- sprintf("No formal annotated object resolved for layer `%s`; run upstream 03/04 outputs before scDesign3 simulation.", layer_id)
  } else {
    gate_status <- "WARN"
    gate_level <- "preflight_ready_not_simulated"
    interpretation_allowed <- "exploratory"
    recommended_resolution <- strsplit(target$resolution_grid[[1]], ",", fixed = TRUE)[[1]][[1]] %||% ""
    reason <- sprintf(
      "Input object resolved at `%s`; target is ready for locked scDesign3 fit/simulate/recluster implementation, but this run only materialized the question gate layer.",
      object_path
    )
  }

  data.frame(
    target_id = target_id,
    target_type = target_type,
    layer_id = layer_id,
    questions_covered = paste(target_questions, collapse = ";"),
    gate_status = gate_status,
    gate_level = gate_level,
    interpretation_allowed = interpretation_allowed,
    recommended_resolution = recommended_resolution,
    input_object = target$input_object[[1]],
    resolved_input_path = object_path,
    primary_metric = target$primary_metric[[1]],
    pass_threshold = target$pass_threshold[[1]],
    warn_threshold = target$warn_threshold[[1]],
    fail_threshold = target$fail_threshold[[1]],
    reason = reason,
    stringsAsFactors = FALSE
  )
}

target_gates <- if (nrow(targets) > 0) {
  dplyr::bind_rows(lapply(seq_len(nrow(targets)), function(idx) target_gate_row_04d(targets[idx, , drop = FALSE])))
} else {
  empty_df_04d(c(
    "target_id", "target_type", "layer_id", "questions_covered",
    "gate_status", "gate_level", "interpretation_allowed",
    "recommended_resolution", "input_object", "resolved_input_path",
    "primary_metric", "pass_threshold", "warn_threshold",
    "fail_threshold", "reason"
  ))
}

question_to_target <- if (nrow(targets) > 0) {
  rows <- list()
  for (idx in seq_len(nrow(targets))) {
    target <- targets[idx, , drop = FALSE]
    qids <- target_questions_04d(target)
    if (length(qids) == 0) {
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(
      question_id = qids,
      target_id = target$target_id[[1]],
      target_type = target$target_type[[1]],
      layer_id = target$layer_id[[1]],
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) > 0) dplyr::bind_rows(rows) else empty_df_04d(c("question_id", "target_id", "target_type", "layer_id"))
} else {
  empty_df_04d(c("question_id", "target_id", "target_type", "layer_id"))
}

question_gate_status <- question_map %>%
  dplyr::left_join(question_to_target, by = "question_id") %>%
  dplyr::left_join(target_gates[, c("target_id", "gate_status", "gate_level", "interpretation_allowed", "reason"), drop = FALSE], by = "target_id", suffix = c("", ".target")) %>%
  dplyr::mutate(
    gate_status = dplyr::case_when(
      scdesign3_role == "not_applicable" ~ "NOT_APPLICABLE",
      scdesign3_role == "derived_from_validated_parent" ~ "DERIVED",
      scdesign3_role == "planned_waiting_input" | enabled == "planned" ~ "PLANNED",
      is.na(gate_status) | !nzchar(gate_status) ~ "PLANNED",
      TRUE ~ gate_status
    ),
    gate_level = dplyr::case_when(
      scdesign3_role == "not_applicable" ~ "context_only",
      scdesign3_role == "derived_from_validated_parent" ~ "derived_parent",
      scdesign3_role == "planned_waiting_input" | enabled == "planned" ~ "planned_input",
      is.na(gate_level) | !nzchar(gate_level) ~ "registered_question",
      TRUE ~ gate_level
    ),
    interpretation_allowed = dplyr::case_when(
      gate_status == "PASS" ~ "core",
      gate_status == "WARN" ~ "exploratory",
      gate_status == "NOT_APPLICABLE" ~ "context",
      gate_status == "DERIVED" ~ "inherit_parent",
      TRUE ~ "no"
    ),
    reason = dplyr::case_when(
      gate_status == "NOT_APPLICABLE" ~ reason,
      gate_status == "DERIVED" ~ paste0("Derived from parent question(s): ", derived_parent_question_id),
      !is.na(reason.target) & nzchar(reason.target) ~ reason.target,
      TRUE ~ reason
    ),
    scdesign3_target_id = target_id,
    question_status = status
  ) %>%
  dplyr::select(
    question_id, section, question_family, question_status,
    scdesign3_role, scdesign3_target_id, gate_status, gate_level,
    affects_modules, interpretation_allowed, reason,
    validation_layer, ground_truth_col, gate_target,
    required_before_core_interpretation, derived_parent_question_id,
    enabled
  )

cluster_metrics <- target_gates %>%
  dplyr::filter(target_type %in% c("cluster_robustness", "communication_robustness")) %>%
  dplyr::transmute(
    target_id, target_type, layer_id,
    primary_metric,
    metric_value = "pending_engine",
    gate_status, gate_level,
    engine_status = "not_run_metadata_gate",
    reason
  )

subtype_recovery_metrics <- empty_df_04d(c(
  "target_id", "layer_id", "truth_label", "best_matching_cluster",
  "jaccard", "precision", "recall", "f1", "status", "reason"
))

composition_recovery_metrics <- target_gates %>%
  dplyr::filter(target_type == "composition_robustness") %>%
  dplyr::transmute(
    target_id, question_id = questions_covered,
    metric = "composition_RMSE", metric_value = "pending_engine",
    status = gate_status, reason
  )

trajectory_recovery_metrics <- target_gates %>%
  dplyr::filter(target_type == "trajectory_robustness") %>%
  dplyr::transmute(
    target_id, question_id = questions_covered,
    metric = "trajectory_order_accuracy", metric_value = "pending_engine",
    status = gate_status, reason
  )

resolution_recommendation <- target_gates %>%
  dplyr::transmute(
    target_id, recommended_resolution, gate_status, reason
  )

simulation_manifest <- if (nrow(targets) > 0) {
  targets %>%
    dplyr::transmute(
      target_id,
      simulation_id = paste0(target_id, "_SIM_REGISTERED"),
      simulation_type = dplyr::case_when(
        target_type %in% c("deconv_validation", "spatial_validation") ~ "synthetic_spots",
        target_type == "composition_robustness" ~ "synthetic_cells_composition",
        target_type == "trajectory_robustness" ~ "synthetic_cells_trajectory",
        target_type == "communication_robustness" ~ "perturbation",
        TRUE ~ "synthetic_cells"
      ),
      engine = "registered_not_executed",
      status = status,
      output_path = "",
      reason = "Simulation design registered; runtime engine not executed in metadata gate materialization."
    )
} else {
  empty_df_04d(c("target_id", "simulation_id", "simulation_type", "engine", "status", "output_path", "reason"))
}

communication_scdesign3_gate <- function() {
  pairs <- read_tsv_optional(communication_pairs_path)
  if (nrow(pairs) == 0) {
    return(empty_df_04d(c(
      "pair_id", "question_id", "sender_label", "receiver_label",
      "sender_recovery_status", "receiver_recovery_status",
      "split_gate_status", "scdesign3_gate_status",
      "allowed_for_cellchat", "allowed_for_nichenet", "reason"
    )))
  }
  for (col in c("pair_id", "source_question_id", "sender", "receiver", "communication_mode")) {
    if (!col %in% colnames(pairs)) pairs[[col]] <- ""
  }
  pairs %>%
    dplyr::left_join(question_gate_status[, c("question_id", "gate_status", "interpretation_allowed", "reason"), drop = FALSE], by = c("source_question_id" = "question_id")) %>%
    dplyr::mutate(
      question_id = source_question_id,
      sender_label = sender,
      receiver_label = receiver,
      sender_recovery_status = gate_status,
      receiver_recovery_status = gate_status,
      split_gate_status = dplyr::if_else(communication_mode == "condition_split", gate_status, "not_split"),
      scdesign3_gate_status = dplyr::coalesce(gate_status, "PLANNED"),
      allowed_for_cellchat = dplyr::if_else(scdesign3_gate_status %in% c("PASS", "WARN", "DERIVED"), "yes", "no"),
      allowed_for_nichenet = dplyr::if_else(scdesign3_gate_status %in% c("PASS", "WARN"), "yes", "no"),
      reason = dplyr::coalesce(reason, "No matching scDesign3 question gate.")
    ) %>%
    dplyr::select(
      pair_id, question_id, sender_label, receiver_label,
      sender_recovery_status, receiver_recovery_status,
      split_gate_status, scdesign3_gate_status,
      allowed_for_cellchat, allowed_for_nichenet, reason
    )
}

communication_gate_df <- communication_scdesign3_gate()

paths <- list(
  cluster_metrics_tsv = file.path(scd_table_dir, "cluster_robustness_metrics.tsv"),
  subtype_recovery_tsv = file.path(scd_table_dir, "subtype_recovery_metrics.tsv"),
  composition_recovery_tsv = file.path(scd_table_dir, "composition_recovery_metrics.tsv"),
  trajectory_recovery_tsv = file.path(scd_table_dir, "trajectory_recovery_metrics.tsv"),
  target_gate_status_tsv = file.path(scd_table_dir, "target_gate_status.tsv"),
  question_gate_status_tsv = file.path(scd_table_dir, "question_gate_status.tsv"),
  question_to_target_map_tsv = file.path(scd_table_dir, "question_to_target_map.tsv"),
  resolution_recommendation_tsv = file.path(scd_table_dir, "resolution_recommendation.tsv"),
  simulation_manifest_tsv = file.path(scd_table_dir, "simulation_manifest.tsv"),
  communication_scdesign3_gate_tsv = file.path(scd_table_dir, "communication_scdesign3_gate.tsv"),
  all_questions_status_tsv = file.path(scd_table_dir, "scdesign3_all_questions_status.tsv"),
  all_questions_status_root_tsv = file.path(cfg$table_dir, "scdesign3_all_questions_status.tsv"),
  report_md = file.path(cfg$subcluster_report_dir, "04d_cluster_robustness.md"),
  all_questions_report_md = file.path(scd_report_dir, "scdesign3_all_questions_validation_report.md")
)

write_tsv_local(cluster_metrics, paths$cluster_metrics_tsv)
write_tsv_local(subtype_recovery_metrics, paths$subtype_recovery_tsv)
write_tsv_local(composition_recovery_metrics, paths$composition_recovery_tsv)
write_tsv_local(trajectory_recovery_metrics, paths$trajectory_recovery_tsv)
write_tsv_local(target_gates, paths$target_gate_status_tsv)
write_tsv_local(question_gate_status, paths$question_gate_status_tsv)
write_tsv_local(question_to_target, paths$question_to_target_map_tsv)
write_tsv_local(resolution_recommendation, paths$resolution_recommendation_tsv)
write_tsv_local(simulation_manifest, paths$simulation_manifest_tsv)
write_tsv_local(communication_gate_df, paths$communication_scdesign3_gate_tsv)
write_tsv_local(question_gate_status, paths$all_questions_status_tsv)
write_tsv_local(question_gate_status, paths$all_questions_status_root_tsv)

section_summary <- question_gate_status %>%
  dplyr::count(section, gate_status, name = "question_n") %>%
  dplyr::arrange(section, gate_status)

role_summary <- question_gate_status %>%
  dplyr::count(scdesign3_role, gate_status, name = "question_n") %>%
  dplyr::arrange(scdesign3_role, gate_status)

target_summary <- target_gates %>%
  dplyr::count(target_type, gate_status, name = "target_n") %>%
  dplyr::arrange(target_type, gate_status)

section_blocks <- list()
for (section in sort(unique(question_gate_status$section))) {
  section_df <- question_gate_status[question_gate_status$section == section, c(
    "question_id", "question_family", "scdesign3_role",
    "gate_status", "interpretation_allowed", "scdesign3_target_id",
    "reason"
  ), drop = FALSE]
  section_blocks[[paste0("Section ", section)]] <- render_markdown_table_local(section_df)
}

report_lines <- build_report_lines_v04(
  title = "04d scDesign3 Question Gate",
  header_bullets = c(
    "scDesign3 is now represented as a full question-driven validation layer.",
    sprintf("question_map_rows: `%s`", nrow(question_map)),
    sprintf("target_rows: `%s`", nrow(targets)),
    "This run materializes gate metadata and preflight status; it does not claim completed scDesign3 simulation metrics."
  ),
  key_files = list(
    scdesign3_question_map = question_map_path,
    scdesign3_targets = targets_path,
    all_questions_status = paths$all_questions_status_root_tsv,
    question_gate_status = paths$question_gate_status_tsv,
    communication_scdesign3_gate = paths$communication_scdesign3_gate_tsv,
    report = paths$all_questions_report_md
  ),
  review_focus = c(
    "Confirm every analysis question appears in scdesign3_all_questions_status.tsv.",
    "Treat WARN rows as preflight-ready but not simulated unless target metrics later show PASS.",
    "Treat PLANNED rows as not available for core interpretation until prerequisite inputs exist."
  ),
  extra_sections = c(
    list(
      "Section Gate Summary" = render_markdown_table_local(section_summary),
      "Role Gate Summary" = render_markdown_table_local(role_summary),
      "Target Gate Summary" = render_markdown_table_local(target_summary)
    ),
    section_blocks
  )
)
ensure_dir(dirname(paths$report_md))
ensure_dir(dirname(paths$all_questions_report_md))
write_markdown_local(report_lines, paths$report_md)
write_markdown_local(report_lines, paths$all_questions_report_md)

if (file.exists(cfg$module_04d_manifest_path)) {
  unlink(cfg$module_04d_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_04d_manifest_path,
  new_outputs = list(
    metrics_tsv = build_output_entry(paths$cluster_metrics_tsv, "tsv", module_name, "scDesign3 cluster/communication target preflight metrics", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_metrics)),
    subtype_recovery_metrics_tsv = build_output_entry(paths$subtype_recovery_tsv, "tsv", module_name, "subtype recovery metrics placeholder schema", base_dir = cfg$project_root, schema = infer_schema_from_df(subtype_recovery_metrics)),
    composition_recovery_metrics_tsv = build_output_entry(paths$composition_recovery_tsv, "tsv", module_name, "composition recovery metrics placeholder schema", base_dir = cfg$project_root, schema = infer_schema_from_df(composition_recovery_metrics)),
    trajectory_recovery_metrics_tsv = build_output_entry(paths$trajectory_recovery_tsv, "tsv", module_name, "trajectory recovery metrics placeholder schema", base_dir = cfg$project_root, schema = infer_schema_from_df(trajectory_recovery_metrics)),
    target_gate_status_tsv = build_output_entry(paths$target_gate_status_tsv, "tsv", module_name, "one row per scDesign3 target gate", base_dir = cfg$project_root, schema = infer_schema_from_df(target_gates)),
    question_gate_status_tsv = build_output_entry(paths$question_gate_status_tsv, "tsv", module_name, "one row per analysis question scDesign3 gate", base_dir = cfg$project_root, schema = infer_schema_from_df(question_gate_status)),
    scdesign3_all_questions_status_tsv = build_output_entry(paths$all_questions_status_root_tsv, "tsv", module_name, "global scDesign3 status table for all questions", base_dir = cfg$project_root, schema = infer_schema_from_df(question_gate_status)),
    communication_scdesign3_gate_tsv = build_output_entry(paths$communication_scdesign3_gate_tsv, "tsv", module_name, "communication pair scDesign3 gate projection", base_dir = cfg$project_root, schema = infer_schema_from_df(communication_gate_df)),
    question_to_target_map_tsv = build_output_entry(paths$question_to_target_map_tsv, "tsv", module_name, "question-to-target expansion", base_dir = cfg$project_root, schema = infer_schema_from_df(question_to_target)),
    resolution_recommendation_tsv = build_output_entry(paths$resolution_recommendation_tsv, "tsv", module_name, "target resolution recommendations or pending reasons", base_dir = cfg$project_root, schema = infer_schema_from_df(resolution_recommendation)),
    simulation_manifest_tsv = build_output_entry(paths$simulation_manifest_tsv, "tsv", module_name, "registered simulation designs and runtime status", base_dir = cfg$project_root, schema = infer_schema_from_df(simulation_manifest)),
    report = build_output_entry(paths$report_md, "md", module_name, "04d scDesign3 question gate report", base_dir = cfg$project_root),
    all_questions_report = build_output_entry(paths$all_questions_report_md, "md", module_name, "A-I scDesign3 all-questions validation report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    scdesign3_question_map = question_map_path,
    scdesign3_targets = targets_path,
    communication_pairs = communication_pairs_path,
    module_04c_manifest = cfg$module_04c_manifest_path,
    module_04b_manifest = cfg$module_04b_manifest_path
  ),
  version = cfg$module_version,
  depends_on = list(
    metadata = list(question_map = question_map_path, targets = targets_path),
    module_04c = cfg$module_04c_manifest_path,
    module_04b = cfg$module_04b_manifest_path
  )
)

message("04d scDesign3 question gate completed. questions: ", nrow(question_gate_status), "; targets: ", nrow(target_gates))
