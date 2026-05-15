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
source_utf8(file.path(.script_dir, "helpers", "scdesign3_engine_helpers.R"))

load_required_packages(c("dplyr", "jsonlite"))

cfg <- get_single_script_config_04()
module_name <- "04f_scdesign3_finalize"
prepare_dirs_04(cfg)

table_dir <- file.path(cfg$table_dir, module_name)
report_dir <- file.path(cfg$results_dir, "reports")
manifest_path <- file.path(cfg$manifest_dir, module_name, "_manifest.json")
ensure_dir(table_dir)
ensure_dir(report_dir)
ensure_dir(dirname(manifest_path))

target_gate_pre_path <- file.path(cfg$table_dir, "04d_cluster_robustness", "target_gate_status.tsv")
question_gate_pre_path <- file.path(cfg$table_dir, "04d_cluster_robustness", "question_gate_status.tsv")
communication_pairs_path <- env_or_default_03("COMMUNICATION_PAIRS_SHEET", file.path(cfg$metadata_dir, "communication_pairs.tsv"))
engine_metrics_path <- file.path(cfg$table_dir, "04e_scdesign3_engine", "target_metrics.tsv")
engine_status_path <- file.path(cfg$table_dir, "04e_scdesign3_engine", "engine_status.tsv")

target_gates_old <- read_tsv_optional(target_gate_pre_path)
question_gates_old <- read_tsv_optional(question_gate_pre_path)
target_metrics <- read_tsv_optional(engine_metrics_path)
engine_status_df <- read_tsv_optional(engine_status_path)

if (nrow(target_gates_old) == 0 || !"target_id" %in% colnames(target_gates_old)) {
  stop(sprintf("Missing 04d target gate status: %s", target_gate_pre_path), call. = FALSE)
}
if (nrow(question_gates_old) == 0 || !"question_id" %in% colnames(question_gates_old)) {
  stop(sprintf("Missing 04d question gate status: %s", question_gate_pre_path), call. = FALSE)
}

for (col in c("target_id", "status", "gate_status", "primary_metric_value", "median_ARI", "median_NMI", "median_max_jaccard", "best_resolution", "n_simulations_done", "reason")) {
  if (!col %in% colnames(target_metrics)) target_metrics[[col]] <- character(nrow(target_metrics))
}
for (col in c("target_id", "fit_ok", "sim_ok", "score_ok", "status", "reason")) {
  if (!col %in% colnames(engine_status_df)) engine_status_df[[col]] <- character(nrow(engine_status_df))
}

engine_projection <- target_metrics %>%
  dplyr::select(
    target_id,
    engine_status = status,
    engine_gate_status = gate_status,
    engine_primary_metric_value = primary_metric_value,
    engine_median_ARI = median_ARI,
    engine_median_NMI = median_NMI,
    engine_median_max_jaccard = median_max_jaccard,
    engine_best_resolution = best_resolution,
    engine_n_simulations_done = n_simulations_done,
    engine_reason = reason
  ) %>%
  dplyr::left_join(
    engine_status_df %>%
      dplyr::select(target_id, fit_ok, sim_ok, score_ok, runtime_status = status, runtime_reason = reason),
    by = "target_id"
  )

for (col in c("gate_status", "gate_level", "interpretation_allowed", "reason")) {
  if (!col %in% colnames(target_gates_old)) target_gates_old[[col]] <- ""
}

target_gates_new <- target_gates_old %>%
  dplyr::mutate(pre_engine_gate_status = gate_status) %>%
  dplyr::left_join(engine_projection, by = "target_id") %>%
  dplyr::mutate(
    engine_applied = dplyr::if_else(engine_status == "ok" & engine_gate_status %in% c("PASS", "WARN", "FAIL"), "yes", "no"),
    gate_status = dplyr::case_when(
      engine_applied == "yes" ~ engine_gate_status,
      TRUE ~ gate_status
    ),
    gate_level = dplyr::case_when(
      engine_applied == "yes" ~ "engine_executed",
      TRUE ~ gate_level
    ),
    interpretation_allowed = dplyr::case_when(
      gate_status == "PASS" ~ "core",
      gate_status == "WARN" ~ "exploratory",
      gate_status == "FAIL" ~ "rejected",
      TRUE ~ interpretation_allowed
    ),
    reason = dplyr::case_when(
      engine_applied == "yes" ~ sprintf(
        "scDesign3 engine executed: primary_metric_value=%s; median_ARI=%s; median_NMI=%s; best_resolution=%s; pre_engine_gate_status=%s; post_engine_gate_status=%s.",
        engine_primary_metric_value,
        engine_median_ARI,
        engine_median_NMI,
        engine_best_resolution,
        pre_engine_gate_status,
        gate_status
      ),
      !is.na(engine_status) & nzchar(engine_status) & engine_status != "ok" ~ sprintf(
        "%s Engine metrics not applied because 04e status=%s: %s",
        reason,
        engine_status,
        dplyr::coalesce(engine_reason, "")
      ),
      TRUE ~ reason
    )
  )

target_projection <- target_gates_new %>%
  dplyr::select(
    scdesign3_target_id = target_id,
    target_gate_status = gate_status,
    target_gate_level = gate_level,
    target_interpretation_allowed = interpretation_allowed,
    target_reason = reason,
    engine_applied,
    target_pre_engine_gate_status = pre_engine_gate_status,
    engine_primary_metric_value,
    engine_median_ARI,
    engine_median_NMI,
    engine_best_resolution,
    engine_n_simulations_done
  )

for (col in c("scdesign3_target_id", "gate_status", "gate_level", "interpretation_allowed", "reason", "scdesign3_role", "enabled")) {
  if (!col %in% colnames(question_gates_old)) question_gates_old[[col]] <- ""
}

question_gates_new <- question_gates_old %>%
  dplyr::mutate(pre_engine_gate_status = gate_status) %>%
  dplyr::left_join(target_projection, by = "scdesign3_target_id") %>%
  dplyr::mutate(
    engine_applied = dplyr::coalesce(engine_applied, "no"),
    gate_status = dplyr::case_when(
      scdesign3_role == "not_applicable" ~ gate_status,
      scdesign3_role == "derived_from_validated_parent" ~ gate_status,
      scdesign3_role == "planned_waiting_input" | enabled == "planned" ~ gate_status,
      engine_applied == "yes" ~ target_gate_status,
      TRUE ~ gate_status
    ),
    gate_level = dplyr::case_when(
      engine_applied == "yes" ~ target_gate_level,
      TRUE ~ gate_level
    ),
    interpretation_allowed = dplyr::case_when(
      gate_status == "PASS" ~ "core",
      gate_status == "WARN" ~ "exploratory",
      gate_status == "FAIL" ~ "rejected",
      gate_status == "NOT_APPLICABLE" ~ "context",
      gate_status == "DERIVED" ~ "inherit_parent",
      TRUE ~ interpretation_allowed
    ),
    reason = dplyr::case_when(
      engine_applied == "yes" & !is.na(target_reason) & nzchar(target_reason) ~ target_reason,
      TRUE ~ reason
    )
  )

drop_cols <- intersect(
  c("target_gate_status", "target_gate_level", "target_interpretation_allowed", "target_reason"),
  colnames(question_gates_new)
)
question_gates_new <- question_gates_new[, setdiff(colnames(question_gates_new), drop_cols), drop = FALSE]

communication_scdesign3_gate <- function(question_gate_status) {
  pairs <- read_tsv_optional(communication_pairs_path)
  if (nrow(pairs) == 0) {
    return(scd_empty_df(c(
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

communication_gate_df <- communication_scdesign3_gate(question_gates_new)

paths <- list(
  target_gate_status_post_engine_tsv = file.path(table_dir, "target_gate_status_post_engine.tsv"),
  question_gate_status_post_engine_tsv = file.path(table_dir, "question_gate_status_post_engine.tsv"),
  communication_scdesign3_gate_post_engine_tsv = file.path(table_dir, "communication_scdesign3_gate_post_engine.tsv"),
  all_questions_status_tsv = file.path(cfg$table_dir, "scdesign3_all_questions_status.tsv"),
  report_md = file.path(report_dir, "04f_scdesign3_finalize.md")
)

write_tsv_local(target_gates_new, paths$target_gate_status_post_engine_tsv)
write_tsv_local(question_gates_new, paths$question_gate_status_post_engine_tsv)
write_tsv_local(communication_gate_df, paths$communication_scdesign3_gate_post_engine_tsv)
write_tsv_local(question_gates_new, paths$all_questions_status_tsv)

target_summary <- target_gates_new %>%
  dplyr::count(target_type, gate_status, engine_applied, name = "target_n") %>%
  dplyr::arrange(target_type, gate_status, engine_applied)

question_summary <- question_gates_new %>%
  dplyr::count(section, gate_status, name = "question_n") %>%
  dplyr::arrange(section, gate_status)

engine_applied_targets <- target_gates_new %>%
  dplyr::filter(engine_applied == "yes") %>%
  dplyr::select(target_id, target_type, pre_engine_gate_status, gate_status, engine_primary_metric_value, engine_median_ARI, engine_median_NMI, engine_best_resolution, engine_n_simulations_done)

report_lines <- build_report_lines_v04(
  title = "04f scDesign3 Finalize",
  header_bullets = c(
    "04f upgrades 04d preflight gates only when 04e produced scoreable engine metrics.",
    sprintf("engine_applied_targets: `%s`", nrow(engine_applied_targets)),
    "Root-level results/tables/scdesign3_all_questions_status.tsv has been rewritten from the post-engine question gate table."
  ),
  key_files = list(
    engine_target_metrics = engine_metrics_path,
    target_gate_status_post_engine = paths$target_gate_status_post_engine_tsv,
    question_gate_status_post_engine = paths$question_gate_status_post_engine_tsv,
    communication_scdesign3_gate_post_engine = paths$communication_scdesign3_gate_post_engine_tsv,
    root_all_questions_status = paths$all_questions_status_tsv
  ),
  review_focus = c(
    "Approve `scdesign3_validated` only after reviewing 04e target metrics and figures.",
    "Rows where engine_applied is `no` still reflect 04d preflight or planned status.",
    "A WARN gate_status can mean either preflight warning or post-engine exploratory support; inspect gate_level to distinguish them.",
    "Downstream modules should consume the root all-question status table or this 04f post-engine table."
  ),
  extra_sections = list(
    "Applied Target Upgrades" = render_markdown_table_local(engine_applied_targets),
    "Target Gate Summary" = render_markdown_table_local(target_summary),
    "Question Gate Summary" = render_markdown_table_local(question_summary)
  )
)
write_markdown_local(report_lines, paths$report_md)

if (file.exists(manifest_path)) {
  unlink(manifest_path)
}
write_manifest_local(
  manifest_path = manifest_path,
  new_outputs = list(
    target_gate_status_post_engine_tsv = build_output_entry(paths$target_gate_status_post_engine_tsv, "tsv", module_name, "one row per scDesign3 target after engine finalization", base_dir = cfg$project_root, schema = infer_schema_from_df(target_gates_new)),
    question_gate_status_post_engine_tsv = build_output_entry(paths$question_gate_status_post_engine_tsv, "tsv", module_name, "one row per analysis question after engine finalization", base_dir = cfg$project_root, schema = infer_schema_from_df(question_gates_new)),
    communication_scdesign3_gate_post_engine_tsv = build_output_entry(paths$communication_scdesign3_gate_post_engine_tsv, "tsv", module_name, "communication pair scDesign3 gate projection after engine finalization", base_dir = cfg$project_root, schema = infer_schema_from_df(communication_gate_df)),
    scdesign3_all_questions_status_tsv = build_output_entry(paths$all_questions_status_tsv, "tsv", module_name, "root all-question scDesign3 status table updated by 04f", base_dir = cfg$project_root, schema = infer_schema_from_df(question_gates_new)),
    report = build_output_entry(paths$report_md, "md", module_name, "04f scDesign3 finalize report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    target_gate_pre_engine = target_gate_pre_path,
    question_gate_pre_engine = question_gate_pre_path,
    engine_target_metrics = engine_metrics_path,
    engine_status = engine_status_path,
    communication_pairs = communication_pairs_path
  ),
  version = cfg$module_version,
  depends_on = list(
    module_04d = cfg$module_04d_manifest_path,
    module_04e = file.path(cfg$manifest_dir, "04e_scdesign3_engine", "_manifest.json")
  )
)

message("04f scDesign3 finalize completed. engine-applied targets: ", nrow(engine_applied_targets), "; questions: ", nrow(question_gates_new))
