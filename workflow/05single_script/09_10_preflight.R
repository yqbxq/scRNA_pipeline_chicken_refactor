#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source(file.path(.script_dir, "helpers", "load_helpers_10.R"), encoding = "UTF-8")

load_required_packages(c("Seurat", "dplyr", "jsonlite"))

cfg <- get_single_script_config_10()
module_name <- "09_10_preflight"
prepare_dirs_10(cfg)

preflight_cols <- c(
  "check_id", "scope", "severity", "status", "expected", "observed",
  "details", "recommended_action"
)

preflight_row <- function(check_id, scope, severity, status, expected = "", observed = "", details = "", recommended_action = "") {
  data.frame(
    check_id = normalize_scalar_value(check_id),
    scope = normalize_scalar_value(scope),
    severity = normalize_scalar_value(severity, "info"),
    status = normalize_scalar_value(status, "pass"),
    expected = normalize_scalar_value(expected),
    observed = normalize_scalar_value(observed),
    details = normalize_scalar_value(details),
    recommended_action = normalize_scalar_value(recommended_action),
    stringsAsFactors = FALSE
  )
}

append_check <- function(rows, ...) {
  rows[[length(rows) + 1L]] <- preflight_row(...)
  rows
}

collapse_values <- function(x, limit = 20L) {
  x <- unique(trimws(as.character(x)))
  x <- x[nzchar(x) & !is.na(x)]
  if (length(x) == 0) {
    return("")
  }
  shown <- utils::head(x, limit)
  suffix <- if (length(x) > limit) sprintf(" ... +%s", length(x) - limit) else ""
  paste0(paste(shown, collapse = ","), suffix)
}

status_overall <- function(checks) {
  if (nrow(checks) == 0) {
    return("pass")
  }
  if (any(checks$status == "fail")) {
    return("fail")
  }
  if (any(checks$status == "warning")) {
    return("warning")
  }
  "pass"
}

sample_tokens_09_10 <- function(value) {
  value <- normalize_scalar_value(value)
  if (!nzchar(value)) {
    return(character(0))
  }
  parts <- trimws(unlist(strsplit(value, "[,;[:space:]]+", perl = TRUE), use.names = FALSE))
  parts[nzchar(parts)]
}

read_velocity_sample_ids_09_10 <- function(cfg) {
  from_env <- sample_tokens_09_10(Sys.getenv("VELOCITY_SAMPLE_IDS", ""))
  if (length(from_env) > 0) {
    return(unique(from_env))
  }
  sheet_path <- if (file.exists(cfg$canonical_sample_sheet)) cfg$canonical_sample_sheet else cfg$sample_sheet
  samples <- read_tsv_optional(sheet_path)
  if (nrow(samples) == 0 || !"sample_id" %in% colnames(samples)) {
    return(character(0))
  }
  if (!"run_velocity" %in% colnames(samples)) {
    samples$run_velocity <- "auto"
  }
  samples <- samples[tolower(normalize_flag(samples$run_velocity, "auto")) != "no", , drop = FALSE]
  unique(vapply(samples$sample_id, normalize_scalar_value, character(1)))
}

existing_file <- function(paths) {
  paths <- paths[nzchar(paths)]
  hit <- paths[file.exists(paths) & file.info(paths)$size > 0]
  if (length(hit) == 0) {
    return("")
  }
  hit[[1]]
}

velocity_loom_candidates_09_10 <- function(cfg, sample_id) {
  c(
    file.path(cfg$velocity_loom_dir, sprintf("%s.loom", sample_id)),
    file.path(cfg$dnbc4tools_out_dir, sample_id, "velocyto", sprintf("%s.loom", sample_id)),
    file.path(cfg$dnbc4tools_out_dir, sample_id, "outs", "velocyto", sprintf("%s.loom", sample_id)),
    file.path(cfg$cellranger_out_dir, sample_id, "velocyto", sprintf("%s.loom", sample_id)),
    file.path(cfg$cellranger_out_dir, sample_id, "outs", "velocyto", sprintf("%s.loom", sample_id))
  )
}

velocity_bam_candidates_09_10 <- function(cfg, sample_id) {
  c(
    file.path(cfg$dnbc4tools_out_dir, sample_id, cfg$velocity_bam_pattern),
    file.path(cfg$dnbc4tools_out_dir, sample_id, "outs", cfg$velocity_bam_pattern),
    file.path(cfg$data_dir, sample_id, cfg$velocity_bam_pattern),
    file.path(cfg$data_dir, sample_id, "outs", cfg$velocity_bam_pattern),
    file.path(cfg$cellranger_out_dir, sample_id, cfg$velocity_bam_pattern),
    file.path(cfg$cellranger_out_dir, sample_id, "outs", cfg$velocity_bam_pattern)
  )
}

velocity_h5_candidates_09_10 <- function(cfg, sample_id) {
  c(
    file.path(cfg$dnbc4tools_out_dir, sample_id, cfg$velocity_h5_pattern),
    file.path(cfg$dnbc4tools_out_dir, sample_id, "outs", cfg$velocity_h5_pattern),
    file.path(cfg$data_dir, sample_id, cfg$velocity_h5_pattern),
    file.path(cfg$data_dir, sample_id, "outs", cfg$velocity_h5_pattern),
    file.path(cfg$cellranger_out_dir, sample_id, cfg$velocity_h5_pattern),
    file.path(cfg$cellranger_out_dir, sample_id, "outs", cfg$velocity_h5_pattern)
  )
}

gate_status_09_10 <- function(cfg, gate_id) {
  gates <- read_tsv_optional(cfg$eda_gate_file)
  if (nrow(gates) == 0 || !"gate_id" %in% colnames(gates) || !"status" %in% colnames(gates)) {
    return("missing")
  }
  hit <- gates[gates$gate_id == gate_id, , drop = FALSE]
  if (nrow(hit) == 0) {
    return("missing")
  }
  tolower(normalize_scalar_value(hit$status[[1]], "pending"))
}

layer_hit_09_10 <- function(layers, layer_scope, cfg) {
  layer_scope <- normalize_scalar_value(layer_scope, cfg$panorama_layer_id)
  aliases <- unique(c(layer_scope, safe_id_09(layer_scope)))
  if (tolower(layer_scope) %in% c("panorama", tolower(cfg$panorama_layer_id))) {
    aliases <- unique(c(aliases, cfg$panorama_layer_id, "panorama"))
  }
  hit <- layers[layers$layer_id %in% aliases, , drop = FALSE]
  if (nrow(hit) == 0) {
    hit <- layers[tolower(layers$layer_id) %in% tolower(aliases), , drop = FALSE]
  }
  hit
}

object_cache <- new.env(parent = emptyenv())

load_preflight_object_09_10 <- function(layer_row) {
  path <- normalize_scalar_value(layer_row$object_rds[[1]])
  if (!nzchar(path) || !file.exists(path)) {
    stop(sprintf("object path is missing: %s", path), call. = FALSE)
  }
  key <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!exists(key, envir = object_cache, inherits = FALSE)) {
    assign(key, load_comm_layer_object_07(layer_row), envir = object_cache)
  }
  get(key, envir = object_cache, inherits = FALSE)
}

get_object_info_09_10 <- function(layers, layer_scope, cfg) {
  hit <- layer_hit_09_10(layers, layer_scope, cfg)
  if (nrow(hit) > 0) {
    obj <- load_preflight_object_09_10(hit[1, , drop = FALSE])
    return(list(ok = TRUE, layer_exists = TRUE, source = hit$source[[1]], layer_id = hit$layer_id[[1]], object = obj))
  }
  fallback <- layer_hit_09_10(layers, cfg$panorama_layer_id, cfg)
  if (nrow(fallback) > 0) {
    obj <- load_preflight_object_09_10(fallback[1, , drop = FALSE])
    return(list(ok = TRUE, layer_exists = FALSE, source = "panorama_fallback", layer_id = fallback$layer_id[[1]], object = obj))
  }
  list(ok = FALSE, layer_exists = FALSE, source = "", layer_id = "", object = NULL)
}

configured_layer_status_09_10 <- function(cfg, layer_scope) {
  cfg_df <- read_tsv_optional(cfg$object_layer_config_file)
  if (nrow(cfg_df) == 0 || !"layer_id" %in% colnames(cfg_df)) {
    return("missing_config")
  }
  if (!"enabled" %in% colnames(cfg_df)) {
    cfg_df$enabled <- ""
  }
  hit <- cfg_df[tolower(cfg_df$layer_id) == tolower(layer_scope), , drop = FALSE]
  if (nrow(hit) == 0) {
    return("missing_layer")
  }
  if (tolower(normalize_scalar_value(hit$enabled[[1]], "yes")) == "no") {
    return("disabled")
  }
  "configured"
}

check_label_value_09_10 <- function(checks, meta, pair_id, label_var, requested_label, role) {
  requested_label <- normalize_scalar_value(requested_label)
  if (!nzchar(requested_label) || requested_label %in% c("*", "auto")) {
    return(checks)
  }
  if (!label_var %in% colnames(meta)) {
    return(append_check(
      checks,
      sprintf("%s_label_column", role),
      pair_id,
      "error",
      "fail",
      label_var,
      collapse_values(colnames(meta)),
      sprintf("%s label column is absent", role),
      "补齐对应 subcluster 注释列，或把 trajectory_pairs.tsv 的 coarse_label_var 指到实际存在的细胞类型列。"
    ))
  }
  labels <- unique(trimws(as.character(meta[[label_var]])))
  if (!requested_label %in% labels) {
    return(append_check(
      checks,
      sprintf("%s_label_value", role),
      pair_id,
      "error",
      "fail",
      requested_label,
      collapse_values(labels),
      sprintf("%s label is absent from %s", role, label_var),
      "确认 root_group/terminal_group 与实际注释标签一致；必要时先完成 04b 子群注释。"
    ))
  }
  append_check(checks, sprintf("%s_label_value", role), pair_id, "info", "pass", requested_label, requested_label)
}

checks <- list()

pairs_raw <- read_tsv_optional(cfg$trajectory_pairs_sheet)
missing_cols <- setdiff(trajectory_pairs_cols_09, colnames(pairs_raw))
extra_cols <- setdiff(colnames(pairs_raw), trajectory_pairs_cols_09)
checks <- append_check(
  checks,
  "trajectory_pairs_schema",
  basename(cfg$trajectory_pairs_sheet),
  if (length(missing_cols) > 0) "error" else "info",
  if (length(missing_cols) > 0) "fail" else "pass",
  paste(trajectory_pairs_cols_09, collapse = ","),
  paste(colnames(pairs_raw), collapse = ","),
  if (length(extra_cols) > 0) sprintf("extra columns: %s", paste(extra_cols, collapse = ",")) else "",
  "保持 trajectory_pairs.tsv 与 helper 中的 17 列 schema 一致。"
)

pairs <- trajectory_read_pairs_09(cfg)
active_pairs <- pairs[pairs$enabled != "no" & pairs$method %in% c("trajectory", "velocity"), , drop = FALSE]
checks <- append_check(
  checks,
  "active_pair_count",
  "trajectory_pairs.tsv",
  if (nrow(active_pairs) == 0) "warning" else "info",
  if (nrow(active_pairs) == 0) "warning" else "pass",
  ">=1 active trajectory or velocity row",
  sprintf("trajectory=%s,velocity=%s", sum(active_pairs$method == "trajectory"), sum(active_pairs$method == "velocity")),
  "",
  "确认 analysis_questions.tsv/Tier 2 metadata 已生成需要的 trajectory/velocity 行。"
)

layers <- tryCatch(communication_layer_status_07(cfg), error = function(e) {
  empty_df_07(c("layer_id", "layer_role", "object_rds", "cluster_column", "source"))
})
checks <- append_check(
  checks,
  "available_layers",
  cfg$layer_status_file,
  if (nrow(layers) == 0) "error" else "info",
  if (nrow(layers) == 0) "fail" else "pass",
  "at least panorama or requested subcluster layer",
  collapse_values(layers$layer_id),
  "",
  "先完成 03d 注释；若 trajectory_pairs.tsv 使用子群 layer_scope，还要完成 04b。"
)

non_panorama <- active_pairs[
  nzchar(active_pairs$layer_scope) &
    !tolower(active_pairs$layer_scope) %in% c("panorama", tolower(cfg$panorama_layer_id)),
  ,
  drop = FALSE
]
if (nrow(non_panorama) > 0) {
  checks <- append_check(
    checks,
    "subcluster_manifest",
    cfg$module_04b_manifest_path,
    if (file.exists(cfg$module_04b_manifest_path)) "info" else "error",
    if (file.exists(cfg$module_04b_manifest_path)) "pass" else "fail",
    "04b manifest present when non-panorama layer_scope is requested",
    if (file.exists(cfg$module_04b_manifest_path)) "present" else "missing",
    sprintf("requested layer_scope: %s", collapse_values(non_panorama$layer_scope)),
    "先完成 04_subcluster，或把 trajectory_pairs.tsv 的 layer_scope 改成已经存在的 layer。"
  )
}

for (gate_id in c("annotation", "subcluster", "trajectory_inputs", "velocity_inputs")) {
  gate_status <- gate_status_09_10(cfg, gate_id)
  approved <- gate_status %in% c("approved", "yes", "true")
  checks <- append_check(
    checks,
    "eda_gate_status",
    gate_id,
    if (approved) "info" else "warning",
    if (approved) "pass" else "warning",
    "approved/yes/true before full stage execution",
    gate_status,
    "",
    "按 EDA gate 设计人工审核后再运行重型阶段；preflight 本身不强制通过 gate。"
  )
}

for (layer_scope in unique(non_panorama$layer_scope)) {
  layer_cfg_status <- configured_layer_status_09_10(cfg, layer_scope)
  checks <- append_check(
    checks,
    "object_layer_config",
    layer_scope,
    if (layer_cfg_status == "configured") "info" else "error",
    if (layer_cfg_status == "configured") "pass" else "fail",
    "configured and enabled in object_layers.tsv",
    layer_cfg_status,
    "",
    "在 config/object_layers.tsv 中配置并启用对应 layer_id，或修改 trajectory_pairs.tsv 的 layer_scope。"
  )
}

velocity_sample_ids <- read_velocity_sample_ids_09_10(cfg)
checks <- append_check(
  checks,
  "velocity_sample_ids",
  "10a",
  if (length(velocity_sample_ids) == 0) "warning" else "info",
  if (length(velocity_sample_ids) == 0) "warning" else "pass",
  "VELOCITY_SAMPLE_IDS or sample sheet sample_id",
  collapse_values(velocity_sample_ids),
  "",
  "如果 velocity 原始输出是 SRR/run 级别，请显式设置 VELOCITY_SAMPLE_IDS 为这些 run id。"
)

if (length(velocity_sample_ids) > 0) {
  for (sample_id in velocity_sample_ids) {
    loom_path <- existing_file(velocity_loom_candidates_09_10(cfg, sample_id))
    bam_path <- existing_file(velocity_bam_candidates_09_10(cfg, sample_id))
    h5_path <- existing_file(velocity_h5_candidates_09_10(cfg, sample_id))
    has_generation_inputs <- nzchar(bam_path) && nzchar(h5_path) && file.exists(cfg$velocity_gtf)
    status <- if (nzchar(loom_path) || has_generation_inputs) "pass" else "fail"
    checks <- append_check(
      checks,
      "velocity_sample_input",
      sample_id,
      if (identical(status, "pass")) "info" else "error",
      status,
      "existing loom, or BAM + filtered_feature_bc_matrix.h5 + velocity GTF",
      sprintf("loom=%s; bam=%s; h5=%s; gtf=%s", display_scalar_value(loom_path), display_scalar_value(bam_path), display_scalar_value(h5_path), display_scalar_value(cfg$velocity_gtf)),
      "",
      "补齐 VELOCITY_SAMPLE_IDS 对应的 Cell Ranger/velocyto 输出，或设置 CELLRANGER_OUT_DIR、VELOCITY_BAM_PATTERN、VELOCITY_GTF。"
    )
  }
}

if (!file.exists(cfg$velocity_gtf)) {
  checks <- append_check(
    checks,
    "velocity_gtf",
    "10a",
    "error",
    "fail",
    "existing VELOCITY_GTF",
    cfg$velocity_gtf,
    "",
    "设置 VELOCITY_GTF/CLEAN_GTF/REFERENCE_GTF 到实际存在的 Ensembl GTF。"
  )
}

for (idx in seq_len(nrow(active_pairs))) {
  pair_row <- active_pairs[idx, , drop = FALSE]
  pair_id <- normalize_scalar_value(pair_row$trajectory_id[[1]])
  layer_scope <- normalize_scalar_value(pair_row$layer_scope[[1]], cfg$panorama_layer_id)
  layer_hit <- layer_hit_09_10(layers, layer_scope, cfg)
  layer_ok <- nrow(layer_hit) > 0
  checks <- append_check(
    checks,
    "layer_scope_resolves",
    pair_id,
    if (layer_ok) "info" else "error",
    if (layer_ok) "pass" else "fail",
    layer_scope,
    collapse_values(layers$layer_id),
    "",
    "确保 trajectory_pairs.tsv 的 layer_scope 与 layer_status.tsv/04b manifest 中的 layer_id 一致。"
  )

  obj_info <- tryCatch(
    get_object_info_09_10(layers, layer_scope, cfg),
    error = function(e) list(ok = FALSE, layer_exists = layer_ok, source = "", layer_id = "", object = NULL, error = conditionMessage(e))
  )
  if (!isTRUE(obj_info$ok)) {
    checks <- append_check(
      checks,
      "layer_object_load",
      pair_id,
      "error",
      "fail",
      layer_scope,
      "",
      obj_info$error %||% "no fallback object available",
      "确认上游 RDS 仍存在且可被当前 R/Seurat 环境读取。"
    )
    next
  }

  meta <- obj_info$object@meta.data
  split_mode <- tolower(normalize_scalar_value(pair_row$split_mode[[1]], "auto"))
  split_var <- normalize_scalar_value(pair_row$condition_split_var[[1]])
  if (!split_mode %in% c("pooled", "force_pooled") && nzchar(split_var)) {
    if (!split_var %in% colnames(meta)) {
      checks <- append_check(
        checks,
        "split_var_exists",
        pair_id,
        "error",
        "fail",
        split_var,
        collapse_values(colnames(meta)),
        sprintf("object source: %s", obj_info$source),
        "把 condition_split_var 改成对象 metadata 中存在的列，或先补齐 metadata。"
      )
    } else {
      expected_values <- split_csv_local(pair_row$condition_split_values[[1]])
      observed_values <- unique(trimws(as.character(meta[[split_var]])))
      observed_values <- observed_values[nzchar(observed_values)]
      missing_values <- setdiff(expected_values, observed_values)
      checks <- append_check(
        checks,
        "split_values_present",
        pair_id,
        if (length(missing_values) == 0) "info" else "error",
        if (length(missing_values) == 0) "pass" else "fail",
        collapse_values(expected_values),
        collapse_values(observed_values),
        if (length(missing_values) > 0) sprintf("missing split values: %s; object source: %s", paste(missing_values, collapse = ","), obj_info$source) else sprintf("object source: %s", obj_info$source),
        "让 condition_split_values 与对象 metadata 的实际取值一致，或补齐 syf/f5 等项目分组列。"
      )
    }
  }

  coarse <- normalize_scalar_value(pair_row$coarse_label_var[[1]])
  fine <- normalize_scalar_value(pair_row$fine_label_var[[1]])
  for (label_var in unique(c(coarse, fine))) {
    if (!nzchar(label_var)) {
      next
    }
    checks <- append_check(
      checks,
      "label_column_exists",
      sprintf("%s:%s", pair_id, label_var),
      if (label_var %in% colnames(meta)) "info" else "error",
      if (label_var %in% colnames(meta)) "pass" else "fail",
      label_var,
      collapse_values(colnames(meta)),
      sprintf("object source: %s", obj_info$source),
      "确认 04b 子群注释已经写入 cell_subtype/cell_type 等列。"
    )
    if (label_var %in% colnames(meta)) {
      labels <- unique(trimws(as.character(meta[[label_var]])))
      labels <- labels[nzchar(labels)]
      unresolved <- length(labels) > 0 && all(grepl("^Uncertain", labels))
      if (unresolved) {
        checks <- append_check(
          checks,
          "label_column_unresolved",
          sprintf("%s:%s", pair_id, label_var),
          "warning",
          "warning",
          "biological labels, not only Uncertain-* cluster labels",
          collapse_values(labels),
          "",
          "补齐 marker panel/人工注释后重跑 annotation/subcluster，否则 root/terminal 生物标签无法解析。"
        )
      }
    }
  }

  if (pair_row$method[[1]] == "trajectory" && coarse %in% colnames(meta)) {
    checks <- check_label_value_09_10(checks, meta, pair_id, coarse, pair_row$root_group[[1]], "root")
    checks <- check_label_value_09_10(checks, meta, pair_id, coarse, pair_row$terminal_group[[1]], "terminal")
  }

  if (pair_row$method[[1]] == "velocity" && length(velocity_sample_ids) > 0) {
    velocity_ids <- velocity_cell_ids_from_meta_10(meta)
    metadata_samples <- unique(sub(":.*$", "", velocity_ids[grepl(":", velocity_ids)]))
    overlap <- intersect(metadata_samples, velocity_sample_ids)
    checks <- append_check(
      checks,
      "velocity_cell_namespace_overlap",
      pair_id,
      if (length(overlap) > 0) "info" else "error",
      if (length(overlap) > 0) "pass" else "fail",
      collapse_values(velocity_sample_ids),
      collapse_values(metadata_samples),
      sprintf("overlap: %s; object source: %s", collapse_values(overlap), obj_info$source),
      "让 Seurat cell_id/sample_id 与 loom sample id 使用同一命名空间；必要时保留 SRR/run 级 sample_id。"
    )
  }
}

checks_df <- if (length(checks) > 0) {
  dplyr::bind_rows(checks)
} else {
  trajectory_empty_df_09(preflight_cols)
}
for (col in preflight_cols) {
  if (!col %in% colnames(checks_df)) {
    checks_df[[col]] <- ""
  }
}
checks_df <- checks_df[, preflight_cols, drop = FALSE]
write_tsv_local(checks_df, cfg$trajectory_velocity_preflight_tsv)

overall <- status_overall(checks_df)
status_counts <- table(factor(checks_df$status, levels = c("pass", "warning", "fail")))
failed <- checks_df[checks_df$status == "fail", , drop = FALSE]
warnings <- checks_df[checks_df$status == "warning", , drop = FALSE]

report_lines <- c(
  "# 09/10 trajectory and velocity preflight",
  "",
  sprintf("- Overall status: `%s`", overall),
  sprintf("- Checks: pass=%s, warning=%s, fail=%s", status_counts[["pass"]], status_counts[["warning"]], status_counts[["fail"]]),
  sprintf("- Active rows: trajectory=%s, velocity=%s", sum(active_pairs$method == "trajectory"), sum(active_pairs$method == "velocity")),
  sprintf("- Available layers: %s", display_scalar_value(collapse_values(layers$layer_id), "none")),
  sprintf("- Velocity samples: %s", display_scalar_value(collapse_values(velocity_sample_ids), "none")),
  "",
  "## Blocking failures",
  ""
)

if (nrow(failed) == 0) {
  report_lines <- c(report_lines, "None.")
} else {
  for (i in seq_len(nrow(failed))) {
    row <- failed[i, , drop = FALSE]
    report_lines <- c(
      report_lines,
      sprintf("- `%s` (%s): %s", row$check_id[[1]], row$scope[[1]], display_scalar_value(row$details[[1]], row$observed[[1]])),
      sprintf("  - Recommended action: %s", display_scalar_value(row$recommended_action[[1]], "review configuration"))
    )
  }
}

report_lines <- c(report_lines, "", "## Warnings", "")
if (nrow(warnings) == 0) {
  report_lines <- c(report_lines, "None.")
} else {
  for (i in seq_len(nrow(warnings))) {
    row <- warnings[i, , drop = FALSE]
    report_lines <- c(
      report_lines,
      sprintf("- `%s` (%s): %s", row$check_id[[1]], row$scope[[1]], display_scalar_value(row$details[[1]], row$observed[[1]])),
      sprintf("  - Recommended action: %s", display_scalar_value(row$recommended_action[[1]], "review configuration"))
    )
  }
}

ensure_dir(dirname(cfg$trajectory_velocity_preflight_report_md))
writeLines(report_lines, cfg$trajectory_velocity_preflight_report_md, useBytes = TRUE)

outputs <- list(
  preflight_checks_tsv = build_output_entry(
    cfg$trajectory_velocity_preflight_tsv,
    "tsv",
    module_name,
    "one row per preflight check",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(checks_df)
  ),
  report_md = build_output_entry(
    cfg$trajectory_velocity_preflight_report_md,
    "md",
    module_name,
    "human-readable 09/10 preflight report",
    base_dir = cfg$project_root
  )
)

velocity_manifest_10(
  cfg,
  cfg$module_09_10_preflight_manifest_path,
  module_name,
  outputs,
  inputs = list(
    trajectory_pairs_sheet = cfg$trajectory_pairs_sheet,
    object_layer_config_file = cfg$object_layer_config_file,
    layer_status_file = cfg$layer_status_file,
    eda_gate_file = cfg$eda_gate_file,
    velocity_gtf = cfg$velocity_gtf,
    velocity_sample_ids = paste(velocity_sample_ids, collapse = ",")
  ),
  depends_on = list(
    `03d_annotate` = cfg$module_03d_manifest_path,
    `04b_subcluster_annotate` = cfg$module_04b_manifest_path
  )
)

message(sprintf("09/10 preflight completed: %s; report: %s", overall, cfg$trajectory_velocity_preflight_report_md))

if (identical(tolower(Sys.getenv("PREFLIGHT_STRICT", "no")), "yes") && identical(overall, "fail")) {
  quit(status = 1L)
}
