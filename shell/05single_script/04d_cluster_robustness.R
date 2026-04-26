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

empty_metrics_04d <- function() {
  data.frame(
    layer_id = character(0),
    n_cells = integer(0),
    n_clusters = integer(0),
    aRI_mean = numeric(0),
    aRI_sd = numeric(0),
    aMI_mean = numeric(0),
    aMI_sd = numeric(0),
    replicates_used = integer(0),
    runtime_sec = numeric(0),
    status = character(0),
    stringsAsFactors = FALSE
  )
}

manifest_annotated_keys_04d <- function(manifest) {
  output_names <- names(manifest$outputs %||% list())
  output_names[startsWith(output_names, "annotated_")]
}

layer_id_from_annotated_key_04d <- function(key) {
  sub("^annotated_", "", key)
}

manifest_04c <- read_manifest_local(cfg$module_04c_manifest_path)
manifest_04b <- read_manifest_local(cfg$module_04b_manifest_path)
annotated_keys <- manifest_annotated_keys_04d(manifest_04b)
layer_ids <- sort(vapply(annotated_keys, layer_id_from_annotated_key_04d, character(1)))

metrics_df <- empty_metrics_04d()
metrics_tsv <- file.path(cfg$subcluster_table_dir, "cluster_robustness_metrics.tsv")
write_tsv_local(metrics_df, metrics_tsv)

report_path <- file.path(cfg$subcluster_report_dir, "04d_cluster_robustness.md")
report_lines <- c(
  "# 04d Cluster Robustness (Placeholder)",
  "",
  "This is a C5 milestone placeholder. Real scDesign3 robustness computation is pending a subsequent milestone.",
  "",
  sprintf("- placeholder_status: `%s`", "pending_scdesign3_implementation"),
  sprintf("- metrics_tsv: `%s`", metrics_tsv),
  sprintf("- source_04c_manifest: `%s`", cfg$module_04c_manifest_path),
  sprintf("- source_04b_manifest: `%s`", cfg$module_04b_manifest_path),
  "",
  "## Expected Layers",
  render_markdown_table_local(data.frame(layer_id = layer_ids, stringsAsFactors = FALSE))
)
ensure_dir(dirname(report_path))
write_markdown_local(report_lines, report_path)

output_entries <- list(
  metrics_tsv = build_output_entry(
    metrics_tsv,
    "tsv",
    module_name,
    "zero-row placeholder robustness metrics; real scDesign3 metrics deferred",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(metrics_df)
  ),
  report = build_output_entry(
    report_path,
    "md",
    module_name,
    "04d cluster robustness placeholder report",
    base_dir = cfg$project_root
  )
)

if (file.exists(cfg$module_04d_manifest_path)) {
  unlink(cfg$module_04d_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_04d_manifest_path,
  new_outputs = output_entries,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_04c_manifest = cfg$module_04c_manifest_path,
    module_04b_manifest = cfg$module_04b_manifest_path,
    subcluster_summary_tsv = resolve_output_local(manifest_04c, "subcluster_summary_tsv")
  ),
  version = cfg$module_version,
  depends_on = list(
    module_04c = cfg$module_04c_manifest_path,
    module_04b = cfg$module_04b_manifest_path
  )
)

message("04d placeholder completed. expected layers: ", length(layer_ids))
