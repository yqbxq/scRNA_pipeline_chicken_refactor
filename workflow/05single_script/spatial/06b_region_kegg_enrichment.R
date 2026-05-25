#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_enrichment_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_06b_region_kegg"

manifest_df <- run_spatial_region_enrichment(cfg, "kegg")
manifest_tsv <- file.path(cfg$spatial_region_kegg_table_dir, "region_kegg_manifest.tsv")
st06_write_tsv(manifest_df, manifest_tsv)

report_path <- file.path(cfg$spatial_enrichment_report_dir, "region_kegg_enrichment.md")
report_lines <- c(
  "# Spatial Region KEGG Enrichment",
  "",
  sprintf("- manifest: `%s`", relative_path_local(manifest_tsv, cfg$project_root)),
  sprintf("- top_n: `%s`", cfg$spatial_enrich_top_n),
  sprintf("- min_genes: `%s`", cfg$spatial_enrich_min_genes),
  sprintf("- qvalue: `%s`", cfg$spatial_enrich_qvalue),
  "",
  "## Status Summary",
  render_markdown_table_local(st06_status_summary(manifest_df)),
  "",
  "## Preview",
  render_markdown_table_local(utils::head(manifest_df, 50))
)
write_markdown_local(report_lines, report_path)

st06_write_manifest_local(
  manifest_path = cfg$module_06b_region_kegg_manifest_path,
  new_outputs = list(
    region_kegg_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "one row per ST05 region gene set and KEGG enrichment result", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    report = build_output_entry(report_path, "md", module_name, "spatial region KEGG enrichment report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    spatial_05_region_marker = cfg$module_05_region_marker_manifest_path,
    spatial_05a_region_pseudobulk = cfg$module_05a_region_pseudobulk_manifest_path,
    spatial_05_spatial_de = cfg$module_05_spatial_de_manifest_path,
    ortholog_manifest = cfg$ortholog_manifest_path
  ),
  version = cfg$module_06_version,
  depends_on = list(spatial_06a_region_go = cfg$module_06a_region_go_manifest_path)
)

message(sprintf("spatial region KEGG enrichment complete: %d manifest rows", nrow(manifest_df)))
