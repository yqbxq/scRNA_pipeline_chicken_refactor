#!/usr/bin/env Rscript
.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})
PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_enrichment_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_08a_decoupler"
prepare_dirs_spatial(cfg)
out_dir <- file.path(cfg$spatial_table_dir, "spatial_regulation", "decoupler")
manifest_path <- file.path(cfg$manifest_dir, module_name, "_manifest.json")
ensure_dir(out_dir)
cmd <- c(
  file.path(cfg$pipeline_root, "workflow", "04python", "spatial_decoupler.py"),
  "run",
  "--st-module", Sys.getenv("SPATIAL_REGULATION_H5AD_MODULE", unset = "spatial_03_region"),
  "--out-dir", out_dir,
  "--results-dir", cfg$results_dir
)
network <- Sys.getenv("SPATIAL_REGULATION_NETWORK_TSV", unset = "")
if (nzchar(network)) cmd <- c(cmd, "--network-tsv", network)
svg_gene_sets <- Sys.getenv("SPATIAL_REGULATION_SVG_GENE_SETS", unset = file.path(cfg$spatial_table_dir, "09_svg", "09c_consensus", "svg_gene_sets.tsv"))
if (file.exists(svg_gene_sets)) cmd <- c(cmd, "--svg-gene-sets", svg_gene_sets)
run <- tryCatch(system2(cfg$py_spatial_bin, args = cmd, stdout = TRUE, stderr = TRUE), error = function(e) structure(conditionMessage(e), status = 127))
manifest_tsv <- file.path(out_dir, "spatial_regulation_h5ad_manifest.tsv")
if (!file.exists(manifest_tsv)) {
  st06_write_tsv(data.frame(status = "failed_sidecar", reason = paste(as.character(run), collapse = " "), input_h5ad_path = "", stringsAsFactors = FALSE), manifest_tsv)
}
df <- st06_read_tsv(manifest_tsv)
st06_write_manifest_local(
  manifest_path,
  list(spatial_regulation_h5ad_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "H5AD-first spatial decoupleR manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(df))),
  module_name,
  cfg$project_root,
  inputs = list(spatial_h5ad_module = Sys.getenv("SPATIAL_REGULATION_H5AD_MODULE", unset = "spatial_03_region"), svg_gene_sets = svg_gene_sets),
  version = cfg$module_07_version
)
message(sprintf("%s complete: %s", module_name, st06_scalar(df$status, "unknown")))
