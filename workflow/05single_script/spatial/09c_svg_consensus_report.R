#!/usr/bin/env Rscript
.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})
PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/manifest_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_enrichment_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_communication_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_svg_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_09c_svg_consensus"
prepare_dirs_spatial(cfg)
svg_root <- file.path(cfg$spatial_table_dir, "09_svg")
out_dir <- file.path(svg_root, "09c_consensus")
ensure_dir(out_dir)
manifest_path <- file.path(cfg$manifest_dir, module_name, "_manifest.json")

method_df <- st09_collect_method_outputs(svg_root)
consensus <- st09_build_consensus(method_df)

method_consistency <- if (nrow(method_df) == 0) {
  data.frame(method = character(), status = character(), status_class = character(), n_genes = integer(), n_sections = integer(), n_conditions = integer(), stringsAsFactors = FALSE)
} else {
  method_df$status_class <- st09_method_status_class(method_df$status)
  keys <- unique(method_df[, c("method", "status", "status_class"), drop = FALSE])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    hit <- method_df[method_df$method == keys$method[[i]] & method_df$status == keys$status[[i]], , drop = FALSE]
    data.frame(
      method = keys$method[[i]],
      status = keys$status[[i]],
      status_class = keys$status_class[[i]],
      n_genes = sum(nzchar(hit$gene_symbol) | nzchar(hit$gene_id)),
      n_sections = length(unique(hit$section_id)),
      n_conditions = length(unique(hit$condition)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

condition_consistency <- if (nrow(consensus) == 0) {
  data.frame(condition = character(), n_high_confidence = integer(), n_stable = integer(), n_single_method = integer(), n_proxy_only = integer(), stringsAsFactors = FALSE)
} else {
  conditions <- unique(consensus$condition)
  do.call(rbind, lapply(conditions, function(cond) {
    hit <- consensus[consensus$condition == cond, , drop = FALSE]
    data.frame(
      condition = cond,
      n_high_confidence = sum(hit$svg_tier == "high_confidence_svg"),
      n_stable = sum(hit$svg_tier == "condition_stable_svg"),
      n_single_method = sum(hit$svg_tier == "single_method_svg"),
      n_proxy_only = sum(hit$svg_tier == "proxy_only_svg"),
      stringsAsFactors = FALSE
    )
  }))
}

evidence_tier <- consensus[, intersect(c("section_id", "condition", "gene_id", "gene_symbol", "svg_tier", "evidence_tier", "support_methods", "status", "reason"), colnames(consensus)), drop = FALSE]
gene_sets <- st09_svg_gene_sets(consensus)
enrichment_handoff <- st09_enrichment_handoff(gene_sets)

candidates <- st09_read_tsv(file.path(cfg$spatial_communication_table_dir, "spatial_communication_consensus.tsv"))
if (nrow(candidates) == 0) candidates <- st09_read_tsv(file.path(cfg$spatial_communication_table_dir, "spatial_comm_candidates.tsv"))
lr_support <- st09_lr_support(consensus, candidates)

scenic_targets <- st09_read_tsv(file.path(cfg$metadata_dir, "scenic_targets.tsv"))
regulation_support <- st09_regulation_support(consensus, scenic_targets)

trajectory_pairs <- st09_read_tsv(file.path(cfg$metadata_dir, "trajectory_pairs.tsv"))
gene_program_targets <- st09_read_tsv(file.path(cfg$metadata_dir, "gene_program_targets.tsv"))
trajectory_support <- st09_trajectory_support(consensus, trajectory_pairs, gene_program_targets)
question_gate_status <- st09_question_gate_status(consensus)

project_manifest <- data.frame(
  status = ifelse(nrow(consensus) > 0 && any(consensus$svg_tier != "unsupported"), "ok", "skipped_no_usable_svg"),
  reason = ifelse(nrow(consensus) > 0 && any(consensus$svg_tier != "unsupported"), "", "No formal, single-method, or proxy SVG evidence was available."),
  method_gene_n = nrow(method_df),
  consensus_gene_n = nrow(consensus),
  high_confidence_gene_n = sum(consensus$svg_tier == "high_confidence_svg"),
  proxy_only_gene_n = sum(consensus$svg_tier == "proxy_only_svg"),
  svg_lr_support_tsv = file.path(out_dir, "svg_lr_support.tsv"),
  svg_regulation_support_tsv = file.path(out_dir, "svg_regulation_support.tsv"),
  svg_enrichment_handoff_tsv = file.path(out_dir, "svg_enrichment_handoff.tsv"),
  svg_trajectory_support_tsv = file.path(out_dir, "svg_trajectory_support.tsv"),
  stringsAsFactors = FALSE
)

paths <- list(
  svg_method_consistency = file.path(out_dir, "svg_method_consistency.tsv"),
  svg_condition_consistency = file.path(out_dir, "svg_condition_consistency.tsv"),
  svg_consensus = file.path(out_dir, "svg_consensus.tsv"),
  svg_evidence_tier = file.path(out_dir, "svg_evidence_tier.tsv"),
  svg_gene_sets = file.path(out_dir, "svg_gene_sets.tsv"),
  svg_lr_support = file.path(out_dir, "svg_lr_support.tsv"),
  svg_regulation_support = file.path(out_dir, "svg_regulation_support.tsv"),
  svg_enrichment_handoff = file.path(out_dir, "svg_enrichment_handoff.tsv"),
  svg_trajectory_support = file.path(out_dir, "svg_trajectory_support.tsv"),
  svg_question_gate_status = file.path(out_dir, "svg_question_gate_status.tsv"),
  svg_project_handoff_manifest = file.path(out_dir, "svg_project_handoff_manifest.tsv"),
  svg_report = file.path(out_dir, "svg_report.md")
)

st09_write_tsv(method_consistency, paths$svg_method_consistency)
st09_write_tsv(condition_consistency, paths$svg_condition_consistency)
st09_write_tsv(consensus, paths$svg_consensus)
st09_write_tsv(evidence_tier, paths$svg_evidence_tier)
st09_write_tsv(gene_sets, paths$svg_gene_sets)
st09_write_tsv(lr_support, paths$svg_lr_support)
st09_write_tsv(regulation_support, paths$svg_regulation_support)
st09_write_tsv(enrichment_handoff, paths$svg_enrichment_handoff)
st09_write_tsv(trajectory_support, paths$svg_trajectory_support)
st09_write_tsv(question_gate_status, paths$svg_question_gate_status)
st09_write_tsv(project_manifest, paths$svg_project_handoff_manifest)

report_lines <- c(
  "# Spatial SVG Consensus",
  "",
  sprintf("- status: `%s`", project_manifest$status[[1]]),
  sprintf("- consensus_gene_n: `%s`", project_manifest$consensus_gene_n[[1]]),
  sprintf("- high_confidence_gene_n: `%s`", project_manifest$high_confidence_gene_n[[1]]),
  sprintf("- proxy_only_gene_n: `%s`", project_manifest$proxy_only_gene_n[[1]]),
  "",
  "## Method Consistency",
  render_markdown_table_local(utils::head(method_consistency, 50)),
  "",
  "## Question Gates",
  render_markdown_table_local(question_gate_status),
  "",
  "SVG evidence is gene-level spatial support. It can annotate communication, regulation, enrichment, and trajectory-in-space interpretation, but it does not by itself upgrade communication to a primary tier or prove pseudotime direction."
)
write_markdown_local(report_lines, paths$svg_report)

outputs <- list(
  svg_method_consistency = build_output_entry(paths$svg_method_consistency, "tsv", module_name, "SVG method status and consistency", base_dir = cfg$project_root, schema = infer_schema_from_df(method_consistency)),
  svg_condition_consistency = build_output_entry(paths$svg_condition_consistency, "tsv", module_name, "SVG condition consistency summary", base_dir = cfg$project_root, schema = infer_schema_from_df(condition_consistency)),
  svg_consensus = build_output_entry(paths$svg_consensus, "tsv", module_name, "per-gene SVG consensus tier", base_dir = cfg$project_root, schema = infer_schema_from_df(consensus)),
  svg_evidence_tier = build_output_entry(paths$svg_evidence_tier, "tsv", module_name, "per-gene SVG evidence tier handoff", base_dir = cfg$project_root, schema = infer_schema_from_df(evidence_tier)),
  svg_gene_sets = build_output_entry(paths$svg_gene_sets, "tsv", module_name, "SVG-derived gene sets for enrichment handoff", base_dir = cfg$project_root, schema = infer_schema_from_df(gene_sets)),
  svg_lr_support = build_output_entry(paths$svg_lr_support, "tsv", module_name, "SVG support annotations for communication candidates", base_dir = cfg$project_root, schema = infer_schema_from_df(lr_support)),
  svg_regulation_support = build_output_entry(paths$svg_regulation_support, "tsv", module_name, "SVG support annotations for regulator targets", base_dir = cfg$project_root, schema = infer_schema_from_df(regulation_support)),
  svg_enrichment_handoff = build_output_entry(paths$svg_enrichment_handoff, "tsv", module_name, "SVG gene set handoff for enrichment", base_dir = cfg$project_root, schema = infer_schema_from_df(enrichment_handoff)),
  svg_trajectory_support = build_output_entry(paths$svg_trajectory_support, "tsv", module_name, "SVG support for trajectory-in-space genes", base_dir = cfg$project_root, schema = infer_schema_from_df(trajectory_support)),
  svg_question_gate_status = build_output_entry(paths$svg_question_gate_status, "tsv", module_name, "I06/I07 exploratory SVG gate status", base_dir = cfg$project_root, schema = infer_schema_from_df(question_gate_status)),
  svg_project_handoff_manifest = build_output_entry(paths$svg_project_handoff_manifest, "tsv", module_name, "single-row ST09 global handoff manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(project_manifest)),
  svg_report = build_output_entry(paths$svg_report, "md", module_name, "ST09 SVG consensus report", base_dir = cfg$project_root)
)

st06_write_manifest_local(
  manifest_path,
  outputs,
  module_name,
  cfg$project_root,
  inputs = list(spatialde2_dir = file.path(svg_root, "spatialde2"), sparkx_dir = file.path(svg_root, "sparkx"), communication_candidates = file.path(cfg$spatial_communication_table_dir, "spatial_comm_candidates.tsv")),
  version = cfg$module_07_version,
  depends_on = list(spatial_09a_spatialde2_svg = file.path(cfg$manifest_dir, "spatial_09a_spatialde2_svg", "_manifest.json"), spatial_09b_sparkx_svg = file.path(cfg$manifest_dir, "spatial_09b_sparkx_svg", "_manifest.json"))
)

message(sprintf("%s complete: %s genes=%d", module_name, project_manifest$status[[1]], nrow(consensus)))
