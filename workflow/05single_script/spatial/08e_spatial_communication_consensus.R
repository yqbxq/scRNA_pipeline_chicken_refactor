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
source(file.path(.script_dir, "helpers", "spatial_deconv_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_communication_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_svg_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_08e_spatial_communication_consensus"
prepare_dirs_spatial(cfg)

candidates <- st07_read_tsv(file.path(cfg$spatial_communication_table_dir, "spatial_comm_candidates.tsv"))
coloc <- st07_read_tsv(file.path(cfg$spatial_communication_table_dir, "lr_colocalization.tsv"))
neighbor <- st07_read_tsv(file.path(cfg$spatial_communication_table_dir, "communication_neighborhood_consistency.tsv"))
commot <- st07_read_tsv(file.path(cfg$spatial_commot_table_dir, "commot_spatial_summary.tsv"))
deconv_status <- st08_deconv_validation_status(cfg)
svg_lr <- st07_read_tsv(Sys.getenv("SVG_LR_SUPPORT_TSV", unset = file.path(cfg$spatial_table_dir, "09_svg", "09c_consensus", "svg_lr_support.tsv")))

merge_optional_08e <- function(left, right, keep_cols) {
  if (nrow(right) == 0 || !"comm_candidate_id" %in% colnames(right)) return(left)
  keep <- intersect(c("comm_candidate_id", keep_cols), colnames(right))
  merge(left, right[, keep, drop = FALSE], by = "comm_candidate_id", all.x = TRUE, sort = FALSE)
}

if (nrow(candidates) == 0) {
  consensus <- st08_empty_consensus()
} else {
  consensus <- candidates
  consensus <- merge_optional_08e(consensus, coloc, c("sender_spatial_abundance", "receiver_spatial_abundance", "sender_receiver_colocalization", "colocalization_support", "ligand_spot_expression_mean", "receptor_spot_expression_mean", "lr_expression_colocalization", "lr_expression_support"))
  consensus <- merge_optional_08e(consensus, neighbor, c("neighborhood_enrichment_score", "co_occurrence_score", "neighborhood_support"))
  if (nrow(commot) > 0) {
    commot$commot_key <- paste(commot$lr_axis_id, commot$section_id, commot$condition_value, sep = "|")
    consensus$commot_key <- paste(consensus$lr_axis_id, consensus$section_id, consensus$condition_value, sep = "|")
    for (col in c("status", "reason", "n_spot_pairs", "distance_threshold")) {
      if (!col %in% colnames(commot)) commot[[col]] <- if (col %in% c("n_spot_pairs", "distance_threshold")) NA_real_ else ""
    }
    commot_keep <- commot[!duplicated(commot$commot_key), c("commot_key", "signal_score", "spatial_support", "status", "reason", "n_spot_pairs", "distance_threshold"), drop = FALSE]
    names(commot_keep) <- c("commot_key", "commot_score", "commot_support", "commot_status", "commot_reason", "commot_n_spot_pairs", "commot_distance_threshold")
    consensus <- merge(consensus, commot_keep, by = "commot_key", all.x = TRUE, sort = FALSE)
    consensus$commot_key <- NULL
  }
  if (nrow(svg_lr) > 0) {
    for (col in c("comm_candidate_id", "lr_axis_id", "condition", "ligand_svg_tier", "receptor_svg_tier", "receiver_target_svg_support_n", "receiver_target_svg_support_genes", "svg_spatial_support_level", "svg_support_reason")) {
      if (!col %in% colnames(svg_lr)) svg_lr[[col]] <- ""
    }
    svg_keep <- svg_lr[!duplicated(svg_lr$comm_candidate_id), c("comm_candidate_id", "ligand_svg_tier", "receptor_svg_tier", "receiver_target_svg_support_n", "receiver_target_svg_support_genes", "svg_spatial_support_level", "svg_support_reason"), drop = FALSE]
    names(svg_keep) <- c("comm_candidate_id", "svg_ligand_support", "svg_receptor_support", "svg_receiver_target_support", "svg_receiver_target_support_genes", "svg_spatial_support_level", "svg_support_reason")
    consensus <- merge(consensus, svg_keep, by = "comm_candidate_id", all.x = TRUE, sort = FALSE)
  }
  numeric_cols <- c("sender_spatial_abundance", "receiver_spatial_abundance", "sender_receiver_colocalization", "ligand_spot_expression_mean", "receptor_spot_expression_mean", "lr_expression_colocalization", "neighborhood_enrichment_score", "co_occurrence_score", "commot_score", "commot_n_spot_pairs", "commot_distance_threshold")
  for (col in setdiff(colnames(st08_empty_consensus()), colnames(consensus))) consensus[[col]] <- if (col %in% numeric_cols) NA_real_ else ""
  if (!"commot_support" %in% colnames(consensus)) consensus$commot_support <- "no"
  consensus$svg_spatial_support_level <- ifelse(nzchar(consensus$svg_spatial_support_level), consensus$svg_spatial_support_level, "none")
  consensus$svg_ligand_support <- ifelse(nzchar(consensus$svg_ligand_support), consensus$svg_ligand_support, "unsupported")
  consensus$svg_receptor_support <- ifelse(nzchar(consensus$svg_receptor_support), consensus$svg_receptor_support, "unsupported")
  consensus$svg_receiver_target_support <- ifelse(nzchar(as.character(consensus$svg_receiver_target_support)), as.character(consensus$svg_receiver_target_support), "0")
  consensus$svg_support_reason <- ifelse(nzchar(consensus$svg_support_reason), consensus$svg_support_reason, "No SVG gene-level support available.")
  consensus$deconv_support_status <- deconv_status
  consensus$deconv_support_level <- st08_deconv_support_level(deconv_status)
  consensus$scrna_support_level <- st08_scrna_support_level(consensus)
  consensus$commot_score <- suppressWarnings(as.numeric(consensus$commot_score))
  scrna_primary <- consensus$scrna_support_level == "primary"
  scrna_supporting <- consensus$scrna_support_level %in% c("primary", "supporting")
  deconv_primary <- consensus$deconv_support_level %in% c("primary", "supporting")
  deconv_supporting <- consensus$deconv_support_level %in% c("primary", "supporting", "smoke_only")
  localized <- tolower(as.character(consensus$colocalization_support)) %in% c("yes", "true", "1", "ok")
  lr_expression <- tolower(as.character(consensus$lr_expression_support)) %in% c("yes", "true", "1", "ok")
  neighborhood <- tolower(as.character(consensus$neighborhood_support)) %in% c("yes", "true", "1", "ok")
  commot_support <- tolower(as.character(consensus$commot_support)) %in% c("yes", "true", "1", "ok")
  spatial_support_n <- as.integer(localized) + as.integer(neighborhood) + as.integer(commot_support | lr_expression)
  consensus$spatial_evidence_tier <- ifelse(
    scrna_primary & deconv_primary & localized & neighborhood & (commot_support | lr_expression),
    "spatial_primary",
    ifelse(scrna_supporting & deconv_supporting & spatial_support_n >= 2 & consensus$deconv_support_level != "smoke_only", "spatial_supporting",
      ifelse(consensus$cellchat_hypothesis == "yes" | consensus$scrna_support_level %in% c("hypothesis_only", "supporting") | consensus$deconv_support_level == "smoke_only", "spatial_hypothesis", "blocked")
    )
  )
  consensus$final_interpretation_level <- ifelse(consensus$spatial_evidence_tier == "spatial_primary", "core",
    ifelse(consensus$spatial_evidence_tier == "spatial_supporting", "exploratory",
      ifelse(consensus$spatial_evidence_tier == "spatial_hypothesis", "hypothesis_only", "blocked")
    )
  )
  consensus$reason <- ifelse(consensus$spatial_evidence_tier == "blocked", "Insufficient scRNA/deconv/spatial support for spatial communication interpretation.",
    ifelse(consensus$deconv_support_level == "smoke_only", "Deconvolution validation is smoke-only; communication remains hypothesis-level.", "")
  )
  strong_svg <- consensus$svg_spatial_support_level == "strong"
  consensus$reason[strong_svg & consensus$spatial_evidence_tier == "spatial_primary"] <- paste(trimws(consensus$reason[strong_svg & consensus$spatial_evidence_tier == "spatial_primary"]), "Core communication with SVG gene-level spatial support.")
  consensus$reason[strong_svg & consensus$spatial_evidence_tier == "spatial_hypothesis"] <- paste(trimws(consensus$reason[strong_svg & consensus$spatial_evidence_tier == "spatial_hypothesis"]), "Hypothesis with SVG-compatible spatial expression pattern.")
  consensus <- consensus[, colnames(st08_empty_consensus()), drop = FALSE]
}

consensus_tsv <- file.path(cfg$spatial_communication_table_dir, "spatial_communication_consensus.tsv")
tier_tsv <- file.path(cfg$spatial_communication_table_dir, "spatial_communication_evidence_tier.tsv")
spatial_write_tsv(consensus, consensus_tsv)
spatial_write_tsv(consensus[, c("comm_candidate_id", "lr_axis_id", "section_id", "scrna_evidence_tier", "spatial_evidence_tier", "final_interpretation_level", "reason"), drop = FALSE], tier_tsv)

st07_write_manifest_local(
  manifest_path = cfg$module_08e_spatial_communication_consensus_manifest_path,
  new_outputs = list(
    spatial_communication_consensus = build_output_entry(consensus_tsv, "tsv", module_name, "integrated scRNA/ST spatial communication evidence table", base_dir = cfg$project_root, schema = infer_schema_from_df(consensus)),
    spatial_communication_evidence_tier = build_output_entry(tier_tsv, "tsv", module_name, "final spatial communication evidence tier per candidate", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(spatial_comm_candidates = file.path(cfg$spatial_communication_table_dir, "spatial_comm_candidates.tsv"), svg_lr_support = file.path(cfg$spatial_table_dir, "09_svg", "09c_consensus", "svg_lr_support.tsv")),
  version = cfg$module_07_version,
  depends_on = list(spatial_08b_spatial_communication_eda = cfg$module_08b_spatial_communication_eda_manifest_path, spatial_08c_lr_colocalization = cfg$module_08c_lr_colocalization_manifest_path, spatial_08d_neighborhood_consistency = cfg$module_08d_neighborhood_consistency_manifest_path)
)

message(sprintf("08e spatial communication consensus completed. rows=%d", nrow(consensus)))
