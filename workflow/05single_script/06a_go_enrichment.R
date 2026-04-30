#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source_utf8 <- function(path) source(path, encoding = "UTF-8")

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_05.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_06.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "comparison_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "deg_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "enrichment_utils.R"))

load_required_packages(c("dplyr", "tibble", "jsonlite", "ggplot2", "AnnotationDbi", "org.Gg.eg.db"))

cfg <- get_single_script_config_06()
module_name <- "06a_go_enrichment"
prepare_dirs_06(cfg)
set.seed(cfg$random_seed)

chicken_org_db <- get_org_db_06("org.Gg.eg.db")
human_enabled <- human_strategy_enabled_06(cfg)
human_org_db <- if (human_enabled) optional_org_db_06("org.Hs.eg.db") else NULL
ortholog_map <- if (human_enabled) load_ortholog_map(cfg) else empty_df_06(c("chicken_symbol", "human_symbol"))
input_grid <- read_enrichment_input_grid_06(cfg)
if (nrow(input_grid) > 0 && "database" %in% colnames(input_grid)) {
  input_grid <- input_grid[toupper(input_grid$database) %in% c("GO", "BOTH", "ALL"), , drop = FALSE]
}

write_go_result <- function(
    layer_id,
    comparison_id,
    cluster_id,
    direction,
    ont,
    source_species,
    deg_source,
    deg_status,
    genes,
    org_db,
    enrichment_source,
    target_id = "",
    formal_status = "",
    result_level = "",
    biological_replicates = "",
    background_genes = character(0),
    background_tsv = "") {
  paths <- enrichment_paths_06(cfg, layer_id, comparison_id, cluster_id, paste0("go_", ont), source_species, direction)
  extra <- list(
    target_id = target_id,
    layer_id = layer_id,
    comparison_id = comparison_id,
    cluster_id = cluster_id,
    gene_direction = direction,
    analysis_type = "go",
    ontology = ont,
    source_species = source_species,
    enrichment_source = enrichment_source,
    deg_source = deg_source,
    deg_inference_status = deg_status,
    formal_status = formal_status,
    result_level = result_level,
    biological_replicates = biological_replicates,
    background_gene_n = length(background_genes)
  )
  background_conversion <- if (!is.null(org_db) && length(background_genes) > 0) {
    convert_symbols_to_entrezid(background_genes, org_db, keytype = "SYMBOL")
  } else {
    list(ids = character(0), mapped_n = 0L)
  }

  if (length(genes) < cfg$enrichment_min_input_genes) {
    result_df <- empty_enrichment_result_06(names(extra))
    counts <- write_enrichment_result_06(result_df, paths, cfg)
    return(build_enrichment_manifest_row_06(
      target_id = target_id,
      layer_id = layer_id,
      comparison_id = comparison_id,
      cluster_id = cluster_id,
      gene_direction = direction,
      analysis_type = "go",
      ontology = ont,
      source_species = source_species,
      deg_source = deg_source,
      deg_inference_status = deg_status,
      formal_status = formal_status,
      result_level = result_level,
      biological_replicates = biological_replicates,
      background_tsv = background_tsv,
      background_gene_n = length(background_genes),
      background_mapped_n = background_conversion$mapped_n %||% 0L,
      status = "too_few_genes",
      reason = sprintf("%s genes before ID conversion", length(genes)),
      input_gene_n = length(genes),
      mapped_gene_n = 0L,
      mapping_rate = NA_real_,
      significant_term_n = counts$significant_n,
      paths = paths
    ))
  }

  if (!clusterprofiler_available_06()) {
    gp_organism <- if (identical(source_species, "human")) "hsapiens" else "ggallus"
    run_res <- run_gprofiler_enrichment_single(
      gene_symbols = genes,
      organism = gp_organism,
      sources = paste0("GO:", ont),
      pvalue_cutoff = cfg$enrichment_pvalue_cutoff,
      background_symbols = background_genes
    )
    result_df <- format_enrichment_result_tsv(run_res$result, extra)
    counts <- write_enrichment_result_06(result_df, paths, cfg)
    dotplot <- ""
    barplot <- ""
    if (nrow(result_df) > 0) {
      plot_title <- sprintf("%s %s %s %s %s %s %s reps=%s", layer_id, comparison_id, cluster_id, direction, ont, source_species, formal_status, biological_replicates)
      dotplot <- plot_enrichment_dotplot(result_df, plot_title, top_n = 20L, output_png = paths$dotplot_png)
      barplot <- plot_enrichment_barplot(result_df, plot_title, top_n = 15L, output_png = paths$barplot_png)
    }
    return(build_enrichment_manifest_row_06(
      target_id = target_id,
      layer_id = layer_id,
      comparison_id = comparison_id,
      cluster_id = cluster_id,
      gene_direction = direction,
      analysis_type = "go",
      ontology = ont,
      source_species = source_species,
      deg_source = deg_source,
      deg_inference_status = deg_status,
      formal_status = formal_status,
      result_level = result_level,
      biological_replicates = biological_replicates,
      background_tsv = background_tsv,
      background_gene_n = length(background_genes),
      background_mapped_n = background_conversion$mapped_n %||% 0L,
      status = run_res$status,
      reason = run_res$reason,
      input_gene_n = length(genes),
      mapped_gene_n = NA_integer_,
      mapping_rate = NA_real_,
      significant_term_n = counts$significant_n,
      paths = paths,
      dotplot_png = normalize_scalar_value(dotplot),
      barplot_png = normalize_scalar_value(barplot)
    ))
  }

  if (is.null(org_db)) {
    result_df <- empty_enrichment_result_06(names(extra))
    counts <- write_enrichment_result_06(result_df, paths, cfg)
    return(build_enrichment_manifest_row_06(
      target_id = target_id,
      layer_id = layer_id,
      comparison_id = comparison_id,
      cluster_id = cluster_id,
      gene_direction = direction,
      analysis_type = "go",
      ontology = ont,
      source_species = source_species,
      deg_source = deg_source,
      deg_inference_status = deg_status,
      formal_status = formal_status,
      result_level = result_level,
      biological_replicates = biological_replicates,
      background_tsv = background_tsv,
      background_gene_n = length(background_genes),
      background_mapped_n = background_conversion$mapped_n %||% 0L,
      status = if (identical(source_species, "human")) "human_orgdb_missing" else "orgdb_missing",
      reason = sprintf("%s OrgDb is not available for clusterProfiler GO backend", source_species),
      input_gene_n = length(genes),
      mapped_gene_n = 0L,
      mapping_rate = NA_real_,
      significant_term_n = counts$significant_n,
      paths = paths
    ))
  }

  conversion <- convert_symbols_to_entrezid(genes, org_db, keytype = "SYMBOL")
  if (length(conversion$ids) < cfg$enrichment_min_input_genes) {
    result_df <- empty_enrichment_result_06(names(extra))
    counts <- write_enrichment_result_06(result_df, paths, cfg)
    return(build_enrichment_manifest_row_06(
      target_id = target_id,
      layer_id = layer_id,
      comparison_id = comparison_id,
      cluster_id = cluster_id,
      gene_direction = direction,
      analysis_type = "go",
      ontology = ont,
      source_species = source_species,
      deg_source = deg_source,
      deg_inference_status = deg_status,
      formal_status = formal_status,
      result_level = result_level,
      biological_replicates = biological_replicates,
      background_tsv = background_tsv,
      background_gene_n = length(background_genes),
      background_mapped_n = background_conversion$mapped_n %||% 0L,
      status = "too_few_mapped_genes",
      reason = sprintf("%s/%s genes mapped to ENTREZID via %s", conversion$mapped_n, conversion$input_n, conversion$keytype),
      input_gene_n = conversion$input_n,
      mapped_gene_n = conversion$mapped_n,
      mapping_rate = conversion$mapping_rate,
      significant_term_n = counts$significant_n,
      paths = paths
    ))
  }

  run_res <- run_go_enrichment_single(
    gene_list = conversion$ids,
    org_db = org_db,
    ont = ont,
    pvalue_cutoff = cfg$enrichment_pvalue_cutoff,
    qvalue_cutoff = cfg$enrichment_qvalue_cutoff,
    min_gs_size = cfg$enrichment_min_gs_size,
    max_gs_size = cfg$enrichment_max_gs_size,
    universe_ids = background_conversion$ids
  )
  result_df <- format_enrichment_result_tsv(run_res$result, extra)
  counts <- write_enrichment_result_06(result_df, paths, cfg)
  dotplot <- ""
  barplot <- ""
  if (nrow(result_df) > 0) {
    plot_title <- sprintf("%s %s %s %s %s %s reps=%s", layer_id, comparison_id, cluster_id, direction, ont, formal_status, biological_replicates)
    dotplot <- plot_enrichment_dotplot(run_res$result, plot_title, top_n = 20L, output_png = paths$dotplot_png)
    barplot <- plot_enrichment_barplot(run_res$result, plot_title, top_n = 15L, output_png = paths$barplot_png)
  }

  build_enrichment_manifest_row_06(
    target_id = target_id,
    layer_id = layer_id,
    comparison_id = comparison_id,
    cluster_id = cluster_id,
    gene_direction = direction,
    analysis_type = "go",
    ontology = ont,
    source_species = source_species,
    deg_source = deg_source,
    deg_inference_status = deg_status,
    formal_status = formal_status,
    result_level = result_level,
    biological_replicates = biological_replicates,
    background_tsv = background_tsv,
    background_gene_n = length(background_genes),
    background_mapped_n = background_conversion$mapped_n %||% 0L,
    status = run_res$status,
    reason = run_res$reason,
    input_gene_n = conversion$input_n,
    mapped_gene_n = conversion$mapped_n,
    mapping_rate = conversion$mapping_rate,
    significant_term_n = counts$significant_n,
    paths = paths,
    dotplot_png = normalize_scalar_value(dotplot),
    barplot_png = normalize_scalar_value(barplot)
  )
}

manifest_rows <- list()

for (idx in seq_len(nrow(input_grid))) {
  item <- input_grid[idx, , drop = FALSE]
  layer_id <- item$layer_id[[1]]
  comparison_id <- item$comparison_id[[1]]
  target_id <- normalize_scalar_value(item$target_id[[1]])
  formal_status <- normalize_scalar_value(item$formal_status[[1]], item$deg_inference_status[[1]])
  result_level <- normalize_scalar_value(item$result_level[[1]])
  biological_replicates <- normalize_scalar_value(item$biological_replicates[[1]], "unknown")
  background_tsv <- normalize_scalar_value(item$background_tsv[[1]])
  background_genes <- read_background_genes_06(background_tsv)
  message("06a GO enrichment: ", layer_id, " / ", comparison_id)

  if (!nzchar(item$deg_tsv[[1]]) || !file.exists(item$deg_tsv[[1]])) {
    manifest_rows[[length(manifest_rows) + 1]] <- build_enrichment_manifest_row_06(
      target_id = target_id,
      layer_id = layer_id,
      comparison_id = comparison_id,
      cluster_id = "",
      gene_direction = "all",
      analysis_type = "go",
      ontology = "",
      source_species = "",
      deg_source = item$deg_source[[1]],
      deg_inference_status = item$deg_inference_status[[1]],
      formal_status = formal_status,
      result_level = result_level,
      biological_replicates = biological_replicates,
      background_tsv = background_tsv,
      background_gene_n = length(background_genes),
      background_mapped_n = 0L,
      status = "missing_deg_source",
      reason = "No targeted gene program table was available from gene_program_registry.tsv",
      input_gene_n = 0L,
      mapped_gene_n = 0L,
      mapping_rate = NA_real_,
      significant_term_n = 0L
    )
    next
  }

  deg_df <- read_tsv_optional(item$deg_tsv[[1]])
  cluster_ids <- cluster_ids_from_deg_table_06(deg_df)
  if (length(cluster_ids) == 0) {
    cluster_ids <- "all"
  }

  for (cluster_id in cluster_ids) {
    gene_sets <- prepare_deg_gene_list_from_df_06(deg_df, cfg, cluster_id)
    for (direction in names(gene_sets)) {
      genes <- gene_sets[[direction]]$genes
      for (ont in c("BP", "CC", "MF")) {
        manifest_rows[[length(manifest_rows) + 1]] <- write_go_result(
          layer_id = layer_id,
          comparison_id = comparison_id,
          cluster_id = cluster_id,
          direction = direction,
          ont = ont,
          source_species = "chicken",
          deg_source = item$deg_source[[1]],
          deg_status = item$deg_inference_status[[1]],
          genes = genes,
          org_db = chicken_org_db,
          enrichment_source = "org.Gg.eg.db",
          target_id = target_id,
          formal_status = formal_status,
          result_level = result_level,
          biological_replicates = biological_replicates,
          background_genes = background_genes,
          background_tsv = background_tsv
        )

        if (human_enabled) {
          human_background <- map_genes_to_human(background_genes, ortholog_map)
          if (is.null(human_org_db) && clusterprofiler_available_06()) {
            paths <- enrichment_paths_06(cfg, layer_id, comparison_id, cluster_id, paste0("go_", ont), "human", direction)
            write_enrichment_result_06(empty_enrichment_result_06(), paths, cfg)
            manifest_rows[[length(manifest_rows) + 1]] <- build_enrichment_manifest_row_06(
              target_id = target_id,
              layer_id = layer_id,
              comparison_id = comparison_id,
              cluster_id = cluster_id,
              gene_direction = direction,
              analysis_type = "go",
              ontology = ont,
              source_species = "human",
              deg_source = item$deg_source[[1]],
              deg_inference_status = item$deg_inference_status[[1]],
              formal_status = formal_status,
              result_level = result_level,
              biological_replicates = biological_replicates,
              background_tsv = background_tsv,
              background_gene_n = length(background_genes),
              background_mapped_n = human_background$mapped_n,
              status = "human_orgdb_missing",
              reason = "org.Hs.eg.db is not installed and clusterProfiler backend is active; human auxiliary GO channel skipped",
              input_gene_n = length(genes),
              mapped_gene_n = 0L,
              mapping_rate = NA_real_,
              significant_term_n = 0L,
              paths = paths
            )
          } else {
            human_map <- map_genes_to_human(genes, ortholog_map)
            manifest_rows[[length(manifest_rows) + 1]] <- write_go_result(
              layer_id = layer_id,
              comparison_id = comparison_id,
              cluster_id = cluster_id,
              direction = direction,
              ont = ont,
              source_species = "human",
              deg_source = item$deg_source[[1]],
              deg_status = item$deg_inference_status[[1]],
              genes = human_map$genes,
              org_db = human_org_db,
              enrichment_source = "ortholog_to_org.Hs.eg.db",
              target_id = target_id,
              formal_status = formal_status,
              result_level = result_level,
              biological_replicates = biological_replicates,
              background_genes = human_background$genes,
              background_tsv = background_tsv
            )
          }
        }
      }
    }
  }
}

manifest_df <- if (length(manifest_rows) > 0) dplyr::bind_rows(manifest_rows) else empty_enrichment_manifest_06()
write_tsv_local(manifest_df, cfg$go_enrichment_manifest_tsv)

report_path <- file.path(cfg$enrichment_report_dir, "06a_go_report.md")
low_mapping <- manifest_df[!is.na(suppressWarnings(as.numeric(manifest_df$mapping_rate))) & suppressWarnings(as.numeric(manifest_df$mapping_rate)) < 0.3, , drop = FALSE]
report_lines <- c(
  "# 06a GO Enrichment",
  "",
  sprintf("- manifest: `%s`", cfg$go_enrichment_manifest_tsv),
  sprintf("- species_strategy: `%s`", cfg$enrichment_species_strategy),
  sprintf("- pvalue_cutoff: `%s`", cfg$enrichment_pvalue_cutoff),
  sprintf("- qvalue_cutoff: `%s`", cfg$enrichment_qvalue_cutoff),
  "",
  "## Status Summary",
  render_markdown_table_local(as.data.frame(table(manifest_df$status), stringsAsFactors = FALSE)),
  "",
  "## Low Mapping Rows",
  if (nrow(low_mapping) == 0) "No rows below 30% ID mapping." else render_markdown_table_local(head(low_mapping[, c("layer_id", "comparison_id", "cluster_id", "gene_direction", "ontology", "source_species", "mapping_rate", "status"), drop = FALSE], 100))
)
write_markdown_local(report_lines, report_path)

if (file.exists(cfg$module_06a_manifest_path)) {
  unlink(cfg$module_06a_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_06a_manifest_path,
  new_outputs = list(
    go_enrichment_manifest_tsv = build_output_entry(cfg$go_enrichment_manifest_tsv, "tsv", module_name, "one row per layer/comparison/cluster/direction/ontology/species GO result", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    report = build_output_entry(report_path, "md", module_name, "GO enrichment report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    enrichment_targets_tsv = cfg$enrichment_targets_sheet,
    gene_program_registry_tsv = cfg$gene_program_registry_tsv,
    ortholog_manifest = cfg$ortholog_manifest_path
  ),
  version = cfg$module_version,
  depends_on = list(
    module_05a = cfg$module_05a_manifest_path,
    module_05b = cfg$module_05b_manifest_path,
    module_05d = cfg$module_05d_manifest_path,
    module_00 = cfg$ortholog_manifest_path
  )
)

message("06a completed. manifest: ", cfg$go_enrichment_manifest_tsv)
