.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})

source_utf8 <- function(path) source(path, encoding = "UTF-8")
source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_05.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_06.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_07.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_mapping_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_08.R"))

load_required_packages(c("dplyr", "tidyr", "tibble", "ggplot2", "jsonlite"))

cfg <- get_single_script_config_08()
module_name <- "08f_regulation_eda"
prepare_dirs_08(cfg)

file_available_08f <- function(path) {
  nzchar(path) && file.exists(path) && !isTRUE(file.info(path)$size == 0)
}

read_csv_optional_08f <- function(path) {
  if (!file_available_08f(path)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

status_row_08f <- function(source_module, input_name, path, required = "no", reason = "") {
  available <- file_available_08f(path)
  data.frame(
    source_module = source_module,
    input_name = input_name,
    path = normalize_path_07(path),
    available = ifelse(available, "yes", "no"),
    status = ifelse(available, "ok", ifelse(identical(required, "yes"), "missing_required", "not_available")),
    reason = ifelse(available, "", reason),
    stringsAsFactors = FALSE
  )
}

clean_tf_name_08f <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("_top[0-9]+$", "", x)
  x <- sub("\\s*\\(.*\\)$", "", x)
  toupper(x)
}

first_existing_col_08f <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) "" else hit[[1]]
}

summarize_scenic_08f <- function(cfg) {
  idx <- read_tsv_optional(cfg$scenic_downstream_index_tsv)
  if (nrow(idx) == 0) {
    return(empty_df_07(c("layer_id", "cell_n", "regulon_n", "top_regulon_n", "status")))
  }
  rows <- lapply(seq_len(nrow(idx)), function(i) {
    row <- idx[i, , drop = FALSE]
    top_path <- normalize_scalar_value(row$top_regulons_by_celltype_csv[[1]])
    top_df <- read_csv_optional_08f(top_path)
    data.frame(
      layer_id = normalize_scalar_value(row$layer_id[[1]]),
      cell_n = suppressWarnings(as.integer(row$cell_n[[1]])),
      regulon_n = suppressWarnings(as.integer(row$regulon_n[[1]])),
      top_regulon_n = if ("regulon" %in% colnames(top_df)) length(unique(top_df$regulon)) else 0L,
      status = "ok",
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

summarize_decoupler_08f <- function(cfg) {
  idx <- read_tsv_optional(cfg$decoupler_index_tsv)
  tf <- read_tsv_optional(cfg$decoupler_tf_activity_tsv)
  pathway <- read_tsv_optional(cfg$decoupler_pathway_activity_tsv)
  if (nrow(idx) == 0) {
    return(empty_df_07(c("layer_id", "group_n", "gene_n", "tf_source_n", "pathway_source_n", "status")))
  }
  rows <- lapply(seq_len(nrow(idx)), function(i) {
    layer_id <- normalize_scalar_value(idx$layer_id[[i]])
    tf_layer <- if (nrow(tf) > 0 && "layer_id" %in% colnames(tf)) tf[tf$layer_id == layer_id & tf$status == "ok", , drop = FALSE] else tf[FALSE, , drop = FALSE]
    pathway_layer <- if (nrow(pathway) > 0 && "layer_id" %in% colnames(pathway)) pathway[pathway$layer_id == layer_id & pathway$status == "ok", , drop = FALSE] else pathway[FALSE, , drop = FALSE]
    data.frame(
      layer_id = layer_id,
      group_n = suppressWarnings(as.integer(idx$group_n[[i]])),
      gene_n = suppressWarnings(as.integer(idx$gene_n[[i]])),
      tf_source_n = length(unique(tf_layer$source_id)),
      pathway_source_n = length(unique(pathway_layer$source_id)),
      status = "ok",
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

build_tf_overlap_08f <- function(cfg) {
  scenic_idx <- read_tsv_optional(cfg$scenic_downstream_index_tsv)
  decoupler_tf <- read_tsv_optional(cfg$decoupler_tf_activity_tsv)
  if (nrow(scenic_idx) == 0 || nrow(decoupler_tf) == 0) {
    return(empty_df_07(c("layer_id", "scenic_tf_n", "decoupler_tf_n", "overlap_n", "jaccard", "overlap_tf")))
  }
  layers <- intersect(unique(scenic_idx$layer_id), unique(decoupler_tf$layer_id))
  rows <- lapply(layers, function(layer_id) {
    scenic_row <- scenic_idx[scenic_idx$layer_id == layer_id, , drop = FALSE][1, , drop = FALSE]
    top_path <- normalize_scalar_value(scenic_row$top_regulons_by_celltype_csv[[1]])
    scenic_top <- read_csv_optional_08f(top_path)
    scenic_col <- first_existing_col_08f(scenic_top, c("regulon", "TF", "source_id"))
    scenic_tf <- if (nzchar(scenic_col)) clean_tf_name_08f(scenic_top[[scenic_col]]) else character(0)
    scenic_tf <- unique(scenic_tf[nzchar(scenic_tf)])

    dec_layer <- decoupler_tf[decoupler_tf$layer_id == layer_id & decoupler_tf$status == "ok", , drop = FALSE]
    if (nrow(dec_layer) > 0) {
      dec_rank <- stats::aggregate(abs(score) ~ source_id, data = dec_layer, FUN = max)
      names(dec_rank) <- c("source_id", "max_abs_score")
      dec_rank <- dec_rank[order(-dec_rank$max_abs_score, dec_rank$source_id), , drop = FALSE]
      dec_tf <- clean_tf_name_08f(head(dec_rank$source_id, cfg$decoupler_top_tf_n))
    } else {
      dec_tf <- character(0)
    }
    dec_tf <- unique(dec_tf[nzchar(dec_tf)])
    overlap <- intersect(scenic_tf, dec_tf)
    union_n <- length(union(scenic_tf, dec_tf))
    data.frame(
      layer_id = layer_id,
      scenic_tf_n = length(scenic_tf),
      decoupler_tf_n = length(dec_tf),
      overlap_n = length(overlap),
      jaccard = if (union_n > 0) length(overlap) / union_n else NA_real_,
      overlap_tf = paste(overlap, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  if (length(rows) == 0) {
    empty_df_07(c("layer_id", "scenic_tf_n", "decoupler_tf_n", "overlap_n", "jaccard", "overlap_tf"))
  } else {
    dplyr::bind_rows(rows)
  }
}

plot_tf_overlap_08f <- function(overlap_df, path) {
  ensure_dir(dirname(path))
  if (nrow(overlap_df) == 0) {
    write_empty_png_07(path, title = "No SCENIC/decoupleR overlap available")
    return(path)
  }
  plot_df <- overlap_df
  plot_df$jaccard[!is.finite(plot_df$jaccard)] <- 0
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = layer_id, y = jaccard)) +
    ggplot2::geom_col(fill = "#2F6CB3", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = overlap_n), vjust = -0.2, size = 3) +
    ggplot2::coord_cartesian(ylim = c(0, max(0.05, max(plot_df$jaccard, na.rm = TRUE) * 1.15))) +
    ggplot2::labs(x = "Layer", y = "Jaccard", title = "SCENIC vs decoupleR TF overlap") +
    ggplot2::theme_classic(base_size = 12)
  ggplot2::ggsave(path, p, width = 7, height = 4, dpi = 300, bg = "white")
  path
}

summarize_resource_fingerprint_08f <- function(cfg) {
  rows <- list()
  if (file_available_08f(cfg$scenic_resource_manifest)) {
    scenic_manifest <- jsonlite::read_json(cfg$scenic_resource_manifest, simplifyVector = FALSE)
    resources <- scenic_manifest$resources %||% list()
    for (resource_name in names(resources)) {
      item <- resources[[resource_name]]
      rows[[length(rows) + 1L]] <- data.frame(
        resource = paste0("scenic_", resource_name),
        source = "cache",
        species_origin = "human",
        mapped_species = "human_symbol_expression",
        row_n_raw = NA_integer_,
        row_n_mapped = NA_integer_,
        package_version = "",
        resource_md5 = item$sha256 %||% "",
        stringsAsFactors = FALSE
      )
    }
  }
  dec <- read_tsv_optional(cfg$decoupler_resource_fingerprint_tsv)
  if (nrow(dec) > 0) {
    rows[[length(rows) + 1L]] <- dec[, intersect(colnames(empty_decoupler_resource_fingerprint_08_local()), colnames(dec)), drop = FALSE]
  }
  if (length(rows) == 0) {
    return(empty_decoupler_resource_fingerprint_08_local())
  }
  out <- dplyr::bind_rows(rows)
  for (col in colnames(empty_decoupler_resource_fingerprint_08_local())) {
    if (!col %in% colnames(out)) {
      out[[col]] <- ""
    }
  }
  out[, colnames(empty_decoupler_resource_fingerprint_08_local()), drop = FALSE]
}

empty_decoupler_resource_fingerprint_08_local <- function() {
  empty_df_07(c(
    "resource", "source", "species_origin", "mapped_species",
    "row_n_raw", "row_n_mapped", "package_version", "resource_md5"
  ))
}

build_ortholog_coverage_08f <- function(cfg) {
  rows <- list()
  scenic_map <- read_tsv_optional(cfg$scenic_export_mapping_summary_tsv)
  if (nrow(scenic_map) > 0) {
    rows[[length(rows) + 1L]] <- data.frame(
      method = "scenic_expression_to_human",
      source_id = scenic_map$layer_id,
      input_n = scenic_map$chicken_gene_n,
      mapped_n = scenic_map$mapped_chicken_gene_n,
      output_n = scenic_map$human_gene_n,
      mapping_rate = scenic_map$mapping_rate,
      stringsAsFactors = FALSE
    )
  }
  dec_map <- read_tsv_optional(cfg$decoupler_network_mapping_summary_tsv)
  if (nrow(dec_map) > 0) {
    rows[[length(rows) + 1L]] <- data.frame(
      method = paste0("decoupler_", dec_map$resource, "_network_to_chicken"),
      source_id = dec_map$source_id,
      input_n = dec_map$human_target_n,
      mapped_n = dec_map$mapped_chicken_target_n,
      output_n = dec_map$mapped_chicken_target_n,
      mapping_rate = dec_map$mapping_rate,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) {
    empty_df_07(c("method", "source_id", "input_n", "mapped_n", "output_n", "mapping_rate"))
  } else {
    dplyr::bind_rows(rows)
  }
}

scenic_ok <- file_available_08f(cfg$scenic_downstream_index_tsv)
decoupler_ok <- file_available_08f(cfg$decoupler_index_tsv)
if (!scenic_ok && !decoupler_ok && !identical(tolower(cfg$allow_regulation_empty_report), "yes")) {
  stop("08f requires at least one completed regulation branch unless ALLOW_REGULATION_EMPTY_REPORT=yes.", call. = FALSE)
}
report_mode <- if (scenic_ok && decoupler_ok) {
  "full"
} else if (scenic_ok) {
  "scenic_only_partial"
} else if (decoupler_ok) {
  "decoupler_only_partial"
} else {
  "empty_allowed"
}

module_status <- dplyr::bind_rows(
  status_row_08f("08a_scenic_export", "manifest", cfg$module_08a_manifest_path, reason = "SCENIC export was not run."),
  status_row_08f("08b_scenic_grn", "manifest", cfg$module_08b_manifest_path, reason = "SCENIC GRN was skipped or failed."),
  status_row_08f("08c_scenic_regulons", "manifest", cfg$module_08c_manifest_path, reason = "SCENIC regulon inference was skipped or failed."),
  status_row_08f("08d_scenic_downstream", "manifest", cfg$module_08d_manifest_path, reason = "SCENIC downstream was skipped or failed."),
  status_row_08f("08e_decoupler", "manifest", cfg$module_08e_manifest_path, reason = "decoupleR was skipped or failed."),
  status_row_08f("05d_deg_eda", "gene_program_registry", cfg$gene_program_registry_tsv, reason = "05 gene_program_registry.tsv not available; upstream DEG linkage omitted."),
  status_row_08f("06c_enrichment_eda", "enrichment_section_summary", cfg$enrichment_section_summary_tsv, reason = "06 enrichment_section_summary.tsv not available; pathway evidence linkage omitted."),
  status_row_08f("07e_communication_eda", "communication_gate_summary", cfg$communication_gate_summary_tsv, reason = "07 communication_gate_summary.tsv not available."),
  status_row_08f("07e_communication_eda", "skipped_low_cells", cfg$communication_skipped_low_cells_tsv, reason = "07 skipped_low_cells.tsv not available."),
  status_row_08f("07e_communication_eda", "communication_fallback_summary", cfg$communication_fallback_summary_tsv, reason = "07 communication_fallback_summary.tsv not available."),
  status_row_08f("07e_communication_eda", "derived_communication_eligibility", cfg$derived_communication_eligibility_tsv, reason = "07 derived_communication_eligibility.tsv not available."),
  status_row_08f("07c_nichenet", "nichenet_index", cfg$nichenet_index_tsv, reason = "07 NicheNet index not available; missing_gene_program linkage omitted.")
)
write_tsv_local(module_status, cfg$regulation_module_status_tsv)

scenic_summary <- summarize_scenic_08f(cfg)
decoupler_summary <- summarize_decoupler_08f(cfg)
resource_fingerprint <- summarize_resource_fingerprint_08f(cfg)
ortholog_coverage <- build_ortholog_coverage_08f(cfg)
tf_overlap <- build_tf_overlap_08f(cfg)
plot_tf_overlap_08f(tf_overlap, cfg$regulation_tf_overlap_png)

triage_rows <- list()
missing_optional <- module_status[module_status$available == "no" & module_status$source_module %in% c("05d_deg_eda", "06c_enrichment_eda", "07e_communication_eda", "07c_nichenet"), , drop = FALSE]
if (nrow(missing_optional) > 0) {
  triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
    sample_id = "__PROJECT__",
    severity = "info",
    signal_id = "regulation_upstream_summary_not_available",
    suspected_issue = "08f upstream evidence section is partial.",
    evidence = paste(missing_optional$path, collapse = "; "),
    recommended_action = "Run or review 05/06/07 summaries before using 08 as final mechanistic interpretation.",
    manual_review_required = "no"
  )
}
if (nrow(tf_overlap) > 0 && any(tf_overlap$overlap_n == 0)) {
  triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
    sample_id = "__PROJECT__",
    severity = "info",
    signal_id = "regulation_tf_overlap_absent",
    suspected_issue = "SCENIC and decoupleR top TFs do not overlap for at least one layer.",
    evidence = paste(tf_overlap$layer_id[tf_overlap$overlap_n == 0], collapse = ","),
    recommended_action = "Treat this as method disagreement, not as proof of failure; inspect target coverage and upstream DEG/enrichment/communication evidence.",
    manual_review_required = "yes"
  )
}
triage <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)
write_tsv_local(resource_fingerprint, cfg$regulation_resource_fingerprint_tsv)
write_tsv_local(ortholog_coverage, cfg$regulation_ortholog_coverage_tsv)
write_tsv_local(scenic_summary, cfg$regulation_scenic_summary_tsv)
write_tsv_local(decoupler_summary, cfg$regulation_decoupler_summary_tsv)
write_tsv_local(tf_overlap, cfg$regulation_tf_overlap_tsv)
write_tsv_local(triage, cfg$regulation_triage_tsv)

gene_program_registry <- read_tsv_optional(cfg$gene_program_registry_tsv)
enrichment_sections <- read_tsv_optional(cfg$enrichment_section_summary_tsv)
communication_gate <- read_tsv_optional(cfg$communication_gate_summary_tsv)
skipped_low_cells <- read_tsv_optional(cfg$communication_skipped_low_cells_tsv)
fallback_summary <- read_tsv_optional(cfg$communication_fallback_summary_tsv)
derived_eligibility <- read_tsv_optional(cfg$derived_communication_eligibility_tsv)
nichenet_index <- read_tsv_optional(cfg$nichenet_index_tsv)
missing_gene_program <- if (nrow(nichenet_index) > 0 && "status" %in% colnames(nichenet_index)) {
  nichenet_index[nichenet_index$status == "missing_gene_program", , drop = FALSE]
} else {
  nichenet_index[FALSE, , drop = FALSE]
}

report_lines <- c(
  "# 08 Regulation EDA",
  "",
  sprintf("- Report mode: `%s`.", report_mode),
  "- SCENIC and decoupleR are complementary evidence streams, not equivalent methods.",
  "- Overlap between the two methods can strengthen confidence, but lack of overlap is not automatically an error.",
  "- SCENIC depends on motif/regulon inference and human cisTarget resources.",
  "- decoupleR depends on curated human prior networks mapped to chicken genes.",
  "- 08 TF/pathway activity must not replace 05 DEG, 06 enrichment, or 07 communication evidence.",
  "",
  "## Module Status",
  render_markdown_table_local(module_status),
  "",
  "## Resource Fingerprint",
  render_markdown_table_local(resource_fingerprint),
  "",
  "## Ortholog Mapping Coverage",
  render_markdown_table_local(head(ortholog_coverage, 80)),
  "",
  "## SCENIC Summary",
  if (nrow(scenic_summary) == 0) "SCENIC outputs are not available in this partial report." else render_markdown_table_local(scenic_summary),
  "",
  "## decoupleR Summary",
  if (nrow(decoupler_summary) == 0) "decoupleR outputs are not available in this partial report." else render_markdown_table_local(decoupler_summary),
  "",
  "## SCENIC vs decoupleR TF Overlap",
  if (nrow(tf_overlap) == 0) "No paired SCENIC/decoupleR TF overlap could be computed." else render_markdown_table_local(tf_overlap),
  "",
  "## Pathway Activity Summary",
  if (nrow(decoupler_summary) == 0) "Pathway activity is unavailable without decoupleR." else render_markdown_table_local(decoupler_summary[, intersect(c("layer_id", "pathway_source_n", "group_n", "status"), colnames(decoupler_summary)), drop = FALSE]),
  "",
  "## Link to 05/06/07 Evidence",
  "05/06/07 inputs are soft optional. Missing files are recorded as `not_available` and do not block this report.",
  "",
  "### 05 Gene Program Registry",
  if (nrow(gene_program_registry) == 0) sprintf("not_available: `%s`", cfg$gene_program_registry_tsv) else render_markdown_table_local(head(gene_program_registry, 40)),
  "",
  "### 06 Enrichment Section Summary",
  if (nrow(enrichment_sections) == 0) sprintf("not_available: `%s`", cfg$enrichment_section_summary_tsv) else render_markdown_table_local(head(enrichment_sections, 60)),
  "",
  "### 07 Communication Gate Summary",
  if (nrow(communication_gate) == 0) sprintf("not_available: `%s`", cfg$communication_gate_summary_tsv) else render_markdown_table_local(head(communication_gate, 60)),
  "",
  "### 07 Skipped Low Cells",
  if (nrow(skipped_low_cells) == 0) sprintf("not_available or empty: `%s`", cfg$communication_skipped_low_cells_tsv) else render_markdown_table_local(head(skipped_low_cells, 60)),
  "",
  "### 07 Fallback / Derived Eligibility",
  if (nrow(fallback_summary) == 0) sprintf("fallback summary not_available or empty: `%s`", cfg$communication_fallback_summary_tsv) else render_markdown_table_local(head(fallback_summary, 40)),
  if (nrow(derived_eligibility) == 0) sprintf("derived eligibility not_available or empty: `%s`", cfg$derived_communication_eligibility_tsv) else render_markdown_table_local(head(derived_eligibility, 40)),
  "",
  "### 07 NicheNet missing_gene_program",
  if (nrow(missing_gene_program) == 0) "No missing_gene_program rows found, or NicheNet index is not available." else render_markdown_table_local(head(missing_gene_program, 60)),
  "",
  "## Triage",
  render_markdown_table_local(triage),
  "",
  "## Recommended Interpretation",
  "- Use 08 as regulation evidence layered on top of 05/06/07, not as a replacement.",
  "- Prefer TF/pathway claims supported by 05 DEG or markers, 06 enrichment, and 07 communication status when available.",
  "- For SCENIC-only or decoupleR-only partial reports, label conclusions as provisional until the missing branch is run or explicitly waived."
)
write_markdown_local(report_lines, cfg$regulation_report_md)

if (file.exists(cfg$module_08f_manifest_path)) {
  unlink(cfg$module_08f_manifest_path)
}
outputs <- list(
  report_md = build_output_entry(cfg$regulation_report_md, "md", module_name, "08 regulation EDA report", base_dir = cfg$project_root),
  module_status_tsv = build_output_entry(cfg$regulation_module_status_tsv, "tsv", module_name, "08 module and upstream input availability", base_dir = cfg$project_root, schema = infer_schema_from_df(module_status)),
  resource_fingerprint_tsv = build_output_entry(cfg$regulation_resource_fingerprint_tsv, "tsv", module_name, "SCENIC and decoupleR resource fingerprints", base_dir = cfg$project_root, schema = infer_schema_from_df(resource_fingerprint)),
  ortholog_mapping_coverage_tsv = build_output_entry(cfg$regulation_ortholog_coverage_tsv, "tsv", module_name, "ortholog mapping coverage by method/resource", base_dir = cfg$project_root, schema = infer_schema_from_df(ortholog_coverage)),
  scenic_summary_tsv = build_output_entry(cfg$regulation_scenic_summary_tsv, "tsv", module_name, "SCENIC regulation summary", base_dir = cfg$project_root, schema = infer_schema_from_df(scenic_summary)),
  decoupler_summary_tsv = build_output_entry(cfg$regulation_decoupler_summary_tsv, "tsv", module_name, "decoupleR regulation summary", base_dir = cfg$project_root, schema = infer_schema_from_df(decoupler_summary)),
  scenic_decoupler_tf_overlap_tsv = build_output_entry(cfg$regulation_tf_overlap_tsv, "tsv", module_name, "SCENIC and decoupleR TF overlap", base_dir = cfg$project_root, schema = infer_schema_from_df(tf_overlap)),
  triage_tsv = build_output_entry(cfg$regulation_triage_tsv, "tsv", module_name, "08 regulation triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage)),
  scenic_decoupler_tf_overlap_png = build_output_entry(cfg$regulation_tf_overlap_png, "png", module_name, "SCENIC and decoupleR TF overlap figure", base_dir = cfg$project_root)
)
write_manifest_local(
  manifest_path = cfg$module_08f_manifest_path,
  new_outputs = outputs,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_08a = cfg$module_08a_manifest_path,
    module_08b = cfg$module_08b_manifest_path,
    module_08c = cfg$module_08c_manifest_path,
    module_08d = cfg$module_08d_manifest_path,
    module_08e = cfg$module_08e_manifest_path,
    gene_program_registry_tsv = cfg$gene_program_registry_tsv,
    enrichment_section_summary_tsv = cfg$enrichment_section_summary_tsv,
    communication_gate_summary_tsv = cfg$communication_gate_summary_tsv,
    communication_skipped_low_cells_tsv = cfg$communication_skipped_low_cells_tsv,
    communication_fallback_summary_tsv = cfg$communication_fallback_summary_tsv,
    derived_communication_eligibility_tsv = cfg$derived_communication_eligibility_tsv,
    nichenet_index_tsv = cfg$nichenet_index_tsv
  ),
  version = cfg$module_version,
  depends_on = list(
    module_08a = cfg$module_08a_manifest_path,
    module_08b = cfg$module_08b_manifest_path,
    module_08c = cfg$module_08c_manifest_path,
    module_08d = cfg$module_08d_manifest_path,
    module_08e = cfg$module_08e_manifest_path,
    module_05d = cfg$module_05d_manifest_path,
    module_06c = cfg$module_06c_manifest_path,
    module_07e = cfg$module_07e_manifest_path
  )
)

message("08f completed. report: ", cfg$regulation_report_md)
