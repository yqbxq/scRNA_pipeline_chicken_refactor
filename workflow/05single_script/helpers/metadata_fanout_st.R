m3_deconv_cols <- c(
  "deconv_id", "source_question_id", "st_scope", "reference_scope",
  "section_filter", "condition_split_var", "condition_split_values",
  "tool", "enabled", "notes"
)

m3_spatial_cols <- c(
  "spatial_pair_id", "source_question_id", "st_scope", "sender", "receiver",
  "contrast_axis", "section_filter", "condition_split_var", "condition_split_values",
  "tool", "enabled", "notes"
)

m3_st_m2_question_ids <- function() {
  c(
    "I01_ST_clustering", "I02_ST_region_marker", "I06_SVG",
    "I07_SVG_by_stage", "I08_deconv_panorama",
    "I12_ST_region_compo", "I13_neighborhood", "I15_ST_enrichment"
  )
}

m3_st_enabled <- function() {
  value <- suppressWarnings(as.integer(Sys.getenv("ST_FANOUT_VERSION", "0")))
  !is.na(value) && value >= 1L
}

m3_st_empty <- function() {
  list(
    comparisons = m3_empty_df(m3_comparison_cols),
    spatial_pairs = m3_empty_df(m3_spatial_cols),
    deconv_pairs = m3_empty_df(m3_deconv_cols),
    enrichment_targets = m3_empty_df(m3_enrichment_cols)
  )
}

m3_st_sections <- function() {
  metadata_dir <- Sys.getenv("METADATA_DIR", unset = "metadata")
  candidates <- unique(c(
    Sys.getenv("SECTION_SHEET", unset = ""),
    file.path(metadata_dir, "sections.tsv")
  ))
  sections <- data.frame(stringsAsFactors = FALSE)
  for (path in candidates[nzchar(candidates)]) {
    sections <- m3_read_tsv(path)
    if (nrow(sections) > 0 && "section_id" %in% colnames(sections)) {
      break
    }
  }
  if (nrow(sections) == 0) {
    sections <- data.frame(
      section_id = c("syf_1", "f5_1"),
      condition = c("syf", "f5"),
      enabled = c("yes", "yes"),
      stringsAsFactors = FALSE
    )
  }
  for (col in c("section_id", "condition", "enabled")) {
    if (!col %in% colnames(sections)) {
      sections[[col]] <- ""
    }
  }
  sections$enabled <- tolower(m3_trim(sections$enabled, "yes"))
  sections <- sections[nzchar(m3_trim(sections$section_id)) & sections$enabled != "no", , drop = FALSE]
  if (nrow(sections) == 0) {
    sections <- data.frame(
      section_id = c("syf_1", "f5_1"),
      condition = c("syf", "f5"),
      enabled = c("yes", "yes"),
      stringsAsFactors = FALSE
    )
  }
  sections
}

m3_st_condition_fields <- function(condition_split) {
  fields <- m3_condition_fields(condition_split)
  if (!nzchar(fields$var)) {
    return(list(var = "", values = ""))
  }
  if (identical(fields$var, "section")) {
    fields$var <- "condition"
  }
  fields
}

m3_st_comparison_row <- function(q, suffix, analysis_mode, group_var,
                                 ident_1, ident_2, composition_group_var = "",
                                 produces_gene_program = "yes",
                                 gene_program_role = "receiver_marker",
                                 notes = "") {
  display_question_id <- m3_question_scalar(q, "display_question_id", q$question_id[[1]])
  output_alias <- m3_question_scalar(q, "output_alias", display_question_id)
  report_title <- m3_question_scalar(q, "report_title", m3_question_scalar(q, "question_zh", display_question_id))
  list(
    comparison_id = paste(q$question_id[[1]], m3_safe_id(suffix), sep = "__"),
    source_question_id = q$question_id[[1]],
    display_question_id = display_question_id,
    output_alias = output_alias,
    report_title = report_title,
    layer_scope = "panorama_st",
    contrast_axis = q$contrast_axis[[1]],
    analysis_mode = analysis_mode,
    analysis_unit = ifelse(identical(analysis_mode, "composition"), "sample_level", "spatial_region"),
    stat_level = ifelse(identical(analysis_mode, "composition"), "composition_formal_if_replicates", "spot_level_exploratory"),
    group_var = group_var,
    ident_1 = ident_1,
    ident_2 = ident_2,
    subset_column = "",
    subset_value = "",
    aggregation_group_var = "",
    composition_group_var = composition_group_var,
    batch_var = "batch",
    enabled = "yes",
    min_biological_replicates = "1",
    force_exploratory = "yes",
    min_cells_per_group = "3",
    logfc_threshold = "0",
    produces_gene_program = produces_gene_program,
    gene_program_role = gene_program_role,
    notes = notes
  )
}

m3_st_spatial_row <- function(q, suffix, sender, receiver, tool,
                              section_filter = "*", split_var = "",
                              split_values = "", notes = "") {
  list(
    spatial_pair_id = paste(q$question_id[[1]], m3_safe_id(suffix), sep = "__"),
    source_question_id = q$question_id[[1]],
    st_scope = "panorama_st",
    sender = sender,
    receiver = receiver,
    contrast_axis = q$contrast_axis[[1]],
    section_filter = section_filter,
    condition_split_var = split_var,
    condition_split_values = split_values,
    tool = tool,
    enabled = "yes",
    notes = notes
  )
}

m3_st_deconv_row <- function(q, suffix, reference_scope, tool,
                             section_filter = "*", split_var = "",
                             split_values = "", notes = "") {
  list(
    deconv_id = paste(q$question_id[[1]], m3_safe_id(suffix), sep = "__"),
    source_question_id = q$question_id[[1]],
    st_scope = "panorama_st",
    reference_scope = reference_scope,
    section_filter = section_filter,
    condition_split_var = split_var,
    condition_split_values = split_values,
    tool = tool,
    enabled = "yes",
    notes = notes
  )
}

m3_st_enrichment_rows <- function(q) {
  comparison_id <- "I02_ST_region_marker__panorama_st"
  rows <- list()
  for (database in c("GO", "KEGG")) {
    rows[[length(rows) + 1L]] <- list(
      target_id = paste(q$question_id[[1]], comparison_id, tolower(database), sep = "__"),
      source_question_id = q$question_id[[1]],
      comparison_id = comparison_id,
      layer_scope = "panorama_st",
      analysis_mode = "subtype_marker",
      gene_program_role = "receiver_marker",
      organism = "chicken_primary",
      database = database,
      enrichment_eligible = "yes",
      enrichment_usage = "identity_baseline_enrichment",
      min_genes = "5",
      enabled = "yes",
      notes = "derived_from_spatial_region_marker_gene_program"
    )
  }
  m3_bind_rows(rows, m3_enrichment_cols)
}

m3_fanout_st <- function(q) {
  if (!m3_st_enabled()) {
    stop(
      sprintf(
        "ST fan-out is disabled (ST_FANOUT_VERSION=0; question_id=%s, contrast_axis=%s). Keep ST_section rows planned or set ST_FANOUT_VERSION=1.",
        q$question_id[[1]],
        q$contrast_axis[[1]]
      ),
      call. = FALSE
    )
  }

  if (!q$question_id[[1]] %in% m3_st_m2_question_ids()) {
    return(m3_st_empty())
  }

  spatial_rows <- list()
  deconv_rows <- list()
  comparison_rows <- list()
  enrichment_targets <- m3_empty_df(m3_enrichment_cols)
  qid <- q$question_id[[1]]
  split <- m3_st_condition_fields(q$condition_split[[1]])
  sections <- m3_st_sections()

  if (identical(qid, "I01_ST_clustering")) {
    spatial_rows[[length(spatial_rows) + 1L]] <- m3_st_spatial_row(
      q, "panorama_st", "*", "*", "spatial_clustering",
      notes = "M2 spatial clustering entry point"
    )
  } else if (identical(qid, "I02_ST_region_marker")) {
    comparison_rows[[length(comparison_rows) + 1L]] <- m3_st_comparison_row(
      q, "panorama_st", "subtype_marker", "spatial_region",
      "__region__", "__rest__",
      notes = "spatial region marker discovery on panorama_st"
    )
  } else if (identical(qid, "I06_SVG")) {
    for (idx in seq_len(nrow(sections))) {
      section_id <- sections$section_id[[idx]]
      spatial_rows[[length(spatial_rows) + 1L]] <- m3_st_spatial_row(
        q, paste(section_id, "sparkx", sep = "__"), "all_spots", "-",
        "sparkx", section_filter = section_id,
        notes = "per-section spatially variable genes"
      )
      spatial_rows[[length(spatial_rows) + 1L]] <- m3_st_spatial_row(
        q, paste(section_id, "spatialde2", sep = "__"), "all_spots", "-",
        "spatialde2", section_filter = section_id,
        notes = "per-section spatially variable genes"
      )
    }
  } else if (identical(qid, "I07_SVG_by_stage")) {
    spatial_rows[[length(spatial_rows) + 1L]] <- m3_st_spatial_row(
      q, "condition_split", "all_spots", "-", "sparkx",
      section_filter = "*", split_var = split$var, split_values = split$values,
      notes = "condition-aware SVG comparison"
    )
  } else if (identical(qid, "I08_deconv_panorama")) {
    deconv_rows[[length(deconv_rows) + 1L]] <- m3_st_deconv_row(
      q, "panorama", "panorama", "rctd,cell2location",
      section_filter = "*", split_var = split$var, split_values = split$values,
      notes = "panorama reference deconvolution into ST sections"
    )
  } else if (identical(qid, "I12_ST_region_compo")) {
    condition <- m3_condition_pair(q$condition_split[[1]])
    if (identical(condition$group_var, "section")) {
      condition$group_var <- "condition"
    }
    comparison_rows[[length(comparison_rows) + 1L]] <- m3_st_comparison_row(
      q, "region_composition", "composition", condition$group_var,
      condition$ident_1, condition$ident_2,
      composition_group_var = "spatial_region",
      produces_gene_program = "no",
      gene_program_role = "none",
      notes = "ST region composition by condition; N=2 remains exploratory"
    )
  } else if (identical(qid, "I13_neighborhood")) {
    spatial_rows[[length(spatial_rows) + 1L]] <- m3_st_spatial_row(
      q, "region_region", "region", "region", "squidpy",
      notes = "region-level spatial neighborhood analysis"
    )
  } else if (identical(qid, "I15_ST_enrichment")) {
    enrichment_targets <- m3_st_enrichment_rows(q)
  }

  list(
    comparisons = m3_bind_rows(comparison_rows, m3_comparison_cols),
    spatial_pairs = m3_bind_rows(spatial_rows, m3_spatial_cols),
    deconv_pairs = m3_bind_rows(deconv_rows, m3_deconv_cols),
    enrichment_targets = enrichment_targets
  )
}
