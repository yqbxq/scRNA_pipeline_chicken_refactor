m3_enrichment_cols <- c(
  "target_id", "source_question_id", "comparison_id", "layer_scope",
  "analysis_mode", "gene_program_role", "organism", "database",
  "enrichment_eligible", "enrichment_usage", "min_genes", "enabled", "notes"
)

m3_gene_program_cols <- c(
  "comparison_id", "source_question_id", "layer_scope", "analysis_mode",
  "gene_program_role", "produces_gene_program", "annotation_only", "qc_only",
  "global_context_only", "nichenet_eligible", "nichenet_usage",
  "enrichment_eligible", "enrichment_usage", "preferred_for_downstream",
  "expected_result_level", "formal_preferred", "formal_status",
  "result_status", "skip_reason", "eligible_reason", "ineligible_reason",
  "notes"
)

m3_fanout_enrichment <- function(questions, comparisons) {
  if (nrow(comparisons) == 0) {
    return(m3_empty_df(m3_enrichment_cols))
  }
  enrichment_questions <- questions[
    questions$status == "active" & grepl("enrichment", questions$tools_to_run, ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  if (nrow(enrichment_questions) == 0) {
    return(m3_empty_df(m3_enrichment_cols))
  }
  source_ids <- enrichment_questions$question_id
  cmp <- comparisons[comparisons$source_question_id %in% source_ids, , drop = FALSE]
  if (nrow(cmp) == 0) {
    return(m3_empty_df(m3_enrichment_cols))
  }

  rows <- list()
  for (idx in seq_len(nrow(cmp))) {
    row <- cmp[idx, , drop = FALSE]
    fields <- m3_gene_program_fields(row)
    if (!identical(fields$produces_gene_program, "yes") ||
        !fields$enrichment_eligible %in% c("yes", "contextual")) {
      next
    }
    for (database in c("GO", "KEGG")) {
      rows[[length(rows) + 1L]] <- list(
        target_id = paste(row$comparison_id[[1]], tolower(database), sep = "__"),
        source_question_id = row$source_question_id[[1]],
        comparison_id = row$comparison_id[[1]],
        layer_scope = row$layer_scope[[1]],
        analysis_mode = row$analysis_mode[[1]],
        gene_program_role = row$gene_program_role[[1]],
        organism = "chicken_primary",
        database = database,
        enrichment_eligible = fields$enrichment_eligible,
        enrichment_usage = fields$enrichment_usage,
        min_genes = "5",
        enabled = "yes",
        notes = "derived_from_05_gene_program"
      )
    }
  }
  m3_bind_rows(rows, m3_enrichment_cols)
}

m3_gene_program_fields <- function(row) {
  mode <- row$analysis_mode[[1]]
  role <- row$gene_program_role[[1]]
  produces <- if ("produces_gene_program" %in% colnames(row)) tolower(trimws(row$produces_gene_program[[1]])) else "yes"
  produces <- ifelse(identical(produces, "yes"), "yes", "no")

  base <- list(
    produces_gene_program = produces,
    annotation_only = "no",
    qc_only = "no",
    global_context_only = "no",
    nichenet_eligible = "no",
    nichenet_usage = "none",
    enrichment_eligible = "no",
    enrichment_usage = "none",
    preferred_for_downstream = "no",
    expected_result_level = "unavailable",
    formal_preferred = "no",
    formal_status = "",
    result_status = "target_defined",
    skip_reason = "",
    eligible_reason = "",
    ineligible_reason = ""
  )

  if (identical(mode, "annotation_cluster_marker") || identical(role, "annotation_marker")) {
    return(utils::modifyList(base, list(
      produces_gene_program = "yes",
      annotation_only = "yes",
      expected_result_level = "annotation_cluster_markers",
      formal_status = "annotation_only",
      result_status = "annotation_only",
      ineligible_reason = "annotation_marker_not_mechanism_gene_program"
    )))
  }
  if (identical(mode, "qc_composition") || identical(role, "qc_only")) {
    return(utils::modifyList(base, list(
      produces_gene_program = "no",
      qc_only = "yes",
      expected_result_level = "qc_composition",
      formal_status = "qc_only",
      result_status = "qc_only",
      skip_reason = "qc_composition_does_not_produce_gene_program",
      ineligible_reason = "qc_composition_does_not_produce_gene_program"
    )))
  }
  if (identical(mode, "composition") || identical(role, "none")) {
    return(utils::modifyList(base, list(
      produces_gene_program = "no",
      expected_result_level = "composition",
      formal_status = "not_gene_program",
      result_status = "not_gene_program",
      skip_reason = "composition_does_not_produce_gene_program",
      ineligible_reason = "composition_does_not_produce_gene_program"
    )))
  }
  if (identical(mode, "global_context") || identical(role, "global_context")) {
    return(utils::modifyList(base, list(
      produces_gene_program = "yes",
      global_context_only = "yes",
      enrichment_eligible = "contextual",
      enrichment_usage = "global_context_enrichment",
      preferred_for_downstream = "contextual",
      expected_result_level = "pseudobulk_formal_or_exploratory_fallback",
      formal_preferred = "yes",
      formal_status = "formal_preferred",
      eligible_reason = "global_context_for_contextual_enrichment_only",
      ineligible_reason = "global_context_not_cell_type_specific_receiver_deg"
    )))
  }
  if (identical(mode, "condition_within_type") && identical(role, "condition_deg")) {
    return(list(
      produces_gene_program = "yes",
      annotation_only = "no",
      qc_only = "no",
      global_context_only = "no",
      nichenet_eligible = "yes",
      nichenet_usage = "receiver_condition_deg",
      enrichment_eligible = "yes",
      enrichment_usage = "mechanism_enrichment",
      preferred_for_downstream = "yes",
      expected_result_level = "pseudobulk_formal_or_exploratory_fallback",
      formal_preferred = "yes",
      formal_status = "formal_preferred",
      result_status = "target_defined",
      skip_reason = "",
      eligible_reason = "condition_within_type_receiver_deg",
      ineligible_reason = ""
    ))
  }
  if (identical(mode, "subtype_marker") && identical(role, "receiver_marker")) {
    return(list(
      produces_gene_program = "yes",
      annotation_only = "no",
      qc_only = "no",
      global_context_only = "no",
      nichenet_eligible = "yes",
      nichenet_usage = "baseline_receiver_marker",
      enrichment_eligible = "yes",
      enrichment_usage = "identity_baseline_enrichment",
      preferred_for_downstream = "baseline_only",
      expected_result_level = "cell_level_exploratory",
      formal_preferred = "no",
      formal_status = "exploratory_marker",
      result_status = "target_defined",
      skip_reason = "",
      eligible_reason = "baseline_receiver_marker_only",
      ineligible_reason = ""
    ))
  }
  if (identical(mode, "subtype_pairwise") && identical(role, "subtype_pairwise_deg")) {
    return(utils::modifyList(base, list(
      produces_gene_program = "yes",
      enrichment_eligible = "yes",
      enrichment_usage = "subtype_pairwise_enrichment",
      preferred_for_downstream = "contextual",
      expected_result_level = "cell_level_exploratory",
      formal_status = "exploratory_marker",
      eligible_reason = "subtype_pairwise_enrichment_only",
      ineligible_reason = "subtype_pairwise_not_default_receiver_deg"
    )))
  }
  utils::modifyList(base, list(
    produces_gene_program = produces,
    expected_result_level = "cell_level_exploratory",
    result_status = "target_defined",
    ineligible_reason = "unsupported_gene_program_role"
  ))
}

m3_gene_program_row_from_target <- function(target) {
  fields <- m3_gene_program_fields(target)
  list(
    comparison_id = target$comparison_id[[1]],
    source_question_id = target$source_question_id[[1]],
    layer_scope = target$layer_scope[[1]],
    analysis_mode = target$analysis_mode[[1]],
    gene_program_role = target$gene_program_role[[1]],
    produces_gene_program = fields$produces_gene_program,
    annotation_only = fields$annotation_only,
    qc_only = fields$qc_only,
    global_context_only = fields$global_context_only,
    nichenet_eligible = fields$nichenet_eligible,
    nichenet_usage = fields$nichenet_usage,
    enrichment_eligible = fields$enrichment_eligible,
    enrichment_usage = fields$enrichment_usage,
    preferred_for_downstream = fields$preferred_for_downstream,
    expected_result_level = fields$expected_result_level,
    formal_preferred = fields$formal_preferred,
    formal_status = fields$formal_status,
    result_status = fields$result_status,
    skip_reason = fields$skip_reason,
    eligible_reason = fields$eligible_reason,
    ineligible_reason = fields$ineligible_reason,
    notes = target$notes[[1]]
  )
}

m3_fanout_gene_program <- function(comparisons, annotation_marker_targets = NULL) {
  if (nrow(comparisons) == 0 && (is.null(annotation_marker_targets) || nrow(annotation_marker_targets) == 0)) {
    return(m3_empty_df(m3_gene_program_cols))
  }

  rows <- list()
  if (nrow(comparisons) > 0) {
    for (idx in seq_len(nrow(comparisons))) {
      row <- comparisons[idx, , drop = FALSE]
      rows[[length(rows) + 1L]] <- m3_gene_program_row_from_target(row)
    }
  }
  if (!is.null(annotation_marker_targets) && nrow(annotation_marker_targets) > 0) {
    for (idx in seq_len(nrow(annotation_marker_targets))) {
      target <- annotation_marker_targets[idx, , drop = FALSE]
      target$comparison_id <- target$target_id[[1]]
      target$produces_gene_program <- "yes"
      rows[[length(rows) + 1L]] <- m3_gene_program_row_from_target(target)
    }
  }
  m3_bind_rows(rows, m3_gene_program_cols)
}
