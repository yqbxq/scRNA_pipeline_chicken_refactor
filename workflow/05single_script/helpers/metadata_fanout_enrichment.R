m3_enrichment_cols <- c(
  "target_id", "source_question_id", "comparison_id", "layer_scope",
  "analysis_mode", "gene_program_role", "organism", "database",
  "min_genes", "enabled", "notes"
)

m3_gene_program_cols <- c(
  "comparison_id", "source_question_id", "layer_scope", "analysis_mode",
  "gene_program_role", "preferred_for_downstream", "expected_result_level",
  "formal_preferred", "nichenet_eligible", "enrichment_eligible", "notes"
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
    produces <- if ("produces_gene_program" %in% colnames(row)) row$produces_gene_program[[1]] else "yes"
    if (!identical(tolower(trimws(produces)), "yes")) {
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
  if (identical(mode, "condition_within_type")) {
    return(list(
      preferred_for_downstream = "yes",
      expected_result_level = "pseudobulk_formal_or_exploratory_fallback",
      formal_preferred = "yes",
      nichenet_eligible = "yes",
      enrichment_eligible = "yes"
    ))
  }
  if (identical(role, "receiver_marker")) {
    return(list(
      preferred_for_downstream = "yes",
      expected_result_level = "cell_level_exploratory",
      formal_preferred = "no",
      nichenet_eligible = "yes",
      enrichment_eligible = "yes"
    ))
  }
  list(
    preferred_for_downstream = "contextual",
    expected_result_level = "cell_level_exploratory",
    formal_preferred = "no",
    nichenet_eligible = "no",
    enrichment_eligible = "yes"
  )
}

m3_fanout_gene_program <- function(comparisons) {
  if (nrow(comparisons) == 0) {
    return(m3_empty_df(m3_gene_program_cols))
  }
  cmp <- comparisons
  if ("produces_gene_program" %in% colnames(cmp)) {
    cmp <- cmp[tolower(trimws(cmp$produces_gene_program)) == "yes", , drop = FALSE]
  }
  if (nrow(cmp) == 0) {
    return(m3_empty_df(m3_gene_program_cols))
  }

  rows <- list()
  for (idx in seq_len(nrow(cmp))) {
    row <- cmp[idx, , drop = FALSE]
    fields <- m3_gene_program_fields(row)
    rows[[length(rows) + 1L]] <- list(
      comparison_id = row$comparison_id[[1]],
      source_question_id = row$source_question_id[[1]],
      layer_scope = row$layer_scope[[1]],
      analysis_mode = row$analysis_mode[[1]],
      gene_program_role = row$gene_program_role[[1]],
      preferred_for_downstream = fields$preferred_for_downstream,
      expected_result_level = fields$expected_result_level,
      formal_preferred = fields$formal_preferred,
      nichenet_eligible = fields$nichenet_eligible,
      enrichment_eligible = fields$enrichment_eligible,
      notes = row$notes[[1]]
    )
  }
  m3_bind_rows(rows, m3_gene_program_cols)
}
