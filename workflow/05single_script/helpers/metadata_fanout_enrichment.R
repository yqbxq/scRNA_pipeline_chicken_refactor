m3_enrichment_cols <- c(
  "target_id", "source_question_id", "comparison_id", "layer_scope",
  "organism", "database", "min_genes", "enabled", "notes"
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
    for (database in c("GO", "KEGG")) {
      rows[[length(rows) + 1L]] <- list(
        target_id = paste(row$comparison_id[[1]], tolower(database), sep = "__"),
        source_question_id = row$source_question_id[[1]],
        comparison_id = row$comparison_id[[1]],
        layer_scope = row$layer_scope[[1]],
        organism = "chicken_primary",
        database = database,
        min_genes = "5",
        enabled = "yes",
        notes = "derived_from_deg_question"
      )
    }
  }
  m3_bind_rows(rows, m3_enrichment_cols)
}
