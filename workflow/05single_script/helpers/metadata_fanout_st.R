m3_deconv_cols <- c(
  "deconv_id", "source_question_id", "st_scope", "reference_scope",
  "section_filter", "condition_split_var", "condition_split_values",
  "tool", "enabled", "notes"
)

m3_spatial_cols <- c(
  "spatial_pair_id", "source_question_id", "st_scope", "sender", "receiver",
  "contrast_axis", "condition_split_var", "condition_split_values",
  "tool", "enabled", "notes"
)

m3_fanout_st <- function(q) {
  stop(
    sprintf(
      "ST fan-out is not implemented yet (question_id=%s, contrast_axis=%s). Keep ST_section rows planned until ST-H.",
      q$question_id,
      q$contrast_axis
    ),
    call. = FALSE
  )
}
