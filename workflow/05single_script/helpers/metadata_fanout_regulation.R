m3_scenic_cols <- c(
  "target_id", "source_question_id", "layer_scope", "cell_subset",
  "contrast_axis", "condition_split_var", "condition_split_values",
  "method", "enabled", "notes"
)

m3_regulation_method <- function(q) {
  mode <- m3_tool_mode(q$tools_to_run, c("scenic", "decoupler"), both_label = "both")
  if (nzchar(mode)) mode else "decoupler"
}

m3_scenic_row <- function(q, suffix, cell_subset, notes = "") {
  split <- m3_condition_fields(q$condition_split)
  list(
    target_id = paste(q$question_id, m3_safe_id(suffix), sep = "__"),
    source_question_id = q$question_id,
    layer_scope = q$scope,
    cell_subset = cell_subset,
    contrast_axis = q$contrast_axis,
    condition_split_var = split$var,
    condition_split_values = split$values,
    method = m3_regulation_method(q),
    enabled = "yes",
    notes = notes
  )
}

m3_fanout_regulation <- function(q) {
  axis <- q$contrast_axis
  rows <- list()

  if (axis == "regulation_per") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    if (length(groups) == 0) groups <- "all_cells"
    for (group in groups) {
      rows[[length(rows) + 1L]] <- m3_scenic_row(q, group, group)
    }
  } else if (axis == "regulation_pair") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    if (length(groups) >= 2) {
      combos <- utils::combn(groups, 2, simplify = FALSE)
      for (combo in combos) {
        rows[[length(rows) + 1L]] <- m3_scenic_row(
          q, paste(combo[[1]], "vs", combo[[2]], sep = "_"),
          m3_collapse(combo),
          notes = paste(combo[[1]], combo[[2]], sep = "_vs_")
        )
      }
    }
  } else if (axis == "regulation_stage") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    if (length(groups) == 0) groups <- "all_cells"
    for (group in groups) {
      rows[[length(rows) + 1L]] <- m3_scenic_row(q, paste(group, "stage", sep = "_"), group)
    }
  }

  m3_bind_rows(rows, m3_scenic_cols)
}
