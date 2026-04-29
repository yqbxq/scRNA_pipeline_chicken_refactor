m3_comparison_cols <- c(
  "comparison_id", "source_question_id", "layer_scope", "contrast_axis",
  "ident_1", "ident_2", "subset_column", "subset_value", "group_var",
  "batch_var", "enabled", "min_biological_replicates", "force_exploratory",
  "min_cells_per_group", "logfc_threshold", "notes"
)

m3_comparison_row <- function(q, suffix, ident_1, ident_2, group_var = "cell_subtype",
                              subset_column = "", subset_value = "", notes = "") {
  list(
    comparison_id = paste(q$question_id, m3_safe_id(suffix), sep = "__"),
    source_question_id = q$question_id,
    layer_scope = q$scope,
    contrast_axis = q$contrast_axis,
    ident_1 = ident_1,
    ident_2 = ident_2,
    subset_column = subset_column,
    subset_value = subset_value,
    group_var = group_var,
    batch_var = "batch",
    enabled = "yes",
    min_biological_replicates = "2",
    force_exploratory = "no",
    min_cells_per_group = "3",
    logfc_threshold = "0",
    notes = notes
  )
}

m3_fanout_comparison <- function(q) {
  axis <- q$contrast_axis
  rows <- list()
  condition_subset <- m3_condition_subset(q$condition_split)

  if (axis == "cluster_marker") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    for (group in groups) {
      rows[[length(rows) + 1L]] <- m3_comparison_row(
        q, paste0(group, "_vs_rest"), group, "__rest__",
        group_var = "cell_subtype",
        notes = "one_vs_rest"
      )
    }
  } else if (axis == "directional_DEG") {
    senders <- m3_flat_groups(q$sender_groups, q$scope)
    receivers <- m3_flat_groups(q$receiver_groups, q$scope)
    for (sender in senders) {
      for (receiver in receivers) {
        rows[[length(rows) + 1L]] <- m3_comparison_row(
          q, paste(sender, "vs", receiver, condition_subset$subset_value, sep = "_"),
          sender, receiver,
          subset_column = condition_subset$subset_column,
          subset_value = condition_subset$subset_value
        )
      }
    }
  } else if (axis == "pairwise") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    if (length(groups) >= 2) {
      combos <- utils::combn(groups, 2, simplify = FALSE)
      for (combo in combos) {
        rows[[length(rows) + 1L]] <- m3_comparison_row(
          q, paste(combo[[1]], "vs", combo[[2]], condition_subset$subset_value, sep = "_"),
          combo[[1]], combo[[2]],
          subset_column = condition_subset$subset_column,
          subset_value = condition_subset$subset_value
        )
      }
    }
  } else if (axis == "contrast_only") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    condition <- m3_condition_pair(q$condition_split)
    if (length(groups) == 0) {
      groups <- "all_cells"
    }
    for (group in groups) {
      subset_column <- if (group %in% c("all_cells")) "" else "cell_subtype"
      subset_value <- if (group %in% c("all_cells")) "" else group
      rows[[length(rows) + 1L]] <- m3_comparison_row(
        q, paste(group, condition$ident_1, "vs", condition$ident_2, sep = "_"),
        condition$ident_1, condition$ident_2,
        group_var = condition$group_var,
        subset_column = subset_column,
        subset_value = subset_value
      )
    }
  } else if (axis == "composition") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    condition <- m3_condition_pair(q$condition_split)
    subset_column <- if (length(groups) == 0 || any(groups %in% c("all_cells"))) "" else "cell_subtype"
    subset_value <- if (nzchar(subset_column)) m3_collapse(groups) else ""
    suffix <- if (nzchar(subset_value)) paste("composition", m3_safe_id(subset_value), sep = "_") else "composition_all"
    rows[[length(rows) + 1L]] <- m3_comparison_row(
      q, suffix,
      condition$ident_1, condition$ident_2,
      group_var = condition$group_var,
      subset_column = subset_column,
      subset_value = subset_value,
      notes = "composition_only"
    )
  }

  m3_bind_rows(rows, m3_comparison_cols)
}
