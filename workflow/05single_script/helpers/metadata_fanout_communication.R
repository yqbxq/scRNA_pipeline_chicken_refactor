m3_communication_cols <- c(
  "pair_id", "source_question_id", "layer_scope", "sender", "receiver",
  "subset_column", "subset_value", "condition_split_var",
  "condition_split_values", "tool", "enabled", "notes"
)

m3_communication_row <- function(q, suffix, sender, receiver, tool, notes = "", pair_id = "") {
  subset <- m3_condition_subset(q$condition_split)
  split <- m3_condition_fields(q$condition_split)
  list(
    pair_id = if (nzchar(pair_id)) pair_id else paste(q$question_id, m3_safe_id(suffix), sep = "__"),
    source_question_id = q$question_id,
    layer_scope = q$scope,
    sender = sender,
    receiver = receiver,
    subset_column = subset$subset_column,
    subset_value = subset$subset_value,
    condition_split_var = split$var,
    condition_split_values = split$values,
    tool = tool,
    enabled = "yes",
    notes = notes
  )
}

m3_comm_tool <- function(q) {
  mode <- m3_tool_mode(q$tools_to_run, c("cellchat", "nichenet"), both_label = "both")
  if (nzchar(mode)) mode else "cellchat"
}

m3_add_differential_comm <- function(rows, q) {
  split <- m3_condition_fields(q$condition_split)
  if (!nzchar(split$var)) {
    return(rows)
  }
  sender <- m3_collapse(m3_flat_groups(q$sender_groups, q$scope))
  receiver <- m3_collapse(m3_flat_groups(q$receiver_groups, q$scope))
  derived_id <- switch(
    q$question_id,
    F02_TC_GC_bidir_split = "F23_TC_GC_diff_overall",
    F04_TC_to_GC_each_split = "F24_TC_GC_diff_subtype_TC_to_GC",
    F06_GC_each_to_TC_split = "F24_TC_GC_diff_subtype_GC_to_TC",
    F08_GC_dev_seq_split = "F25_GC_internal_diff",
    F15_panorama_screen_split = "F26_panorama_screen_diff",
    ""
  )
  rows[[length(rows) + 1L]] <- m3_communication_row(
    q, "diff_summary", sender, receiver, "differential",
    notes = "derived_differential_communication",
    pair_id = derived_id
  )
  rows
}

m3_fanout_communication <- function(q) {
  axis <- q$contrast_axis
  tool <- m3_comm_tool(q)
  rows <- list()

  if (axis == "bidirectional") {
    senders <- m3_flat_groups(q$sender_groups, q$scope)
    receivers <- m3_flat_groups(q$receiver_groups, q$scope)
    for (sender in senders) {
      for (receiver in receivers) {
        rows[[length(rows) + 1L]] <- m3_communication_row(q, paste(sender, "to", receiver, sep = "_"), sender, receiver, tool)
        rows[[length(rows) + 1L]] <- m3_communication_row(q, paste(receiver, "to", sender, sep = "_"), receiver, sender, tool)
      }
    }
  } else if (axis == "directional") {
    sender_sets <- m3_group_sets(q$sender_groups, q$scope)
    receiver_sets <- m3_group_sets(q$receiver_groups, q$scope)
    for (sender_set in sender_sets) {
      for (receiver_set in receiver_sets) {
        sender <- m3_collapse(sender_set)
        receiver <- m3_collapse(receiver_set)
        rows[[length(rows) + 1L]] <- m3_communication_row(q, paste(sender, "to", receiver, sep = "_"), sender, receiver, tool)
      }
    }
  } else if (axis == "sequential") {
    seq_groups <- m3_group_sequence(q$sender_groups, q$scope)
    if (length(seq_groups) >= 2) {
      for (idx in seq_len(length(seq_groups) - 1L)) {
        sender <- seq_groups[[idx]]
        receiver <- seq_groups[[idx + 1L]]
        rows[[length(rows) + 1L]] <- m3_communication_row(q, paste(sender, "to", receiver, sep = "_"), sender, receiver, tool)
      }
    }
  } else if (axis == "symmetric") {
    rows[[length(rows) + 1L]] <- m3_communication_row(q, "all_by_all", "*", "*", tool)
  } else if (axis == "pairwise_comm") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    for (sender in groups) {
      for (receiver in groups) {
        if (!identical(sender, receiver)) {
          rows[[length(rows) + 1L]] <- m3_communication_row(q, paste(sender, "to", receiver, sep = "_"), sender, receiver, tool)
        }
      }
    }
  }

  rows <- m3_add_differential_comm(rows, q)
  m3_bind_rows(rows, m3_communication_cols)
}
