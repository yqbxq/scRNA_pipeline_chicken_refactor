m3_trajectory_cols <- c(
  "trajectory_id", "source_question_id", "layer_scope", "root_group",
  "terminal_group", "condition_split_var", "condition_split_values",
  "method", "tools_to_run", "methods_extra", "outlier_qc_policy",
  "regress_cell_cycle", "coarse_label_var", "fine_label_var",
  "split_mode", "enabled", "notes"
)

m3_trajectory_default_tools <- function(method) {
  if (identical(method, "velocity")) {
    return("scvelo_dynamical,scvelo_stochastic,velocyto,cellrank")
  }
  "slingshot,monocle3,paga_dpt,tradeseq,palantir"
}

m3_trajectory_coarse_label_var <- function(scope) {
  if (identical(scope, "panorama")) {
    return("cell_type")
  }
  "cell_subtype"
}

m3_trajectory_notes <- function(q, notes, split) {
  row_notes <- m3_trim(notes)
  if (!nzchar(row_notes)) {
    row_notes <- m3_trim(q$notes)
  }
  if (!nzchar(row_notes) && nzchar(split$var)) {
    row_notes <- paste0("split=", split$var, ":", split$values)
  }
  if (!nzchar(row_notes)) "-" else row_notes
}

m3_trajectory_row <- function(q, suffix, root_group, terminal_group, method, notes = "") {
  split <- m3_condition_fields(q$condition_split)
  list(
    trajectory_id = paste(q$question_id, m3_safe_id(suffix), sep = "__"),
    source_question_id = q$question_id,
    layer_scope = q$scope,
    root_group = root_group,
    terminal_group = terminal_group,
    condition_split_var = split$var,
    condition_split_values = split$values,
    method = method,
    tools_to_run = m3_trajectory_default_tools(method),
    methods_extra = "",
    outlier_qc_policy = "standard",
    regress_cell_cycle = "auto",
    coarse_label_var = m3_trajectory_coarse_label_var(q$scope),
    fine_label_var = "seurat_clusters",
    split_mode = "auto",
    enabled = "yes",
    notes = m3_trajectory_notes(q, notes, split)
  )
}

m3_fanout_trajectory <- function(q) {
  rows <- list()
  if (identical(q$contrast_axis, "lineage")) {
    root <- m3_flat_groups(q$sender_groups, q$scope)
    terminal <- m3_flat_groups(q$receiver_groups, q$scope)
    root_group <- if (length(root) > 0) root[[1]] else "auto"
    terminal_group <- if (length(terminal) > 0) terminal[[1]] else "auto"
    rows[[1L]] <- m3_trajectory_row(
      q, paste(root_group, "to", terminal_group, sep = "_"),
      root_group,
      terminal_group,
      "trajectory"
    )
  } else if (identical(q$contrast_axis, "velocity")) {
    rows[[1L]] <- m3_trajectory_row(q, "velocity", "*", "*", "velocity")
  }
  m3_bind_rows(rows, m3_trajectory_cols)
}
