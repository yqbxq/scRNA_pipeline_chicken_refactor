empty_communication_pairs_07 <- function() {
  empty_df_07(c(
    "pair_id", "source_question_id", "layer_scope", "sender", "receiver",
    "subset_column", "subset_value", "condition_split_var", "condition_split_values",
    "tool", "communication_mode", "receiver_gene_program_source",
    "activation_policy", "min_sender_cells", "min_receiver_cells",
    "min_cells_per_condition", "fallback_pair_id", "derived_from_pair_id",
    "run_baseline_if_split_fails", "requires_all_derived_inputs_pass",
    "baseline_marker_comparison_id", "receiver_deg_comparison_id",
    "direction_filter", "requires_cell_subtype", "notes", "enabled"
  ))
}

read_communication_pairs <- function(cfg) {
  df <- read_tsv_optional(cfg$communication_pairs_sheet)
  expected_cols <- colnames(empty_communication_pairs_07())
  if (nrow(df) == 0) {
    return(empty_communication_pairs_07())
  }
  for (col in expected_cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  df <- df[, expected_cols, drop = FALSE]
  for (col in expected_cols) {
    df[[col]] <- vapply(df[[col]], normalize_scalar_value, character(1))
  }
  df$tool[!nzchar(df$tool)] <- "both"
  df$tool <- tolower(df$tool)
  df$communication_mode[!nzchar(df$communication_mode)] <- "baseline"
  df$activation_policy[!nzchar(df$activation_policy)] <- ifelse(df$communication_mode == "differential_summary", "derived_from_split", "always")
  df$activation_policy <- tolower(df$activation_policy)
  df$run_baseline_if_split_fails[!nzchar(df$run_baseline_if_split_fails)] <- "no"
  df$run_baseline_if_split_fails <- tolower(df$run_baseline_if_split_fails)
  df$requires_all_derived_inputs_pass[!nzchar(df$requires_all_derived_inputs_pass)] <- "no"
  df$requires_all_derived_inputs_pass <- tolower(df$requires_all_derived_inputs_pass)
  df$receiver_gene_program_source[!nzchar(df$receiver_gene_program_source)] <- "none"
  df$receiver_gene_program_source <- tolower(df$receiver_gene_program_source)
  df$direction_filter[!nzchar(df$direction_filter)] <- "yes"
  df$direction_filter <- tolower(df$direction_filter)
  df$requires_cell_subtype[!nzchar(df$requires_cell_subtype)] <- "no"
  df$requires_cell_subtype <- tolower(df$requires_cell_subtype)
  df$enabled[!nzchar(df$enabled)] <- "yes"
  df$enabled <- tolower(df$enabled)
  df <- df[nzchar(df$pair_id), , drop = FALSE]
  df
}

communication_pair_enabled <- function(pair_row) {
  enabled <- tolower(normalize_scalar_value(pair_row$enabled[[1]], "yes"))
  enabled %in% c("yes", "true", "1", "on")
}

communication_pair_tool_enabled <- function(pair_row, tool) {
  requested <- tolower(normalize_scalar_value(pair_row$tool[[1]], "both"))
  if (identical(requested, "both")) {
    return(TRUE)
  }
  if (identical(tool, "cellchat")) {
    return(requested %in% c("cellchat", "cellchat_only"))
  }
  if (identical(tool, "nichenet")) {
    return(requested %in% c("nichenet", "nichenet_only"))
  }
  FALSE
}

communication_pair_direction_filter <- function(pair_row) {
  flag <- tolower(normalize_scalar_value(pair_row$direction_filter[[1]], "yes"))
  flag %in% c("yes", "true", "1", "on")
}

communication_pair_requires_cell_subtype <- function(pair_row) {
  flag <- tolower(normalize_scalar_value(pair_row$requires_cell_subtype[[1]], "no"))
  flag %in% c("yes", "true", "1", "on")
}

communication_pair_activation_policy <- function(pair_row) {
  tolower(normalize_scalar_value(pair_row$activation_policy[[1]], "always"))
}

communication_pair_auto_gate <- function(pair_row) {
  identical(communication_pair_activation_policy(pair_row), "auto_if_min_cells")
}

communication_pair_derived <- function(pair_row) {
  identical(communication_pair_activation_policy(pair_row), "derived_from_split")
}

communication_pair_run_baseline_if_split_fails <- function(pair_row) {
  flag <- tolower(normalize_scalar_value(pair_row$run_baseline_if_split_fails[[1]], "no"))
  flag %in% c("yes", "true", "1", "on")
}

communication_pair_requires_all_derived_inputs_pass <- function(pair_row) {
  flag <- tolower(normalize_scalar_value(pair_row$requires_all_derived_inputs_pass[[1]], "no"))
  flag %in% c("yes", "true", "1", "on")
}

communication_pair_min_int <- function(pair_row, col) {
  if (!col %in% colnames(pair_row)) {
    return(NA_integer_)
  }
  suppressWarnings(as.integer(normalize_scalar_value(pair_row[[col]][[1]])))
}

communication_gate_counts_07 <- function(obj, cell_type_col, roles) {
  meta <- obj@meta.data
  sender_n <- sum(as.character(meta[[cell_type_col]]) %in% roles$sender, na.rm = TRUE)
  receiver_n <- sum(as.character(meta[[cell_type_col]]) %in% roles$receiver, na.rm = TRUE)
  list(
    sender_n = as.integer(sender_n),
    receiver_n = as.integer(receiver_n),
    condition_pair_cell_n = as.integer(sender_n + receiver_n)
  )
}

evaluate_communication_min_cell_gate_07 <- function(pair_row, sender_n, receiver_n, condition_pair_cell_n) {
  policy <- communication_pair_activation_policy(pair_row)
  min_sender <- communication_pair_min_int(pair_row, "min_sender_cells")
  min_receiver <- communication_pair_min_int(pair_row, "min_receiver_cells")
  min_condition <- communication_pair_min_int(pair_row, "min_cells_per_condition")
  if (!identical(policy, "auto_if_min_cells")) {
    return(list(
      gate_status = ifelse(identical(policy, "derived_from_split"), "derived_not_executable", "not_applicable"),
      pass = TRUE,
      reason = "",
      min_sender_cells = min_sender,
      min_receiver_cells = min_receiver,
      min_cells_per_condition = min_condition
    ))
  }
  missing_threshold <- any(is.na(c(min_sender, min_receiver, min_condition))) ||
    any(c(min_sender, min_receiver, min_condition) <= 0)
  if (missing_threshold) {
    return(list(
      gate_status = "fail",
      pass = FALSE,
      reason = "auto_if_min_cells requires positive min_sender_cells, min_receiver_cells, and min_cells_per_condition",
      min_sender_cells = min_sender,
      min_receiver_cells = min_receiver,
      min_cells_per_condition = min_condition
    ))
  }
  failures <- character(0)
  if (sender_n < min_sender) failures <- c(failures, sprintf("sender_n=%s < min_sender_cells=%s", sender_n, min_sender))
  if (receiver_n < min_receiver) failures <- c(failures, sprintf("receiver_n=%s < min_receiver_cells=%s", receiver_n, min_receiver))
  if (condition_pair_cell_n < min_condition) failures <- c(failures, sprintf("condition_pair_cell_n=%s < min_cells_per_condition=%s", condition_pair_cell_n, min_condition))
  list(
    gate_status = ifelse(length(failures) == 0, "pass", "fail"),
    pass = length(failures) == 0,
    reason = paste(failures, collapse = "; "),
    min_sender_cells = min_sender,
    min_receiver_cells = min_receiver,
    min_cells_per_condition = min_condition
  )
}

pair_scope_applies_to_layer_07 <- function(pair_row, layer_id) {
  scope <- normalize_scalar_value(pair_row$layer_scope[[1]], "*")
  if (!nzchar(scope) || identical(scope, "*")) {
    return(TRUE)
  }
  layer_id %in% split_csv_local(scope)
}

communication_pairs_for_layer <- function(pairs, layer_id, tool) {
  if (nrow(pairs) == 0) {
    return(pairs)
  }
  keep <- vapply(seq_len(nrow(pairs)), function(i) {
    row <- pairs[i, , drop = FALSE]
    communication_pair_enabled(row) &&
      communication_pair_tool_enabled(row, tool) &&
      pair_scope_applies_to_layer_07(row, layer_id)
  }, logical(1))
  pairs[keep, , drop = FALSE]
}

resolve_role_pattern_07 <- function(pattern, all_cell_types) {
  pattern <- normalize_scalar_value(pattern, "*")
  if (!nzchar(pattern) || identical(pattern, "*")) {
    return(all_cell_types)
  }
  tokens <- split_csv_local(pattern)
  if (length(tokens) == 0) {
    tokens <- pattern
  }
  out <- character(0)
  for (token in tokens) {
    token <- normalize_scalar_value(token)
    if (!nzchar(token) || identical(token, "*")) {
      out <- c(out, all_cell_types)
    } else if (endsWith(token, "_*")) {
      prefix <- sub("_\\*$", "", token)
      out <- c(out, all_cell_types[startsWith(all_cell_types, prefix)])
    } else if (grepl("\\*$", token)) {
      prefix <- sub("\\*$", "", token)
      out <- c(out, all_cell_types[startsWith(all_cell_types, prefix)])
    } else {
      out <- c(out, intersect(token, all_cell_types))
    }
  }
  sort(unique(out[nzchar(out)]))
}

resolve_sender_receiver_sets <- function(pair_row, seu, cell_type_col) {
  all_cell_types <- sort(unique(as.character(seu@meta.data[[cell_type_col]])))
  all_cell_types <- all_cell_types[nzchar(all_cell_types)]
  sender_spec <- normalize_scalar_value(pair_row$sender[[1]], "*")
  receiver_spec <- normalize_scalar_value(pair_row$receiver[[1]], "*")
  sender <- resolve_role_pattern_07(sender_spec, all_cell_types)
  receiver <- resolve_role_pattern_07(receiver_spec, all_cell_types)
  requested <- unique(c(split_csv_local(sender_spec), split_csv_local(receiver_spec)))
  requested <- requested[nzchar(requested) & !requested %in% c("*")]
  concrete_requested <- requested[!grepl("\\*$", requested)]
  missing <- setdiff(concrete_requested, all_cell_types)
  list(
    sender = sender,
    receiver = receiver,
    sender_label = paste(sender, collapse = ","),
    receiver_label = paste(receiver, collapse = ","),
    missing = missing
  )
}

subset_by_pair_filter_07 <- function(seu, pair_row) {
  subset_column <- normalize_scalar_value(pair_row$subset_column[[1]])
  subset_value <- normalize_scalar_value(pair_row$subset_value[[1]])
  if (!nzchar(subset_column) || subset_column %in% c("-", "NA", "na")) {
    return(list(object = seu, status = "ok", reason = "", subset_n = ncol(seu)))
  }
  if (!subset_column %in% colnames(seu@meta.data)) {
    return(list(object = NULL, status = "subset_column_missing", reason = sprintf("subset_column missing: %s", subset_column), subset_n = 0L))
  }
  values <- split_csv_local(subset_value)
  if (length(values) == 0 || subset_value %in% c("-", "NA", "na")) {
    return(list(object = seu, status = "ok", reason = "", subset_n = ncol(seu)))
  }
  keep <- as.character(seu@meta.data[[subset_column]]) %in% values
  if (!any(keep)) {
    return(list(object = NULL, status = "empty_subset", reason = sprintf("0 cells match %s in %s", subset_value, subset_column), subset_n = 0L))
  }
  obj <- subset(seu, cells = rownames(seu@meta.data)[keep])
  if (exists("maybe_join_layers", mode = "function")) {
    obj <- maybe_join_layers(obj)
  }
  list(object = obj, status = "ok", reason = "", subset_n = ncol(obj))
}

split_by_condition <- function(seu, pair_row) {
  split_var <- normalize_scalar_value(pair_row$condition_split_var[[1]])
  split_values <- split_csv_local(pair_row$condition_split_values[[1]])
  if (!nzchar(split_var) || split_var %in% c("-", "NA", "na")) {
    return(list(all = list(object = seu, condition_value = "all", status = "ok", reason = "")))
  }
  if (!split_var %in% colnames(seu@meta.data)) {
    return(list(missing = list(object = NULL, condition_value = "", status = "condition_split_var_missing", reason = sprintf("condition_split_var missing: %s", split_var))))
  }
  if (length(split_values) == 0) {
    split_values <- sort(unique(as.character(seu@meta.data[[split_var]])))
    split_values <- split_values[nzchar(split_values)]
  }
  out <- list()
  for (value in split_values) {
    keep <- as.character(seu@meta.data[[split_var]]) == value
    if (!any(keep)) {
      out[[value]] <- list(object = NULL, condition_value = value, status = "empty_condition", reason = sprintf("0 cells for %s=%s", split_var, value))
      next
    }
    obj <- subset(seu, cells = rownames(seu@meta.data)[keep])
    if (exists("maybe_join_layers", mode = "function")) {
      obj <- maybe_join_layers(obj)
    }
    out[[value]] <- list(object = obj, condition_value = value, status = "ok", reason = "")
  }
  out
}

pair_condition_task_id_07 <- function(pair_id, condition_value) {
  paste(safe_id_07(pair_id), safe_id_07(condition_value), sep = "_")
}
