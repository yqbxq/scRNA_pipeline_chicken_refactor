empty_communication_pairs_07 <- function() {
  empty_df_07(c(
    "pair_id", "source_question_id", "layer_scope", "sender", "receiver",
    "subset_column", "subset_value", "condition_split_var", "condition_split_values",
    "tool", "communication_mode", "receiver_gene_program_source",
    "baseline_marker_comparison_id", "receiver_deg_comparison_id",
    "direction_filter", "requires_cell_subtype", "enabled", "notes"
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
