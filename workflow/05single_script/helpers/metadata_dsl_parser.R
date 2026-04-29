m3_null_or <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

`%||%` <- m3_null_or

m3_trim <- function(x, default = "") {
  x <- as.character(x %||% default)
  x[is.na(x)] <- default
  trimws(x)
}

m3_empty_df <- function(cols) {
  out <- as.data.frame(setNames(rep(list(character(0)), length(cols)), cols), stringsAsFactors = FALSE)
  out[, cols, drop = FALSE]
}

m3_safe_id <- function(x) {
  x <- m3_trim(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
  x <- gsub("^_+|_+$", "", x, perl = TRUE)
  if (!nzchar(x)) "na" else x
}

m3_collapse <- function(x) {
  x <- unique(m3_trim(x))
  x <- x[nzchar(x)]
  paste(x, collapse = ",")
}

m3_read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- sub("^﻿", "", lines)
  keep <- nzchar(trimws(lines)) & !startsWith(trimws(lines), "#")
  if (!any(keep)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  con <- textConnection(lines[keep])
  on.exit(close(con), add = TRUE)
  read.delim(con, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "", fill = TRUE)
}

m3_write_generated_tsv <- function(df, path, cols) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  for (col in cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- character(nrow(df))
    }
  }
  df <- df[, cols, drop = FALSE]
  if ("notes" %in% colnames(df) && nrow(df) > 0) {
    df$notes[!nzchar(df$notes)] <- "-"
  }
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines("# AUTO-GENERATED. Edit metadata/analysis_questions.tsv instead.", con, useBytes = TRUE)
  utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
}

m3_normalize_expr <- function(x) {
  x <- m3_trim(x)
  x <- gsub("→", "->", x, fixed = TRUE)
  gsub("\\s+", "", x, perl = TRUE)
}

m3_scope_groups <- function(scope) {
  switch(
    scope,
    panorama = c("GC", "TC"),
    GC_subcluster = c("pGC", "eGC", "rgGC", "lGC"),
    TC_subcluster = c("TC_*"),
    ST_section = c("regions"),
    character(0)
  )
}

m3_expand_special_group <- function(token, scope) {
  token <- m3_trim(token)
  if (!nzchar(token) || identical(token, "-")) {
    return(character(0))
  }
  if (identical(token, "*")) {
    return(m3_scope_groups(scope))
  }
  if (identical(token, "GC_subtypes")) {
    return(c("pGC", "eGC", "rgGC", "lGC"))
  }
  if (identical(token, "TC_subtypes")) {
    return("TC_*")
  }
  if (identical(token, "regions")) {
    return("regions")
  }
  token
}

m3_parse_group_set <- function(expr, scope) {
  expr <- m3_normalize_expr(expr)
  if (!nzchar(expr) || identical(expr, "-")) {
    return(character(0))
  }
  if (startsWith(expr, "[") && endsWith(expr, "]")) {
    expr <- substring(expr, 2, nchar(expr) - 1)
  }
  tokens <- unlist(strsplit(expr, ",", fixed = TRUE), use.names = FALSE)
  tokens <- unlist(lapply(tokens, m3_expand_special_group, scope = scope), use.names = FALSE)
  unique(tokens[nzchar(tokens)])
}

m3_parse_group_expr <- function(expr, scope) {
  expr <- m3_normalize_expr(expr)
  if (!nzchar(expr) || identical(expr, "-")) {
    return(list(kind = "empty", sets = list(), sequence = character(0), raw = expr))
  }
  if (grepl("->", expr, fixed = TRUE)) {
    parts <- unlist(strsplit(expr, "->", fixed = TRUE), use.names = FALSE)
    seq_tokens <- unlist(lapply(parts, function(x) m3_parse_group_set(x, scope)), use.names = FALSE)
    return(list(kind = "sequence", sets = as.list(seq_tokens), sequence = seq_tokens, raw = expr))
  }
  if (grepl("|", expr, fixed = TRUE)) {
    parts <- unlist(strsplit(expr, "|", fixed = TRUE), use.names = FALSE)
    sets <- lapply(parts, m3_parse_group_set, scope = scope)
    return(list(kind = "fanout", sets = sets, sequence = character(0), raw = expr))
  }
  list(kind = "set", sets = list(m3_parse_group_set(expr, scope)), sequence = character(0), raw = expr)
}

m3_flat_groups <- function(expr, scope) {
  parsed <- m3_parse_group_expr(expr, scope)
  unique(unlist(parsed$sets, use.names = FALSE))
}

m3_group_sets <- function(expr, scope) {
  m3_parse_group_expr(expr, scope)$sets
}

m3_group_sequence <- function(expr, scope) {
  parsed <- m3_parse_group_expr(expr, scope)
  if (length(parsed$sequence) > 0) parsed$sequence else m3_flat_groups(expr, scope)
}

m3_parse_condition_split <- function(x) {
  x <- m3_trim(x)
  if (!nzchar(x) || identical(x, "-")) {
    return(list(type = "none", column = "", values = character(0)))
  }
  if (x %in% c("syf_only", "f5_only")) {
    return(list(type = "subset", column = "group", values = sub("_only$", "", x)))
  }
  parts <- strsplit(x, ":", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    return(list(type = "split", column = parts[[1]], values = strsplit(parts[[2]], ",", fixed = TRUE)[[1]]))
  }
  list(type = "invalid", column = "", values = character(0))
}

m3_condition_subset <- function(condition_split) {
  split <- m3_parse_condition_split(condition_split)
  if (identical(split$type, "subset")) {
    return(list(subset_column = split$column, subset_value = m3_collapse(split$values)))
  }
  list(subset_column = "", subset_value = "")
}

m3_condition_pair <- function(condition_split) {
  split <- m3_parse_condition_split(condition_split)
  if (identical(split$type, "split") && length(split$values) >= 2) {
    return(list(group_var = split$column, ident_1 = split$values[[1]], ident_2 = split$values[[2]]))
  }
  list(group_var = "group", ident_1 = "syf", ident_2 = "f5")
}

m3_condition_fields <- function(condition_split) {
  split <- m3_parse_condition_split(condition_split)
  if (identical(split$type, "split")) {
    return(list(var = split$column, values = m3_collapse(split$values)))
  }
  list(var = "", values = "")
}

m3_tool_mode <- function(tools_to_run, candidates, both_label = "both") {
  tools <- tolower(unlist(strsplit(m3_trim(tools_to_run), "+", fixed = TRUE), use.names = FALSE))
  tools <- intersect(tools, candidates)
  if (length(tools) == 0) {
    return("")
  }
  if (length(tools) > 1) {
    return(both_label)
  }
  tools[[1]]
}

m3_bind_rows <- function(rows, cols) {
  if (length(rows) == 0) {
    return(m3_empty_df(cols))
  }
  df <- do.call(rbind, lapply(rows, function(row) {
    row <- as.list(row)
    for (col in cols) {
      if (is.null(row[[col]])) row[[col]] <- ""
    }
    as.data.frame(row[cols], stringsAsFactors = FALSE)
  }))
  rownames(df) <- NULL
  df
}

m3_assert_unique <- function(df, id_col, table_name) {
  if (nrow(df) == 0) {
    return(invisible(TRUE))
  }
  duplicated_ids <- unique(df[[id_col]][duplicated(df[[id_col]])])
  duplicated_ids <- duplicated_ids[nzchar(duplicated_ids)]
  if (length(duplicated_ids) > 0) {
    stop(sprintf("%s has duplicated %s values: %s", table_name, id_col, paste(duplicated_ids, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}
