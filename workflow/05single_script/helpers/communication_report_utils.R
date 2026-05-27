if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
}

comm_report_scalar <- function(x, default = "") {
  if (exists("normalize_scalar_value", mode = "function")) {
    return(normalize_scalar_value(x, default))
  }
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(default)
  }
  trimws(as.character(x))
}

comm_report_truthy <- function(x) {
  tolower(comm_report_scalar(x, "false")) %in% c("1", "true", "yes", "y", "ok")
}

communication_report_tier_rank <- function(tier) {
  ranks <- c(blocked = 0L, candidate = 1L, exploratory = 2L, primary = 3L)
  tier <- tolower(as.character(tier))
  tier[is.na(tier) | !nzchar(tier)] <- "blocked"
  unname(ranks[tier] %||% rep(NA_integer_, length(tier)))
}

communication_report_tier_colors <- function() {
  c(primary = "#2E7D32", exploratory = "#F9A825", candidate = "#EF6C00", blocked = "#9E9E9E")
}

communication_report_required_columns <- function(df, cols, default = "") {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (col in cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- if (nrow(df) == 0L) character(0) else rep(default, nrow(df))
    }
  }
  df
}

communication_report_normalize_consensus <- function(df) {
  cols <- c(
    "lr_axis_id", "condition_value", "pair_id", "layer_id", "source", "target",
    "ligand", "receptor", "cellchat_hit", "liana_consensus_hit", "nichenet_hit",
    "commot_spatial_hit", "evidence_method_n", "evidence_tier", "can_be_primary",
    "evidence_reason", "nichenet_n_targets_in_receiver_de"
  )
  df <- communication_report_required_columns(df, cols)
  if (nrow(df) == 0) {
    return(df)
  }
  for (col in c("cellchat_hit", "liana_consensus_hit", "nichenet_hit", "commot_spatial_hit")) {
    df[[col]] <- vapply(df[[col]], comm_report_truthy, logical(1))
  }
  df$evidence_method_n <- suppressWarnings(as.integer(df$evidence_method_n))
  missing_n <- is.na(df$evidence_method_n)
  df$evidence_method_n[missing_n] <- as.integer(df$cellchat_hit[missing_n]) +
    as.integer(df$liana_consensus_hit[missing_n]) +
    as.integer(df$nichenet_hit[missing_n]) +
    as.integer(df$commot_spatial_hit[missing_n])
  df$evidence_tier <- tolower(as.character(df$evidence_tier))
  df$evidence_tier[is.na(df$evidence_tier) | !nzchar(df$evidence_tier)] <- "blocked"
  df$methods_hit_pattern <- communication_report_methods_pattern(df)
  df
}

communication_report_methods_pattern <- function(df) {
  if (nrow(df) == 0) {
    return(character(0))
  }
  vapply(seq_len(nrow(df)), function(i) {
    hits <- character(0)
    if (isTRUE(df$cellchat_hit[[i]])) hits <- c(hits, "cellchat")
    if (isTRUE(df$liana_consensus_hit[[i]])) hits <- c(hits, "liana")
    if (isTRUE(df$nichenet_hit[[i]])) hits <- c(hits, "nichenet")
    if (isTRUE(df$commot_spatial_hit[[i]])) hits <- c(hits, "commot")
    if (length(hits) == 0) "none" else paste(hits, collapse = "+")
  }, character(1))
}

communication_report_pair_requirements <- function(pairs) {
  cols <- c(
    "pair_id", "methods_required", "min_methods_agreed", "require_downstream_de",
    "require_spatial_support", "evidence_tier_required", "cellchat_min_samples_consistent"
  )
  pairs <- communication_report_required_columns(pairs, cols)
  if (nrow(pairs) == 0) {
    return(pairs[, cols, drop = FALSE])
  }
  pairs$methods_required[!nzchar(pairs$methods_required)] <- "cellchat,liana,nichenet"
  pairs$min_methods_agreed[!nzchar(pairs$min_methods_agreed)] <- "2"
  pairs$require_downstream_de[!nzchar(pairs$require_downstream_de)] <- "auto"
  pairs$require_spatial_support[!nzchar(pairs$require_spatial_support)] <- "auto"
  pairs$evidence_tier_required[!nzchar(pairs$evidence_tier_required)] <- "exploratory"
  pairs[, cols, drop = FALSE]
}

communication_report_filter_by_tier <- function(consensus_df, pairs, filter_tier = "auto") {
  consensus_df <- communication_report_normalize_consensus(consensus_df)
  if (nrow(consensus_df) == 0) {
    return(consensus_df)
  }
  filter_tier <- tolower(comm_report_scalar(filter_tier, "auto"))
  if (!identical(filter_tier, "auto")) {
    required <- rep(filter_tier, nrow(consensus_df))
  } else {
    pair_req <- communication_report_pair_requirements(pairs)
    required <- rep("exploratory", nrow(consensus_df))
    if (nrow(pair_req) > 0 && "pair_id" %in% colnames(consensus_df)) {
      idx <- match(consensus_df$pair_id, pair_req$pair_id)
      has <- !is.na(idx)
      required[has] <- pair_req$evidence_tier_required[idx[has]]
    }
  }
  required[!nzchar(required)] <- "exploratory"
  keep <- communication_report_tier_rank(consensus_df$evidence_tier) >= communication_report_tier_rank(required)
  consensus_df[keep %in% TRUE, , drop = FALSE]
}

communication_report_primary_axes <- function(df, top_n = 30L) {
  df <- communication_report_normalize_consensus(df)
  if (nrow(df) == 0) {
    return(df)
  }
  order_idx <- order(
    -communication_report_tier_rank(df$evidence_tier),
    -suppressWarnings(as.integer(df$evidence_method_n)),
    df$pair_id,
    df$lr_axis_id
  )
  df <- df[order_idx, , drop = FALSE]
  top_n <- suppressWarnings(as.integer(top_n))
  if (is.na(top_n) || top_n <= 0L) top_n <- 30L
  df[seq_len(min(nrow(df), top_n)), , drop = FALSE]
}

communication_report_md_table <- function(df, max_rows = 50L) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) == 0 || ncol(df) == 0) {
    return(c("_No rows._", ""))
  }
  df <- df[seq_len(min(nrow(df), max_rows)), , drop = FALSE]
  df[] <- lapply(df, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    gsub("\\|", "/", x)
  })
  header <- paste0("| ", paste(colnames(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  c(header, sep, rows, "")
}

communication_report_write_basic_html <- function(md_path, html_path) {
  if (!file.exists(md_path)) {
    return(invisible(FALSE))
  }
  lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  escape <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }
  html <- c(
    "<!doctype html>",
    "<html><head><meta charset=\"utf-8\"><title>Communication EDA</title>",
    "<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:32px auto;line-height:1.45} table{border-collapse:collapse} td,th{border:1px solid #ddd;padding:4px 8px} code,pre{background:#f6f8fa}</style>",
    "</head><body><pre>",
    escape(lines),
    "</pre></body></html>"
  )
  if (exists("ensure_dir", mode = "function")) ensure_dir(dirname(html_path)) else dir.create(dirname(html_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, html_path, useBytes = TRUE)
  invisible(TRUE)
}

communication_report_panels_to_tsv <- function(panels) {
  rows <- lapply(names(panels), function(name) {
    panel <- panels[[name]]
    data.frame(
      panel_id = name,
      panel_name = comm_report_scalar(panel$name %||% name, name),
      status = comm_report_scalar(panel$status %||% "ok", "ok"),
      n_rows = as.integer(nrow(panel$tsv %||% data.frame())),
      plot_paths = paste(unlist(panel$plots %||% list(), use.names = FALSE), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  if (length(rows) == 0) {
    return(data.frame(panel_id = character(), panel_name = character(), status = character(), n_rows = integer(), plot_paths = character()))
  }
  do.call(rbind, rows)
}

communication_report_assemble <- function(panels, consensus_df, filtered_df, pairs, cfg) {
  consensus_df <- communication_report_normalize_consensus(consensus_df)
  filtered_df <- communication_report_normalize_consensus(filtered_df)
  pair_req <- communication_report_pair_requirements(pairs)
  lines <- c(
    "# Communication EDA Report",
    "",
    "CellChat evidence is hypothesis-only in this pipeline. A CellChat-only axis remains candidate and cannot be promoted to primary without independent support.",
    "",
    sprintf("- Total axes: %d", nrow(consensus_df)),
    sprintf("- Filtered axes: %d", nrow(filtered_df)),
    sprintf("- Report filter: %s", comm_report_scalar(cfg$communication_report_filter_tier %||% "auto", "auto")),
    ""
  )
  for (name in names(panels)) {
    lines <- c(lines, panels[[name]]$md_lines %||% character(0), "")
  }
  lines <- c(
    lines,
    "## Pair Requirements",
    "",
    communication_report_md_table(pair_req, max_rows = 80L),
    "## Filtered Axis Preview",
    "",
    communication_report_md_table(
      communication_report_primary_axes(filtered_df, cfg$communication_report_top_n_primary %||% 30L)[
        , intersect(c("pair_id", "condition_value", "lr_axis_id", "methods_hit_pattern", "evidence_method_n", "evidence_tier", "can_be_primary"), colnames(filtered_df)), drop = FALSE
      ],
      max_rows = cfg$communication_report_top_n_primary %||% 30L
    )
  )
  lines
}

communication_report_json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", as.character(x))
  x <- gsub("\"", "\\\\\"", x)
  x <- gsub("\n", "\\\\n", x, fixed = TRUE)
  paste0("\"", x, "\"")
}

communication_report_json_value <- function(x, indent = 0L) {
  pad <- paste(rep(" ", indent), collapse = "")
  if (is.null(x)) return("null")
  if (is.list(x) && is.null(names(x))) {
    inner <- vapply(x, communication_report_json_value, character(1), indent = indent + 2L)
    return(paste0("[", paste(inner, collapse = ", "), "]"))
  }
  if (is.list(x)) {
    nms <- names(x)
    inner <- vapply(nms, function(nm) {
      paste0(pad, "  ", communication_report_json_escape(nm), ": ", communication_report_json_value(x[[nm]], indent + 2L))
    }, character(1))
    return(paste0("{\n", paste(inner, collapse = ",\n"), "\n", pad, "}"))
  }
  if (is.logical(x)) return(ifelse(isTRUE(x), "true", "false"))
  if (is.numeric(x) && length(x) == 1L && !is.na(x)) return(as.character(x))
  communication_report_json_escape(x)
}

communication_report_write_manifest <- function(manifest_path, outputs, module_name, base_dir, inputs, version, depends_on = list()) {
  manifest <- list(
    module = module_name,
    version = version,
    timestamp = if (exists("manifest_timestamp_now_local", mode = "function")) manifest_timestamp_now_local() else format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    base_dir = base_dir,
    inputs = inputs,
    outputs = outputs,
    depends_on = depends_on
  )
  if (exists("ensure_dir", mode = "function")) ensure_dir(dirname(manifest_path)) else dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
  } else {
    writeLines(communication_report_json_value(manifest), manifest_path, useBytes = TRUE)
  }
}
