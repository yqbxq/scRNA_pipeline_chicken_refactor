#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

env_or <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) default else value
}

project_root <- env_or("PROJECT_ROOT", getwd())
metadata_dir <- env_or("METADATA_DIR", file.path(project_root, "metadata"))
results_dir <- env_or("RESULTS_DIR", file.path(project_root, "results"))
table_dir <- env_or("TABLE_DIR", file.path(results_dir, "tables"))
question_path <- env_or("ANALYSIS_QUESTIONS_FILE", file.path(metadata_dir, "analysis_questions.tsv"))
report_path <- env_or("METADATA_VALIDATION_REPORT", file.path(results_dir, "00_validation", "metadata_validation.json"))
layer_status_file <- env_or("LAYER_STATUS_FILE", file.path(table_dir, "layer_status.tsv"))
panorama_layer_id <- env_or("PANORAMA_LAYER_ID", "panorama")

records <- list()

json_escape <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("\\\\", "\\\\\\\\", x, perl = TRUE)
  x <- gsub("\"", "\\\\\"", x, perl = TRUE)
  x <- gsub("\n", "\\\\n", x, perl = TRUE)
  x <- gsub("\r", "\\\\r", x, perl = TRUE)
  x <- gsub("\t", "\\\\t", x, perl = TRUE)
  x
}

json_value <- function(x) {
  if (is.logical(x)) {
    return(ifelse(isTRUE(x), "true", "false"))
  }
  if (is.numeric(x) && length(x) == 1 && !is.na(x)) {
    return(as.character(x))
  }
  sprintf("\"%s\"", json_escape(x))
}

json_object <- function(x) {
  fields <- names(x)
  paste0(
    "{",
    paste(sprintf("\"%s\":%s", fields, vapply(x, json_value, character(1))), collapse = ","),
    "}"
  )
}

add_record <- function(level, check_id, status, message, question_id = "", fix = "") {
  records[[length(records) + 1L]] <<- list(
    level = level,
    check_id = check_id,
    status = status,
    question_id = question_id,
    message = message,
    fix = fix
  )
}

pass <- function(check_id, message, question_id = "") add_record("PASS", check_id, "pass", message, question_id)
warn <- function(check_id, message, question_id = "", fix = "") add_record("WARN", check_id, "warn", message, question_id, fix)
fail <- function(check_id, message, question_id = "", fix = "") add_record("FAIL", check_id, "fail", message, question_id, fix)
skip <- function(check_id, message, question_id = "") add_record("SKIP", check_id, "skip", message, question_id)

read_tsv <- function(path) {
  if (!nzchar(path) || !file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- sub("^\ufeff", "", lines)
  keep <- nzchar(trimws(lines)) & !startsWith(trimws(lines), "#")
  if (!any(keep)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  con <- textConnection(lines[keep])
  on.exit(close(con), add = TRUE)
  read.delim(con, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "", fill = TRUE)
}

trim <- function(x) trimws(as.character(x %||% ""))

split_list <- function(x, sep = ";") {
  x <- trim(x)
  if (!nzchar(x)) {
    return(character(0))
  }
  out <- trimws(unlist(strsplit(x, sep, fixed = TRUE), use.names = FALSE))
  out[nzchar(out)]
}

split_dependency_ids <- function(x) {
  x <- trim(x)
  if (!nzchar(x)) {
    return(character(0))
  }
  out <- trimws(unlist(strsplit(x, "[;,]", perl = TRUE), use.names = FALSE))
  out[nzchar(out)]
}

required_cols <- c(
  "question_id", "question_zh", "scope", "sender_groups", "receiver_groups",
  "condition_split", "contrast_axis", "tools_to_run", "priority", "status",
  "depends_on", "notes"
)

allowed_scopes <- c("panorama", "GC_subcluster", "TC_subcluster", "ST_section")
allowed_priorities <- c("P0", "P1", "P2", "P3")
allowed_status <- c("active", "planned")
axis_targets <- c(
  cluster_marker = "comparisons.tsv",
  directional_DEG = "comparisons.tsv",
  pairwise = "comparisons.tsv",
  contrast_only = "comparisons.tsv",
  composition = "comparisons.tsv",
  bidirectional = "communication_pairs.tsv",
  directional = "communication_pairs.tsv",
  sequential = "communication_pairs.tsv",
  symmetric = "communication_pairs.tsv",
  pairwise_comm = "communication_pairs.tsv",
  regulation_per = "scenic_targets.tsv",
  regulation_pair = "scenic_targets.tsv",
  regulation_stage = "scenic_targets.tsv",
  lineage = "trajectory_pairs.tsv",
  velocity = "trajectory_pairs.tsv",
  spatial_clustering = "spatial_pairs.tsv",
  spatial_neighbor = "spatial_pairs.tsv",
  spatial_overlay = "spatial_pairs.tsv",
  deconv_pair = "deconv_pairs.tsv",
  SVG = "spatial_pairs.tsv",
  SVG_diff = "spatial_pairs.tsv",
  deconv_validation = "deconv_pairs.tsv",
  enrichment_target = "enrichment_targets.tsv"
)

axis_tool_allow <- list(
  cluster_marker = c("deg", "enrichment"),
  directional_DEG = c("deg", "enrichment"),
  pairwise = c("deg", "enrichment"),
  contrast_only = c("deg", "enrichment"),
  composition = c("composition_test"),
  bidirectional = c("cellchat", "nichenet"),
  directional = c("cellchat", "nichenet"),
  sequential = c("cellchat", "nichenet"),
  symmetric = c("cellchat", "nichenet"),
  pairwise_comm = c("cellchat", "nichenet"),
  regulation_per = c("scenic", "decoupler"),
  regulation_pair = c("scenic", "decoupler"),
  regulation_stage = c("scenic", "decoupler"),
  lineage = c("trajectory", "velocity"),
  velocity = c("velocity"),
  enrichment_target = c("deg", "enrichment")
)

reserved_group_tokens <- c(
  "-", "*", "", "all_cells", "GC_subtypes", "TC_subtypes", "regions", "all_spots",
  "auto", "region", "synthetic", "panorama_ref", "GC_sub_ref", "TC_sub_ref",
  "H01_results", "F02_results", "I13_results", "I19_results", "I19+I20"
)

normalize_dsl <- function(x) {
  x <- trim(x)
  x <- gsub("→", "->", x, fixed = TRUE)
  gsub("\\s+", "", x, perl = TRUE)
}

valid_dsl <- function(x) {
  x <- normalize_dsl(x)
  if (!nzchar(x) || x %in% c("-", "*")) {
    return(TRUE)
  }
  if (x %in% c("(sequential)", "(fromI05+I06)")) {
    return(TRUE)
  }
  left_n <- gregexpr("\\[", x, perl = TRUE)[[1]]
  right_n <- gregexpr("\\]", x, perl = TRUE)[[1]]
  left_count <- if (identical(left_n, -1L)) 0L else length(left_n)
  right_count <- if (identical(right_n, -1L)) 0L else length(right_n)
  if (left_count != right_count) {
    return(FALSE)
  }
  if (grepl("-->", x, fixed = TRUE)) {
    return(FALSE)
  }
  if (grepl(">", gsub("->", "", x, fixed = TRUE), fixed = TRUE)) {
    return(FALSE)
  }
  grepl("^[A-Za-z0-9_*,+\\[\\]|>.-]+$", x, perl = TRUE)
}

extract_group_tokens <- function(x) {
  x <- normalize_dsl(x)
  if (!nzchar(x) || x %in% c("-", "*")) {
    return(character(0))
  }
  if (x %in% c("(sequential)", "(fromI05+I06)")) {
    return(character(0))
  }
  x <- gsub("\\[|\\]", "", x, perl = TRUE)
  pieces <- trimws(unlist(strsplit(x, "\\||->|,", perl = TRUE), use.names = FALSE))
  pieces <- pieces[nzchar(pieces)]
  pieces <- pieces[!pieces %in% reserved_group_tokens]
  pieces <- pieces[!grepl("[*]$", pieces)]
  unique(pieces)
}

parse_condition_split <- function(x) {
  x <- trim(x)
  if (!nzchar(x) || x == "-") {
    return(list(type = "none", column = "", values = character(0), ok = TRUE))
  }
  if (x %in% c("syf_only", "f5_only")) {
    return(list(type = "shortcut", column = "group_id", values = sub("_only$", "", x), ok = TRUE))
  }
  if (!grepl("^[A-Za-z0-9_.:-]+:[A-Za-z0-9_.-]+(,[A-Za-z0-9_.-]+)*$", x, perl = TRUE)) {
    return(list(type = "invalid", column = "", values = character(0), ok = FALSE))
  }
  parts <- strsplit(x, ":", fixed = TRUE)[[1]]
  list(type = "column", column = parts[[1]], values = strsplit(parts[[2]], ",", fixed = TRUE)[[1]], ok = TRUE)
}

tools_for_axis_ok <- function(axis, tools_to_run) {
  allow <- axis_tool_allow[[axis]]
  if (is.null(allow)) {
    return(TRUE)
  }
  tools <- tolower(split_list(tools_to_run, sep = "+"))
  any(tools %in% allow)
}

find_dependency_cycle <- function(df) {
  ids <- df$question_id
  resolve_dep <- function(dep) {
    dep <- trim(dep)
    if (!nzchar(dep)) return("")
    if (dep %in% ids) return(dep)
    hit <- ids[startsWith(ids, paste0(dep, "_"))]
    if (length(hit) == 1L) hit[[1]] else ""
  }
  graph <- setNames(vector("list", length(ids)), ids)
  for (idx in seq_len(nrow(df))) {
    deps <- vapply(split_dependency_ids(df$depends_on[[idx]]), resolve_dep, character(1))
    graph[[df$question_id[[idx]]]] <- deps[nzchar(deps)]
  }

  state <- setNames(rep(0L, length(ids)), ids)
  stack <- character(0)
  cycle <- character(0)

  visit <- function(node) {
    if (length(cycle) > 0) {
      return(invisible(FALSE))
    }
    state[[node]] <<- 1L
    stack <<- c(stack, node)
    for (dep in graph[[node]]) {
      if (state[[dep]] == 0L) {
        visit(dep)
      } else if (state[[dep]] == 1L) {
        start <- match(dep, stack)
        cycle <<- c(stack[start:length(stack)], dep)
        return(invisible(FALSE))
      }
    }
    stack <<- head(stack, -1L)
    state[[node]] <<- 2L
    invisible(TRUE)
  }

  for (id in ids) {
    if (state[[id]] == 0L) {
      visit(id)
    }
    if (length(cycle) > 0) {
      break
    }
  }
  cycle
}

scope_rds_candidates <- function(scope) {
  env_name <- paste0("METADATA_VALIDATION_RDS_", toupper(gsub("[^A-Za-z0-9]", "_", scope)))
  env_hit <- env_or(env_name, "")
  candidates <- character(0)
  if (nzchar(env_hit)) {
    candidates <- c(candidates, env_hit)
  }
  if (scope == "panorama") {
    candidates <- c(candidates, env_or("ANNOTATION_HUB_PATH", file.path(results_dir, "checkpoints", "03_after_annotation.rds")))
  }
  status_df <- read_tsv(layer_status_file)
  if (nrow(status_df) > 0 && all(c("layer_id", "annotated_rds") %in% colnames(status_df))) {
    hit <- status_df[status_df$layer_id == scope, "annotated_rds", drop = TRUE]
    hit <- hit[nzchar(hit)]
    candidates <- c(candidates, hit)
  }
  candidates <- c(candidates, file.path(results_dir, "04_subcluster", scope, "annotated.rds"))
  unique(candidates[nzchar(candidates)])
}

metadata_from_rds <- function(path) {
  obj <- readRDS(path)
  meta <- tryCatch(obj@meta.data, error = function(e) NULL)
  if (is.null(meta) && is.data.frame(obj)) {
    meta <- obj
  }
  if (is.null(meta) && is.list(obj) && is.data.frame(obj$meta.data)) {
    meta <- obj$meta.data
  }
  if (!is.data.frame(meta)) {
    stop("RDS does not expose a data.frame-like metadata table", call. = FALSE)
  }
  meta
}

load_scope_metadata <- function(scopes) {
  out <- list()
  for (scope in scopes) {
    candidates <- scope_rds_candidates(scope)
    hit <- candidates[file.exists(candidates) & file.info(candidates)$size > 0]
    if (length(hit) == 0) {
      out[[scope]] <- list(path = "", meta = NULL, reason = "no annotated RDS found")
      next
    }
    meta <- tryCatch(metadata_from_rds(hit[[1]]), error = function(e) e)
    if (inherits(meta, "error")) {
      out[[scope]] <- list(path = hit[[1]], meta = NULL, reason = conditionMessage(meta))
    } else {
      out[[scope]] <- list(path = hit[[1]], meta = meta, reason = "")
    }
  }
  out
}

questions <- read_tsv(question_path)
if (nrow(questions) == 0) {
  fail("questions.exists", sprintf("analysis questions table is missing or empty: %s", question_path))
} else {
  pass("questions.exists", sprintf("read %s rows from %s", nrow(questions), question_path))
}

missing_cols <- setdiff(required_cols, colnames(questions))
if (length(missing_cols) > 0) {
  fail("questions.required_columns", paste("missing columns:", paste(missing_cols, collapse = ", ")))
} else {
  pass("questions.required_columns", "all required columns are present")
}

if (nrow(questions) > 0 && length(missing_cols) == 0) {
  expected_materialized <- 80L
  if (nrow(questions) == expected_materialized) {
    pass("questions.row_count", "80 materialized rows present; F23-F26 are intentionally derived")
  } else {
    fail(
      "questions.row_count",
      sprintf("expected 80 materialized rows after derived-row removal, found %s", nrow(questions)),
      fix = "Check analysis_questions_FULL.md and keep F23-F26 as generator-derived rows."
    )
  }

  duplicated_ids <- unique(questions$question_id[duplicated(questions$question_id)])
  if (length(duplicated_ids) > 0) {
    fail("questions.unique_id", paste("duplicated question_id:", paste(duplicated_ids, collapse = ", ")))
  } else {
    pass("questions.unique_id", "question_id values are unique")
  }

  enum_checks <- list(
    scope = allowed_scopes,
    priority = allowed_priorities,
    status = allowed_status,
    contrast_axis = names(axis_targets)
  )
  for (col in names(enum_checks)) {
    bad <- unique(questions[[col]][!questions[[col]] %in% enum_checks[[col]]])
    bad <- bad[nzchar(bad)]
    if (length(bad) > 0) {
      fail(
        paste0("questions.enum.", col),
        sprintf("%s has unsupported values: %s", col, paste(bad, collapse = ", ")),
        fix = sprintf("Allowed values: %s", paste(enum_checks[[col]], collapse = ", "))
      )
    } else {
      pass(paste0("questions.enum.", col), sprintf("%s values are valid", col))
    }
  }

  for (idx in seq_len(nrow(questions))) {
    row <- questions[idx, , drop = FALSE]
    qid <- row$question_id[[1]]
    for (dsl_col in c("sender_groups", "receiver_groups")) {
      if (!valid_dsl(row[[dsl_col]][[1]])) {
        fail(
          "questions.dsl_syntax",
          sprintf("%s has invalid DSL: %s", dsl_col, row[[dsl_col]][[1]]),
          qid,
          "Use *, -, [A], [A,B], [A]|[B], [A]->[B], or prefix_*."
        )
      }
    }
    split <- parse_condition_split(row$condition_split[[1]])
    if (!isTRUE(split$ok)) {
      fail(
        "questions.condition_split_syntax",
        sprintf("invalid condition_split: %s", row$condition_split[[1]]),
        qid,
        "Use -, syf_only, f5_only, or column:value1,value2."
      )
    }
    if (!grepl("^[A-Za-z0-9_.+-]+$", row$tools_to_run[[1]], perl = TRUE)) {
      fail("questions.tools_syntax", sprintf("invalid tools_to_run: %s", row$tools_to_run[[1]]), qid)
    }
    if (identical(row$status[[1]], "active") && !tools_for_axis_ok(row$contrast_axis[[1]], row$tools_to_run[[1]])) {
      fail(
        "questions.axis_tools",
        sprintf("tools_to_run=%s is inconsistent with contrast_axis=%s", row$tools_to_run[[1]], row$contrast_axis[[1]]),
        qid,
        sprintf("Use one of: %s", paste(axis_tool_allow[[row$contrast_axis[[1]]]], collapse = ", "))
      )
    }
  }

  st_rows <- questions[questions$scope == "ST_section", , drop = FALSE]
  if (nrow(st_rows) > 0 && any(st_rows$status != "planned")) {
    fail("questions.st_planned", "all ST_section rows must remain planned until ST-H")
  } else {
    pass("questions.st_planned", sprintf("all %s ST_section rows are planned", nrow(st_rows)))
  }

  active_rows <- questions[questions$status == "active", , drop = FALSE]
  unmapped_axis <- unique(active_rows$contrast_axis[!active_rows$contrast_axis %in% names(axis_targets)])
  if (length(unmapped_axis) > 0) {
    fail("questions.coverage_mapping", paste("active axes without Tier 2 mapping:", paste(unmapped_axis, collapse = ", ")))
  } else {
    pass("questions.coverage_mapping", sprintf("%s active rows have known Tier 2 target mappings", nrow(active_rows)))
  }

  dep_ids <- unique(unlist(lapply(questions$depends_on[nzchar(questions$depends_on)], split_dependency_ids), use.names = FALSE))
  dep_ids <- trimws(dep_ids[nzchar(trimws(dep_ids))])
  dep_known <- vapply(dep_ids, function(dep) {
    dep %in% questions$question_id || any(startsWith(questions$question_id, paste0(dep, "_")))
  }, logical(1))
  missing_dep_ids <- dep_ids[!dep_known]
  if (length(missing_dep_ids) > 0) {
    fail("questions.depends_on", paste("unknown depends_on IDs:", paste(missing_dep_ids, collapse = ", ")))
  } else {
    pass("questions.depends_on", "depends_on references are valid or empty")
  }

  dep_cycle <- find_dependency_cycle(questions)
  if (length(dep_cycle) > 0) {
    fail("questions.depends_on_cycle", paste("depends_on cycle detected:", paste(dep_cycle, collapse = " -> ")))
  } else {
    pass("questions.depends_on_cycle", "depends_on graph is acyclic")
  }

  scope_meta <- load_scope_metadata(setdiff(unique(active_rows$scope), "ST_section"))
  for (scope in names(scope_meta)) {
    item <- scope_meta[[scope]]
    if (is.null(item$meta)) {
      warn(
        "semantic.scope_metadata_skipped",
        sprintf("semantic checks skipped for %s: %s", scope, item$reason),
        fix = "Set METADATA_VALIDATION_RDS_<SCOPE> or provide layer_status.tsv annotated_rds paths."
      )
    } else {
      pass("semantic.scope_metadata_loaded", sprintf("loaded metadata for %s from %s", scope, item$path))
    }
  }

  for (idx in seq_len(nrow(active_rows))) {
    row <- active_rows[idx, , drop = FALSE]
    qid <- row$question_id[[1]]
    item <- scope_meta[[row$scope[[1]]]]
    if (is.null(item) || is.null(item$meta)) {
      next
    }
    meta <- item$meta
    split <- parse_condition_split(row$condition_split[[1]])
    if (nzchar(split$column) && !split$column %in% colnames(meta)) {
      alternatives <- c("group_id", "group", "condition", "section")
      if (!any(alternatives %in% colnames(meta))) {
        fail(
          "semantic.condition_column",
          sprintf("metadata for %s lacks condition_split column %s", row$scope[[1]], split$column),
          qid,
          "Add the split column to object metadata or adjust condition_split."
        )
      } else {
        warn(
          "semantic.condition_column",
          sprintf("metadata for %s lacks %s; available fallback-like columns: %s", row$scope[[1]], split$column, paste(intersect(alternatives, colnames(meta)), collapse = ", ")),
          qid
        )
      }
    }

    if ("cell_subtype" %in% colnames(meta)) {
      observed <- unique(trimws(as.character(meta$cell_subtype)))
      observed <- observed[nzchar(observed)]
      tokens <- unique(c(extract_group_tokens(row$sender_groups[[1]]), extract_group_tokens(row$receiver_groups[[1]])))
      missing_tokens <- setdiff(tokens, observed)
      if (length(missing_tokens) > 0) {
        fail(
          "semantic.cell_subtype_exists",
          sprintf("cell_subtype values not found in %s metadata: %s", row$scope[[1]], paste(missing_tokens, collapse = ", ")),
          qid,
          "Fix group spelling or rerun annotation/backfill."
        )
      }
    } else {
      warn("semantic.cell_subtype_column", sprintf("metadata for %s lacks cell_subtype column", row$scope[[1]]), qid)
    }
  }
}

generated_tables <- list(
  comparisons = list(path = file.path(metadata_dir, "comparisons.tsv"), id = "comparison_id"),
  communication_pairs = list(path = file.path(metadata_dir, "communication_pairs.tsv"), id = "pair_id"),
  trajectory_pairs = list(path = file.path(metadata_dir, "trajectory_pairs.tsv"), id = "trajectory_id"),
  scenic_targets = list(path = file.path(metadata_dir, "scenic_targets.tsv"), id = "target_id"),
  enrichment_targets = list(path = file.path(metadata_dir, "enrichment_targets.tsv"), id = "target_id"),
  gene_program_targets = list(path = file.path(metadata_dir, "gene_program_targets.tsv"), id = "comparison_id"),
  deconv_pairs = list(path = file.path(metadata_dir, "deconv_pairs.tsv"), id = "deconv_id"),
  spatial_pairs = list(path = file.path(metadata_dir, "spatial_pairs.tsv"), id = "spatial_pair_id")
)

loaded_generated <- list()
for (name in names(generated_tables)) {
  spec <- generated_tables[[name]]
  df <- read_tsv(spec$path)
  loaded_generated[[name]] <- df
  if (nrow(df) == 0) {
    if (spec$id %in% colnames(df)) {
      pass("tier2.empty_table", sprintf("%s is empty with a valid header", basename(spec$path)))
    } else {
      warn("tier2.empty_table", sprintf("%s has no data rows or header", basename(spec$path)))
    }
    next
  }
  if (!spec$id %in% colnames(df)) {
    fail("tier2.id_column", sprintf("%s missing ID column %s", basename(spec$path), spec$id))
    next
  }
  duplicated_ids <- unique(df[[spec$id]][duplicated(df[[spec$id]])])
  duplicated_ids <- duplicated_ids[nzchar(duplicated_ids)]
  if (length(duplicated_ids) > 0) {
    fail("tier2.unique_id", sprintf("%s duplicated IDs: %s", basename(spec$path), paste(duplicated_ids, collapse = ", ")))
  } else {
    pass("tier2.unique_id", sprintf("%s ID values are unique", basename(spec$path)))
  }
}

if (exists("questions") && nrow(questions) > 0 && length(missing_cols) == 0) {
  comparisons <- loaded_generated$comparisons
  nichenet_rows <- questions[questions$status == "active" & grepl("nichenet", questions$tools_to_run, ignore.case = TRUE), , drop = FALSE]
  comparisons_is_m3_generated <- nrow(comparisons) > 0 &&
    all(c("source_question_id", "contrast_axis") %in% colnames(comparisons))
  if (nrow(nichenet_rows) > 0 && !comparisons_is_m3_generated) {
    if (nrow(comparisons) == 0) {
      pass("semantic.nichenet_deg_consistency", "comparisons.tsv is empty, so NicheNet/DEG consistency is deferred until M3")
    } else {
      pass("semantic.nichenet_deg_consistency", "comparisons.tsv is not M3-generated, so NicheNet/DEG consistency is deferred until M3")
    }
  } else if (nrow(nichenet_rows) > 0) {
    searchable <- paste(apply(comparisons, 1, paste, collapse = " "), collapse = "\n")
    missing_nichenet_support <- 0L
    for (idx in seq_len(nrow(nichenet_rows))) {
      row <- nichenet_rows[idx, , drop = FALSE]
      receivers <- extract_group_tokens(row$receiver_groups[[1]])
      if (length(receivers) == 0) {
        next
      }
      hits <- vapply(receivers, function(token) grepl(token, searchable, fixed = TRUE), logical(1))
      if (!any(hits)) {
        missing_nichenet_support <- missing_nichenet_support + 1L
        fail(
          "semantic.nichenet_deg_consistency",
          sprintf("no comparison row appears to support NicheNet receiver(s): %s", paste(receivers, collapse = ", ")),
          row$question_id[[1]],
          "Ensure M3 fan-out creates matching receiver DEG comparisons."
        )
      }
    }
    if (missing_nichenet_support == 0L) {
      pass("semantic.nichenet_deg_consistency", sprintf("%s active NicheNet rows have generated comparison support", nrow(nichenet_rows)))
    }
  }
}

level_counts <- table(factor(vapply(records, `[[`, character(1), "level"), levels = c("PASS", "WARN", "FAIL", "SKIP")))
summary <- list(
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  question_file = question_path,
  report_file = report_path,
  pass = unname(level_counts[["PASS"]]),
  warn = unname(level_counts[["WARN"]]),
  fail = unname(level_counts[["FAIL"]]),
  skip = unname(level_counts[["SKIP"]])
)

dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
json <- paste0(
  "{\n",
  "  \"summary\": ", json_object(summary), ",\n",
  "  \"records\": [\n    ",
  paste(vapply(records, json_object, character(1)), collapse = ",\n    "),
  "\n  ]\n",
  "}\n"
)
writeLines(json, report_path, useBytes = TRUE)

cat("metadata validation report\n")
cat(sprintf("PASS=%s WARN=%s FAIL=%s SKIP=%s\n", summary$pass, summary$warn, summary$fail, summary$skip))
for (rec in records) {
  cat(sprintf("[%s] %s", rec$level, rec$check_id))
  if (nzchar(rec$question_id)) {
    cat(sprintf(" %s", rec$question_id))
  }
  cat(sprintf(" - %s\n", rec$message))
  if (nzchar(rec$fix)) {
    cat(sprintf("      fix: %s\n", rec$fix))
  }
}
cat(sprintf("JSON=%s\n", report_path))

if (summary$fail > 0) {
  quit(status = 1)
}
