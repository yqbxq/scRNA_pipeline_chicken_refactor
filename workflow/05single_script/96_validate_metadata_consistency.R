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

semantic_object_problem <- function(check_id, message, question_id = "", fix = "") {
  hard <- tolower(env_or("METADATA_VALIDATION_SEMANTIC_HARD", "no")) %in% c("yes", "true", "1", "on")
  if (hard) {
    fail(check_id, message, question_id, fix)
  } else {
    warn(check_id, message, question_id, fix)
  }
}

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
  "depends_on", "notes", "activation_policy", "min_sender_cells",
  "min_receiver_cells", "min_cells_per_condition", "fallback_pair_id",
  "derived_from_pair_id", "run_baseline_if_split_fails"
)

allowed_scopes <- c("panorama", "GC_subcluster", "TC_subcluster", "ST_section")
allowed_priorities <- c("P0", "P1", "P2", "P3")
allowed_status <- c("active", "planned")
allowed_activation_policy <- c("always", "auto_if_min_cells", "derived_from_split")
axis_targets <- c(
  cluster_marker = "annotation_marker_targets.tsv",
  identity_marker = "comparisons.tsv",
  directional_DEG = "comparisons.tsv",
  pairwise = "comparisons.tsv",
  contrast_only = "comparisons.tsv",
  global_stage_context = "comparisons.tsv",
  composition = "comparisons.tsv",
  qc_composition = "comparisons.tsv",
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
  cluster_marker = c("annotation", "deg", "enrichment"),
  identity_marker = c("deg", "enrichment"),
  directional_DEG = c("deg", "enrichment"),
  pairwise = c("deg", "enrichment"),
  contrast_only = c("deg", "enrichment"),
  global_stage_context = c("deg", "enrichment"),
  composition = c("composition_test"),
  qc_composition = c("composition_qc", "composition_test"),
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
  spatial_clustering = c("spatial_clustering"),
  spatial_neighbor = c("squidpy"),
  spatial_overlay = c("rctd", "stlearn", "custom", "visualization"),
  deconv_pair = c("rctd", "cell2location", "card", "seurat_transfer"),
  SVG = c("spark", "sparkx", "spatialde", "spatialde2"),
  SVG_diff = c("spark", "sparkx", "spatialde", "spatialde2"),
  deconv_validation = c("scdesign3", "rctd", "cell2location"),
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

positive_int_value <- function(x) {
  value <- suppressWarnings(as.integer(trim(x)))
  if (is.na(value) || value <= 0L) NA_integer_ else value
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
st_fanout_version <- suppressWarnings(as.integer(env_or("ST_FANOUT_VERSION", "0")))
if (is.na(st_fanout_version)) {
  st_fanout_version <- 0L
}
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

optional_question_cols <- c("display_question_id", "output_alias", "report_title")
for (col in optional_question_cols) {
  if (!col %in% colnames(questions)) {
    questions[[col]] <- character(nrow(questions))
  }
}

if (nrow(questions) > 0 && length(missing_cols) == 0) {
  expected_materialized <- 83L
  if (nrow(questions) == expected_materialized) {
    pass("questions.row_count", "83 materialized rows present; F23-F26 are intentionally derived")
  } else {
    warn(
      "questions.row_count",
      sprintf("expected 83 materialized rows after derived-row removal, found %s", nrow(questions)),
      fix = "Check analysis_questions_FULL.md and keep F23-F26 as generator-derived rows."
    )
  }

  duplicated_ids <- unique(questions$question_id[duplicated(questions$question_id)])
  if (length(duplicated_ids) > 0) {
    fail("questions.unique_id", paste("duplicated question_id:", paste(duplicated_ids, collapse = ", ")))
  } else {
    pass("questions.unique_id", "question_id values are unique")
  }

  e03_question <- questions[questions$question_id == "E03_layer_compo", , drop = FALSE]
  if (nrow(e03_question) != 1L) {
    fail("questions.E03_qc_alias", "E03_layer_compo must exist exactly once")
  } else {
    e03_alias_bad <- trim(e03_question$display_question_id[[1]]) != "E03_scRNA_GC_TC_capture_balance" ||
      trim(e03_question$output_alias[[1]]) != "E03_scRNA_GC_TC_capture_balance" ||
      !grepl("captured-cell balance QC", trim(e03_question$report_title[[1]]), fixed = TRUE) ||
      e03_question$question_class[[1]] != "qc" ||
      e03_question$contrast_axis[[1]] != "qc_composition" ||
      e03_question$tools_to_run[[1]] != "composition_qc"
    if (e03_alias_bad) {
      fail(
        "questions.E03_qc_alias",
        "E03 must be qc_composition with display/output alias E03_scRNA_GC_TC_capture_balance and captured-cell QC title",
        "E03_layer_compo"
      )
    } else {
      pass("questions.E03_qc_alias", "E03 has the required captured-cell balance QC alias and Tier1 mode")
    }
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
    policy <- tolower(trim(row$activation_policy[[1]]))
    if (!nzchar(policy)) {
      policy <- "always"
    }
    if (!policy %in% allowed_activation_policy) {
      fail(
        "questions.activation_policy",
        sprintf("%s has unsupported activation_policy=%s", qid, policy),
        qid,
        sprintf("Allowed values: %s", paste(allowed_activation_policy, collapse = ", "))
      )
    }
    if (identical(policy, "auto_if_min_cells")) {
      min_values <- c(
        min_sender_cells = positive_int_value(row$min_sender_cells[[1]]),
        min_receiver_cells = positive_int_value(row$min_receiver_cells[[1]]),
        min_cells_per_condition = positive_int_value(row$min_cells_per_condition[[1]])
      )
      if (any(is.na(min_values))) {
        fail(
          "questions.activation_policy.min_cells",
          sprintf("%s auto_if_min_cells requires positive integer min_* fields", qid),
          qid
        )
      }
      run_baseline <- tolower(trim(row$run_baseline_if_split_fails[[1]])) %in% c("yes", "true", "1", "on")
      if (run_baseline && !nzchar(trim(row$fallback_pair_id[[1]]))) {
        fail(
          "questions.activation_policy.fallback",
          sprintf("%s run_baseline_if_split_fails=yes requires fallback_pair_id", qid),
          qid
        )
      }
    }
  }

  st_rows <- questions[questions$scope == "ST_section", , drop = FALSE]
  allowed_st_m2_active <- c(
    "I01_ST_clustering", "I02_ST_region_marker", "I06_SVG",
    "I07_SVG_by_stage", "I08_deconv_panorama",
    "I12_ST_region_compo", "I13_neighborhood", "I15_ST_enrichment"
  )
  if (st_fanout_version < 1L) {
    if (nrow(st_rows) > 0 && any(st_rows$status != "planned")) {
      fail("questions.st_planned", "all ST_section rows must remain planned when ST_FANOUT_VERSION=0")
    } else {
      pass("questions.st_planned", sprintf("all %s ST_section rows are planned", nrow(st_rows)))
    }
  } else {
    active_st_rows <- st_rows[st_rows$status == "active", , drop = FALSE]
    disallowed_active <- setdiff(active_st_rows$question_id, allowed_st_m2_active)
    if (length(disallowed_active) > 0) {
      fail(
        "questions.st_planned",
        sprintf("ST_FANOUT_VERSION=1 only allows M2 ST active IDs; disallowed active IDs: %s", paste(disallowed_active, collapse = ", "))
      )
    } else {
      pass(
        "questions.st_planned",
        sprintf("ST_FANOUT_VERSION=1 allows M2 ST fan-out while non-M2 ST rows remain planned (%s ST rows total)", nrow(st_rows))
      )
    }
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
  annotation_marker_targets = list(path = file.path(metadata_dir, "annotation_marker_targets.tsv"), id = "target_id"),
  communication_pairs = list(path = file.path(metadata_dir, "communication_pairs.tsv"), id = "pair_id"),
  trajectory_pairs = list(path = file.path(metadata_dir, "trajectory_pairs.tsv"), id = "trajectory_id"),
  scenic_targets = list(path = file.path(metadata_dir, "scenic_targets.tsv"), id = "target_id"),
  enrichment_targets = list(path = file.path(metadata_dir, "enrichment_targets.tsv"), id = "target_id"),
  gene_program_targets = list(path = file.path(metadata_dir, "gene_program_targets.tsv"), id = "comparison_id"),
  deconv_pairs = list(path = file.path(metadata_dir, "deconv_pairs.tsv"), id = "deconv_id"),
  spatial_pairs = list(path = file.path(metadata_dir, "spatial_pairs.tsv"), id = "spatial_pair_id"),
  scdesign3_question_map = list(path = file.path(metadata_dir, "scdesign3_question_map.tsv"), id = "question_id"),
  scdesign3_targets = list(path = file.path(metadata_dir, "scdesign3_targets.tsv"), id = "target_id"),
  scdesign3_simulation_designs = list(path = file.path(metadata_dir, "scdesign3_simulation_designs.tsv"), id = "simulation_design_id"),
  scdesign3_thresholds = list(path = file.path(metadata_dir, "scdesign3_thresholds.tsv"), id = "threshold_id")
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

require_generated_cols <- function(table_name, df, cols) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    fail(
      paste0("tier2.schema.", table_name),
      sprintf("%s missing required columns: %s", table_name, paste(missing, collapse = ", "))
    )
    return(FALSE)
  }
  pass(paste0("tier2.schema.", table_name), sprintf("%s contains required schema columns", table_name))
  TRUE
}

split_tokens_any <- function(x) {
  x <- trim(x)
  if (!nzchar(x) || x %in% c("-", "*")) {
    return(if (identical(x, "*")) "*" else character(0))
  }
  out <- trimws(unlist(strsplit(x, "[,;+]", perl = TRUE), use.names = FALSE))
  out[nzchar(out)]
}

known_spatial_filter_tokens <- function() {
  sections <- read_tsv(env_or("SECTION_SHEET", file.path(metadata_dir, "sections.tsv")))
  samples <- read_tsv(env_or("CANONICAL_SAMPLE_SHEET", file.path(metadata_dir, "samples.canonical.tsv")))
  if (nrow(samples) == 0) {
    samples <- read_tsv(file.path(metadata_dir, "samples.tsv"))
  }
  section_ids <- character(0)
  conditions <- character(0)
  if (nrow(sections) > 0) {
    if ("section_id" %in% colnames(sections)) {
      section_ids <- c(section_ids, trim(sections$section_id))
    }
    if ("condition" %in% colnames(sections)) {
      conditions <- c(conditions, trim(sections$condition))
    }
  }
  if (nrow(samples) > 0) {
    if ("section_id" %in% colnames(samples)) {
      section_ids <- c(section_ids, trim(samples$section_id))
    }
    if ("condition" %in% colnames(samples)) {
      conditions <- c(conditions, trim(samples$condition))
    }
    if ("group_id" %in% colnames(samples)) {
      conditions <- c(conditions, trim(samples$group_id))
    }
  }
  unique(c("*", "syf_only", "f5_only", "syf_1", "f5_1", section_ids[nzchar(section_ids)], conditions[nzchar(conditions)]))
}

split_comm_tokens <- function(x) {
  x <- trim(x)
  if (!nzchar(x) || x %in% c("-", "*")) {
    return(character(0))
  }
  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE))
  out <- out[nzchar(out) & !out %in% reserved_group_tokens]
  out <- out[!grepl("[*]$", out)]
  unique(out)
}

allowed_analysis_modes <- c(
  "annotation_cluster_marker", "subtype_marker", "subtype_pairwise",
  "condition_within_type", "composition", "qc_composition", "global_context"
)
allowed_gene_program_roles <- c(
  "annotation_marker", "receiver_marker", "subtype_pairwise_deg",
  "condition_deg", "global_context", "qc_only", "none"
)
allowed_yes_no_contextual <- c("yes", "no", "contextual")
allowed_nichenet_usage <- c("none", "baseline_receiver_marker", "receiver_condition_deg")
allowed_enrichment_usage <- c(
  "none", "identity_baseline_enrichment", "subtype_pairwise_enrichment",
  "mechanism_enrichment", "global_context_enrichment"
)

if (!is.null(loaded_generated$scdesign3_question_map)) {
  scdesign3_question_map <- loaded_generated$scdesign3_question_map
  scdesign3_map_schema <- c(
    "question_id", "section", "question_family", "status",
    "scdesign3_role", "validation_layer", "validation_object",
    "ground_truth_col", "simulation_type", "simulation_unit",
    "primary_metric", "secondary_metrics", "gate_target", "affects_modules",
    "required_before_core_interpretation", "derived_parent_question_id",
    "enabled", "reason", "materialized", "question_class", "scope",
    "contrast_axis", "tools_to_run", "question_zh"
  )
  if (require_generated_cols("scdesign3_question_map.tsv", scdesign3_question_map, scdesign3_map_schema) && nrow(scdesign3_question_map) > 0) {
    required_derived <- c(
      "F23_TC_GC_diff_overall",
      "F24_TC_GC_diff_subtype",
      "F25_GC_internal_diff",
      "F26_panorama_screen_diff"
    )
    missing_materialized <- if (exists("questions") && nrow(questions) > 0 && "question_id" %in% colnames(questions)) {
      setdiff(questions$question_id, scdesign3_question_map$question_id)
    } else {
      character(0)
    }
    missing_derived <- setdiff(required_derived, scdesign3_question_map$question_id)
    if (length(missing_materialized) > 0) {
      fail("tier2.scdesign3_question_map.materialized_coverage", sprintf("materialized questions missing from scDesign3 map: %s", paste(missing_materialized, collapse = ", ")))
    } else {
      pass("tier2.scdesign3_question_map.materialized_coverage", sprintf("all %s materialized questions are present in scDesign3 map", if (exists("questions")) nrow(questions) else 0L))
    }
    if (length(missing_derived) > 0) {
      fail("tier2.scdesign3_question_map.derived_coverage", sprintf("derived F23-F26 questions missing from scDesign3 map: %s", paste(missing_derived, collapse = ", ")))
    } else {
      pass("tier2.scdesign3_question_map.derived_coverage", "F23-F26 derived communication questions are present in scDesign3 map")
    }

    allowed_roles <- c(
      "direct_validate", "upstream_gate", "downstream_consistency",
      "not_applicable", "planned_waiting_input",
      "derived_from_validated_parent"
    )
    bad_roles <- unique(scdesign3_question_map$scdesign3_role[!scdesign3_question_map$scdesign3_role %in% allowed_roles])
    bad_roles <- bad_roles[nzchar(bad_roles)]
    if (length(bad_roles) > 0) {
      fail("tier2.scdesign3_question_map.roles", sprintf("unsupported scdesign3_role values: %s", paste(bad_roles, collapse = ", ")))
    } else {
      pass("tier2.scdesign3_question_map.roles", "scDesign3 roles use the approved role vocabulary")
    }

    bad_enabled <- unique(scdesign3_question_map$enabled[!scdesign3_question_map$enabled %in% c("yes", "planned", "no")])
    bad_enabled <- bad_enabled[nzchar(bad_enabled)]
    if (length(bad_enabled) > 0) {
      fail("tier2.scdesign3_question_map.enabled", sprintf("enabled must be yes/planned/no: %s", paste(bad_enabled, collapse = ", ")))
    } else {
      pass("tier2.scdesign3_question_map.enabled", "scDesign3 enabled values are valid")
    }

    derived_bad <- scdesign3_question_map[
      scdesign3_question_map$question_id %in% required_derived &
        (scdesign3_question_map$scdesign3_role != "derived_from_validated_parent" |
          !nzchar(scdesign3_question_map$derived_parent_question_id)),
      ,
      drop = FALSE
    ]
    if (nrow(derived_bad) > 0) {
      fail("tier2.scdesign3_question_map.derived_role", sprintf("derived rows must inherit parent gates: %s", paste(derived_bad$question_id, collapse = ", ")))
    } else {
      pass("tier2.scdesign3_question_map.derived_role", "derived scDesign3 rows inherit parent gate status")
    }
  }
}

if (!is.null(loaded_generated$scdesign3_targets)) {
  scdesign3_targets <- loaded_generated$scdesign3_targets
  scdesign3_target_schema <- c(
    "target_id", "target_type", "layer_id", "input_object", "truth_col",
    "questions_covered", "n_simulations", "resolution_grid",
    "mixture_design", "primary_metric", "pass_threshold", "warn_threshold",
    "fail_threshold", "output_dir", "status"
  )
  if (require_generated_cols("scdesign3_targets.tsv", scdesign3_targets, scdesign3_target_schema)) {
    allowed_target_types <- c(
      "cluster_robustness", "composition_robustness", "trajectory_robustness",
      "communication_robustness", "deconv_validation", "spatial_validation"
    )
    bad_target_types <- unique(scdesign3_targets$target_type[!scdesign3_targets$target_type %in% allowed_target_types])
    bad_target_types <- bad_target_types[nzchar(bad_target_types)]
    if (length(bad_target_types) > 0) {
      fail("tier2.scdesign3_targets.target_type", sprintf("unsupported target_type values: %s", paste(bad_target_types, collapse = ", ")))
    } else {
      pass("tier2.scdesign3_targets.target_type", "scDesign3 target types are valid")
    }

    if (!is.null(loaded_generated$scdesign3_question_map) && nrow(scdesign3_targets) > 0) {
      covered <- unique(unlist(lapply(scdesign3_targets$questions_covered, split_dependency_ids), use.names = FALSE))
      covered <- covered[nzchar(covered)]
      target_required <- loaded_generated$scdesign3_question_map[
        !loaded_generated$scdesign3_question_map$scdesign3_role %in% c("not_applicable", "derived_from_validated_parent") &
          nzchar(loaded_generated$scdesign3_question_map$validation_object),
        ,
        drop = FALSE
      ]
      missing_target_coverage <- setdiff(target_required$question_id, covered)
      if (length(missing_target_coverage) > 0) {
        fail("tier2.scdesign3_targets.question_coverage", sprintf("scDesign3 target-required questions missing from targets: %s", paste(missing_target_coverage, collapse = ", ")))
      } else {
        pass("tier2.scdesign3_targets.question_coverage", sprintf("%s target-required scDesign3 questions are covered by targets", nrow(target_required)))
      }
    }
  }
}

if (!is.null(loaded_generated$scdesign3_simulation_designs)) {
  require_generated_cols(
    "scdesign3_simulation_designs.tsv",
    loaded_generated$scdesign3_simulation_designs,
    c(
      "simulation_design_id", "target_id", "simulation_type",
      "simulation_unit", "n_simulations", "resolution_grid",
      "mixture_design", "max_cells_per_label", "n_hvg", "status", "notes"
    )
  )
}

if (!is.null(loaded_generated$scdesign3_thresholds)) {
  require_generated_cols(
    "scdesign3_thresholds.tsv",
    loaded_generated$scdesign3_thresholds,
    c(
      "threshold_id", "target_type", "primary_metric", "pass_threshold",
      "warn_threshold", "fail_threshold", "notes"
    )
  )
}

if (!is.null(loaded_generated$spatial_pairs)) {
  spatial_pairs <- loaded_generated$spatial_pairs
  spatial_schema <- c(
    "spatial_pair_id", "source_question_id", "st_scope", "sender", "receiver",
    "contrast_axis", "section_filter", "condition_split_var",
    "condition_split_values", "tool", "enabled", "notes"
  )
  if (require_generated_cols("spatial_pairs.tsv", spatial_pairs, spatial_schema) && nrow(spatial_pairs) > 0) {
    allowed_st_tools <- c(
      "spatial_clustering", "squidpy", "rctd", "cell2location", "card",
      "seurat_transfer", "sparkx", "spark", "spatialde2", "spatialde",
      "stagate", "scdesign3", "stlearn", "decoupler", "visualization", "custom"
    )
    allowed_filters <- known_spatial_filter_tokens()
    spatial_fail_n <- 0L
    active_spatial <- spatial_pairs[tolower(trim(spatial_pairs$enabled)) != "no", , drop = FALSE]
    for (idx in seq_len(nrow(active_spatial))) {
      row <- active_spatial[idx, , drop = FALSE]
      tools <- tolower(split_tokens_any(row$tool[[1]]))
      bad_tools <- setdiff(tools, allowed_st_tools)
      if (length(tools) == 0 || length(bad_tools) > 0) {
        fail("tier2.spatial_pairs.tool", sprintf("%s has unsupported tool tokens: %s", row$spatial_pair_id[[1]], paste(bad_tools %||% tools, collapse = ", ")))
        spatial_fail_n <- spatial_fail_n + 1L
      }
      filters <- split_tokens_any(row$section_filter[[1]])
      bad_filters <- setdiff(filters, allowed_filters)
      if (length(filters) == 0 || length(bad_filters) > 0) {
        fail("tier2.spatial_pairs.section_filter", sprintf("%s has unsupported section_filter tokens: %s", row$spatial_pair_id[[1]], paste(bad_filters %||% filters, collapse = ", ")))
        spatial_fail_n <- spatial_fail_n + 1L
      }
    }
    if (spatial_fail_n == 0L) {
      pass("tier2.spatial_pairs.contract", sprintf("%s enabled spatial rows use valid tools and section filters", nrow(active_spatial)))
    }
  }
}

if (!is.null(loaded_generated$deconv_pairs)) {
  deconv_pairs <- loaded_generated$deconv_pairs
  deconv_schema <- c(
    "deconv_id", "source_question_id", "st_scope", "reference_scope",
    "section_filter", "condition_split_var", "condition_split_values",
    "tool", "enabled", "notes"
  )
  if (require_generated_cols("deconv_pairs.tsv", deconv_pairs, deconv_schema) && nrow(deconv_pairs) > 0) {
    allowed_deconv_tools <- c("rctd", "cell2location", "card", "seurat_transfer", "scdesign3", "custom")
    allowed_reference_scopes <- c("panorama", "GC_subcluster", "TC_subcluster", "synthetic_panorama", "*")
    allowed_filters <- known_spatial_filter_tokens()
    deconv_fail_n <- 0L
    active_deconv <- deconv_pairs[tolower(trim(deconv_pairs$enabled)) != "no", , drop = FALSE]
    for (idx in seq_len(nrow(active_deconv))) {
      row <- active_deconv[idx, , drop = FALSE]
      tools <- tolower(split_tokens_any(row$tool[[1]]))
      bad_tools <- setdiff(tools, allowed_deconv_tools)
      if (length(tools) == 0 || length(bad_tools) > 0) {
        fail("tier2.deconv_pairs.tool", sprintf("%s has unsupported tool tokens: %s", row$deconv_id[[1]], paste(bad_tools %||% tools, collapse = ", ")))
        deconv_fail_n <- deconv_fail_n + 1L
      }
      scopes <- split_tokens_any(row$reference_scope[[1]])
      bad_scopes <- setdiff(scopes, allowed_reference_scopes)
      if (length(scopes) == 0 || length(bad_scopes) > 0) {
        fail("tier2.deconv_pairs.reference_scope", sprintf("%s has unsupported reference_scope tokens: %s", row$deconv_id[[1]], paste(bad_scopes %||% scopes, collapse = ", ")))
        deconv_fail_n <- deconv_fail_n + 1L
      }
      filters <- split_tokens_any(row$section_filter[[1]])
      bad_filters <- setdiff(filters, allowed_filters)
      if (length(filters) == 0 || length(bad_filters) > 0) {
        fail("tier2.deconv_pairs.section_filter", sprintf("%s has unsupported section_filter tokens: %s", row$deconv_id[[1]], paste(bad_filters %||% filters, collapse = ", ")))
        deconv_fail_n <- deconv_fail_n + 1L
      }
    }
    if (deconv_fail_n == 0L) {
      pass("tier2.deconv_pairs.contract", sprintf("%s enabled deconv rows use valid tools, references, and section filters", nrow(active_deconv)))
    }
  }
}

if (!is.null(loaded_generated$trajectory_pairs)) {
  trajectory_pairs <- loaded_generated$trajectory_pairs
  trajectory_schema <- c(
    "trajectory_id", "source_question_id", "layer_scope", "root_group",
    "terminal_group", "condition_split_var", "condition_split_values",
    "method", "tools_to_run", "methods_extra", "outlier_qc_policy",
    "regress_cell_cycle", "coarse_label_var", "fine_label_var",
    "split_mode", "enabled", "notes"
  )
  if (require_generated_cols("trajectory_pairs.tsv", trajectory_pairs, trajectory_schema)) {
    observed <- colnames(trajectory_pairs)[seq_along(trajectory_schema)]
    if (identical(observed, trajectory_schema)) {
      pass("tier2.trajectory_pairs.column_order", "trajectory_pairs.tsv uses the v4 17-column order")
    } else {
      fail(
        "tier2.trajectory_pairs.column_order",
        sprintf("trajectory_pairs.tsv first %s columns do not match v4 order", length(trajectory_schema))
      )
    }
  }

  if (nrow(trajectory_pairs) > 0 && all(trajectory_schema %in% colnames(trajectory_pairs))) {
    trajectory_pairs$method <- tolower(vapply(trajectory_pairs$method, trim, character(1)))
    trajectory_pairs$enabled <- tolower(vapply(trajectory_pairs$enabled, trim, character(1)))
    active_trajectory_pairs <- trajectory_pairs[trajectory_pairs$enabled != "no", , drop = FALSE]
    allowed_methods <- c("trajectory", "velocity")
    bad_methods <- unique(active_trajectory_pairs$method[!active_trajectory_pairs$method %in% allowed_methods])
    bad_methods <- bad_methods[nzchar(bad_methods)]
    if (length(bad_methods) > 0) {
      fail("tier2.trajectory_pairs.method", sprintf("unsupported method values: %s", paste(bad_methods, collapse = ", ")))
    } else {
      pass("tier2.trajectory_pairs.method", "trajectory pair methods are valid")
    }

    enum_checks <- list(
      outlier_qc_policy = c("strict", "standard", "off"),
      regress_cell_cycle = c("yes", "no", "auto"),
      split_mode = c("auto", "force_split", "force_pooled", "pooled")
    )
    for (col in names(enum_checks)) {
      values <- tolower(vapply(active_trajectory_pairs[[col]], trim, character(1)))
      values[!nzchar(values)] <- switch(col, outlier_qc_policy = "standard", regress_cell_cycle = "auto", split_mode = "auto")
      bad <- unique(values[!values %in% enum_checks[[col]]])
      bad <- bad[nzchar(bad)]
      if (length(bad) > 0) {
        fail(paste0("tier2.trajectory_pairs.", col), sprintf("unsupported %s values: %s", col, paste(bad, collapse = ", ")))
      } else {
        pass(paste0("tier2.trajectory_pairs.", col), sprintf("%s values are valid", col))
      }
    }

    split_method_tokens <- function(x) {
      parts <- trimws(unlist(strsplit(trim(x), "[,;[:space:]]+", perl = TRUE), use.names = FALSE))
      parts[nzchar(parts)]
    }
    allowed_tools <- list(
      trajectory = c("slingshot", "monocle3", "monocle2", "paga-dpt", "tradeseq", "palantir"),
      velocity = c("scvelo-dynamical", "scvelo-stochastic", "velocyto", "cellrank")
    )
    tool_fail_n <- 0L
    extra_fail_n <- 0L
    for (idx in seq_len(nrow(active_trajectory_pairs))) {
      row <- active_trajectory_pairs[idx, , drop = FALSE]
      method <- row$method[[1]]
      tools <- tolower(gsub("_", "-", split_method_tokens(row$tools_to_run[[1]]), fixed = TRUE))
      allowed <- allowed_tools[[method]] %||% character(0)
      bad_tools <- setdiff(tools, allowed)
      if (length(tools) == 0) {
        fail("tier2.trajectory_pairs.tools_to_run", sprintf("%s has empty tools_to_run", row$trajectory_id[[1]]))
        tool_fail_n <- tool_fail_n + 1L
      } else if (length(bad_tools) > 0) {
        fail("tier2.trajectory_pairs.tools_to_run", sprintf("%s has unsupported tools_to_run tokens: %s", row$trajectory_id[[1]], paste(bad_tools, collapse = ", ")))
        tool_fail_n <- tool_fail_n + length(bad_tools)
      }
      extra_tokens <- split_method_tokens(row$methods_extra[[1]])
      bad_extra <- extra_tokens[!grepl("^(\\+|-|no[-_])?[A-Za-z0-9_.-]+$", extra_tokens, perl = TRUE)]
      if (length(bad_extra) > 0) {
        fail("tier2.trajectory_pairs.methods_extra", sprintf("%s has invalid methods_extra tokens: %s", row$trajectory_id[[1]], paste(bad_extra, collapse = ", ")))
        extra_fail_n <- extra_fail_n + length(bad_extra)
      }
    }
    if (tool_fail_n == 0L) {
      pass("tier2.trajectory_pairs.tools_to_run", "tools_to_run tokens are valid for trajectory/velocity rows")
    }
    if (extra_fail_n == 0L) {
      pass("tier2.trajectory_pairs.methods_extra", "methods_extra +/- tokens are syntactically valid")
    }

    missing_labels <- active_trajectory_pairs[
      !nzchar(trim(active_trajectory_pairs$coarse_label_var)) |
        !nzchar(trim(active_trajectory_pairs$fine_label_var)),
      ,
      drop = FALSE
    ]
    if (nrow(missing_labels) > 0) {
      fail("tier2.trajectory_pairs.label_vars", sprintf("missing label vars: %s", paste(missing_labels$trajectory_id, collapse = ", ")))
    } else {
      pass("tier2.trajectory_pairs.label_vars", "coarse_label_var and fine_label_var are populated")
    }
  }
}

if (!is.null(loaded_generated$annotation_marker_targets)) {
  annotation_marker_targets <- loaded_generated$annotation_marker_targets
  annotation_marker_schema <- c(
    "target_id", "source_question_id", "layer_scope", "object_layer",
    "cluster_column", "annotation_label_column", "group_var",
    "ident_1", "ident_2", "analysis_mode", "gene_program_role",
    "output_dir", "annotation_only", "enabled", "notes"
  )
  if (require_generated_cols("annotation_marker_targets.tsv", annotation_marker_targets, annotation_marker_schema) && nrow(annotation_marker_targets) > 0) {
    bad_mode <- annotation_marker_targets[annotation_marker_targets$analysis_mode != "annotation_cluster_marker", , drop = FALSE]
    bad_role <- annotation_marker_targets[annotation_marker_targets$gene_program_role != "annotation_marker", , drop = FALSE]
    bad_flag <- annotation_marker_targets[tolower(annotation_marker_targets$annotation_only) != "yes", , drop = FALSE]
    if (nrow(bad_mode) > 0) {
      fail("tier2.annotation_marker_targets.analysis_mode", sprintf("annotation targets must use annotation_cluster_marker: %s", paste(bad_mode$target_id, collapse = ", ")))
    }
    if (nrow(bad_role) > 0) {
      fail("tier2.annotation_marker_targets.gene_program_role", sprintf("annotation targets must use annotation_marker: %s", paste(bad_role$target_id, collapse = ", ")))
    }
    if (nrow(bad_flag) > 0) {
      fail("tier2.annotation_marker_targets.annotation_only", sprintf("annotation targets must set annotation_only=yes: %s", paste(bad_flag$target_id, collapse = ", ")))
    }
    if (nrow(bad_mode) == 0 && nrow(bad_role) == 0 && nrow(bad_flag) == 0) {
      pass("tier2.annotation_marker_targets.contract", "annotation marker targets are annotation-only 03/04 inputs")
    }
  }
}

if (!is.null(loaded_generated$comparisons)) {
  comparisons <- loaded_generated$comparisons
  comparison_schema <- c(
    "comparison_id", "source_question_id", "layer_scope", "contrast_axis",
    "display_question_id", "output_alias", "report_title",
    "analysis_mode", "analysis_unit", "stat_level", "group_var",
    "ident_1", "ident_2", "subset_column", "subset_value",
    "aggregation_group_var", "composition_group_var", "batch_var",
    "enabled", "min_biological_replicates", "force_exploratory",
    "min_cells_per_group", "logfc_threshold", "produces_gene_program",
    "gene_program_role", "notes"
  )
  if (require_generated_cols("comparisons.tsv", comparisons, comparison_schema) && nrow(comparisons) > 0) {
    bad_modes <- unique(comparisons$analysis_mode[!comparisons$analysis_mode %in% allowed_analysis_modes])
    bad_modes <- bad_modes[nzchar(bad_modes)]
    if (length(bad_modes) > 0) {
      fail("tier2.comparisons.analysis_mode", sprintf("unsupported analysis_mode values: %s", paste(bad_modes, collapse = ", ")))
    } else {
      pass("tier2.comparisons.analysis_mode", "analysis_mode values are valid")
    }

    bad_roles <- unique(comparisons$gene_program_role[!comparisons$gene_program_role %in% allowed_gene_program_roles])
    bad_roles <- bad_roles[nzchar(bad_roles)]
    if (length(bad_roles) > 0) {
      fail("tier2.comparisons.gene_program_role", sprintf("unsupported gene_program_role values: %s", paste(bad_roles, collapse = ", ")))
    } else {
      pass("tier2.comparisons.gene_program_role", "gene_program_role values are valid")
    }

    rest_bad <- comparisons[
      (comparisons$ident_1 == "__rest__" | comparisons$ident_2 == "__rest__") &
        !comparisons$analysis_mode %in% c("annotation_cluster_marker", "subtype_marker"),
      ,
      drop = FALSE
    ]
    if (nrow(rest_bad) > 0) {
      fail(
        "tier2.comparisons.rest_scope",
        sprintf("__rest__ used outside marker modes: %s", paste(rest_bad$comparison_id, collapse = ", ")),
        fix = "__rest__ is only valid for annotation_cluster_marker or subtype_marker."
      )
    } else {
      pass("tier2.comparisons.rest_scope", "__rest__ only appears in marker modes")
    }

    condition_bad <- comparisons[
      comparisons$analysis_mode == "condition_within_type" &
        ((!nzchar(comparisons$subset_column) & comparisons$aggregation_group_var != "all_cells") | comparisons$group_var != "group_id" |
           !comparisons$ident_1 %in% c("syf", "f5") | !comparisons$ident_2 %in% c("syf", "f5")),
      ,
      drop = FALSE
    ]
    if (nrow(condition_bad) > 0) {
      fail(
        "tier2.comparisons.condition_mode",
        sprintf("condition_within_type rows have invalid group/subset/idents: %s", paste(condition_bad$comparison_id, collapse = ", ")),
        fix = "Use group_var=group_id, syf/f5 idents, and a non-empty subset_column unless the row is an explicit all-cell condition DEG."
      )
    } else {
      pass("tier2.comparisons.condition_mode", "condition_within_type rows have syf/f5 group_id semantics")
    }

    composition_bad <- comparisons[
      comparisons$analysis_mode %in% c("composition", "qc_composition") &
        (!nzchar(comparisons$composition_group_var) | comparisons$produces_gene_program == "yes"),
      ,
      drop = FALSE
    ]
    if (nrow(composition_bad) > 0) {
      fail(
        "tier2.comparisons.composition_mode",
        sprintf("composition/qc_composition rows missing composition_group_var or producing gene programs: %s", paste(composition_bad$comparison_id, collapse = ", ")),
        fix = "composition and qc_composition rows must set composition_group_var and produces_gene_program=no."
      )
    } else {
      pass("tier2.comparisons.composition_mode", "composition and qc_composition rows are sample-level targets, not gene programs")
    }

    qc_bad <- comparisons[
      comparisons$analysis_mode == "qc_composition" &
        (comparisons$gene_program_role != "qc_only" | comparisons$produces_gene_program == "yes"),
      ,
      drop = FALSE
    ]
    if (nrow(qc_bad) > 0) {
      fail("tier2.comparisons.qc_composition", sprintf("qc_composition rows must be qc_only and produce no gene program: %s", paste(qc_bad$comparison_id, collapse = ", ")))
    } else {
      pass("tier2.comparisons.qc_composition", "qc_composition rows are qc_only")
    }

    e03_cmp <- comparisons[comparisons$source_question_id == "E03_layer_compo", , drop = FALSE]
    if (nrow(e03_cmp) != 1L) {
      fail("tier2.comparisons.E03_qc_alias", sprintf("E03 must fan out to exactly one qc_composition row, found %s", nrow(e03_cmp)))
    } else {
      e03_cmp_bad <- e03_cmp$analysis_mode[[1]] != "qc_composition" ||
        e03_cmp$gene_program_role[[1]] != "qc_only" ||
        e03_cmp$produces_gene_program[[1]] != "no" ||
        e03_cmp$display_question_id[[1]] != "E03_scRNA_GC_TC_capture_balance" ||
        e03_cmp$output_alias[[1]] != "E03_scRNA_GC_TC_capture_balance" ||
        !grepl("captured-cell balance QC", e03_cmp$report_title[[1]], fixed = TRUE)
      if (e03_cmp_bad) {
        fail("tier2.comparisons.E03_qc_alias", "E03 comparison row must carry qc_only semantics and captured-cell balance alias/title")
      } else {
        pass("tier2.comparisons.E03_qc_alias", "E03 comparison row carries qc_only alias/title metadata")
      }
    }

    global_bad <- comparisons[
      comparisons$analysis_mode == "global_context" & comparisons$gene_program_role != "global_context",
      ,
      drop = FALSE
    ]
    if (nrow(global_bad) > 0) {
      fail("tier2.comparisons.global_context", sprintf("global_context rows must use gene_program_role=global_context: %s", paste(global_bad$comparison_id, collapse = ", ")))
    } else {
      pass("tier2.comparisons.global_context", "global_context rows use global_context role")
    }

    if (exists("scope_meta")) {
      for (idx in seq_len(nrow(comparisons))) {
        row <- comparisons[idx, , drop = FALSE]
        scopes <- split_comm_tokens(row$layer_scope[[1]])
        if (length(scopes) == 0) scopes <- row$layer_scope[[1]]
        for (scope in scopes) {
          item <- scope_meta[[scope]]
          if (is.null(item) || is.null(item$meta)) next
          meta <- item$meta
          group_var <- row$group_var[[1]]
          if (group_var %in% c("cell_type", "cell_subtype")) {
            if (!group_var %in% colnames(meta)) {
              semantic_object_problem(
                "semantic.comparison_group_var_column",
                sprintf("%s uses group_var=%s but %s metadata lacks that column", row$comparison_id[[1]], group_var, scope)
              )
              next
            }
            observed <- unique(trimws(as.character(meta[[group_var]])))
            observed <- observed[nzchar(observed)]
            tokens <- setdiff(c(row$ident_1[[1]], row$ident_2[[1]]), "__rest__")
            tokens <- tokens[nzchar(tokens)]
            missing_tokens <- setdiff(tokens, observed)
            if (length(missing_tokens) > 0) {
              semantic_object_problem(
                "semantic.comparison_group_var_values",
                sprintf("%s group_var=%s values absent in %s metadata: %s", row$comparison_id[[1]], group_var, scope, paste(missing_tokens, collapse = ", "))
              )
            }
          }
        }
      }
    }
  }
}

if (!is.null(loaded_generated$gene_program_targets)) {
  gene_program_targets <- loaded_generated$gene_program_targets
  gene_program_schema <- c(
    "comparison_id", "source_question_id", "layer_scope", "analysis_mode",
    "gene_program_role", "produces_gene_program", "annotation_only", "qc_only",
    "global_context_only", "nichenet_eligible", "nichenet_usage",
    "enrichment_eligible", "enrichment_usage", "preferred_for_downstream",
    "expected_result_level", "formal_preferred", "formal_status",
    "result_status", "skip_reason", "eligible_reason", "ineligible_reason",
    "notes"
  )
  if (require_generated_cols("gene_program_targets.tsv", gene_program_targets, gene_program_schema) && nrow(gene_program_targets) > 0) {
    gp_fail_n <- 0L
    bad_modes <- unique(gene_program_targets$analysis_mode[!gene_program_targets$analysis_mode %in% allowed_analysis_modes])
    bad_modes <- bad_modes[nzchar(bad_modes)]
    if (length(bad_modes) > 0) {
      fail("tier2.gene_program_targets.analysis_mode", sprintf("unsupported analysis_mode values: %s", paste(bad_modes, collapse = ", ")))
      gp_fail_n <- gp_fail_n + length(bad_modes)
    }
    bad_roles <- unique(gene_program_targets$gene_program_role[!gene_program_targets$gene_program_role %in% allowed_gene_program_roles])
    bad_roles <- bad_roles[nzchar(bad_roles)]
    if (length(bad_roles) > 0) {
      fail("tier2.gene_program_targets.gene_program_role", sprintf("unsupported gene_program_role values: %s", paste(bad_roles, collapse = ", ")))
      gp_fail_n <- gp_fail_n + length(bad_roles)
    }
    bad_nichenet <- unique(gene_program_targets$nichenet_eligible[!gene_program_targets$nichenet_eligible %in% c("yes", "no")])
    bad_nichenet <- bad_nichenet[nzchar(bad_nichenet)]
    if (length(bad_nichenet) > 0) {
      fail("tier2.gene_program_targets.nichenet_eligible", sprintf("nichenet_eligible must be yes/no: %s", paste(bad_nichenet, collapse = ", ")))
      gp_fail_n <- gp_fail_n + length(bad_nichenet)
    }
    bad_enrichment <- unique(gene_program_targets$enrichment_eligible[!gene_program_targets$enrichment_eligible %in% allowed_yes_no_contextual])
    bad_enrichment <- bad_enrichment[nzchar(bad_enrichment)]
    if (length(bad_enrichment) > 0) {
      fail("tier2.gene_program_targets.enrichment_eligible", sprintf("enrichment_eligible must be yes/no/contextual: %s", paste(bad_enrichment, collapse = ", ")))
      gp_fail_n <- gp_fail_n + length(bad_enrichment)
    }
    bad_usage <- gene_program_targets[
      !gene_program_targets$nichenet_usage %in% allowed_nichenet_usage |
        !gene_program_targets$enrichment_usage %in% allowed_enrichment_usage,
      ,
      drop = FALSE
    ]
    if (nrow(bad_usage) > 0) {
      fail("tier2.gene_program_targets.usage", sprintf("unsupported nichenet/enrichment usage in: %s", paste(bad_usage$comparison_id, collapse = ", ")))
      gp_fail_n <- gp_fail_n + nrow(bad_usage)
    }
    annotation_bad <- gene_program_targets[
      gene_program_targets$analysis_mode == "annotation_cluster_marker" &
        (gene_program_targets$gene_program_role != "annotation_marker" |
           gene_program_targets$annotation_only != "yes" |
           gene_program_targets$nichenet_eligible != "no" |
           gene_program_targets$enrichment_eligible != "no"),
      ,
      drop = FALSE
    ]
    if (nrow(annotation_bad) > 0) {
      fail("tier2.gene_program_targets.annotation_marker", sprintf("annotation_cluster_marker must be annotation_only and ineligible for NicheNet/mechanism enrichment: %s", paste(annotation_bad$comparison_id, collapse = ", ")))
      gp_fail_n <- gp_fail_n + nrow(annotation_bad)
    }
    qc_bad <- gene_program_targets[
      gene_program_targets$analysis_mode == "qc_composition" &
        (gene_program_targets$produces_gene_program == "yes" |
           gene_program_targets$gene_program_role != "qc_only" |
           gene_program_targets$qc_only != "yes" |
           gene_program_targets$nichenet_eligible != "no" |
           gene_program_targets$enrichment_eligible != "no"),
      ,
      drop = FALSE
    ]
    if (nrow(qc_bad) > 0) {
      fail("tier2.gene_program_targets.qc_composition", sprintf("qc_composition must be qc_only, no gene program, no NicheNet/enrichment: %s", paste(qc_bad$comparison_id, collapse = ", ")))
      gp_fail_n <- gp_fail_n + nrow(qc_bad)
    }
    global_bad <- gene_program_targets[
      (gene_program_targets$analysis_mode == "global_context" |
         gene_program_targets$gene_program_role == "global_context") &
        (gene_program_targets$analysis_mode != "global_context" |
           gene_program_targets$gene_program_role != "global_context" |
           gene_program_targets$global_context_only != "yes" |
           gene_program_targets$nichenet_eligible != "no" |
           gene_program_targets$nichenet_usage != "none" |
           gene_program_targets$preferred_for_downstream != "contextual" |
           !(gene_program_targets$enrichment_eligible %in% c("contextual", "no"))),
      ,
      drop = FALSE
    ]
    if (nrow(global_bad) > 0) {
      fail("tier2.gene_program_targets.global_context", sprintf("global_context rows must be contextual-only and NicheNet-ineligible: %s", paste(global_bad$comparison_id, collapse = ", ")))
      gp_fail_n <- gp_fail_n + nrow(global_bad)
    }
    e03_bad <- gene_program_targets[
      gene_program_targets$source_question_id == "E03_layer_compo" &
        (gene_program_targets$analysis_mode != "qc_composition" |
           gene_program_targets$gene_program_role != "qc_only" |
           gene_program_targets$produces_gene_program != "no" |
           gene_program_targets$qc_only != "yes" |
           gene_program_targets$nichenet_eligible != "no" |
           gene_program_targets$enrichment_eligible != "no" |
           gene_program_targets$preferred_for_downstream != "no" |
           gene_program_targets$result_status != "qc_only" |
           gene_program_targets$skip_reason != "qc_composition_does_not_produce_gene_program" |
           gene_program_targets$ineligible_reason != "qc_composition_does_not_produce_gene_program"),
      ,
      drop = FALSE
    ]
    if (nrow(e03_bad) > 0) {
      fail("tier2.gene_program_targets.E03_qc_only", "E03 must be qc_only, not produce gene programs, and remain ineligible for 06/07")
      gp_fail_n <- gp_fail_n + nrow(e03_bad)
    }
    d05_bad <- gene_program_targets[
      gene_program_targets$source_question_id == "D05_panorama_global_stage" &
        (gene_program_targets$analysis_mode != "global_context" |
           gene_program_targets$gene_program_role != "global_context" |
           gene_program_targets$global_context_only != "yes" |
           gene_program_targets$nichenet_eligible != "no" |
           gene_program_targets$nichenet_usage != "none" |
           gene_program_targets$enrichment_eligible != "contextual" |
           gene_program_targets$enrichment_usage != "global_context_enrichment" |
           gene_program_targets$preferred_for_downstream != "contextual"),
      ,
      drop = FALSE
    ]
    if (nrow(d05_bad) > 0) {
      fail("tier2.gene_program_targets.D05_global_context", "D05 must be contextual global_context and NicheNet-ineligible")
      gp_fail_n <- gp_fail_n + nrow(d05_bad)
    }
    if (gp_fail_n == 0L) {
      pass("tier2.gene_program_targets.contract", "gene program role and eligibility matrix is valid")
    }
  }
}

if (!is.null(loaded_generated$enrichment_targets)) {
  enrichment_targets <- loaded_generated$enrichment_targets
  if (nrow(enrichment_targets) > 0 && "source_question_id" %in% colnames(enrichment_targets)) {
    e03_enrichment <- enrichment_targets[enrichment_targets$source_question_id == "E03_layer_compo", , drop = FALSE]
    if (nrow(e03_enrichment) > 0) {
      fail("tier2.enrichment_targets.E03_absent", "E03 qc_composition must not enter enrichment_targets.tsv")
    } else {
      pass("tier2.enrichment_targets.E03_absent", "E03 qc_composition is absent from enrichment targets")
    }
    d05_enrichment <- enrichment_targets[enrichment_targets$source_question_id == "D05_panorama_global_stage", , drop = FALSE]
    if (nrow(d05_enrichment) == 0) {
      fail("tier2.enrichment_targets.D05_contextual", "D05 must have contextual enrichment targets")
    } else if (any(d05_enrichment$analysis_mode != "global_context" |
        d05_enrichment$gene_program_role != "global_context" |
        d05_enrichment$enrichment_eligible != "contextual" |
        d05_enrichment$enrichment_usage != "global_context_enrichment")) {
      fail("tier2.enrichment_targets.D05_contextual", "D05 enrichment targets must be global_context_enrichment and contextual")
    } else {
      pass("tier2.enrichment_targets.D05_contextual", "D05 enrichment targets are contextual global context targets")
    }
    if (!is.null(loaded_generated$gene_program_targets)) {
      gp_targets <- loaded_generated$gene_program_targets
      for (col in c("comparison_id", "layer_scope", "produces_gene_program", "enrichment_eligible")) {
        if (!col %in% colnames(gp_targets)) gp_targets[[col]] <- character(nrow(gp_targets))
      }
      for (col in c("target_id", "comparison_id", "layer_scope", "enrichment_eligible", "enabled")) {
        if (!col %in% colnames(enrichment_targets)) enrichment_targets[[col]] <- character(nrow(enrichment_targets))
      }
      enabled_enrichment <- enrichment_targets[
        tolower(trim(enrichment_targets$enabled)) %in% c("", "yes", "true", "1", "on"),
        ,
        drop = FALSE
      ]
      enrichment_fail_n <- 0L
      for (idx in seq_len(nrow(enabled_enrichment))) {
        row <- enabled_enrichment[idx, , drop = FALSE]
        hit <- gp_targets[
          gp_targets$comparison_id == row$comparison_id[[1]] &
            gp_targets$layer_scope == row$layer_scope[[1]],
          ,
          drop = FALSE
        ]
        if (nrow(hit) == 0) {
          fail("tier2.enrichment_targets.registry_join", sprintf("%s comparison_id/layer_scope is absent from gene_program_targets.tsv: %s/%s", row$target_id[[1]], row$comparison_id[[1]], row$layer_scope[[1]]))
          enrichment_fail_n <- enrichment_fail_n + 1L
          next
        }
        eligible <- trim(row$enrichment_eligible[[1]])
        if (eligible %in% c("yes", "contextual") && !identical(hit$produces_gene_program[[1]], "yes")) {
          fail("tier2.enrichment_targets.registry_join", sprintf("%s is enrichment eligible but registry target does not produce a gene program: %s", row$target_id[[1]], row$comparison_id[[1]]))
          enrichment_fail_n <- enrichment_fail_n + 1L
        }
      }
      if (enrichment_fail_n == 0L) {
        pass("tier2.enrichment_targets.registry_join", sprintf("%s enabled enrichment targets join gene_program_targets.tsv", nrow(enabled_enrichment)))
      }
    }
  }
}

if (!is.null(loaded_generated$communication_pairs)) {
  communication_pairs <- loaded_generated$communication_pairs
  communication_schema <- c(
    "pair_id", "source_question_id", "layer_scope", "sender", "receiver",
    "condition_split_var", "condition_split_values", "tool", "communication_mode",
    "activation_policy", "min_sender_cells", "min_receiver_cells",
    "min_cells_per_condition", "fallback_pair_id", "derived_from_pair_id",
    "run_baseline_if_split_fails", "requires_all_derived_inputs_pass",
    "receiver_gene_program_source", "baseline_marker_comparison_id",
    "receiver_deg_comparison_id", "direction_filter", "requires_cell_subtype",
    "notes", "enabled"
  )
  if (require_generated_cols("communication_pairs.tsv", communication_pairs, communication_schema) && nrow(communication_pairs) > 0) {
    communication_pairs$activation_policy <- tolower(vapply(communication_pairs$activation_policy, trim, character(1)))
    communication_pairs$activation_policy[!nzchar(communication_pairs$activation_policy)] <- "always"
    communication_pair_ids <- unique(communication_pairs$pair_id)
    policy_fail_n <- 0L
    bad_policy <- unique(communication_pairs$activation_policy[!communication_pairs$activation_policy %in% allowed_activation_policy])
    bad_policy <- bad_policy[nzchar(bad_policy)]
    if (length(bad_policy) > 0) {
      fail("tier2.communication.activation_policy", sprintf("unsupported activation_policy values: %s", paste(bad_policy, collapse = ", ")))
      policy_fail_n <- policy_fail_n + length(bad_policy)
    }
    auto_pairs <- communication_pairs[communication_pairs$activation_policy == "auto_if_min_cells", , drop = FALSE]
    if (nrow(auto_pairs) > 0) {
      for (idx in seq_len(nrow(auto_pairs))) {
        row <- auto_pairs[idx, , drop = FALSE]
        pair_id <- row$pair_id[[1]]
        min_values <- c(
          min_sender_cells = positive_int_value(row$min_sender_cells[[1]]),
          min_receiver_cells = positive_int_value(row$min_receiver_cells[[1]]),
          min_cells_per_condition = positive_int_value(row$min_cells_per_condition[[1]])
        )
        if (any(is.na(min_values))) {
          fail("tier2.communication.auto_min_cells", sprintf("%s auto_if_min_cells requires positive integer min_* fields", pair_id))
          policy_fail_n <- policy_fail_n + 1L
        }
        run_baseline <- tolower(trim(row$run_baseline_if_split_fails[[1]])) %in% c("yes", "true", "1", "on")
        fallback_id <- trim(row$fallback_pair_id[[1]])
        if (run_baseline && !nzchar(fallback_id)) {
          fail("tier2.communication.fallback_pair_id", sprintf("%s run_baseline_if_split_fails=yes requires fallback_pair_id", pair_id))
          policy_fail_n <- policy_fail_n + 1L
        }
        if (nzchar(fallback_id)) {
          fallback <- communication_pairs[communication_pairs$pair_id == fallback_id, , drop = FALSE]
          if (nrow(fallback) == 0) {
            fail("tier2.communication.fallback_pair_id", sprintf("%s fallback_pair_id not found: %s", pair_id, fallback_id))
            policy_fail_n <- policy_fail_n + 1L
          } else if (!identical(fallback$activation_policy[[1]], "always") || !identical(fallback$communication_mode[[1]], "baseline")) {
            fail("tier2.communication.fallback_pair_id", sprintf("%s fallback_pair_id=%s must point to activation_policy=always baseline", pair_id, fallback_id))
            policy_fail_n <- policy_fail_n + 1L
          }
        }
      }
    }
    derived_pairs <- communication_pairs[communication_pairs$activation_policy == "derived_from_split", , drop = FALSE]
    if (nrow(derived_pairs) > 0) {
      executable_derived <- derived_pairs[tolower(derived_pairs$tool) %in% c("both", "cellchat", "cellchat_only", "nichenet", "nichenet_only"), , drop = FALSE]
      if (nrow(executable_derived) > 0) {
        fail("tier2.communication.derived_not_executable", sprintf("derived rows cannot be executable by 07a/07c: %s", paste(executable_derived$pair_id, collapse = ", ")))
        policy_fail_n <- policy_fail_n + nrow(executable_derived)
      }
      for (idx in seq_len(nrow(derived_pairs))) {
        row <- derived_pairs[idx, , drop = FALSE]
        source_ids <- split_comm_tokens(row$derived_from_pair_id[[1]])
        missing_source <- source_ids[!source_ids %in% communication_pair_ids]
        if (length(missing_source) > 0) {
          fail("tier2.communication.derived_from_pair_id", sprintf("%s derived_from_pair_id not found: %s", row$pair_id[[1]], paste(missing_source, collapse = ", ")))
          policy_fail_n <- policy_fail_n + length(missing_source)
        }
      }
    }
    f08_rows <- communication_pairs[
      communication_pairs$source_question_id == "F08_GC_dev_seq_split" &
        startsWith(communication_pairs$pair_id, "F08_GC_dev_seq_split__"),
      ,
      drop = FALSE
    ]
    if (nrow(f08_rows) > 0 && any(f08_rows$activation_policy != "auto_if_min_cells")) {
      fail("tier2.communication.A1_F08_policy", "F08 expanded rows must use activation_policy=auto_if_min_cells")
      policy_fail_n <- policy_fail_n + 1L
    }
    f17_rows <- communication_pairs[communication_pairs$source_question_id == "F17_GC_screen_split", , drop = FALSE]
    if (nrow(f17_rows) > 0 && any(f17_rows$activation_policy != "auto_if_min_cells")) {
      fail("tier2.communication.A1_F17_policy", "F17 expanded rows must use activation_policy=auto_if_min_cells")
      policy_fail_n <- policy_fail_n + 1L
    }
    f25 <- communication_pairs[communication_pairs$pair_id == "F25_GC_internal_diff", , drop = FALSE]
    if (nrow(f25) == 1L) {
      f25_sources <- split_comm_tokens(f25$derived_from_pair_id[[1]])
      if (length(f25_sources) == 0 || any(!startsWith(f25_sources, "F08_GC_dev_seq_split__"))) {
        fail("tier2.communication.F25_sources", "F25_GC_internal_diff must derive from expanded F08 split rows")
        policy_fail_n <- policy_fail_n + 1L
      }
      if (!tolower(trim(f25$requires_all_derived_inputs_pass[[1]])) %in% c("yes", "true", "1", "on")) {
        fail("tier2.communication.F25_requires_all", "F25_GC_internal_diff must require all derived inputs to pass")
        policy_fail_n <- policy_fail_n + 1L
      }
    }
    if (policy_fail_n == 0L) {
      pass("tier2.communication.activation_policy", "communication activation policies, fallbacks, and derived inputs are valid")
    }

    gene_targets <- loaded_generated$gene_program_targets
    if (!is.null(gene_targets)) {
      for (col in c("comparison_id", "analysis_mode", "gene_program_role", "nichenet_usage")) {
        if (!col %in% colnames(gene_targets)) {
          gene_targets[[col]] <- character(nrow(gene_targets))
        }
      }
    }
    gene_target_ids <- if (!is.null(gene_targets) && "comparison_id" %in% colnames(gene_targets)) unique(gene_targets$comparison_id) else character(0)
    comparison_ids <- if (exists("comparisons") && "comparison_id" %in% colnames(comparisons)) unique(comparisons$comparison_id) else character(0)
    nichenet_requested <- tolower(communication_pairs$tool) %in% c("both", "nichenet", "nichenet_only")
    nichenet_pairs <- communication_pairs[nichenet_requested & communication_pairs$enabled != "no", , drop = FALSE]
    comm_fail_n <- 0L
    e03_comparison_ids <- if (exists("comparisons") && all(c("comparison_id", "source_question_id") %in% colnames(comparisons))) {
      comparisons$comparison_id[comparisons$source_question_id == "E03_layer_compo"]
    } else {
      character(0)
    }
    if (length(e03_comparison_ids) > 0) {
      e03_comm_refs <- communication_pairs[
        communication_pairs$baseline_marker_comparison_id %in% e03_comparison_ids |
          communication_pairs$receiver_deg_comparison_id %in% e03_comparison_ids,
        ,
        drop = FALSE
      ]
      if (nrow(e03_comm_refs) > 0) {
        fail("tier2.communication.E03_absent", sprintf("E03 qc_composition must not be referenced by 07 communication rows: %s", paste(e03_comm_refs$pair_id, collapse = ", ")))
        comm_fail_n <- comm_fail_n + nrow(e03_comm_refs)
      }
    }
    global_context_ids <- if (!is.null(gene_targets) && all(c("comparison_id", "analysis_mode", "gene_program_role") %in% colnames(gene_targets))) {
      gene_targets$comparison_id[
        gene_targets$analysis_mode == "global_context" |
          gene_targets$gene_program_role == "global_context"
      ]
    } else {
      character(0)
    }
    if (length(global_context_ids) > 0) {
      global_receiver_ref <- vapply(communication_pairs$receiver_deg_comparison_id, function(x) {
        any(split_comm_tokens(x) %in% global_context_ids)
      }, logical(1))
      global_comm_refs <- communication_pairs[global_receiver_ref, , drop = FALSE]
      if (nrow(global_comm_refs) > 0) {
        fail("tier2.communication.global_context_receiver_deg", sprintf("receiver_deg_comparison_id must not point to global_context rows: %s", paste(global_comm_refs$pair_id, collapse = ", ")))
        comm_fail_n <- comm_fail_n + nrow(global_comm_refs)
      }
    }
    for (idx in seq_len(nrow(nichenet_pairs))) {
      row <- nichenet_pairs[idx, , drop = FALSE]
      pair_id <- row$pair_id[[1]]
      source <- tolower(row$receiver_gene_program_source[[1]])
      if (!source %in% c("condition_deg", "receiver_marker", "none")) {
        fail("tier2.communication.gene_program_source", sprintf("%s has unsupported receiver_gene_program_source=%s", pair_id, source))
        comm_fail_n <- comm_fail_n + 1L
      }
      if (identical(source, "none")) {
        fail("tier2.communication.gene_program_source", sprintf("%s is NicheNet-capable but receiver_gene_program_source=none", pair_id))
        comm_fail_n <- comm_fail_n + 1L
      }
      mode <- tolower(trim(row$communication_mode[[1]]))
      split_requested <- nzchar(trim(row$condition_split_var[[1]])) || grepl("split", mode, fixed = TRUE)
      if (split_requested && !identical(source, "condition_deg")) {
        fail("tier2.communication.gene_program_source", sprintf("%s split NicheNet rows must use receiver_gene_program_source=condition_deg", pair_id))
        comm_fail_n <- comm_fail_n + 1L
      }
      if (!split_requested && identical(mode, "baseline") && !identical(source, "receiver_marker")) {
        fail("tier2.communication.gene_program_source", sprintf("%s baseline NicheNet rows must use receiver_gene_program_source=receiver_marker", pair_id))
        comm_fail_n <- comm_fail_n + 1L
      }
      baseline_id <- trim(row$baseline_marker_comparison_id[[1]])
      receiver_id <- trim(row$receiver_deg_comparison_id[[1]])
      if (source != "none" && !nzchar(baseline_id)) {
        fail("tier2.communication.baseline_marker_id", sprintf("%s missing baseline_marker_comparison_id", pair_id))
        comm_fail_n <- comm_fail_n + 1L
      }
      baseline_ids <- split_comm_tokens(baseline_id)
      missing_baseline <- baseline_ids[!baseline_ids %in% comparison_ids | !baseline_ids %in% gene_target_ids]
      if (length(missing_baseline) > 0) {
        fail("tier2.communication.baseline_marker_id", sprintf("%s baseline_marker_comparison_id not found in comparisons/gene_program_targets: %s", pair_id, paste(missing_baseline, collapse = ", ")))
        comm_fail_n <- comm_fail_n + length(missing_baseline)
      }
      if (length(baseline_ids) > 0 && !is.null(gene_targets) && nrow(gene_targets) > 0) {
        baseline_rows <- gene_targets[gene_targets$comparison_id %in% baseline_ids, , drop = FALSE]
        bad_baseline <- baseline_rows[
          baseline_rows$analysis_mode != "subtype_marker" |
            baseline_rows$gene_program_role != "receiver_marker" |
            baseline_rows$nichenet_usage != "baseline_receiver_marker",
          ,
          drop = FALSE
        ]
        if (nrow(bad_baseline) > 0) {
          fail("tier2.communication.baseline_marker_role", sprintf("%s baseline_marker_comparison_id must point to subtype_marker + receiver_marker baseline usage: %s", pair_id, paste(bad_baseline$comparison_id, collapse = ", ")))
          comm_fail_n <- comm_fail_n + nrow(bad_baseline)
        }
      }
      if (identical(source, "condition_deg") && !nzchar(receiver_id)) {
        fail("tier2.communication.receiver_deg_id", sprintf("%s source=condition_deg but receiver_deg_comparison_id is empty", pair_id))
        comm_fail_n <- comm_fail_n + 1L
      }
      receiver_ids <- split_comm_tokens(receiver_id)
      missing_receiver <- receiver_ids[!receiver_ids %in% comparison_ids | !receiver_ids %in% gene_target_ids]
      if (length(missing_receiver) > 0) {
        fail("tier2.communication.receiver_deg_id", sprintf("%s receiver_deg_comparison_id not found in comparisons/gene_program_targets: %s", pair_id, paste(missing_receiver, collapse = ", ")))
        comm_fail_n <- comm_fail_n + length(missing_receiver)
      }
      if (length(receiver_ids) > 0 && !is.null(gene_targets) && nrow(gene_targets) > 0) {
        receiver_rows <- gene_targets[gene_targets$comparison_id %in% receiver_ids, , drop = FALSE]
        bad_receiver <- receiver_rows[
          receiver_rows$analysis_mode != "condition_within_type" |
            receiver_rows$gene_program_role != "condition_deg" |
            receiver_rows$nichenet_usage != "receiver_condition_deg",
          ,
          drop = FALSE
        ]
        if (nrow(bad_receiver) > 0) {
          fail("tier2.communication.receiver_deg_role", sprintf("%s receiver_deg_comparison_id must point to condition_within_type + condition_deg receiver usage: %s", pair_id, paste(bad_receiver$comparison_id, collapse = ", ")))
          comm_fail_n <- comm_fail_n + nrow(bad_receiver)
        }
      }
    }
    if (comm_fail_n == 0L) {
      pass("tier2.communication.nichenet_gene_program_ids", sprintf("%s NicheNet-capable communication rows have explicit gene-program IDs", nrow(nichenet_pairs)))
    }

    strict_pairs <- communication_pairs[tolower(communication_pairs$requires_cell_subtype) %in% c("yes", "true", "1", "on"), , drop = FALSE]
    if (nrow(strict_pairs) > 0 && exists("scope_meta")) {
      for (idx in seq_len(nrow(strict_pairs))) {
        row <- strict_pairs[idx, , drop = FALSE]
        scopes <- split_comm_tokens(row$layer_scope[[1]])
        if (length(scopes) == 0) scopes <- row$layer_scope[[1]]
        tokens <- unique(c(split_comm_tokens(row$sender[[1]]), split_comm_tokens(row$receiver[[1]])))
        for (scope in scopes) {
          item <- scope_meta[[scope]]
          if (is.null(item) || is.null(item$meta)) next
          meta <- item$meta
          if (!"cell_subtype" %in% colnames(meta)) {
            fail("semantic.communication_requires_cell_subtype", sprintf("%s requires_cell_subtype=yes but %s metadata lacks cell_subtype", row$pair_id[[1]], scope))
            next
          }
          observed <- unique(trimws(as.character(meta$cell_subtype)))
          observed <- observed[nzchar(observed)]
          missing_tokens <- setdiff(tokens, observed)
          if (length(missing_tokens) > 0) {
            fail("semantic.communication_requires_cell_subtype", sprintf("%s requires_cell_subtype=yes but strict roles are absent from %s cell_subtype: %s", row$pair_id[[1]], scope, paste(missing_tokens, collapse = ", ")))
          }
        }
      }
    }
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
