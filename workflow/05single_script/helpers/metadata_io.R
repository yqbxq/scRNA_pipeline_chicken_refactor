normalize_flag <- function(x, default = "yes") {
  if (length(x) == 0) {
    return(character(0))
  }
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x[x == ""] <- default
  tolower(x)
}

normalize_scalar_value <- function(x, default = "") {
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(default)
  }
  trimws(as.character(x))
}

display_scalar_value <- function(x, default = "NA") {
  value <- normalize_scalar_value(x)
  if (nzchar(value)) {
    return(value)
  }
  default
}

collapse_unique_values <- function(x, sep = ", ") {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return("")
  }
  paste(unique(x), collapse = sep)
}

read_tsv_optional <- function(path) {
  if (!nzchar(path) || !file.exists(path) || isTRUE(file.info(path)$size == 0)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  lines <- sub("^\ufeff", "", lines)
  trimmed <- trimws(lines)
  keep <- nzchar(trimmed) & !startsWith(trimmed, "#")
  if (!any(keep)) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  con <- textConnection(lines[keep])
  on.exit(close(con), add = TRUE)
  read.delim(
    con,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = "",
    fill = TRUE
  )
}

write_tsv_local <- function(df, path) {
  ensure_dir(dirname(path))
  write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

write_csv_local <- function(df, path) {
  ensure_dir(dirname(path))
  write.csv(df, path, row.names = FALSE)
}

read_identifier_list <- function(path) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path) || !file.exists(path) || isTRUE(file.info(path)$size == 0)) {
    return(character(0))
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) == 0) {
    return(character(0))
  }

  lines <- sub("^\ufeff", "", lines)
  trimmed <- trimws(lines)
  trimmed <- trimmed[nzchar(trimmed) & !startsWith(trimmed, "#")]
  if (length(trimmed) == 0) {
    return(character(0))
  }

  tokens <- unlist(strsplit(trimmed, "[,\t ]+", perl = TRUE), use.names = FALSE)
  tokens <- trimws(tokens)
  tokens <- tokens[nzchar(tokens)]
  tokens <- tokens[!tolower(tokens) %in% c("gene_id", "gene_name", "feature", "feature_name")]
  unique(tokens)
}

active_sample_sheet_path_local <- function(cfg) {
  if (nzchar(cfg$canonical_sample_sheet) && file.exists(cfg$canonical_sample_sheet)) {
    return(cfg$canonical_sample_sheet)
  }
  cfg$sample_sheet
}

read_sample_sheet_local <- function(cfg) {
  df <- read_tsv_optional(active_sample_sheet_path_local(cfg))
  if (nrow(df) == 0) {
    stop(sprintf("缺少样本表或样本表为空: %s", cfg$sample_sheet), call. = FALSE)
  }

  expected_cols <- c(
    "sample_id", "condition", "cell_type_broad", "biological_replicate", "technical_replicate",
    "batch", "input_mode", "input_source", "source_path",
    "platform", "gene_id_type", "reference_version", "group_id", "timepoint",
    "tissue", "chemistry", "run_main", "run_velocity", "run_scenic"
  )
  for (col in expected_cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  df
}

main_sample_sheet_local <- function(cfg) {
  df <- read_sample_sheet_local(cfg)
  df[normalize_flag(df$run_main, "yes") != "no", , drop = FALSE]
}

read_input_inventory_local <- function(cfg) {
  read_tsv_optional(cfg$input_inventory_file)
}

read_branch_readiness_local <- function(cfg) {
  read_tsv_optional(cfg$branch_readiness_file)
}

require_columns_local <- function(df, required_cols, table_name, path) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "%s 缺少必需字段: %s。请在 upstream intake 层适配字段 contract: %s",
        table_name,
        paste(missing_cols, collapse = ", "),
        path
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

require_input_inventory_contract_local <- function(cfg, inventory_df) {
  if (!file.exists(cfg$input_inventory_file) || nrow(inventory_df) == 0) {
    stop(sprintf("缺少正式 intake 产物 input_inventory.tsv: %s", cfg$input_inventory_file), call. = FALSE)
  }
  require_columns_local(
    inventory_df,
    c(
      "sample_id",
      "filtered_matrix_dir",
      "raw_matrix_dir",
      "metrics_path",
      "platform_resolved",
      "gene_id_type_resolved",
      "feature_name_profile",
      "reference_version"
    ),
    "input_inventory.tsv",
    cfg$input_inventory_file
  )
}

require_branch_readiness_contract_local <- function(cfg, branch_readiness_df) {
  if (!file.exists(cfg$branch_readiness_file) || nrow(branch_readiness_df) == 0) {
    stop(sprintf("缺少正式 intake 产物 branch_readiness.tsv: %s", cfg$branch_readiness_file), call. = FALSE)
  }
  require_columns_local(branch_readiness_df, c("sample_id", "notes"), "branch_readiness.tsv", cfg$branch_readiness_file)
  if (!any(branch_readiness_df$sample_id == "__PROJECT__")) {
    stop("branch_readiness.tsv 缺少 sample_id == '__PROJECT__' 的项目级 readiness 行。", call. = FALSE)
  }
}

read_qc_threshold_overrides_local <- function(cfg) {
  df <- read_tsv_optional(cfg$qc_threshold_file)
  if (nrow(df) == 0) {
    return(df)
  }

  required <- c("sample_id", "qc_min_nfeature", "qc_min_ncount", "qc_min_log10umi", "qc_max_mito_pct")
  for (col in required) {
    if (!col %in% colnames(df)) {
      df[[col]] <- NA_character_
    }
  }
  for (col in setdiff(required, "sample_id")) {
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }
  df
}

is_mex_matrix_dir_local <- function(path) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path) || !dir.exists(path)) {
    return(FALSE)
  }
  matrix_exists <- file.exists(file.path(path, "matrix.mtx.gz")) || file.exists(file.path(path, "matrix.mtx"))
  features_exists <- any(file.exists(file.path(path, c("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"))))
  barcodes_exists <- file.exists(file.path(path, "barcodes.tsv.gz")) || file.exists(file.path(path, "barcodes.tsv"))
  matrix_exists && features_exists && barcodes_exists
}

resolve_inventory_row_local <- function(cfg, sample_id, inventory_df = NULL) {
  if (is.null(inventory_df)) {
    inventory_df <- read_input_inventory_local(cfg)
  }
  if (nrow(inventory_df) == 0 || !"sample_id" %in% colnames(inventory_df)) {
    return(NULL)
  }
  hit <- inventory_df[inventory_df$sample_id == sample_id, , drop = FALSE]
  if (nrow(hit) == 0) {
    return(NULL)
  }
  hit[1, , drop = FALSE]
}

first_existing_path <- function(paths) {
  paths <- paths[nzchar(paths)]
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return("")
  }
  existing[1]
}

pick_expression_matrix <- function(x) {
  if (!is.list(x)) {
    return(x)
  }
  if ("Gene Expression" %in% names(x)) {
    return(x[["Gene Expression"]])
  }
  x[[1]]
}

resolve_metrics_summary_path_local <- function(cfg, sample_id, inventory_df = NULL) {
  inventory_row <- resolve_inventory_row_local(cfg, sample_id, inventory_df = inventory_df)
  candidates <- character(0)

  if (!is.null(inventory_row) && "metrics_path" %in% colnames(inventory_row)) {
    metrics_path <- normalize_scalar_value(inventory_row$metrics_path[1])
    if (nzchar(metrics_path) && file.exists(metrics_path)) {
      return(metrics_path)
    }
  }

  if (!is.null(inventory_row) && "cellranger_sample_dir" %in% colnames(inventory_row)) {
    sample_dir <- normalize_scalar_value(inventory_row$cellranger_sample_dir[1])
    if (nzchar(sample_dir)) {
      candidates <- c(
        candidates,
        file.path(sample_dir, "outs", "metrics_summary.csv"),
        file.path(sample_dir, "metrics_summary.csv"),
        file.path(sample_dir, "metrics_summary.xls")
      )
    }
  }

  if (nzchar(cfg$cellranger_out_dir)) {
    candidates <- c(
      candidates,
      file.path(cfg$cellranger_out_dir, sample_id, "outs", "metrics_summary.csv"),
      file.path(cfg$cellranger_out_dir, sample_id, "metrics_summary.xls")
    )
  }

  first_existing_path(candidates)
}

canonicalize_metric_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

extract_metric_value <- function(metrics_df, candidates) {
  if (is.null(metrics_df) || nrow(metrics_df) == 0) {
    return(NA_character_)
  }
  metric_map <- setNames(colnames(metrics_df), canonicalize_metric_name(colnames(metrics_df)))
  for (candidate in candidates) {
    key <- canonicalize_metric_name(candidate)
    if (key %in% names(metric_map)) {
      value <- metrics_df[[metric_map[[key]]]][1]
      if (length(value) > 0 && !is.na(value)) {
        return(as.character(value))
      }
    }
  }
  NA_character_
}

as_numeric_metric <- function(value) {
  if (length(value) == 0 || is.na(value) || !nzchar(value)) {
    return(NA_real_)
  }
  as.numeric(gsub("[^0-9.]+", "", value))
}

read_metrics_table_local <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    return(NULL)
  }

  readers <- list(
    function() read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    function() read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    function() utils::read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  )

  df <- NULL
  for (reader in readers) {
    df <- tryCatch(reader(), error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) {
      break
    }
  }
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  df
}

read_cellranger_metrics_local <- function(path) {
  df <- read_metrics_table_local(path)
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  list(
    estimated_cells = as_numeric_metric(extract_metric_value(df, c("Estimated Number of Cells"))),
    mean_reads_per_cell = as_numeric_metric(extract_metric_value(df, c("Mean Reads per Cell"))),
    median_genes_per_cell = as_numeric_metric(extract_metric_value(df, c("Median Genes per Cell", "Median Genes per Cell Associated with Cell"))),
    sequencing_saturation = as_numeric_metric(extract_metric_value(df, c("Sequencing Saturation"))),
    fraction_reads_in_cells = as_numeric_metric(extract_metric_value(df, c("Fraction Reads in Cells", "Reads Mapped Confidently to Transcriptome")))
  )
}

get_sample_qc_thresholds_local <- function(cfg, sample_id, overrides = NULL) {
  thresholds <- list(
    qc_min_nfeature = cfg$qc_min_nfeature,
    qc_min_ncount = cfg$qc_min_ncount,
    qc_min_log10umi = cfg$qc_min_log10umi,
    qc_max_mito_pct = cfg$qc_max_mito_pct
  )

  if (is.null(overrides)) {
    overrides <- read_qc_threshold_overrides_local(cfg)
  }
  if (nrow(overrides) == 0 || !"sample_id" %in% colnames(overrides)) {
    return(thresholds)
  }

  hit <- overrides[overrides$sample_id == sample_id, , drop = FALSE]
  if (nrow(hit) == 0) {
    return(thresholds)
  }
  for (name in names(thresholds)) {
    override_value <- suppressWarnings(as.numeric(hit[[name]][1]))
    if (!is.na(override_value)) {
      thresholds[[name]] <- override_value
    }
  }
  thresholds
}
