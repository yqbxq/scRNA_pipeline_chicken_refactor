read_manifest_local <- function(manifest_path) {
  if (!file.exists(manifest_path)) {
    stop(sprintf("缺少 manifest: %s", manifest_path), call. = FALSE)
  }
  jsonlite::read_json(manifest_path, simplifyVector = FALSE)
}

resolve_output_local <- function(manifest, key) {
  entry <- manifest$outputs[[key]]
  if (is.null(entry) || is.null(entry$path)) {
    stop(sprintf("manifest 缺少输出键: %s", key), call. = FALSE)
  }
  raw_path <- as.character(entry$path)
  if (grepl("^/", raw_path)) {
    normalizePath(raw_path, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(as.character(manifest$base_dir), raw_path), winslash = "/", mustWork = FALSE)
  }
}

relative_path_local <- function(path, base_dir) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  normalized_base <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
  prefix <- paste0(normalized_base, "/")
  if (startsWith(normalized_path, prefix)) {
    substring(normalized_path, nchar(prefix) + 1L)
  } else {
    normalized_path
  }
}

infer_schema_from_df <- function(df) {
  schema_type <- function(column) {
    if (is.logical(column)) {
      return("logical")
    }
    if (is.integer(column) || is.numeric(column)) {
      return("numeric")
    }
    "character"
  }

  schema <- lapply(df, schema_type)
  names(schema) <- names(df)
  schema
}

build_output_entry <- function(path, type, produced_by, row_semantics, base_dir, schema = NULL) {
  entry <- list(
    path = relative_path_local(path, base_dir),
    type = type,
    produced_by = produced_by,
    row_semantics = row_semantics
  )
  if (!is.null(schema)) {
    entry$schema <- schema
  }
  entry
}

write_manifest_local <- function(manifest_path, new_outputs, module_name = NULL, base_dir = NULL, inputs = NULL, version = "1.0", depends_on = list()) {
  if (file.exists(manifest_path)) {
    existing <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
    existing$outputs <- modifyList(existing$outputs %||% list(), new_outputs)
    existing$timestamp <- timestamp_now()
    existing$module <- existing$module %||% module_name
    existing$version <- existing$version %||% version
    existing$base_dir <- existing$base_dir %||% base_dir
    existing$inputs <- existing$inputs %||% inputs
    existing$depends_on <- existing$depends_on %||% depends_on
    manifest <- existing
  } else {
    if (is.null(module_name) || is.null(base_dir) || is.null(inputs)) {
      stop("manifest 不存在且缺少初始化参数", call. = FALSE)
    }
    manifest <- list(
      module = module_name,
      version = version,
      timestamp = timestamp_now(),
      base_dir = base_dir,
      inputs = inputs,
      outputs = new_outputs,
      depends_on = depends_on
    )
  }

  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}
