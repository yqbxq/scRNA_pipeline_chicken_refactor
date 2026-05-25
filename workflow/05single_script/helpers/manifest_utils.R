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

manifest_timestamp_now_local <- function() {
  if (exists("timestamp_now", mode = "function")) {
    return(timestamp_now())
  }
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

compute_file_sha256 <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!file.exists(path) || dir.exists(path)) {
    return(if (dir.exists(path)) "dir" else "missing")
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    con <- file(path, "rb")
    on.exit(close(con), add = TRUE)
    return(as.character(openssl::sha256(con)))
  }
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    out <- system2(sha256sum, path, stdout = TRUE)
    return(strsplit(out[[1]], "[[:space:]]+")[[1]][[1]])
  }
  shasum <- Sys.which("shasum")
  if (nzchar(shasum)) {
    out <- system2(shasum, c("-a", "256", path), stdout = TRUE)
    return(strsplit(out[[1]], "[[:space:]]+")[[1]][[1]])
  }
  stop("missing sha256 implementation for R manifest fingerprints", call. = FALSE)
}

compute_string_sha256_local <- function(value) {
  if (requireNamespace("openssl", quietly = TRUE)) {
    return(as.character(openssl::sha256(charToRaw(value))))
  }
  tmp <- tempfile("checkpoint_hash_")
  on.exit(unlink(tmp), add = TRUE)
  writeBin(charToRaw(value), tmp)
  compute_file_sha256(tmp)
}

expand_manifest_inputs_R <- function(path, depth = 0L, visited = character()) {
  if (depth > 5L) {
    stop(sprintf("manifest nesting exceeds 5: %s", path), call. = FALSE)
  }
  if (!file.exists(path) || dir.exists(path) || !grepl("\\.json$", path, ignore.case = TRUE)) {
    return(path)
  }
  manifest <- tryCatch(jsonlite::read_json(path, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(manifest) || is.null(manifest$outputs)) {
    return(path)
  }
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (normalized %in% visited) {
    return(character())
  }
  visited <- c(visited, normalized)
  base_dir <- manifest$base_dir %||% dirname(path)
  out <- character()
  for (entry in manifest$outputs) {
    raw_path <- entry$path %||% ""
    if (!nzchar(raw_path)) next
    next_path <- if (grepl("^/", raw_path)) raw_path else file.path(base_dir, raw_path)
    out <- c(out, expand_manifest_inputs_R(next_path, depth + 1L, visited))
  }
  out
}

expand_input_paths_R <- function(input_paths) {
  out <- character()
  for (path in input_paths) {
    if (!nzchar(path)) next
    expanded <- expand_manifest_inputs_R(path)
    for (item in expanded) {
      if (dir.exists(item)) {
        files <- list.files(item, recursive = TRUE, all.files = FALSE, full.names = TRUE, no.. = TRUE)
        files <- files[file.exists(files) & !dir.exists(files)]
        out <- c(out, files)
      } else {
        out <- c(out, item)
      }
    }
  }
  out
}

compute_input_hash_R <- function(input_paths) {
  items <- expand_input_paths_R(input_paths)
  if (length(items) == 0L) {
    return(compute_string_sha256_local(""))
  }
  rows <- vapply(sort(unique(items)), function(path) {
    sprintf("%s\t%s", path, compute_file_sha256(path))
  }, character(1))
  compute_string_sha256_local(paste0(paste(rows, collapse = "\n"), "\n"))
}

compute_script_hash_R <- function(script_path) {
  files <- script_path
  if (file.exists(script_path)) {
    content <- paste(readLines(script_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    if (grepl("\\.R$", script_path)) {
      matches <- gregexpr("(source|source_utf8)\\s*\\(\\s*[\"']([^\"']+)[\"']\\s*\\)", content, perl = TRUE)
      hits <- regmatches(content, matches)[[1]]
      if (length(hits) > 0L && hits[[1]] != "-1") {
        helpers <- sub(".*[\"']([^\"']+)[\"'].*", "\\1", hits)
        files <- c(files, file.path(dirname(script_path), helpers))
      }
    } else if (grepl("\\.py$", script_path)) {
      matches <- gregexpr("(from\\s+helpers\\.([A-Za-z0-9_]+)\\s+import|import\\s+helpers\\.([A-Za-z0-9_]+))", content, perl = TRUE)
      hits <- regmatches(content, matches)[[1]]
      if (length(hits) > 0L && hits[[1]] != "-1") {
        mods <- sub(".*helpers\\.([A-Za-z0-9_]+).*", "\\1", hits)
        files <- c(files, file.path(dirname(script_path), "helpers", paste0(mods, ".py")))
      }
    }
  }
  rows <- vapply(sort(unique(normalizePath(files, winslash = "/", mustWork = FALSE))), function(path) {
    sprintf("%s\t%s", path, compute_file_sha256(path))
  }, character(1))
  compute_string_sha256_local(paste0(paste(rows, collapse = "\n"), "\n"))
}

compute_params_hash_R <- function(params_names, params_env = Sys.getenv) {
  params_names <- sort(unique(params_names[nzchar(params_names)]))
  if (length(params_names) == 0L) {
    return("")
  }
  values <- vapply(params_names, function(name) {
    value <- params_env(name, unset = "__UNSET__")
    sprintf("%s=%s", name, value)
  }, character(1))
  compute_string_sha256_local(paste0(paste(values, collapse = "\n"), "\n"))
}

build_fingerprints <- function(input_paths, script_path, params_names = character(), params_env = Sys.getenv, checkpoint_mode = "fingerprint") {
  params_names <- sort(unique(params_names[nzchar(params_names)]))
  list(
    algorithm = "sha256",
    computed_at = manifest_timestamp_now_local(),
    input_hash = compute_input_hash_R(input_paths),
    script_hash = compute_script_hash_R(script_path),
    params_hash = compute_params_hash_R(params_names, params_env),
    params_names = as.list(params_names),
    checkpoint_mode = checkpoint_mode
  )
}

write_manifest_local <- function(manifest_path, new_outputs, module_name = NULL, base_dir = NULL, inputs = NULL, version = "1.0", depends_on = list(), fingerprints = NULL) {
  if (file.exists(manifest_path)) {
    existing <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
    existing$outputs <- modifyList(existing$outputs %||% list(), new_outputs)
    existing$timestamp <- manifest_timestamp_now_local()
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
      timestamp = manifest_timestamp_now_local(),
      base_dir = base_dir,
      inputs = inputs,
      outputs = new_outputs,
      depends_on = depends_on
    )
  }

  if (!is.null(fingerprints)) {
    manifest$fingerprints <- fingerprints
  }

  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}
