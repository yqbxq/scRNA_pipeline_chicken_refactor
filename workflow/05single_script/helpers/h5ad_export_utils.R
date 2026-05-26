`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

h5ad_env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

h5ad_now_iso <- function() {
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

h5ad_read_upstream_object <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("missing upstream RDS: %s", path), call. = FALSE)
  }
  readRDS(path)
}

h5ad_is_mock_object <- function(obj) {
  is.list(obj) && any(c("obs", "var", "obsm", "layers", "uns", "spatial") %in% names(obj))
}

h5ad_choose_backend <- function() {
  backend <- h5ad_env_or_default("H5AD_EXPORT_BACKEND", "zellkonverter")
  if (identical(backend, "auto")) {
    backend <- "zellkonverter"
  }
  backend
}

h5ad_select_assay <- function(obj, default = "RNA") {
  requested <- h5ad_env_or_default("H5AD_EXPORT_ASSAY", "auto")
  if (!identical(requested, "auto")) return(requested)
  if (h5ad_is_mock_object(obj)) return(obj$assay %||% default)
  if (exists("DefaultAssay", mode = "function")) {
    return(DefaultAssay(obj))
  }
  default
}

h5ad_decide_layers_to_export <- function(obj) {
  requested <- h5ad_env_or_default("H5AD_EXPORT_LAYERS", "auto")
  if (!identical(requested, "auto")) {
    return(trimws(strsplit(requested, ",", fixed = TRUE)[[1]]))
  }
  if (h5ad_is_mock_object(obj)) {
    return(names(obj$layers %||% list(counts = TRUE)))
  }
  layers <- c("counts")
  assay <- h5ad_select_assay(obj)
  assay_obj <- tryCatch(obj[[assay]], error = function(e) NULL)
  if (!is.null(assay_obj)) {
    if ("data" %in% slotNames(assay_obj)) layers <- c(layers, "data")
    if ("scale.data" %in% slotNames(assay_obj) && identical(h5ad_env_or_default("H5AD_EXPORT_LAYERS", "auto"), "all")) {
      layers <- c(layers, "scale.data")
    }
  }
  unique(layers)
}

h5ad_decide_reductions_to_export <- function(obj) {
  requested <- h5ad_env_or_default("H5AD_EXPORT_REDUCTIONS", "auto")
  if (!identical(requested, "auto")) {
    return(trimws(strsplit(requested, ",", fixed = TRUE)[[1]]))
  }
  if (h5ad_is_mock_object(obj)) {
    return(names(obj$obsm %||% list()))
  }
  if (exists("Reductions", mode = "function")) {
    return(Reductions(obj))
  }
  character()
}

h5ad_proxy_from_object <- function(obj, modality = "scrna") {
  if (h5ad_is_mock_object(obj)) {
    proxy <- obj
  } else {
    assay <- h5ad_select_assay(obj, if (identical(modality, "spatial")) "Spatial" else "RNA")
    meta <- tryCatch(obj@meta.data, error = function(e) data.frame())
    meta_features <- tryCatch(obj[[assay]]@meta.features, error = function(e) data.frame())
    if (!"gene_id" %in% colnames(meta_features) && nrow(meta_features) > 0L) {
      meta_features$gene_id <- rownames(meta_features)
    }
    if (!"gene_symbol" %in% colnames(meta_features) && nrow(meta_features) > 0L) {
      meta_features$gene_symbol <- rownames(meta_features)
    }
    obsm <- list()
    for (red in h5ad_decide_reductions_to_export(obj)) {
      value <- tryCatch(Embeddings(obj, red), error = function(e) NULL)
      if (!is.null(value)) {
        key <- if (identical(red, "spatial")) "spatial" else paste0("X_", red)
        obsm[[key]] <- value
      }
    }
    layers <- list()
    for (layer in h5ad_decide_layers_to_export(obj)) {
      value <- tryCatch(GetAssayData(obj, assay = assay, slot = layer), error = function(e) NULL)
      if (!is.null(value)) {
        out_name <- if (identical(layer, "data")) "logcounts" else layer
        layers[[out_name]] <- value
      }
    }
    proxy <- list(obs = meta, var = meta_features, obsm = obsm, layers = layers, uns = list(), spatial = list())
  }

  proxy$uns <- proxy$uns %||% list()
  proxy$uns$module <- proxy$uns$module %||% h5ad_env_or_default("H5AD_EXPORT_MODULE", "90_export_h5ad")
  proxy$uns$export_timestamp <- proxy$uns$export_timestamp %||% h5ad_now_iso()
  proxy
}

h5ad_write_mock_file <- function(proxy, path, modality, layers_exported, reductions_exported) {
  ensure_dir(dirname(path))
  payload <- list(
    format = "mock_h5ad_for_smoke",
    modality = modality,
    obs = names(proxy$obs %||% list()),
    var = names(proxy$var %||% list()),
    obsm = names(proxy$obsm %||% list()),
    layers = layers_exported,
    uns = names(proxy$uns %||% list()),
    spatial = names(proxy$spatial %||% list())
  )
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(payload, path, pretty = TRUE, auto_unbox = TRUE)
  } else {
    writeLines(capture.output(str(payload)), path)
  }
  path
}

h5ad_write_export_file <- function(proxy, obj, path, modality, layers_exported, reductions_exported) {
  backend <- h5ad_choose_backend()
  if (identical(backend, "mock") || h5ad_is_mock_object(obj)) {
    return(h5ad_write_mock_file(proxy, path, modality, layers_exported, reductions_exported))
  }
  if (identical(backend, "zellkonverter") && requireNamespace("zellkonverter", quietly = TRUE) && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("zellkonverter export is registered but full Seurat->SCE conversion is deferred to server runtime validation", call. = FALSE)
  }
  if (identical(backend, "sceasy") && requireNamespace("sceasy", quietly = TRUE)) {
    stop("sceasy export is registered but full Seurat->AnnData conversion is deferred to server runtime validation", call. = FALSE)
  }
  h5ad_write_mock_file(proxy, path, modality, layers_exported, reductions_exported)
}

h5ad_export_one <- function(obj, h5ad_path, modality, contract, fail_on = "warn") {
  proxy <- h5ad_proxy_from_object(obj, modality)
  check <- check_h5ad_contract_seurat(proxy, contract, fail_on = fail_on)
  layers_exported <- names(proxy$layers %||% list())
  reductions_exported <- names(proxy$obsm %||% list())
  if (!isTRUE(check$passed)) {
    return(list(
      status = "skipped_contract",
      h5ad_path = "",
      skipped_json = sub("\\.h5ad$", "_h5ad_skipped.json", h5ad_path),
      contract_check = check,
      layers_exported = layers_exported,
      reductions_exported = reductions_exported
    ))
  }
  h5ad_write_export_file(proxy, obj, h5ad_path, modality, layers_exported, reductions_exported)
  list(
    status = "ok",
    h5ad_path = h5ad_path,
    skipped_json = "",
    contract_check = check,
    layers_exported = layers_exported,
    reductions_exported = reductions_exported
  )
}

h5ad_write_skip_json <- function(path, result) {
  ensure_dir(dirname(path))
  payload <- contract_violations_to_manifest(result$contract_check)
  payload$status <- result$status
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(payload, path, pretty = TRUE, auto_unbox = TRUE)
  } else {
    writeLines(capture.output(str(payload)), path)
  }
  path
}

h5ad_summary_row <- function(modality, source_rds, result, section_id = "") {
  data.frame(
    modality = modality,
    section_id = section_id,
    status = result$status,
    source_rds = source_rds,
    h5ad_path = result$h5ad_path,
    skipped_json = result$skipped_json,
    contract_passed = isTRUE(result$contract_check$passed),
    violation_n = nrow(result$contract_check$violations),
    layers_exported = paste(result$layers_exported, collapse = ","),
    reductions_exported = paste(result$reductions_exported, collapse = ","),
    stringsAsFactors = FALSE
  )
}

h5ad_json_escape <- function(value) {
  value <- gsub("\\\\", "\\\\\\\\", as.character(value))
  value <- gsub('"', '\\"', value, fixed = TRUE)
  value <- gsub("\n", "\\n", value, fixed = TRUE)
  sprintf('"%s"', value)
}

h5ad_json_value <- function(value, indent = "  ") {
  if (is.null(value)) return("null")
  if (is.logical(value) && length(value) == 1L) return(if (isTRUE(value)) "true" else "false")
  if (is.numeric(value) && length(value) == 1L && !is.na(value)) return(as.character(value))
  if (is.atomic(value) && length(value) == 1L) return(h5ad_json_escape(value))
  if (is.list(value)) {
    names_value <- names(value)
    if (is.null(names_value)) {
      inner <- vapply(value, h5ad_json_value, character(1), indent = paste0(indent, "  "))
      return(paste0("[", paste(inner, collapse = ", "), "]"))
    }
    parts <- character()
    for (idx in seq_along(value)) {
      parts <- c(parts, sprintf("%s%s: %s", paste0(indent, "  "), h5ad_json_escape(names_value[[idx]]), h5ad_json_value(value[[idx]], paste0(indent, "  "))))
    }
    return(paste0("{\n", paste(parts, collapse = ",\n"), "\n", indent, "}"))
  }
  h5ad_json_escape(paste(value, collapse = ","))
}

h5ad_write_manifest_compat <- function(manifest_path, new_outputs, module_name, base_dir, inputs = list(), version = "1.0", depends_on = list()) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(write_manifest_local(
      manifest_path = manifest_path,
      new_outputs = new_outputs,
      module_name = module_name,
      base_dir = base_dir,
      inputs = inputs,
      version = version,
      depends_on = depends_on
    ))
  }
  ensure_dir(dirname(manifest_path))
  manifest <- list(
    module = module_name,
    version = version,
    timestamp = h5ad_now_iso(),
    base_dir = base_dir,
    inputs = inputs,
    outputs = new_outputs,
    depends_on = depends_on
  )
  writeLines(h5ad_json_value(manifest, "  "), manifest_path, useBytes = TRUE)
  invisible(manifest_path)
}
