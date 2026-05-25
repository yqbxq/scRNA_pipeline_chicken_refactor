`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

h5ad_contract_repo_template_path <- function() {
  workflow_root <- Sys.getenv("WORKFLOW_ROOT", unset = "")
  if (nzchar(workflow_root)) {
    return(file.path(workflow_root, "..", "metadata", "h5ad_export_contract.tsv.template"))
  }
  file.path(getwd(), "metadata", "h5ad_export_contract.tsv.template")
}

load_h5ad_contract <- function(contract_path = Sys.getenv("H5AD_CONTRACT_FILE", "metadata/h5ad_export_contract.tsv"), modality = c("scrna", "spatial")) {
  modality <- match.arg(modality)
  if (!file.exists(contract_path)) {
    contract_path <- h5ad_contract_repo_template_path()
  }
  if (!file.exists(contract_path)) {
    stop(sprintf("missing H5AD contract: %s", contract_path), call. = FALSE)
  }
  contract <- read.delim(contract_path, sep = "\t", comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)
  validate_h5ad_contract_schema(contract)
  contract[contract$modality %in% c(modality, "both"), , drop = FALSE]
}

validate_h5ad_contract_schema <- function(contract) {
  required_cols <- c("field_path", "field_type", "required", "dtype", "gate", "modality", "semantics", "notes")
  missing_cols <- setdiff(required_cols, colnames(contract))
  if (length(missing_cols) > 0L) {
    stop(sprintf("H5AD contract missing columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  allowed_type <- c("scalar", "matrix", "dict", "dataframe")
  allowed_required <- c("yes", "no", "conditional")
  allowed_gate <- c("primary", "exploratory", "conditional")
  allowed_modality <- c("scrna", "spatial", "both")
  bad_type <- unique(contract$field_type[!contract$field_type %in% allowed_type])
  bad_required <- unique(contract$required[!contract$required %in% allowed_required])
  bad_gate <- unique(contract$gate[!contract$gate %in% allowed_gate])
  bad_modality <- unique(contract$modality[!contract$modality %in% allowed_modality])
  if (length(bad_type) > 0L) stop(sprintf("invalid field_type: %s", paste(bad_type, collapse = ", ")), call. = FALSE)
  if (length(bad_required) > 0L) stop(sprintf("invalid required: %s", paste(bad_required, collapse = ", ")), call. = FALSE)
  if (length(bad_gate) > 0L) stop(sprintf("invalid gate: %s", paste(bad_gate, collapse = ", ")), call. = FALSE)
  if (length(bad_modality) > 0L) stop(sprintf("invalid modality: %s", paste(bad_modality, collapse = ", ")), call. = FALSE)

  core <- c("obs.sample_id", "obs.condition", "var.gene_id", "layers.counts")
  missing_core <- setdiff(core, contract$field_path)
  if (length(missing_core) > 0L) {
    stop(sprintf("H5AD contract removed core fields: %s", paste(missing_core, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

is_mock_h5ad_object <- function(obj) {
  is.list(obj) && any(c("obs", "var", "obsm", "layers", "uns", "spatial") %in% names(obj))
}

contract_location_field <- function(field_path) {
  parts <- strsplit(field_path, ".", fixed = TRUE)[[1]]
  list(location = parts[[1]], field = paste(parts[-1], collapse = "."))
}

mock_contract_value <- function(obj, field_path) {
  parsed <- contract_location_field(field_path)
  location <- parsed$location
  field <- parsed$field
  container <- obj[[location]]
  if (is.null(container)) {
    return(list(found = FALSE, value = NULL))
  }
  if (is.data.frame(container)) {
    if (field %in% colnames(container)) return(list(found = TRUE, value = container[[field]]))
    return(list(found = FALSE, value = NULL))
  }
  if (is.list(container)) {
    parts <- strsplit(field, ".", fixed = TRUE)[[1]]
    cursor <- container
    for (part in parts) {
      if (!is.list(cursor) || is.null(cursor[[part]])) {
        return(list(found = FALSE, value = NULL))
      }
      cursor <- cursor[[part]]
    }
    return(list(found = TRUE, value = cursor))
  }
  list(found = FALSE, value = NULL)
}

seurat_contract_value <- function(seurat_obj, field_path) {
  parsed <- contract_location_field(field_path)
  location <- parsed$location
  field <- parsed$field
  value <- NULL
  found <- FALSE

  if (location == "obs" && !is.null(seurat_obj@meta.data) && field %in% colnames(seurat_obj@meta.data)) {
    found <- TRUE
    value <- seurat_obj@meta.data[[field]]
  } else if (location == "var") {
    assay_name <- if (exists("DefaultAssay", mode = "function")) DefaultAssay(seurat_obj) else names(seurat_obj@assays)[[1]]
    meta_features <- seurat_obj@assays[[assay_name]]@meta.features
    if (!is.null(meta_features) && field %in% colnames(meta_features)) {
      found <- TRUE
      value <- meta_features[[field]]
    }
  } else if (location == "obsm") {
    reductions <- if (exists("Reductions", mode = "function")) Reductions(seurat_obj) else names(seurat_obj@reductions)
    reduction_name <- sub("^X_", "", field)
    if (field %in% reductions || reduction_name %in% reductions) {
      found <- TRUE
      value <- matrix(0, nrow = ncol(seurat_obj), ncol = 1)
    }
  } else if (location == "layers") {
    assay_name <- if (exists("DefaultAssay", mode = "function")) DefaultAssay(seurat_obj) else names(seurat_obj@assays)[[1]]
    assay <- seurat_obj@assays[[assay_name]]
    if (field %in% c("counts", "data", "scale.data", "logcounts", "SCT")) {
      found <- TRUE
      value <- tryCatch({
        if (field %in% slotNames(assay)) slot(assay, field) else matrix(0, nrow = 1, ncol = 1)
      }, error = function(e) matrix(0, nrow = 1, ncol = 1))
    }
  } else if (location == "uns" && !is.null(seurat_obj@misc[[field]])) {
    found <- TRUE
    value <- seurat_obj@misc[[field]]
  } else if (location == "spatial" && length(seurat_obj@images) > 0L) {
    found <- TRUE
    value <- seurat_obj@images
  }

  list(found = found, value = value)
}

detect_h5ad_contract_dtype <- function(value) {
  if (is.null(value)) return("missing")
  if (is.data.frame(value)) return("dataframe")
  if (is.list(value) && !is.data.frame(value)) return("dict")
  if (is.matrix(value) || inherits(value, "Matrix")) {
    if (is.integer(value)) return("int64")
    if (is.numeric(value)) return("float32")
    if (is.character(value)) return("str")
    return(class(value)[[1]])
  }
  if (is.integer(value)) return("int64")
  if (is.numeric(value)) return("float32")
  if (is.character(value) || is.factor(value)) return("str")
  if (is.logical(value)) return("bool")
  class(value)[[1]]
}

dtype_compatible_h5ad_contract <- function(actual, expected) {
  expected <- expected %||% ""
  if (!nzchar(expected) || expected == "-") return(TRUE)
  if (actual == expected) return(TRUE)
  if (expected %in% c("float32", "float64") && actual %in% c("float32", "float64", "int64")) return(TRUE)
  if (expected %in% c("int64", "int32") && actual %in% c("int64", "int32")) return(TRUE)
  if (expected == "str" && actual %in% c("str", "factor")) return(TRUE)
  FALSE
}

contract_violation_row <- function(field, severity, reason, expected_dtype = "", actual_dtype = "") {
  data.frame(
    field = field,
    severity = severity,
    reason = reason,
    expected_dtype = expected_dtype,
    actual_dtype = actual_dtype,
    stringsAsFactors = FALSE
  )
}

check_h5ad_contract_seurat <- function(seurat_obj, contract, fail_on = c("error", "warn", "none")) {
  fail_on <- match.arg(fail_on)
  validate_h5ad_contract_schema(contract)
  violations <- list()
  is_mock <- is_mock_h5ad_object(seurat_obj)

  for (i in seq_len(nrow(contract))) {
    row <- contract[i, , drop = FALSE]
    field_path <- row$field_path[[1]]
    required <- row$required[[1]]
    gate <- row$gate[[1]]
    severity <- if (gate == "primary") "error" else "warn"
    value_info <- if (is_mock) mock_contract_value(seurat_obj, field_path) else seurat_contract_value(seurat_obj, field_path)

    if (!isTRUE(value_info$found)) {
      if (required == "yes") {
        parsed <- contract_location_field(field_path)
        violations[[length(violations) + 1L]] <- contract_violation_row(
          field_path,
          severity,
          sprintf("required=yes but field not found in %s", parsed$location),
          row$dtype[[1]],
          "missing"
        )
      }
      next
    }

    actual_dtype <- detect_h5ad_contract_dtype(value_info$value)
    if (!dtype_compatible_h5ad_contract(actual_dtype, row$dtype[[1]])) {
      violations[[length(violations) + 1L]] <- contract_violation_row(
        field_path,
        severity,
        sprintf("expected dtype=%s got %s", row$dtype[[1]], actual_dtype),
        row$dtype[[1]],
        actual_dtype
      )
    }
  }

  vdf <- if (length(violations) > 0L) do.call(rbind, violations) else data.frame(
    field = character(), severity = character(), reason = character(),
    expected_dtype = character(), actual_dtype = character(), stringsAsFactors = FALSE
  )
  passed <- !any(vdf$severity == "error")
  result <- list(
    passed = passed,
    violations = vdf,
    summary = data.frame(
      passed = passed,
      violation_n = nrow(vdf),
      error_n = sum(vdf$severity == "error"),
      warn_n = sum(vdf$severity == "warn"),
      stringsAsFactors = FALSE
    ),
    contract = contract
  )

  if (!passed && fail_on == "error") {
    stop(sprintf(
      "H5AD contract violations: %s",
      paste(vdf$field[vdf$severity == "error"], collapse = ", ")
    ), call. = FALSE)
  }
  if (nrow(vdf) > 0L && fail_on == "warn") {
    warning(sprintf("H5AD contract violations: %d field(s)", nrow(vdf)), call. = FALSE)
  }
  result
}

contract_violations_to_manifest <- function(check_result) {
  violations <- check_result$violations %||% data.frame()
  list(
    contract_passed = isTRUE(check_result$passed),
    contract_violations_count = nrow(violations),
    contract_violations_summary = lapply(seq_len(nrow(violations)), function(i) as.list(violations[i, , drop = FALSE]))
  )
}
