inventory_gate_flag <- function(name, default = "yes") {
  value <- tolower(trimws(Sys.getenv(name, unset = default)))
  value %in% c("yes", "true", "1", "on")
}

inventory_gate_csv <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  out <- trimws(unlist(strsplit(value, "[,;]+", perl = TRUE), use.names = FALSE))
  out[nzchar(out)]
}

read_inventory_gate_eligibility <- function(manifest_path, required = NULL) {
  required <- if (is.null(required)) inventory_gate_flag("INVENTORY_GATE_REQUIRED", "yes") else isTRUE(required)
  path <- ""
  if (nzchar(manifest_path %||% "") && file.exists(manifest_path)) {
    manifest <- tryCatch(read_manifest_local(manifest_path), error = function(e) NULL)
    if (!is.null(manifest) && !is.null(manifest$outputs$cluster_eligibility_tsv)) {
      path <- tryCatch(resolve_output_local(manifest, "cluster_eligibility_tsv"), error = function(e) "")
    }
  }
  if (!nzchar(path) || !file.exists(path)) {
    if (required) {
      stop(sprintf("inventory gate requires cluster_eligibility_tsv from manifest: %s", manifest_path), call. = FALSE)
    }
    return(list(path = path, eligibility = data.frame(stringsAsFactors = FALSE), status = "missing"))
  }
  list(path = path, eligibility = read_tsv_optional(path), status = "ok")
}

filter_eligibility_for_communication <- function(eligibility_df, comparison_id = NULL, allowed_tiers = NULL) {
  if (is.null(allowed_tiers)) {
    allowed_tiers <- inventory_gate_csv("INVENTORY_GATE_ALLOWED_TIERS_COMMUNICATION", "primary,exploratory,primary_merged")
  }
  if (nrow(eligibility_df) == 0 || !"evidence_tier" %in% colnames(eligibility_df)) {
    return(eligibility_df[FALSE, , drop = FALSE])
  }
  df <- eligibility_df
  if (!is.null(comparison_id) && nzchar(comparison_id) && "comparison_id" %in% colnames(df)) {
    hit <- df[as.character(df$comparison_id) == comparison_id, , drop = FALSE]
    if (nrow(hit) > 0) {
      df <- hit
    }
  }
  df[as.character(df$evidence_tier) %in% allowed_tiers, , drop = FALSE]
}

inventory_gate_log_for_labels <- function(labels, eligibility_df, layer_id = "", pair_id = "", comparison_id = NULL, allowed_tiers = NULL, strict = NULL) {
  labels <- sort(unique(as.character(labels)))
  labels <- labels[nzchar(labels)]
  if (is.null(allowed_tiers)) {
    allowed_tiers <- inventory_gate_csv("INVENTORY_GATE_ALLOWED_TIERS_COMMUNICATION", "primary,exploratory,primary_merged")
  }
  strict <- if (is.null(strict)) inventory_gate_flag("INVENTORY_GATE_STRICT_MODE", "no") else isTRUE(strict)
  if (length(labels) == 0) {
    return(data.frame(layer_id = character(), pair_id = character(), cluster_id = character(), passed_inventory_gate = logical(), evidence_tiers = character(), reason_if_filtered_out = character(), stringsAsFactors = FALSE))
  }
  rows <- lapply(labels, function(label) {
    df <- eligibility_df
    if ("layer_id" %in% colnames(df) && nzchar(layer_id)) {
      layer_hit <- df[as.character(df$layer_id) == layer_id, , drop = FALSE]
      if (nrow(layer_hit) > 0) {
        df <- layer_hit
      }
    }
    if ("cluster_id" %in% colnames(df)) {
      df <- df[as.character(df$cluster_id) == label, , drop = FALSE]
    } else {
      df <- df[FALSE, , drop = FALSE]
    }
    if (!is.null(comparison_id) && nzchar(comparison_id) && "comparison_id" %in% colnames(df)) {
      cmp_hit <- df[as.character(df$comparison_id) == comparison_id, , drop = FALSE]
      if (nrow(cmp_hit) > 0) {
        df <- cmp_hit
      }
    }
    tiers <- if (nrow(df) > 0 && "evidence_tier" %in% colnames(df)) unique(as.character(df$evidence_tier)) else character(0)
    passed <- if (length(tiers) == 0) {
      FALSE
    } else if (strict) {
      all(tiers %in% allowed_tiers)
    } else {
      any(tiers %in% allowed_tiers)
    }
    data.frame(
      layer_id = layer_id,
      pair_id = pair_id,
      cluster_id = label,
      passed_inventory_gate = passed,
      evidence_tiers = paste(tiers, collapse = ","),
      reason_if_filtered_out = if (passed) "passed" else if (length(tiers) == 0) "no_matching_eligibility" else paste(tiers, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

subset_by_inventory_gate <- function(seu, cell_type_col, eligibility_df, layer_id = "", pair_id = "", comparison_id = NULL) {
  meta <- seu@meta.data
  labels <- as.character(meta[[cell_type_col]])
  log_df <- inventory_gate_log_for_labels(labels, eligibility_df, layer_id, pair_id, comparison_id)
  eligible <- log_df$cluster_id[isTRUE(nrow(log_df) > 0) & log_df$passed_inventory_gate]
  keep <- rownames(meta)[labels %in% eligible]
  if (length(keep) == 0) {
    fail <- log_df[!log_df$passed_inventory_gate, , drop = FALSE]
    reason <- if (nrow(fail) == 0) "no labels available for inventory gate" else paste(unique(fail$reason_if_filtered_out), collapse = ";")
    return(list(object = NULL, status = "skipped_inventory_gate", reason = reason, log = log_df, eligible_clusters = eligible))
  }
  list(object = subset(seu, cells = keep), status = "ok", reason = "", log = log_df, eligible_clusters = eligible)
}
