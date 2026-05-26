evidence_norm_05 <- function(value, default = "") {
  if (is.null(value) || length(value) == 0) {
    return(default)
  }
  value <- trimws(as.character(value[[1]]))
  if (is.na(value) || !nzchar(value)) default else value
}

manifest_output_optional_05 <- function(manifest_path, key) {
  manifest_path <- evidence_norm_05(manifest_path)
  if (!nzchar(manifest_path) || !file.exists(manifest_path)) {
    return("")
  }
  manifest <- tryCatch(read_manifest_local(manifest_path), error = function(e) NULL)
  if (is.null(manifest) || is.null(manifest$outputs[[key]])) {
    return("")
  }
  tryCatch(resolve_output_local(manifest, key), error = function(e) "")
}

read_cluster_eligibility_05 <- function(manifest_path, key = "cluster_eligibility_tsv") {
  path <- manifest_output_optional_05(manifest_path, key)
  if (!nzchar(path) || !file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")
}

evidence_empty_decision_05 <- function(status = "missing_cluster_eligibility") {
  list(
    evidence_tier = status,
    recommended_action = "skip",
    reason = "cluster_eligibility_tsv is missing or empty",
    warning_banner = "WARNING: 04c cluster eligibility is unavailable; no cell-level fallback will be run.",
    parent_cluster_if_merged = "",
    evidence_row_n = 0L,
    evidence_tier_summary = ""
  )
}

evidence_tier_priority_05 <- function(tier) {
  tier <- as.character(tier)
  priority <- c(
    skip = 60L,
    candidate_only = 50L,
    module_score_only = 40L,
    merge_to_parent = 30L,
    exploratory = 20L,
    primary = 10L
  )
  out <- unname(priority[tier])
  out[is.na(out)] <- 55L
  out
}

evidence_decision_05 <- function(eligibility_df, layer_id, comparison_id, group_var = "", cluster_id = "") {
  if (is.null(eligibility_df) || nrow(eligibility_df) == 0) {
    return(evidence_empty_decision_05())
  }
  required <- c("comparison_id", "evidence_tier", "recommended_action")
  if (length(setdiff(required, colnames(eligibility_df))) > 0) {
    return(evidence_empty_decision_05("malformed_cluster_eligibility"))
  }

  rows <- eligibility_df[as.character(eligibility_df$comparison_id) == comparison_id, , drop = FALSE]
  if ("layer_id" %in% colnames(rows) && nzchar(layer_id)) {
    rows <- rows[as.character(rows$layer_id) == layer_id, , drop = FALSE]
  }
  if ("group_var" %in% colnames(rows) && nzchar(group_var)) {
    rows <- rows[as.character(rows$group_var) == group_var, , drop = FALSE]
  }
  if ("cluster_id" %in% colnames(rows) && nzchar(cluster_id)) {
    rows <- rows[as.character(rows$cluster_id) == cluster_id, , drop = FALSE]
  }
  if (nrow(rows) == 0) {
    return(evidence_empty_decision_05("no_matching_cluster_eligibility"))
  }

  tier_counts <- as.data.frame(table(evidence_tier = rows$evidence_tier), stringsAsFactors = FALSE)
  tier_summary <- paste(sprintf("%s=%s", tier_counts$evidence_tier, tier_counts$Freq), collapse = ";")
  tier <- as.character(rows$evidence_tier[order(evidence_tier_priority_05(rows$evidence_tier), decreasing = TRUE)][[1]])
  selected <- rows[as.character(rows$evidence_tier) == tier, , drop = FALSE][1, , drop = FALSE]
  parent <- if ("parent_cluster_id" %in% colnames(selected)) evidence_norm_05(selected$parent_cluster_id) else ""
  reason <- if ("reason" %in% colnames(selected)) evidence_norm_05(selected$reason) else tier_summary
  action <- evidence_norm_05(selected$recommended_action, "skip")

  warning <- switch(
    tier,
    primary = "",
    exploratory = "WARNING: evidence_tier=exploratory; pseudobulk output must be reported as exploratory.",
    module_score_only = "WARNING: evidence_tier=module_score_only; formal DE is not run.",
    merge_to_parent = sprintf("WARNING: evidence_tier=merge_to_parent; target parent=%s.", parent),
    candidate_only = "WARNING: evidence_tier=candidate_only; report as hypothesis only, no formal DE.",
    skip = "WARNING: evidence_tier=skip; comparison is not biologically comparable.",
    sprintf("WARNING: unknown evidence_tier=%s; formal DE is not run.", tier)
  )

  list(
    evidence_tier = tier,
    recommended_action = action,
    reason = reason,
    warning_banner = warning,
    parent_cluster_if_merged = parent,
    evidence_row_n = nrow(rows),
    evidence_tier_summary = tier_summary
  )
}

evidence_allows_formal_05 <- function(decision) {
  decision$evidence_tier %in% c("primary", "exploratory")
}

evidence_manifest_columns_05 <- function(df) {
  defaults <- list(
    evidence_tier = "",
    recommended_action = "",
    evidence_tier_summary = "",
    warning_banner = "",
    parent_cluster_if_merged = "",
    cluster_eligibility_tsv = ""
  )
  for (name in names(defaults)) {
    if (!name %in% colnames(df)) {
      df[[name]] <- defaults[[name]]
    }
  }
  df
}
