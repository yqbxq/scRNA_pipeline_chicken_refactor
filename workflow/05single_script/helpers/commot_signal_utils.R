if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
}

commot_scalar <- function(x, default = "") {
  if (exists("normalize_scalar_value", mode = "function")) {
    return(normalize_scalar_value(x, default))
  }
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(as.character(x)))) default else trimws(as.character(x))
}

commot_truthy <- function(x) {
  tolower(commot_scalar(x, "false")) %in% c("1", "true", "yes", "y", "ok")
}

empty_commot_lr_candidates <- function() {
  data.frame(
    lr_axis_id = character(),
    ligand = character(),
    receptor = character(),
    pathway = character(),
    sender = character(),
    receiver = character(),
    pair_id = character(),
    condition_value = character(),
    stringsAsFactors = FALSE
  )
}

prepare_commot_lr_candidates <- function(consensus_df) {
  if (is.null(consensus_df) || nrow(consensus_df) == 0) {
    return(empty_commot_lr_candidates())
  }
  df <- as.data.frame(consensus_df, stringsAsFactors = FALSE)
  for (col in c("lr_axis_id", "ligand", "receptor", "ligand_human", "receptor_human", "source", "target", "pair_id", "condition_value", "evidence_tier")) {
    if (!col %in% colnames(df)) df[[col]] <- ""
  }
  if ("can_be_primary" %in% colnames(df)) {
    keep <- commot_truthy(df$can_be_primary) | df$evidence_tier %in% c("primary", "exploratory", "candidate")
    df <- df[keep %in% TRUE, , drop = FALSE]
  }
  out <- data.frame(
    lr_axis_id = df$lr_axis_id,
    ligand = ifelse(nzchar(df$ligand_human), df$ligand_human, df$ligand),
    receptor = ifelse(nzchar(df$receptor_human), df$receptor_human, df$receptor),
    pathway = if ("pathway_name" %in% colnames(df)) df$pathway_name else "user",
    sender = df$source,
    receiver = df$target,
    pair_id = df$pair_id,
    condition_value = df$condition_value,
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$lr_axis_id) & nzchar(out$ligand) & nzchar(out$receptor), , drop = FALSE]
  out[!duplicated(out[, c("lr_axis_id", "condition_value"), drop = FALSE]), , drop = FALSE]
}

empty_commot_spatial_summary <- function() {
  data.frame(
    lr_axis_id = character(),
    section_id = character(),
    pair_id = character(),
    condition_value = character(),
    ligand = character(),
    receptor = character(),
    sender = character(),
    receiver = character(),
    signal_score = numeric(),
    spatial_support = character(),
    status = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

summarize_commot_signal <- function(commot_df, signal_threshold = 0) {
  if (is.null(commot_df) || nrow(commot_df) == 0) {
    return(empty_commot_spatial_summary())
  }
  df <- as.data.frame(commot_df, stringsAsFactors = FALSE)
  for (col in colnames(empty_commot_spatial_summary())) {
    if (!col %in% colnames(df)) df[[col]] <- ""
  }
  df$signal_score <- suppressWarnings(as.numeric(df$signal_score))
  df$spatial_support <- ifelse(!is.na(df$signal_score) & df$signal_score > signal_threshold & df$status %in% c("ok", "ok_commot_run"), "yes", "no")
  df[, colnames(empty_commot_spatial_summary()), drop = FALSE]
}

attach_commot_spatial_summary_to_consensus <- function(consensus_df, commot_summary, signal_threshold = 0) {
  if (is.null(consensus_df) || nrow(consensus_df) == 0) {
    return(consensus_df)
  }
  if (is.null(commot_summary) || nrow(commot_summary) == 0) {
    consensus_df$commot_spatial_hit <- FALSE
    return(consensus_df)
  }
  summary <- summarize_commot_signal(commot_summary, signal_threshold)
  supported <- unique(summary$lr_axis_id[summary$spatial_support == "yes"])
  consensus_df$commot_spatial_hit <- consensus_df$lr_axis_id %in% supported
  consensus_df
}
