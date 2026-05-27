if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
}

standardize_lr_axis_id <- function(ligand, receptor, sender_cell_type, receiver_cell_type) {
  ligand_complex <- sort_lr_complex_07(ligand)
  receptor_complex <- sort_lr_complex_07(receptor)
  sender <- normalize_lr_axis_part_07(sender_cell_type)
  receiver <- normalize_lr_axis_part_07(receiver_cell_type)
  sprintf("%s|%s|%s->%s", ligand_complex, receptor_complex, sender, receiver)
}

normalize_lr_axis_part_07 <- function(value) {
  value <- trimws(as.character(value))
  value[is.na(value) | value %in% c("", "-", "NA", "NaN", "NULL")] <- ""
  value
}

sort_lr_complex_07 <- function(name, sep = "[_+/]") {
  value <- normalize_lr_axis_part_07(name)
  vapply(value, function(item) {
    if (!nzchar(item)) {
      return("")
    }
    parts <- unlist(strsplit(item, sep, perl = TRUE), use.names = FALSE)
    parts <- trimws(parts[nzchar(trimws(parts))])
    paste(sort(parts), collapse = "|")
  }, character(1), USE.NAMES = FALSE)
}

standardize_lr_axis_id_in_df <- function(df,
                                         ligand_col = "ligand_human",
                                         receptor_col = "receptor_human",
                                         source_col = "source",
                                         target_col = "target",
                                         out_col = "lr_axis_id",
                                         ortholog_lut = NULL) {
  if (!is.null(ortholog_lut)) {
    if (!ligand_col %in% colnames(df) && "ligand" %in% colnames(df)) {
      df[[ligand_col]] <- ortholog_chicken_to_human_complex_vec(df$ligand, ortholog_lut)
    }
    if (!receptor_col %in% colnames(df) && "receptor" %in% colnames(df)) {
      df[[receptor_col]] <- ortholog_chicken_to_human_complex_vec(df$receptor, ortholog_lut)
    }
  }
  missing_cols <- setdiff(c(ligand_col, receptor_col, source_col, target_col), colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Cannot build lr_axis_id; missing columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  df[[out_col]] <- standardize_lr_axis_id(df[[ligand_col]], df[[receptor_col]], df[[source_col]], df[[target_col]])
  df
}

as_logical_flag_07d <- function(value) {
  value <- tolower(trimws(as.character(value)))
  value %in% c("1", "true", "yes", "y", "ok")
}

first_nonempty_07d <- function(...) {
  values <- list(...)
  n <- max(vapply(values, length, integer(1)), 0L)
  if (n == 0L) {
    return(character(0))
  }
  out <- rep("", n)
  for (value in values) {
    value <- rep_len(as.character(value), n)
    value[is.na(value)] <- ""
    take <- !nzchar(out) & nzchar(value)
    out[take] <- value[take]
  }
  out
}

empty_method_consensus_07d <- function() {
  data.frame(
    lr_axis_id = character(),
    condition_value = character(),
    pair_id = character(),
    layer_id = character(),
    source = character(),
    target = character(),
    ligand = character(),
    receptor = character(),
    ligand_human = character(),
    receptor_human = character(),
    cellchat_hit = logical(),
    cellchat_method_evidence_class = character(),
    cellchat_can_be_primary = character(),
    liana_consensus_hit = logical(),
    liana_n_methods_agreed = numeric(),
    liana_consensus_score = numeric(),
    nichenet_hit = logical(),
    nichenet_method_evidence_class = character(),
    nichenet_n_targets_in_receiver_de = numeric(),
    commot_spatial_hit = logical(),
    evidence_method_n = numeric(),
    evidence_tier = character(),
    evidence_reason = character(),
    can_be_primary = character(),
    stringsAsFactors = FALSE
  )
}

method_prefixed_df_07d <- function(df, method) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(lr_axis_id = character(), condition_value = character(), stringsAsFactors = FALSE))
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (col in c("lr_axis_id", "condition_value")) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  if (!"lr_axis_id" %in% colnames(df) || any(!nzchar(df$lr_axis_id))) {
    required <- c("ligand_human", "receptor_human", "source", "target")
    if (all(required %in% colnames(df))) {
      missing <- !nzchar(as.character(df$lr_axis_id))
      df$lr_axis_id[missing] <- standardize_lr_axis_id(
        df$ligand_human[missing],
        df$receptor_human[missing],
        df$source[missing],
        df$target[missing]
      )
    }
  }
  df <- df[nzchar(as.character(df$lr_axis_id)), , drop = FALSE]
  if (nrow(df) == 0) {
    return(data.frame(lr_axis_id = character(), condition_value = character(), stringsAsFactors = FALSE))
  }
  common_cols <- intersect(
    c("pair_id", "layer_id", "source", "target", "ligand", "receptor", "ligand_human", "receptor_human"),
    colnames(df)
  )
  keep <- unique(c("lr_axis_id", "condition_value", common_cols))
  if (identical(method, "cellchat")) {
    df$cellchat_hit <- TRUE
    if (!"method_evidence_class" %in% colnames(df)) df$method_evidence_class <- "hypothesis_only"
    if (!"can_be_primary" %in% colnames(df)) df$can_be_primary <- "no"
    keep <- unique(c(keep, intersect(c("prob", "pval", "pathway_name", "method_evidence_class", "can_be_primary"), colnames(df)), "cellchat_hit"))
  } else if (identical(method, "liana")) {
    if (!"liana_consensus_hit" %in% colnames(df)) df$liana_consensus_hit <- FALSE
    keep <- unique(c(keep, intersect(c("liana_consensus_hit", "n_methods_agreed", "liana_consensus_score", "status", "reason"), colnames(df))))
  } else if (identical(method, "nichenet")) {
    if (!"nichenet_hit" %in% colnames(df)) df$nichenet_hit <- as_logical_flag_07d(df$success %||% TRUE)
    keep <- unique(c(keep, intersect(c("nichenet_hit", "status", "reason", "method_evidence_class", "n_targets_in_receiver_de"), colnames(df))))
  }
  out <- df[, keep, drop = FALSE]
  rename <- setdiff(colnames(out), c("lr_axis_id", "condition_value"))
  colnames(out)[match(rename, colnames(out))] <- paste(method, rename, sep = "_")
  out[!duplicated(out[, c("lr_axis_id", "condition_value"), drop = FALSE]), , drop = FALSE]
}

coalesce_consensus_columns_07d <- function(df) {
  if (nrow(df) == 0) {
    return(empty_method_consensus_07d())
  }
  for (col in c("cellchat_cellchat_hit", "liana_liana_consensus_hit", "nichenet_nichenet_hit")) {
    if (!col %in% colnames(df)) df[[col]] <- FALSE
  }
  for (col in c("cellchat_method_evidence_class", "cellchat_can_be_primary", "liana_n_methods_agreed", "liana_liana_consensus_score", "nichenet_method_evidence_class", "nichenet_n_targets_in_receiver_de")) {
    if (!col %in% colnames(df)) df[[col]] <- ""
  }
  out <- df
  out$pair_id <- first_nonempty_07d(out$cellchat_pair_id %||% "", out$nichenet_pair_id %||% "", out$liana_pair_id %||% "")
  out$layer_id <- first_nonempty_07d(out$cellchat_layer_id %||% "", out$nichenet_layer_id %||% "", out$liana_layer_id %||% "")
  out$source <- first_nonempty_07d(out$cellchat_source %||% "", out$liana_source %||% "", out$nichenet_source %||% "")
  out$target <- first_nonempty_07d(out$cellchat_target %||% "", out$liana_target %||% "", out$nichenet_target %||% "")
  out$ligand <- first_nonempty_07d(out$cellchat_ligand %||% "", out$liana_ligand %||% "", out$nichenet_ligand %||% "")
  out$receptor <- first_nonempty_07d(out$cellchat_receptor %||% "", out$liana_receptor %||% "", out$nichenet_receptor %||% "")
  out$ligand_human <- first_nonempty_07d(out$cellchat_ligand_human %||% "", out$liana_ligand_human %||% "", out$nichenet_ligand_human %||% "")
  out$receptor_human <- first_nonempty_07d(out$cellchat_receptor_human %||% "", out$liana_receptor_human %||% "", out$nichenet_receptor_human %||% "")
  out$cellchat_hit <- as_logical_flag_07d(out$cellchat_cellchat_hit)
  out$liana_consensus_hit <- as_logical_flag_07d(out$liana_liana_consensus_hit)
  out$nichenet_hit <- as_logical_flag_07d(out$nichenet_nichenet_hit)
  out$cellchat_method_evidence_class <- out$cellchat_method_evidence_class %||% ""
  out$cellchat_can_be_primary <- out$cellchat_can_be_primary %||% ""
  out$liana_n_methods_agreed <- suppressWarnings(as.numeric(out$liana_n_methods_agreed %||% NA_real_))
  out$liana_consensus_score <- suppressWarnings(as.numeric(out$liana_liana_consensus_score %||% NA_real_))
  out$nichenet_method_evidence_class <- out$nichenet_method_evidence_class %||% ""
  out$nichenet_n_targets_in_receiver_de <- suppressWarnings(as.numeric(out$nichenet_n_targets_in_receiver_de %||% NA_real_))
  out$commot_spatial_hit <- FALSE
  out
}

full_outer_join_lr_tables <- function(cellchat_lr, liana_lr, multinichenet_lr) {
  method_tables <- list(
    method_prefixed_df_07d(cellchat_lr, "cellchat"),
    method_prefixed_df_07d(liana_lr, "liana"),
    method_prefixed_df_07d(multinichenet_lr, "nichenet")
  )
  method_tables <- method_tables[vapply(method_tables, nrow, integer(1)) > 0]
  if (length(method_tables) == 0) {
    return(empty_method_consensus_07d())
  }
  merged <- Reduce(
    function(left, right) merge(left, right, by = c("lr_axis_id", "condition_value"), all = TRUE, sort = FALSE),
    method_tables
  )
  coalesce_consensus_columns_07d(merged)
}

assign_evidence_tier_07d <- function(consensus_df, env_thresholds = NULL) {
  if (is.null(consensus_df) || nrow(consensus_df) == 0) {
    return(empty_method_consensus_07d())
  }
  df <- as.data.frame(consensus_df, stringsAsFactors = FALSE)
  env_thresholds <- env_thresholds %||% list()
  min_methods <- suppressWarnings(as.integer(env_thresholds$min_methods_for_primary %||% Sys.getenv("CONSENSUS_MIN_METHODS_FOR_PRIMARY", "2")))
  if (is.na(min_methods) || min_methods < 1L) min_methods <- 2L
  require_nichenet <- tolower(as.character(env_thresholds$require_nichenet_for_primary %||% Sys.getenv("CONSENSUS_REQUIRE_NICHENET_FOR_PRIMARY", "yes"))) %in% c("yes", "true", "1", "on")

  df$cellchat_hit <- as_logical_flag_07d(df$cellchat_hit %||% FALSE)
  df$liana_consensus_hit <- as_logical_flag_07d(df$liana_consensus_hit %||% FALSE)
  df$nichenet_hit <- as_logical_flag_07d(df$nichenet_hit %||% FALSE)
  df$commot_spatial_hit <- as_logical_flag_07d(df$commot_spatial_hit %||% FALSE)
  df$evidence_method_n <- as.integer(df$cellchat_hit) + as.integer(df$liana_consensus_hit) +
    as.integer(df$nichenet_hit) + as.integer(df$commot_spatial_hit)

  primary_ok <- df$liana_consensus_hit & df$evidence_method_n >= min_methods
  if (isTRUE(require_nichenet)) {
    primary_ok <- primary_ok & df$nichenet_hit
  }
  df$evidence_tier <- ifelse(
    primary_ok,
    "primary",
    ifelse(
      df$liana_consensus_hit | df$nichenet_hit | (df$cellchat_hit & df$evidence_method_n >= 2L),
      "exploratory",
      ifelse(df$cellchat_hit, "candidate", "blocked")
    )
  )
  cellchat_only <- df$cellchat_hit & !df$liana_consensus_hit & !df$nichenet_hit & !df$commot_spatial_hit
  df$evidence_tier[cellchat_only] <- "candidate"
  df$can_be_primary <- ifelse(df$evidence_tier == "primary", "yes", "no")
  df$evidence_reason <- ifelse(
    df$evidence_tier == "primary",
    "LIANA consensus plus downstream validation met primary thresholds",
    ifelse(
      df$evidence_tier == "exploratory",
      "At least one non-CellChat evidence layer supports the axis",
      ifelse(df$evidence_tier == "candidate", "CellChat-only or hypothesis-level evidence; not primary", "No supporting communication evidence")
    )
  )
  df
}

attach_commot_spatial <- function(consensus_df, commot_df) {
  if (is.null(commot_df) || nrow(commot_df) == 0) {
    consensus_df$commot_spatial_hit <- FALSE
    return(consensus_df)
  }
  stop("attach_commot_spatial: COMMOT integration is reserved for P-R03-G", call. = FALSE)
}
