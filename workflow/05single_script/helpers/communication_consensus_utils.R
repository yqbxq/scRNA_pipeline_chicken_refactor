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

full_outer_join_lr_tables <- function(cellchat_lr, liana_lr, multinichenet_lr) {
  stop("full_outer_join_lr_tables: not yet implemented (P-R03-E)", call. = FALSE)
}

assign_evidence_tier_07d <- function(consensus_df, env_thresholds = NULL) {
  stop("assign_evidence_tier_07d: not yet implemented (P-R03-E)", call. = FALSE)
}

attach_commot_spatial <- function(consensus_df, commot_df) {
  stop("attach_commot_spatial: not yet implemented (P-R03-G)", call. = FALSE)
}
