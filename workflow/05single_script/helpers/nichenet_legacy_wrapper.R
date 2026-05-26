legacy_nichenet_mode_label_07c <- function() {
  "nichenet_legacy"
}

attach_receiver_de_target_count_07c <- function(lr_table, ligand_target_matrix = NULL, receiver_de_genes = character(), top_n = 20L) {
  if (nrow(lr_table) == 0 || is.null(ligand_target_matrix) || length(receiver_de_genes) == 0) {
    lr_table$n_targets_in_receiver_de <- NA_integer_
    return(lr_table)
  }
  ligand_col <- intersect(c("ligand", "test_ligand", "ligand_human"), colnames(lr_table))[1]
  if (is.na(ligand_col)) {
    lr_table$n_targets_in_receiver_de <- NA_integer_
    return(lr_table)
  }
  lr_table$n_targets_in_receiver_de <- vapply(lr_table[[ligand_col]], function(ligand) {
    if (!ligand %in% rownames(ligand_target_matrix)) {
      return(NA_integer_)
    }
    scores <- ligand_target_matrix[ligand, ]
    top_targets <- names(sort(scores, decreasing = TRUE))[seq_len(min(top_n, length(scores)))]
    sum(top_targets %in% receiver_de_genes)
  }, integer(1))
  lr_table
}
