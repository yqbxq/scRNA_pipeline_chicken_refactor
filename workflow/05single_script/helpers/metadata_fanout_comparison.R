m3_comparison_cols <- c(
  "comparison_id", "source_question_id", "layer_scope", "contrast_axis",
  "analysis_mode", "analysis_unit", "stat_level",
  "group_var", "ident_1", "ident_2", "subset_column", "subset_value",
  "aggregation_group_var", "composition_group_var",
  "batch_var", "enabled", "min_biological_replicates", "force_exploratory",
  "min_cells_per_group", "logfc_threshold",
  "produces_gene_program", "gene_program_role", "notes"
)

m3_comparison_row <- function(q, suffix, ident_1, ident_2, group_var = "cell_subtype",
                              subset_column = "", subset_value = "",
                              analysis_mode = "subtype_pairwise",
                              analysis_unit = "whole_layer",
                              stat_level = "cell_level_exploratory",
                              aggregation_group_var = "",
                              composition_group_var = "",
                              produces_gene_program = "yes",
                              gene_program_role = "subtype_pairwise_deg",
                              notes = "") {
  list(
    comparison_id = paste(q$question_id, m3_safe_id(suffix), sep = "__"),
    source_question_id = q$question_id,
    layer_scope = q$scope,
    contrast_axis = q$contrast_axis,
    analysis_mode = analysis_mode,
    analysis_unit = analysis_unit,
    stat_level = stat_level,
    aggregation_group_var = aggregation_group_var,
    composition_group_var = composition_group_var,
    produces_gene_program = produces_gene_program,
    gene_program_role = gene_program_role,
    group_var = group_var,
    ident_1 = ident_1,
    ident_2 = ident_2,
    subset_column = subset_column,
    subset_value = subset_value,
    batch_var = "batch",
    enabled = "yes",
    min_biological_replicates = "2",
    force_exploratory = "no",
    min_cells_per_group = "3",
    logfc_threshold = "0",
    notes = notes
  )
}

m3_uses_cell_type <- function(groups) {
  groups <- unique(m3_trim(groups))
  groups <- groups[nzchar(groups)]
  length(groups) > 0 && all(groups %in% c("GC", "TC"))
}

m3_identity_group_var <- function(groups) {
  if (m3_uses_cell_type(groups)) "cell_type" else "cell_subtype"
}

m3_marker_note <- function(group_var) {
  if (identical(group_var, "cell_type")) {
    "broad cell-type marker; use cell_type, not cell_subtype"
  } else {
    "GC subtype one-vs-rest marker"
  }
}

m3_pairwise_note <- function(qid, group_var, subset_value = "") {
  if (identical(group_var, "cell_type")) {
    if (nzchar(subset_value)) {
      return(sprintf("broad GC-vs-TC identity comparison within %s; use cell_type", subset_value))
    }
    return("broad GC-vs-TC identity comparison; use cell_type")
  }
  if (grepl("_syf$", qid)) return("GC subtype pairwise within syf")
  if (grepl("_f5$", qid)) return("GC subtype pairwise within f5")
  "GC subtype pairwise identity comparison"
}

m3_stage_subset <- function(group) {
  if (group %in% c("all_cells")) {
    return(list(subset_column = "", subset_value = "", analysis_unit = "whole_layer", aggregation_group_var = "all_cells", notes = "global stage signature across all cells"))
  }
  if (group %in% c("GC", "TC")) {
    return(list(
      subset_column = "cell_type",
      subset_value = group,
      analysis_unit = "within_cell_type",
      aggregation_group_var = "cell_type",
      notes = sprintf("formal DEG priority; %s cells only", group)
    ))
  }
  list(
    subset_column = "cell_subtype",
    subset_value = group,
    analysis_unit = "within_cell_subtype",
    aggregation_group_var = "cell_subtype",
    notes = sprintf("formal DEG priority; %s only", group)
  )
}

m3_composition_fields <- function(q, groups) {
  if (identical(q$question_id[[1]], "E01_GC_subtype_compo")) {
    return(list(
      suffix = "composition_GC_subtypes",
      subset_column = "cell_subtype",
      subset_value = "pGC,eGC,rgGC,lGC",
      composition_group_var = "cell_subtype",
      notes = "GC subtype composition; sample-level proportions"
    ))
  }
  if (identical(q$question_id[[1]], "E03_layer_compo")) {
    return(list(
      suffix = "composition_GC_TC",
      subset_column = "cell_type",
      subset_value = "GC,TC",
      composition_group_var = "cell_type",
      notes = "broad GC-vs-TC composition; sample-level proportions"
    ))
  }
  subset_column <- if (length(groups) == 0 || any(groups %in% c("all_cells"))) "" else m3_identity_group_var(groups)
  subset_value <- if (nzchar(subset_column)) m3_collapse(groups) else ""
  list(
    suffix = if (nzchar(subset_value)) paste("composition", m3_safe_id(subset_value), sep = "_") else "composition_all",
    subset_column = subset_column,
    subset_value = subset_value,
    composition_group_var = if (m3_uses_cell_type(groups)) "cell_type" else "cell_subtype",
    notes = "composition_only"
  )
}

m3_fanout_comparison <- function(q) {
  axis <- q$contrast_axis
  rows <- list()
  condition_subset <- m3_condition_subset(q$condition_split)

  if (axis == "cluster_marker") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    group_var <- m3_identity_group_var(groups)
    for (group in groups) {
      rows[[length(rows) + 1L]] <- m3_comparison_row(
        q, paste0(group, "_vs_rest"), group, "__rest__",
        group_var = group_var,
        analysis_mode = "subtype_marker",
        analysis_unit = "whole_layer",
        stat_level = "cell_level_exploratory",
        aggregation_group_var = "",
        composition_group_var = "",
        produces_gene_program = "yes",
        gene_program_role = "receiver_marker",
        notes = m3_marker_note(group_var)
      )
    }
  } else if (axis == "directional_DEG") {
    senders <- m3_flat_groups(q$sender_groups, q$scope)
    receivers <- m3_flat_groups(q$receiver_groups, q$scope)
    group_var <- m3_identity_group_var(c(senders, receivers))
    for (sender in senders) {
      for (receiver in receivers) {
        rows[[length(rows) + 1L]] <- m3_comparison_row(
          q, paste(sender, "vs", receiver, condition_subset$subset_value, sep = "_"),
          sender, receiver,
          group_var = group_var,
          subset_column = condition_subset$subset_column,
          subset_value = condition_subset$subset_value,
          analysis_mode = "subtype_pairwise",
          analysis_unit = "whole_layer",
          stat_level = "cell_level_exploratory",
          aggregation_group_var = "",
          composition_group_var = "",
          produces_gene_program = "yes",
          gene_program_role = "subtype_pairwise_deg",
          notes = m3_pairwise_note(q$question_id[[1]], group_var, condition_subset$subset_value)
        )
      }
    }
  } else if (axis == "pairwise") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    group_var <- m3_identity_group_var(groups)
    if (length(groups) >= 2) {
      combos <- utils::combn(groups, 2, simplify = FALSE)
      for (combo in combos) {
        rows[[length(rows) + 1L]] <- m3_comparison_row(
          q, paste(combo[[1]], "vs", combo[[2]], condition_subset$subset_value, sep = "_"),
          combo[[1]], combo[[2]],
          group_var = group_var,
          subset_column = condition_subset$subset_column,
          subset_value = condition_subset$subset_value,
          analysis_mode = "subtype_pairwise",
          analysis_unit = "whole_layer",
          stat_level = "cell_level_exploratory",
          aggregation_group_var = "",
          composition_group_var = "",
          produces_gene_program = "yes",
          gene_program_role = "subtype_pairwise_deg",
          notes = m3_pairwise_note(q$question_id[[1]], group_var, condition_subset$subset_value)
        )
      }
    }
  } else if (axis == "contrast_only") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    condition <- m3_condition_pair(q$condition_split)
    if (length(groups) == 0) {
      groups <- "all_cells"
    }
    for (group in groups) {
      stage <- m3_stage_subset(group)
      rows[[length(rows) + 1L]] <- m3_comparison_row(
        q, paste(group, condition$ident_1, "vs", condition$ident_2, sep = "_"),
        condition$ident_1, condition$ident_2,
        group_var = condition$group_var,
        subset_column = stage$subset_column,
        subset_value = stage$subset_value,
        analysis_mode = "condition_within_type",
        analysis_unit = stage$analysis_unit,
        stat_level = "pseudobulk_formal_if_replicates",
        aggregation_group_var = stage$aggregation_group_var,
        composition_group_var = "",
        produces_gene_program = "yes",
        gene_program_role = "condition_deg",
        notes = stage$notes
      )
    }
  } else if (axis == "composition") {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    condition <- m3_condition_pair(q$condition_split)
    fields <- m3_composition_fields(q, groups)
    rows[[length(rows) + 1L]] <- m3_comparison_row(
      q, fields$suffix,
      condition$ident_1, condition$ident_2,
      group_var = condition$group_var,
      subset_column = fields$subset_column,
      subset_value = fields$subset_value,
      analysis_mode = "composition",
      analysis_unit = "sample_level",
      stat_level = "composition_formal_if_replicates",
      aggregation_group_var = "",
      composition_group_var = fields$composition_group_var,
      produces_gene_program = "no",
      gene_program_role = "none",
      notes = fields$notes
    )
  }

  m3_bind_rows(rows, m3_comparison_cols)
}
