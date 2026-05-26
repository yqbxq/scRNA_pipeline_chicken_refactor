m3_comparison_cols <- c(
  "comparison_id", "source_question_id", "display_question_id",
  "output_alias", "report_title", "layer_scope", "contrast_axis",
  "analysis_mode", "analysis_unit", "stat_level",
  "group_var", "ident_1", "ident_2", "subset_column", "subset_value",
  "aggregation_group_var", "composition_group_var",
  "batch_var", "enabled", "min_biological_replicates", "force_exploratory",
  "min_cells_per_group", "logfc_threshold",
  "produces_gene_program", "gene_program_role", "notes",
  "min_cells_override", "min_samples_override", "total_umi_override",
  "single_sample_frac_override"
)

m3_annotation_marker_cols <- c(
  "target_id", "source_question_id", "layer_scope", "object_layer",
  "cluster_column", "annotation_label_column", "group_var",
  "ident_1", "ident_2", "analysis_mode", "gene_program_role",
  "output_dir", "annotation_only", "enabled", "notes"
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
  comparison_id <- paste(q$question_id, m3_safe_id(suffix), sep = "__")
  display_question_id <- m3_question_scalar(q, "display_question_id", q$question_id[[1]])
  output_alias <- m3_question_scalar(q, "output_alias", display_question_id)
  report_title <- m3_question_scalar(q, "report_title", m3_question_scalar(q, "question_zh", display_question_id))
  list(
    comparison_id = comparison_id,
    source_question_id = q$question_id,
    display_question_id = display_question_id,
    output_alias = output_alias,
    report_title = report_title,
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
    notes = notes,
    min_cells_override = "-",
    min_samples_override = "-",
    total_umi_override = "-",
    single_sample_frac_override = "-"
  )
}

m3_question_scalar <- function(q, col, default = "") {
  if (!col %in% colnames(q)) {
    return(default)
  }
  value <- m3_trim(q[[col]][[1]], default)
  if (!nzchar(value)) default else value
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
    "post-annotation broad cell-type identity marker; not raw cluster annotation evidence"
  } else {
    "post-annotation subtype identity marker; not raw cluster annotation evidence"
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
      notes = "QC-only broad GC-vs-TC capture balance; no biological abundance conclusion"
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

m3_annotation_marker_target_row <- function(target_id, source_question_id, layer_scope, object_layer,
                                            cluster_column, annotation_label_column,
                                            output_dir, notes) {
  list(
    target_id = target_id,
    source_question_id = source_question_id,
    layer_scope = layer_scope,
    object_layer = object_layer,
    cluster_column = cluster_column,
    annotation_label_column = annotation_label_column,
    group_var = cluster_column,
    ident_1 = "__cluster__",
    ident_2 = "__rest__",
    analysis_mode = "annotation_cluster_marker",
    gene_program_role = "annotation_marker",
    output_dir = output_dir,
    annotation_only = "yes",
    enabled = "yes",
    notes = notes
  )
}

m3_fanout_annotation_marker_targets <- function(questions) {
  if (nrow(questions) == 0) {
    return(m3_empty_df(m3_annotation_marker_cols))
  }
  active_ids <- questions$question_id[questions$status == "active"]
  rows <- list()
  if ("A01_panorama_marker" %in% active_ids) {
    rows[[length(rows) + 1L]] <- m3_annotation_marker_target_row(
      target_id = "ANN01_panorama_cluster_marker",
      source_question_id = "A01_panorama_marker",
      layer_scope = "panorama",
      object_layer = "panorama",
      cluster_column = "panorama_cluster",
      annotation_label_column = "cell_type",
      output_dir = "annotation/layers/panorama",
      notes = "03d raw panorama cluster marker evidence for annotation only"
    )
  }
  if ("A02_GC_subtype_marker" %in% active_ids) {
    rows[[length(rows) + 1L]] <- m3_annotation_marker_target_row(
      target_id = "ANN02_GC_subcluster_cluster_marker",
      source_question_id = "A02_GC_subtype_marker",
      layer_scope = "GC_subcluster",
      object_layer = "GC_subcluster",
      cluster_column = "GC_subcluster_cluster",
      annotation_label_column = "cell_subtype",
      output_dir = "annotation/layers/GC_subcluster",
      notes = "04b raw GC subcluster marker evidence for annotation only"
    )
  }
  m3_bind_rows(rows, m3_annotation_marker_cols)
}

m3_fanout_comparison <- function(q) {
  axis <- q$contrast_axis
  rows <- list()
  condition_subset <- m3_condition_subset(q$condition_split)

  if (axis %in% c("identity_marker", "cluster_marker")) {
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
  } else if (axis %in% c("contrast_only", "global_stage_context")) {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    condition <- m3_condition_pair(q$condition_split)
    if (length(groups) == 0) {
      groups <- "all_cells"
    }
    for (group in groups) {
      stage <- m3_stage_subset(group)
      is_global_context <- identical(axis, "global_stage_context") || identical(q$question_class[[1]], "global_context")
      rows[[length(rows) + 1L]] <- m3_comparison_row(
        q, paste(group, condition$ident_1, "vs", condition$ident_2, sep = "_"),
        condition$ident_1, condition$ident_2,
        group_var = condition$group_var,
        subset_column = stage$subset_column,
        subset_value = stage$subset_value,
        analysis_mode = ifelse(is_global_context, "global_context", "condition_within_type"),
        analysis_unit = stage$analysis_unit,
        stat_level = "pseudobulk_formal_if_replicates",
        aggregation_group_var = stage$aggregation_group_var,
        composition_group_var = "",
        produces_gene_program = "yes",
        gene_program_role = ifelse(is_global_context, "global_context", "condition_deg"),
        notes = stage$notes
      )
    }
  } else if (axis %in% c("composition", "qc_composition")) {
    groups <- m3_flat_groups(q$sender_groups, q$scope)
    condition <- m3_condition_pair(q$condition_split)
    fields <- m3_composition_fields(q, groups)
    is_qc <- identical(axis, "qc_composition") || identical(q$question_class[[1]], "qc")
    rows[[length(rows) + 1L]] <- m3_comparison_row(
      q, fields$suffix,
      condition$ident_1, condition$ident_2,
      group_var = condition$group_var,
      subset_column = fields$subset_column,
      subset_value = fields$subset_value,
      analysis_mode = ifelse(is_qc, "qc_composition", "composition"),
      analysis_unit = "sample_level",
      stat_level = "composition_formal_if_replicates",
      aggregation_group_var = "",
      composition_group_var = fields$composition_group_var,
      produces_gene_program = "no",
      gene_program_role = ifelse(is_qc, "qc_only", "none"),
      notes = fields$notes
    )
  }

  m3_bind_rows(rows, m3_comparison_cols)
}
