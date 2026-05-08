m3_communication_cols <- c(
  "pair_id", "source_question_id", "layer_scope", "sender", "receiver",
  "condition_split_var", "condition_split_values", "tool", "communication_mode",
  "activation_policy", "min_sender_cells", "min_receiver_cells",
  "min_cells_per_condition", "fallback_pair_id", "derived_from_pair_id",
  "run_baseline_if_split_fails", "requires_all_derived_inputs_pass",
  "receiver_gene_program_source", "baseline_marker_comparison_id",
  "receiver_deg_comparison_id", "direction_filter", "requires_cell_subtype",
  "notes", "enabled"
)

m3_gc_subtypes <- c("pGC", "eGC", "rgGC", "lGC")

m3_comm_tool <- function(q) {
  mode <- m3_tool_mode(q$tools_to_run, c("cellchat", "nichenet"), both_label = "both")
  if (nzchar(mode)) mode else "cellchat"
}

m3_comm_enabled <- function(q) {
  if (identical(q$status[[1]], "active")) "yes" else "no"
}

m3_comm_mode <- function(q, tool = "") {
  if (identical(tool, "differential")) {
    return("differential_summary")
  }
  split <- m3_condition_fields(q$condition_split)
  if (nzchar(split$var)) "condition_split" else "baseline"
}

m3_comm_question_field <- function(q, col, default = "") {
  if (col %in% colnames(q)) {
    return(m3_trim(q[[col]][[1]], default))
  }
  default
}

m3_comm_activation_policy <- function(q, tool = "", default = NULL) {
  if (!is.null(default) && nzchar(default)) {
    return(default)
  }
  if (identical(tool, "differential")) {
    return("derived_from_split")
  }
  policy <- m3_comm_question_field(q, "activation_policy", "always")
  if (nzchar(policy)) policy else "always"
}

m3_comm_resolve_fallback_pair_id <- function(q, pair_id, fallback_pair_id) {
  fallback_pair_id <- m3_trim(fallback_pair_id)
  if (!nzchar(fallback_pair_id)) {
    return("")
  }
  if (grepl("__", fallback_pair_id, fixed = TRUE)) {
    return(fallback_pair_id)
  }
  qid <- q$question_id[[1]]
  if (startsWith(pair_id, paste0(qid, "__"))) {
    suffix <- sub(paste0("^", qid, "__"), "", pair_id)
    return(paste(fallback_pair_id, suffix, sep = "__"))
  }
  fallback_pair_id
}

m3_comm_requires_cell_subtype <- function(q, sender, receiver) {
  if (q$question_id[[1]] %in% c("F14_panorama_screen_baseline", "F15_panorama_screen_split", "F18_TC_screen_baseline", "F19_TC_screen_split")) {
    return("no")
  }
  if (q$scope[[1]] %in% c("GC_subcluster")) {
    return("yes")
  }
  tokens <- unlist(strsplit(paste(sender, receiver, sep = ","), ",", fixed = TRUE), use.names = FALSE)
  if (any(tokens %in% c(m3_gc_subtypes, "TC_*"))) "yes" else "no"
}

m3_comm_direction_filter <- function(q, tool, pair_id = "") {
  if (tool %in% c("differential")) {
    return(ifelse(pair_id %in% c("F26_panorama_screen_diff"), "no", "yes"))
  }
  if (q$question_id[[1]] %in% c(
    "F09_GC_dev_skip", "F10_GC_feedback", "F13_GC_subtype_pairwise",
    "F14_panorama_screen_baseline", "F15_panorama_screen_split",
    "F16_GC_screen_baseline", "F17_GC_screen_split",
    "F18_TC_screen_baseline", "F19_TC_screen_split"
  )) {
    return("no")
  }
  "yes"
}

m3_marker_for_receiver <- function(receiver) {
  receiver <- m3_trim(receiver)
  if (!nzchar(receiver) || receiver %in% c("*", "TC_*")) {
    return("")
  }
  tokens <- unlist(strsplit(receiver, ",", fixed = TRUE), use.names = FALSE)
  tokens <- m3_trim(tokens)
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) == 1L && tokens %in% m3_gc_subtypes) {
    return(sprintf("A02_GC_subtype_marker__%s_vs_rest", tokens))
  }
  if (length(tokens) > 1L && all(tokens %in% m3_gc_subtypes)) {
    return(paste(sprintf("A02_GC_subtype_marker__%s_vs_rest", tokens), collapse = ","))
  }
  if (identical(receiver, "TC")) {
    return("A01_panorama_marker__TC_vs_rest")
  }
  if (identical(receiver, "GC") || all(tokens %in% m3_gc_subtypes)) {
    return("A01_panorama_marker__GC_vs_rest")
  }
  ""
}

m3_condition_deg_for_receiver <- function(receiver) {
  receiver <- m3_trim(receiver)
  tokens <- unlist(strsplit(receiver, ",", fixed = TRUE), use.names = FALSE)
  tokens <- m3_trim(tokens)
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) == 1L && tokens %in% m3_gc_subtypes) {
    return(sprintf("D03_GC_subtype_stage_DEG__%s_syf_vs_f5", tokens))
  }
  if (length(tokens) > 1L && all(tokens %in% m3_gc_subtypes)) {
    return(paste(sprintf("D03_GC_subtype_stage_DEG__%s_syf_vs_f5", tokens), collapse = ","))
  }
  if (identical(receiver, "TC")) {
    return("D02_TC_stage_DEG__TC_syf_vs_f5")
  }
  if (identical(receiver, "GC") || all(tokens %in% m3_gc_subtypes)) {
    return("D01_GC_stage_DEG__GC_syf_vs_f5")
  }
  ""
}

m3_comm_row <- function(q, pair_id, sender, receiver, tool, notes = "",
                        baseline_marker_comparison_id = NULL,
                        receiver_deg_comparison_id = NULL,
                        receiver_gene_program_source = NULL,
                        direction_filter = NULL,
                        requires_cell_subtype = NULL,
                        activation_policy = NULL,
                        fallback_pair_id = NULL,
                        derived_from_pair_id = NULL,
                        requires_all_derived_inputs_pass = NULL) {
  split <- m3_condition_fields(q$condition_split)
  mode <- m3_comm_mode(q, tool)
  policy <- m3_comm_activation_policy(q, tool, activation_policy)
  fallback <- fallback_pair_id %||% m3_comm_question_field(q, "fallback_pair_id")
  fallback <- m3_comm_resolve_fallback_pair_id(q, pair_id, fallback)
  run_baseline_if_split_fails <- m3_comm_question_field(q, "run_baseline_if_split_fails")
  if (!nzchar(run_baseline_if_split_fails) && nzchar(fallback)) {
    run_baseline_if_split_fails <- "yes"
  }
  if (identical(policy, "derived_from_split")) {
    run_baseline_if_split_fails <- "no"
  }
  min_sender <- if (identical(policy, "derived_from_split")) "" else m3_comm_question_field(q, "min_sender_cells")
  min_receiver <- if (identical(policy, "derived_from_split")) "" else m3_comm_question_field(q, "min_receiver_cells")
  min_condition <- if (identical(policy, "derived_from_split")) "" else m3_comm_question_field(q, "min_cells_per_condition")
  baseline <- baseline_marker_comparison_id %||% m3_marker_for_receiver(receiver)
  deg <- receiver_deg_comparison_id %||% if (identical(mode, "condition_split")) m3_condition_deg_for_receiver(receiver) else ""
  source <- receiver_gene_program_source %||% if (identical(mode, "condition_split") && nzchar(deg)) {
    "condition_deg"
  } else if (nzchar(baseline)) {
    "receiver_marker"
  } else {
    "none"
  }
  list(
    pair_id = pair_id,
    source_question_id = q$question_id,
    layer_scope = q$scope,
    sender = sender,
    receiver = receiver,
    condition_split_var = split$var,
    condition_split_values = split$values,
    tool = tool,
    communication_mode = mode,
    activation_policy = policy,
    min_sender_cells = min_sender,
    min_receiver_cells = min_receiver,
    min_cells_per_condition = min_condition,
    fallback_pair_id = fallback,
    derived_from_pair_id = derived_from_pair_id %||% m3_comm_question_field(q, "derived_from_pair_id"),
    run_baseline_if_split_fails = run_baseline_if_split_fails,
    requires_all_derived_inputs_pass = requires_all_derived_inputs_pass %||% "no",
    receiver_gene_program_source = source,
    baseline_marker_comparison_id = baseline,
    receiver_deg_comparison_id = deg,
    direction_filter = direction_filter %||% m3_comm_direction_filter(q, tool, pair_id),
    requires_cell_subtype = requires_cell_subtype %||% m3_comm_requires_cell_subtype(q, sender, receiver),
    enabled = m3_comm_enabled(q),
    notes = notes
  )
}

m3_add_differential_comm <- function(rows, q) {
  split <- m3_condition_fields(q$condition_split)
  if (!nzchar(split$var)) {
    return(rows)
  }
  derived <- switch(
    q$question_id[[1]],
    F02_TC_GC_bidir_split = list(id = "F23_TC_GC_diff_overall", sender = "TC", receiver = m3_collapse(m3_gc_subtypes), notes = "derived F23; compare TC-GC communication between syf and f5", direction = "yes"),
    F04_TC_to_GC_each_split = list(id = "F24_TC_to_GC_subtype_diff", sender = "TC", receiver = m3_collapse(m3_gc_subtypes), notes = "derived F24; compare TC-to-GC-subtype communication between syf and f5", direction = "yes"),
    F08_GC_dev_seq_split = list(id = "F25_GC_internal_diff", sender = m3_collapse(m3_gc_subtypes), receiver = m3_collapse(m3_gc_subtypes), notes = "derived F25; compare GC internal developmental communication between syf and f5", direction = "yes"),
    F15_panorama_screen_split = list(id = "F26_panorama_screen_diff", sender = "TC,pGC,eGC,rgGC,lGC", receiver = "TC,pGC,eGC,rgGC,lGC", notes = "derived F26; compare full panorama communication network between syf and f5", direction = "no"),
    NULL
  )
  if (is.null(derived)) {
    return(rows)
  }
  source_pair_ids <- vapply(rows, function(row) m3_trim(row$pair_id), character(1))
  source_pair_ids <- source_pair_ids[nzchar(source_pair_ids)]
  rows[[length(rows) + 1L]] <- m3_comm_row(
    q, derived$id, derived$sender, derived$receiver, "differential",
    notes = derived$notes,
    baseline_marker_comparison_id = "",
    receiver_deg_comparison_id = "",
    receiver_gene_program_source = "none",
    direction_filter = derived$direction,
    requires_cell_subtype = "yes",
    activation_policy = "derived_from_split",
    fallback_pair_id = "",
    derived_from_pair_id = paste(source_pair_ids, collapse = ","),
    requires_all_derived_inputs_pass = ifelse(identical(derived$id, "F25_GC_internal_diff"), "yes", "no")
  )
  rows
}

m3_fanout_communication <- function(q) {
  qid <- q$question_id[[1]]
  tool <- m3_comm_tool(q)
  rows <- list()

  if (qid %in% c("F01_TC_GC_bidir_baseline", "F02_TC_GC_bidir_split")) {
    rows[[length(rows) + 1L]] <- m3_comm_row(q, paste(qid, "TC_to_GC", sep = "__"), "TC", m3_collapse(m3_gc_subtypes), tool, notes = ifelse(grepl("split", qid), "TC->GC split; GC represented by GC subtypes via M5", "TC->GC baseline; GC represented by GC subtypes via M5"), baseline_marker_comparison_id = "A01_panorama_marker__GC_vs_rest", receiver_deg_comparison_id = ifelse(grepl("split", qid), "D01_GC_stage_DEG__GC_syf_vs_f5", ""))
    rows[[length(rows) + 1L]] <- m3_comm_row(q, paste(qid, "GC_to_TC", sep = "__"), m3_collapse(m3_gc_subtypes), "TC", tool, notes = ifelse(grepl("split", qid), "GC->TC split; GC represented by GC subtypes via M5", "GC->TC baseline; GC represented by GC subtypes via M5"), baseline_marker_comparison_id = "A01_panorama_marker__TC_vs_rest")
  } else if (qid %in% c("F03_TC_to_GC_each_baseline", "F04_TC_to_GC_each_split")) {
    for (receiver in m3_gc_subtypes) {
      rows[[length(rows) + 1L]] <- m3_comm_row(q, paste(qid, paste("TC", "to", receiver, sep = "_"), sep = "__"), "TC", receiver, tool, notes = ifelse(receiver == "pGC" && qid == "F04_TC_to_GC_each_split", "核心", ""))
    }
  } else if (qid %in% c("F05_GC_each_to_TC_baseline", "F06_GC_each_to_TC_split")) {
    for (sender in m3_gc_subtypes) {
      rows[[length(rows) + 1L]] <- m3_comm_row(q, paste(qid, paste(sender, "to", "TC", sep = "_"), sep = "__"), sender, "TC", tool, notes = ifelse(sender == "pGC" && qid == "F05_GC_each_to_TC_baseline", "反向", ""))
    }
  } else if (qid %in% c("F07_GC_dev_seq_baseline", "F08_GC_dev_seq_split")) {
    for (idx in seq_len(length(m3_gc_subtypes) - 1L)) {
      sender <- m3_gc_subtypes[[idx]]
      receiver <- m3_gc_subtypes[[idx + 1L]]
      rows[[length(rows) + 1L]] <- m3_comm_row(q, paste(qid, paste(sender, "to", receiver, sep = "_"), sep = "__"), sender, receiver, tool, notes = ifelse(qid == "F07_GC_dev_seq_baseline", sprintf("发育流seg%s", idx), ""))
    }
  } else if (qid %in% c("F09_GC_dev_skip", "F10_GC_feedback")) {
    sender_sets <- m3_group_sets(q$sender_groups, q$scope)
    receiver_sets <- m3_group_sets(q$receiver_groups, q$scope)
    for (sender_set in sender_sets) {
      for (receiver_set in receiver_sets) {
        sender <- m3_collapse(sender_set)
        receiver <- m3_collapse(receiver_set)
        rows[[length(rows) + 1L]] <- m3_comm_row(q, paste(qid, paste(sender, "to", receiver, sep = "_"), sep = "__"), sender, receiver, tool, notes = ifelse(qid == "F09_GC_dev_skip", "跨段", "反馈"))
      }
    }
  } else if (qid == "F13_GC_subtype_pairwise") {
    rows[[length(rows) + 1L]] <- m3_comm_row(q, "F13_GC_subtype_pairwise__GC_pairwise", m3_collapse(m3_gc_subtypes), "*", tool, notes = "12对", baseline_marker_comparison_id = "", receiver_deg_comparison_id = "", receiver_gene_program_source = "none", direction_filter = "no")
  } else if (qid %in% c("F14_panorama_screen_baseline", "F15_panorama_screen_split")) {
    rows[[length(rows) + 1L]] <- m3_comm_row(q, paste(qid, "panorama_full", sep = "__"), "*", "*", tool, notes = ifelse(qid == "F14_panorama_screen_baseline", "无偏筛查", "重连图谱"), baseline_marker_comparison_id = "", receiver_deg_comparison_id = "", receiver_gene_program_source = "none", direction_filter = "no", requires_cell_subtype = "no")
  } else if (qid %in% c("F16_GC_screen_baseline", "F17_GC_screen_split")) {
    rows[[length(rows) + 1L]] <- m3_comm_row(q, paste(qid, "GC_full", sep = "__"), "*", "*", tool, notes = ifelse(qid == "F17_GC_screen_split", "细胞数风险", ""), baseline_marker_comparison_id = "", receiver_deg_comparison_id = "", receiver_gene_program_source = "none", direction_filter = "no", requires_cell_subtype = "yes")
  } else if (qid %in% c("F18_TC_screen_baseline", "F19_TC_screen_split")) {
    rows[[length(rows) + 1L]] <- m3_comm_row(q, paste(qid, "TC_full", sep = "__"), "*", "*", tool, notes = "等TC", baseline_marker_comparison_id = "", receiver_deg_comparison_id = "", receiver_gene_program_source = "none", direction_filter = "no", requires_cell_subtype = "no")
  } else if (qid == "F20_TC_sub_to_GC_baseline") {
    rows[[length(rows) + 1L]] <- m3_comm_row(q, "F20_TC_sub_to_GC_baseline__TCsub_to_GC", "TC_*", m3_collapse(m3_gc_subtypes), tool, notes = "等TC")
  } else if (qid == "F21_TC_sub_to_GC_split") {
    rows[[length(rows) + 1L]] <- m3_comm_row(q, "F21_TC_sub_to_GC_split__TCsub_to_GC", "TC_*", m3_collapse(m3_gc_subtypes), tool, notes = "等TC")
  } else if (qid == "F22_TC_sub_to_GC_each") {
    for (receiver in m3_gc_subtypes) {
      rows[[length(rows) + 1L]] <- m3_comm_row(q, paste(qid, paste("TCsub", "to", receiver, sep = "_"), sep = "__"), "TC_*", receiver, tool, notes = "等TC")
    }
  }

  rows <- m3_add_differential_comm(rows, q)
  m3_bind_rows(rows, m3_communication_cols)
}
