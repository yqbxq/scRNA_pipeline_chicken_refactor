m3_scdesign3_question_map_cols <- c(
  "question_id", "section", "question_family", "status",
  "scdesign3_role", "validation_layer", "validation_object",
  "ground_truth_col", "simulation_type", "simulation_unit",
  "primary_metric", "secondary_metrics", "gate_target", "affects_modules",
  "required_before_core_interpretation", "derived_parent_question_id",
  "enabled", "reason", "materialized", "question_class", "scope",
  "contrast_axis", "tools_to_run", "question_zh"
)

m3_scdesign3_target_cols <- c(
  "target_id", "target_type", "layer_id", "input_object", "truth_col",
  "questions_covered", "n_simulations", "resolution_grid",
  "mixture_design", "primary_metric", "pass_threshold", "warn_threshold",
  "fail_threshold", "max_cells_per_label", "n_hvg", "n_pcs",
  "output_dir", "status"
)

m3_scdesign3_simulation_design_cols <- c(
  "simulation_design_id", "target_id", "simulation_type",
  "simulation_unit", "n_simulations", "resolution_grid", "mixture_design",
  "max_cells_per_label", "n_hvg", "n_pcs", "status", "notes"
)

m3_scdesign3_threshold_cols <- c(
  "threshold_id", "target_type", "primary_metric", "pass_threshold",
  "warn_threshold", "fail_threshold", "notes"
)

m3_scdesign3_simulation_override_cols <- c(
  "override_id", "target_id", "target_type", "n_simulations",
  "resolution_grid", "max_cells_per_label", "n_hvg", "n_pcs",
  "status", "notes"
)

m3_scdesign3_derived_questions <- data.frame(
  question_id = c(
    "F23_TC_GC_diff_overall",
    "F24_TC_GC_diff_subtype",
    "F25_GC_internal_diff",
    "F26_panorama_screen_diff"
  ),
  question_class = "communication",
  question_zh = c(
    "derived TC-GC overall differential communication",
    "derived TC-GC subtype differential communication",
    "derived GC internal differential communication",
    "derived panorama screen differential communication"
  ),
  scope = "panorama",
  sender_groups = "*",
  receiver_groups = "*",
  condition_split = "group_id:syf,f5",
  contrast_axis = "derived_communication",
  tools_to_run = "cellchat+nichenet",
  priority = "P1",
  status = "derived",
  depends_on = c("F02", "F04;F06", "F08;F17", "F15"),
  notes = "Generated scDesign3 derived question; not materialized in analysis_questions.tsv.",
  activation_policy = "derived_from_split",
  min_sender_cells = "",
  min_receiver_cells = "",
  min_cells_per_condition = "",
  fallback_pair_id = "",
  derived_from_pair_id = c("F02", "F04;F06", "F08;F17", "F15"),
  display_question_id = "",
  output_alias = "",
  report_title = "",
  run_baseline_if_split_fails = "no",
  materialized = "no",
  stringsAsFactors = FALSE
)

m3_scdesign3_question_code <- function(question_id) {
  sub("_.*$", "", m3_trim(question_id))
}

m3_scdesign3_section <- function(question_id) {
  substr(m3_trim(question_id), 1L, 1L)
}

m3_scdesign3_question_family <- function(q) {
  class <- m3_trim(q$question_class[[1]])
  axis <- m3_trim(q$contrast_axis[[1]])
  section <- m3_scdesign3_section(q$question_id[[1]])
  if (class %in% c("marker", "deg", "composition", "communication", "regulation", "trajectory", "spatial", "qc", "global_context")) {
    return(switch(
      class,
      marker = "marker",
      deg = "DEG",
      composition = "composition",
      communication = "communication",
      regulation = "regulation",
      trajectory = "trajectory",
      spatial = "spatial",
      qc = "QC",
      global_context = "global_context",
      class
    ))
  }
  switch(
    section,
    A = "marker",
    B = "DEG",
    C = "DEG",
    D = "DEG",
    E = "composition",
    F = "communication",
    G = "regulation",
    H = "trajectory",
    I = "spatial",
    ifelse(nzchar(axis), axis, "unknown")
  )
}

m3_scdesign3_explicit_role <- function(code) {
  direct <- c(
    "A01", "A02", "B02", "C03", "C04", "D03", "E01",
    "F07", "F08", "F13", "F16", "F17",
    "H01", "H02", "I08", "I09", "I11"
  )
  upstream <- c(
    "A04", "B01", "C01", "C02", "D01", "D02", "D06", "E04",
    "F01", "F02", "F03", "F04", "F05", "F06", "F09", "F10",
    "F14", "F15", "G01", "G02", "G04", "G06", "G07", "G08",
    "I01", "I02", "I03", "I04", "I05", "I12", "I15", "I16",
    "I17", "I21"
  )
  downstream <- c("H05", "H06", "I06", "I07", "I13", "I14", "I18", "I19", "I20")
  planned <- c(
    "A03", "B03", "C05", "C06", "D04", "E02",
    "F11", "F12", "F18", "F19", "F20", "F21", "F22",
    "G03", "G05", "G09", "H03", "H04", "H07", "H08", "I10"
  )
  not_applicable <- c("D05", "E03")
  derived <- c("F23", "F24", "F25", "F26")

  if (code %in% derived) return("derived_from_validated_parent")
  if (code %in% not_applicable) return("not_applicable")
  if (code %in% planned) return("planned_waiting_input")
  if (code %in% direct) return("direct_validate")
  if (code %in% upstream) return("upstream_gate")
  if (code %in% downstream) return("downstream_consistency")
  ""
}

m3_scdesign3_infer_role <- function(q, family) {
  code <- m3_scdesign3_question_code(q$question_id[[1]])
  role <- m3_scdesign3_explicit_role(code)
  if (nzchar(role)) {
    return(role)
  }
  if (identical(m3_trim(q$status[[1]]), "derived")) {
    return("derived_from_validated_parent")
  }
  if (identical(family, "QC") || identical(family, "global_context")) {
    return("not_applicable")
  }
  if (identical(q$scope[[1]], "ST_section")) {
    return("upstream_gate")
  }
  if (identical(family, "marker")) {
    return("direct_validate")
  }
  if (identical(family, "communication")) {
    return("upstream_gate")
  }
  if (identical(family, "trajectory")) {
    return("downstream_consistency")
  }
  "upstream_gate"
}

m3_scdesign3_validation_layer <- function(q) {
  scope <- m3_trim(q$scope[[1]])
  code <- m3_scdesign3_question_code(q$question_id[[1]])
  if (scope == "ST_section" || startsWith(code, "I")) return("ST")
  if (scope == "GC_subcluster") return("GC_subcluster")
  if (scope == "TC_subcluster") return("TC_subcluster")
  if (code %in% c("A02", "B02", "C03", "C04", "D03", "E01", "F07", "F08", "F09", "F10", "F13", "F16", "F17", "G02", "G04", "G08", "H01", "H02", "H05", "H06")) {
    return("GC_subcluster")
  }
  if (code %in% c("A03", "B03", "C05", "C06", "D04", "E02", "F18", "F19", "F20", "F21", "F22", "G03", "G05", "G09", "H03", "H04", "H07", "H08")) {
    return("TC_subcluster")
  }
  "panorama"
}

m3_scdesign3_ground_truth_col <- function(layer) {
  switch(
    layer,
    ST = "region_label",
    GC_subcluster = "cell_subtype",
    TC_subcluster = "cell_subtype",
    panorama = "cell_type",
    "cell_type"
  )
}

m3_scdesign3_target_type <- function(role, family, layer, code) {
  if (role %in% c("not_applicable", "derived_from_validated_parent")) {
    return("")
  }
  if (layer == "ST") {
    if (code %in% c("I08", "I09", "I10", "I11")) {
      return("deconv_validation")
    }
    return("spatial_validation")
  }
  if (identical(family, "composition")) {
    return("composition_robustness")
  }
  if (identical(family, "trajectory")) {
    return("trajectory_robustness")
  }
  if (identical(family, "communication")) {
    return("communication_robustness")
  }
  "cluster_robustness"
}

m3_scdesign3_simulation_type <- function(target_type, family, layer) {
  if (target_type == "deconv_validation") return("synthetic_spots")
  if (target_type == "spatial_validation") return("synthetic_spots")
  if (target_type == "composition_robustness") return("synthetic_cells_composition")
  if (target_type == "trajectory_robustness") return("synthetic_cells_trajectory")
  if (target_type == "communication_robustness") return("perturbation")
  "synthetic_cells"
}

m3_scdesign3_primary_metric <- function(target_type) {
  switch(
    target_type,
    deconv_validation = "RMSE",
    spatial_validation = "consistency_score",
    composition_robustness = "composition_RMSE",
    trajectory_robustness = "trajectory_order_accuracy",
    communication_robustness = "sender_receiver_recovery",
    cluster_robustness = "ARI",
    "consistency_score"
  )
}

m3_scdesign3_secondary_metrics <- function(target_type) {
  switch(
    target_type,
    deconv_validation = "correlation;MAE;dominant_cell_type_accuracy",
    spatial_validation = "region_recovery;neighborhood_consistency;gradient_correlation",
    composition_robustness = "composition_MAE;dominant_subtype_accuracy",
    trajectory_robustness = "root_recovery;terminal_recovery;pseudotime_correlation;branch_consistency",
    communication_robustness = "sender_recovery;receiver_recovery;split_gate_status",
    cluster_robustness = "NMI;Jaccard;purity",
    "NMI;Jaccard"
  )
}

m3_scdesign3_gate_target <- function(family, target_type) {
  if (target_type %in% c("deconv_validation", "spatial_validation")) return("spatial_support_validation")
  if (family == "communication") return("communication_interpretation")
  if (family == "regulation") return("regulation_interpretation")
  if (family == "trajectory") return("trajectory_interpretation")
  if (family %in% c("DEG", "marker", "composition")) return("cell_identity_interpretation")
  "context"
}

m3_scdesign3_affects_modules <- function(section) {
  switch(
    section,
    A = "03;04;05;06;07",
    B = "05;06;07",
    C = "05;06;07",
    D = "05;06;07;08",
    E = "05;06",
    F = "07",
    G = "08",
    H = "09;10",
    I = "ST;05;06;07;08;09",
    ""
  )
}

m3_scdesign3_target_id <- function(target_type, layer, family, code) {
  if (!nzchar(target_type)) return("")
  if (target_type == "deconv_validation") {
    if (code == "I08") return("SCD_DECONV_PANORAMA")
    if (code == "I09") return("SCD_DECONV_GC_SUBTYPE")
    if (code == "I10") return("SCD_DECONV_TC_SUBTYPE")
    return("SCD_DECONV_VALIDATION")
  }
  if (target_type == "spatial_validation") return("SCD_SPATIAL_SUPPORT")
  if (target_type == "composition_robustness") {
    if (layer == "GC_subcluster") return("SCD_COMPOSITION_GC_SUBTYPE")
    if (layer == "TC_subcluster") return("SCD_COMPOSITION_TC_SUBTYPE")
    return("SCD_COMPOSITION_PANORAMA")
  }
  if (target_type == "trajectory_robustness") {
    if (layer == "GC_subcluster") return("SCD_TRAJECTORY_GC")
    if (layer == "TC_subcluster") return("SCD_TRAJECTORY_TC")
    return("SCD_TRAJECTORY_PANORAMA")
  }
  if (target_type == "communication_robustness") {
    if (layer == "GC_subcluster") return("SCD_COMM_GC_INTERNAL")
    if (layer == "TC_subcluster") return("SCD_COMM_TC_INTERNAL")
    return("SCD_COMM_PANORAMA")
  }
  if (layer == "GC_subcluster") return("SCD_CLUSTER_GC_SUBTYPE")
  if (layer == "TC_subcluster") return("SCD_CLUSTER_TC_SUBTYPE")
  "SCD_CLUSTER_PANORAMA"
}

m3_scdesign3_input_object <- function(layer) {
  switch(
    layer,
    panorama = "manifest:03d_annotated_object",
    GC_subcluster = "layer_status:GC_subcluster",
    TC_subcluster = "layer_status:TC_subcluster",
    ST = "ST-E2:synthetic_spots",
    ""
  )
}

m3_scdesign3_thresholds <- function(target_type, metric) {
  if (target_type == "cluster_robustness") {
    return(c(pass = "ARI>=0.80", warn = "ARI>=0.60;NMI>=0.60;min_Jaccard>=0.45", fail = "ARI<0.60"))
  }
  if (target_type == "composition_robustness") {
    return(c(pass = "RMSE<=0.10;MAE<=0.10", warn = "RMSE<=0.20", fail = "RMSE>0.20"))
  }
  if (target_type == "trajectory_robustness") {
    return(c(pass = "order_accuracy>=0.80;pseudotime_correlation>=0.70", warn = "order_accuracy>=0.60", fail = "order_accuracy<0.60"))
  }
  if (target_type == "communication_robustness") {
    return(c(pass = "sender_receiver_recovery=PASS;split_gate!=FAIL", warn = "sender_receiver_recovery=WARN", fail = "sender_receiver_recovery=FAIL"))
  }
  if (target_type == "deconv_validation") {
    return(c(pass = "RMSE<=0.10;correlation>=0.80", warn = "RMSE<=0.20;correlation>=0.60", fail = "RMSE>0.20"))
  }
  c(pass = "consistency_score>=0.80", warn = "consistency_score>=0.60", fail = "consistency_score<0.60")
}

m3_scdesign3_add_derived_questions <- function(questions) {
  if (!"materialized" %in% colnames(questions)) {
    questions$materialized <- "yes"
  }
  missing_cols <- setdiff(colnames(questions), colnames(m3_scdesign3_derived_questions))
  derived <- m3_scdesign3_derived_questions
  for (col in missing_cols) {
    derived[[col]] <- ""
  }
  extra_cols <- setdiff(colnames(derived), colnames(questions))
  for (col in extra_cols) {
    questions[[col]] <- ""
  }
  derived <- derived[, colnames(questions), drop = FALSE]
  existing <- unique(questions$question_id)
  rbind(questions, derived[!derived$question_id %in% existing, , drop = FALSE])
}

m3_build_scdesign3_question_map <- function(questions) {
  if (nrow(questions) == 0) {
    return(m3_empty_df(m3_scdesign3_question_map_cols))
  }
  questions <- m3_scdesign3_add_derived_questions(questions)
  rows <- lapply(seq_len(nrow(questions)), function(idx) {
    q <- questions[idx, , drop = FALSE]
    qid <- m3_trim(q$question_id[[1]])
    code <- m3_scdesign3_question_code(qid)
    section <- m3_scdesign3_section(qid)
    family <- m3_scdesign3_question_family(q)
    role <- m3_scdesign3_infer_role(q, family)
    layer <- m3_scdesign3_validation_layer(q)
    target_type <- m3_scdesign3_target_type(role, family, layer, code)
    primary_metric <- m3_scdesign3_primary_metric(target_type)
    q_status <- m3_trim(q$status[[1]], "active")
    enabled <- if (role == "not_applicable") {
      "no"
    } else if (role == "planned_waiting_input" || q_status == "planned") {
      "planned"
    } else {
      "yes"
    }
    reason <- switch(
      role,
      direct_validate = "scDesign3 directly validates the question object or core assumption.",
      upstream_gate = "Question interpretation depends on an upstream scDesign3 identity, subtype, region, or reference gate.",
      downstream_consistency = "scDesign3 provides a downstream consistency check, not primary discovery evidence.",
      not_applicable = "QC or global context question; no scDesign3 biological gate required.",
      planned_waiting_input = "Waiting for TC subtype, ST, velocity, or split-specific prerequisite input.",
      derived_from_validated_parent = "Derived communication question; inherits validation from parent split questions.",
      "Default scDesign3 gate assignment generated from question metadata."
    )
    data.frame(
      question_id = qid,
      section = section,
      question_family = family,
      status = q_status,
      scdesign3_role = role,
      validation_layer = layer,
      validation_object = m3_scdesign3_target_id(target_type, layer, family, code),
      ground_truth_col = m3_scdesign3_ground_truth_col(layer),
      simulation_type = m3_scdesign3_simulation_type(target_type, family, layer),
      simulation_unit = ifelse(layer == "ST", "spot", "cell"),
      primary_metric = primary_metric,
      secondary_metrics = m3_scdesign3_secondary_metrics(target_type),
      gate_target = m3_scdesign3_gate_target(family, target_type),
      affects_modules = m3_scdesign3_affects_modules(section),
      required_before_core_interpretation = ifelse(role %in% c("direct_validate", "upstream_gate", "downstream_consistency"), "yes", "no"),
      derived_parent_question_id = ifelse(role == "derived_from_validated_parent", m3_trim(q$derived_from_pair_id[[1]], m3_trim(q$depends_on[[1]])), ""),
      enabled = enabled,
      reason = reason,
      materialized = m3_trim(q$materialized[[1]], "yes"),
      question_class = m3_trim(q$question_class[[1]]),
      scope = m3_trim(q$scope[[1]]),
      contrast_axis = m3_trim(q$contrast_axis[[1]]),
      tools_to_run = m3_trim(q$tools_to_run[[1]]),
      question_zh = m3_trim(q$question_zh[[1]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[, m3_scdesign3_question_map_cols, drop = FALSE]
}

m3_scdesign3_thresholds_from_table <- function(threshold_table, target_type, metric) {
  if (!is.null(threshold_table) && nrow(threshold_table) > 0 && "target_type" %in% colnames(threshold_table)) {
    hit <- threshold_table[threshold_table$target_type == target_type, , drop = FALSE]
    if (nrow(hit) > 0) {
      return(c(
        pass = m3_trim(hit$pass_threshold[[1]]),
        warn = m3_trim(hit$warn_threshold[[1]]),
        fail = m3_trim(hit$fail_threshold[[1]])
      ))
    }
  }
  m3_scdesign3_thresholds(target_type, metric)
}

m3_scdesign3_default_simulation_settings <- function(target_type, layer_id) {
  list(
    n_simulations = "5",
    resolution_grid = ifelse(layer_id == "ST", "", "0.2,0.4,0.6,0.8,1.0,1.2"),
    max_cells_per_label = "2000",
    n_hvg = "2000",
    n_pcs = "30"
  )
}

m3_normalize_scdesign3_simulation_overrides <- function(overrides = NULL) {
  if (is.null(overrides) || nrow(overrides) == 0) {
    return(m3_empty_df(m3_scdesign3_simulation_override_cols))
  }
  for (col in m3_scdesign3_simulation_override_cols) {
    if (!col %in% colnames(overrides)) {
      overrides[[col]] <- ""
    }
  }
  overrides <- overrides[, m3_scdesign3_simulation_override_cols, drop = FALSE]
  overrides$status <- ifelse(nzchar(m3_trim(overrides$status)), m3_trim(overrides$status), "active")
  overrides
}

m3_scdesign3_simulation_override_for_target <- function(overrides, target_id, target_type) {
  overrides <- m3_normalize_scdesign3_simulation_overrides(overrides)
  if (nrow(overrides) == 0) {
    return(NULL)
  }
  active <- overrides[!tolower(overrides$status) %in% c("no", "false", "disabled", "inactive"), , drop = FALSE]
  if (nrow(active) == 0) {
    return(NULL)
  }
  target_hit <- active[nzchar(active$target_id) & active$target_id == target_id, , drop = FALSE]
  if (nrow(target_hit) > 0) {
    return(target_hit[1L, , drop = FALSE])
  }
  type_hit <- active[!nzchar(active$target_id) & active$target_type == target_type, , drop = FALSE]
  if (nrow(type_hit) > 0) {
    return(type_hit[1L, , drop = FALSE])
  }
  NULL
}

m3_scdesign3_override_value <- function(override, col, default = "") {
  if (is.null(override) || !col %in% colnames(override)) {
    return(default)
  }
  value <- m3_trim(override[[col]][[1]])
  if (nzchar(value)) value else default
}

m3_build_scdesign3_targets <- function(question_map, thresholds_override = NULL, simulation_overrides = NULL) {
  if (nrow(question_map) == 0) {
    return(m3_empty_df(m3_scdesign3_target_cols))
  }
  runnable <- question_map[
    !question_map$scdesign3_role %in% c("not_applicable", "derived_from_validated_parent") &
      nzchar(question_map$validation_object),
    ,
    drop = FALSE
  ]
  if (nrow(runnable) == 0) {
    return(m3_empty_df(m3_scdesign3_target_cols))
  }
  threshold_table <- m3_build_scdesign3_thresholds(thresholds_override)
  rows <- list()
  for (target_id in sort(unique(runnable$validation_object))) {
    hit <- runnable[runnable$validation_object == target_id, , drop = FALSE]
    target_type <- m3_scdesign3_target_type(hit$scdesign3_role[[1]], hit$question_family[[1]], hit$validation_layer[[1]], m3_scdesign3_question_code(hit$question_id[[1]]))
    metric <- m3_scdesign3_primary_metric(target_type)
    thresholds <- m3_scdesign3_thresholds_from_table(threshold_table, target_type, metric)
    sim_defaults <- m3_scdesign3_default_simulation_settings(target_type, hit$validation_layer[[1]])
    sim_override <- m3_scdesign3_simulation_override_for_target(simulation_overrides, target_id, target_type)
    status <- if (all(hit$enabled == "planned")) "planned" else "active"
    rows[[length(rows) + 1L]] <- data.frame(
      target_id = target_id,
      target_type = target_type,
      layer_id = hit$validation_layer[[1]],
      input_object = m3_scdesign3_input_object(hit$validation_layer[[1]]),
      truth_col = hit$ground_truth_col[[1]],
      questions_covered = paste(sort(unique(hit$question_id)), collapse = ";"),
      n_simulations = m3_scdesign3_override_value(sim_override, "n_simulations", sim_defaults$n_simulations),
      resolution_grid = m3_scdesign3_override_value(sim_override, "resolution_grid", sim_defaults$resolution_grid),
      mixture_design = ifelse(target_type == "composition_robustness", "balanced;observed;perturbed", "observed_balanced"),
      primary_metric = metric,
      pass_threshold = thresholds[["pass"]],
      warn_threshold = thresholds[["warn"]],
      fail_threshold = thresholds[["fail"]],
      max_cells_per_label = m3_scdesign3_override_value(sim_override, "max_cells_per_label", sim_defaults$max_cells_per_label),
      n_hvg = m3_scdesign3_override_value(sim_override, "n_hvg", sim_defaults$n_hvg),
      n_pcs = m3_scdesign3_override_value(sim_override, "n_pcs", sim_defaults$n_pcs),
      output_dir = file.path("results", "tables", "04d_cluster_robustness", target_id),
      status = status,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out[, m3_scdesign3_target_cols, drop = FALSE]
}

m3_build_scdesign3_simulation_designs <- function(targets) {
  if (nrow(targets) == 0) {
    return(m3_empty_df(m3_scdesign3_simulation_design_cols))
  }
  for (col in c("max_cells_per_label", "n_hvg", "n_pcs")) {
    if (!col %in% colnames(targets)) {
      targets[[col]] <- ""
    }
  }
  rows <- lapply(seq_len(nrow(targets)), function(idx) {
    target <- targets[idx, , drop = FALSE]
    defaults <- m3_scdesign3_default_simulation_settings(target$target_type[[1]], target$layer_id[[1]])
    data.frame(
      simulation_design_id = paste0(target$target_id[[1]], "_DESIGN"),
      target_id = target$target_id[[1]],
      simulation_type = switch(
        target$target_type[[1]],
        deconv_validation = "synthetic_spots",
        spatial_validation = "synthetic_spots",
        composition_robustness = "synthetic_cells_composition",
        trajectory_robustness = "synthetic_cells_trajectory",
        communication_robustness = "perturbation",
        "synthetic_cells"
      ),
      simulation_unit = ifelse(target$layer_id[[1]] == "ST", "spot", "cell"),
      n_simulations = target$n_simulations[[1]],
      resolution_grid = target$resolution_grid[[1]],
      mixture_design = target$mixture_design[[1]],
      max_cells_per_label = m3_trim(target$max_cells_per_label[[1]], defaults$max_cells_per_label),
      n_hvg = m3_trim(target$n_hvg[[1]], defaults$n_hvg),
      n_pcs = m3_trim(target$n_pcs[[1]], defaults$n_pcs),
      status = target$status[[1]],
      notes = "Generated from scdesign3_targets.tsv; concrete engines may refine parameters per target.",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[, m3_scdesign3_simulation_design_cols, drop = FALSE]
}

m3_build_scdesign3_thresholds <- function(thresholds_override = NULL) {
  target_types <- c(
    "cluster_robustness", "composition_robustness", "trajectory_robustness",
    "communication_robustness", "deconv_validation", "spatial_validation"
  )
  rows <- lapply(target_types, function(target_type) {
    metric <- m3_scdesign3_primary_metric(target_type)
    thresholds <- m3_scdesign3_thresholds(target_type, metric)
    data.frame(
      threshold_id = paste0("SCD_THRESHOLD_", toupper(target_type)),
      target_type = target_type,
      primary_metric = metric,
      pass_threshold = thresholds[["pass"]],
      warn_threshold = thresholds[["warn"]],
      fail_threshold = thresholds[["fail"]],
      notes = "Default metadata threshold; project config may override during runtime.",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (!is.null(thresholds_override) && nrow(thresholds_override) > 0) {
    extra_cols <- setdiff(colnames(thresholds_override), m3_scdesign3_threshold_cols)
    for (col in extra_cols) {
      if (!col %in% colnames(out)) {
        out[[col]] <- ""
      }
    }
    for (col in m3_scdesign3_threshold_cols) {
      if (!col %in% colnames(thresholds_override)) {
        thresholds_override[[col]] <- ""
      }
    }
    thresholds_override <- thresholds_override[, c(m3_scdesign3_threshold_cols, extra_cols), drop = FALSE]
    for (idx in seq_len(nrow(thresholds_override))) {
      target_type <- m3_trim(thresholds_override$target_type[[idx]])
      if (!nzchar(target_type)) {
        next
      }
      hit <- which(out$target_type == target_type)
      if (length(hit) == 0) {
        for (col in colnames(out)) {
          if (!col %in% colnames(thresholds_override)) {
            thresholds_override[[col]] <- ""
          }
        }
        out <- rbind(out, thresholds_override[idx, colnames(out), drop = FALSE])
      } else {
        for (col in c("primary_metric", "pass_threshold", "warn_threshold", "fail_threshold", "notes", extra_cols)) {
          value <- m3_trim(thresholds_override[[col]][[idx]])
          if (nzchar(value)) {
            out[[col]][[hit[[1]]]] <- value
          }
        }
        override_id <- m3_trim(thresholds_override$threshold_id[[idx]])
        if (nzchar(override_id)) {
          out$threshold_id[[hit[[1]]]] <- override_id
        }
      }
    }
  }
  out[, c(m3_scdesign3_threshold_cols, setdiff(colnames(out), m3_scdesign3_threshold_cols)), drop = FALSE]
}
