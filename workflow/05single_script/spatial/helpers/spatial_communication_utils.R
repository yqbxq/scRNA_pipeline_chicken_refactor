st08_empty_candidates <- function() {
  data.frame(
    comm_candidate_id = character(),
    source_question_id = character(),
    lr_axis_id = character(),
    ligand = character(),
    receptor = character(),
    sender_cell_type = character(),
    receiver_cell_type = character(),
    stage = character(),
    section_id = character(),
    region_id = character(),
    niche_id = character(),
    pair_id = character(),
    condition_value = character(),
    scrna_evidence_tier = character(),
    liana_support = character(),
    multinichenet_support = character(),
    receiver_deg_support = character(),
    cellchat_hypothesis = character(),
    status = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

st08_empty_colocalization <- function() {
  data.frame(
    comm_candidate_id = character(),
    lr_axis_id = character(),
    section_id = character(),
    sender_cell_type = character(),
    receiver_cell_type = character(),
    sender_spatial_abundance = numeric(),
    receiver_spatial_abundance = numeric(),
    sender_receiver_colocalization = numeric(),
    colocalization_support = character(),
    ligand_spot_expression_mean = numeric(),
    receptor_spot_expression_mean = numeric(),
    lr_expression_colocalization = numeric(),
    lr_expression_support = character(),
    status = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

st08_empty_neighborhood <- function() {
  data.frame(
    comm_candidate_id = character(),
    lr_axis_id = character(),
    section_id = character(),
    sender_cell_type = character(),
    receiver_cell_type = character(),
    neighborhood_enrichment_score = numeric(),
    co_occurrence_score = numeric(),
    neighborhood_support = character(),
    status = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

st08_empty_consensus <- function() {
  data.frame(
    comm_candidate_id = character(),
    source_question_id = character(),
    lr_axis_id = character(),
    ligand = character(),
    receptor = character(),
    sender_cell_type = character(),
    receiver_cell_type = character(),
    stage = character(),
    section_id = character(),
    region_id = character(),
    niche_id = character(),
    pair_id = character(),
    condition_value = character(),
    scrna_evidence_tier = character(),
    liana_support = character(),
    multinichenet_support = character(),
    receiver_deg_support = character(),
    cellchat_hypothesis = character(),
    deconv_support_status = character(),
    deconv_support_level = character(),
    scrna_support_level = character(),
    sender_spatial_abundance = numeric(),
    receiver_spatial_abundance = numeric(),
    sender_receiver_colocalization = numeric(),
    colocalization_support = character(),
    ligand_spot_expression_mean = numeric(),
    receptor_spot_expression_mean = numeric(),
    lr_expression_colocalization = numeric(),
    lr_expression_support = character(),
    neighborhood_enrichment_score = numeric(),
    co_occurrence_score = numeric(),
    neighborhood_support = character(),
    commot_score = numeric(),
    commot_support = character(),
    commot_status = character(),
    commot_reason = character(),
    commot_n_spot_pairs = numeric(),
    commot_distance_threshold = numeric(),
    svg_spatial_support_level = character(),
    svg_ligand_support = character(),
    svg_receptor_support = character(),
    svg_receiver_target_support = character(),
    svg_support_reason = character(),
    spatial_evidence_tier = character(),
    final_interpretation_level = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

st08_scalar <- function(x, default = "") {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]])) return(default)
  value <- trimws(as.character(x[[1]]))
  if (!nzchar(value)) default else value
}

st08_bool_text <- function(x) {
  ifelse(tolower(as.character(x)) %in% c("yes", "true", "1", "on", "ok", "primary"), "yes", "no")
}

st08_first_col <- function(df, candidates, default = "") {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) rep(default, nrow(df)) else as.character(df[[hit]])
}

st08_make_candidate_id <- function(lr_axis_id, section_id, condition_value) {
  raw <- paste(lr_axis_id, section_id, condition_value, sep = "|")
  gsub("[^A-Za-z0-9_.-]+", "_", raw)
}

st08_normalize_candidates <- function(consensus_df, sections = data.frame(stringsAsFactors = FALSE)) {
  if (is.null(consensus_df) || nrow(consensus_df) == 0) return(st08_empty_candidates())
  df <- as.data.frame(consensus_df, stringsAsFactors = FALSE)
  for (col in c("lr_axis_id", "ligand", "receptor", "source", "target", "pair_id", "condition_value", "evidence_tier", "can_be_primary")) {
    if (!col %in% colnames(df)) df[[col]] <- ""
  }
  df$sender_cell_type <- st08_first_col(df, c("sender_cell_type", "source", "sender"))
  df$receiver_cell_type <- st08_first_col(df, c("receiver_cell_type", "target", "receiver"))
  df$scrna_evidence_tier <- st08_first_col(df, c("scrna_evidence_tier", "evidence_tier"), "candidate")
  df$cellchat_hypothesis <- st08_bool_text(st08_first_col(df, c("cellchat_hit", "cellchat_hypothesis"), "no"))
  df$liana_support <- st08_bool_text(st08_first_col(df, c("liana_consensus_hit", "liana_support"), "no"))
  df$multinichenet_support <- st08_bool_text(st08_first_col(df, c("multinichenet_hit", "nichenet_hit", "multinichenet_support"), "no"))
  df$receiver_deg_support <- st08_bool_text(suppressWarnings(as.numeric(st08_first_col(df, c("nichenet_n_targets_in_receiver_de", "receiver_deg_target_n"), "0"))) > 0)
  df <- df[df$scrna_evidence_tier != "blocked" & nzchar(df$lr_axis_id) & nzchar(df$ligand) & nzchar(df$receptor), , drop = FALSE]
  if (nrow(df) == 0) return(st08_empty_candidates())

  section_ids <- "all"
  if (nrow(sections) > 0 && "section_id" %in% colnames(sections)) {
    section_ids <- unique(as.character(sections$section_id[nzchar(as.character(sections$section_id))]))
    if (length(section_ids) == 0) section_ids <- "all"
  }
  rows <- lapply(section_ids, function(section_id) {
    out <- data.frame(
      comm_candidate_id = st08_make_candidate_id(df$lr_axis_id, section_id, df$condition_value),
      source_question_id = st08_first_col(df, c("source_question_id", "question_id"), ""),
      lr_axis_id = df$lr_axis_id,
      ligand = df$ligand,
      receptor = df$receptor,
      sender_cell_type = df$sender_cell_type,
      receiver_cell_type = df$receiver_cell_type,
      stage = st08_first_col(df, c("stage", "condition_value"), ""),
      section_id = section_id,
      region_id = st08_first_col(df, c("region_id", "layer_id"), ""),
      niche_id = st08_first_col(df, c("niche_id"), ""),
      pair_id = df$pair_id,
      condition_value = df$condition_value,
      scrna_evidence_tier = df$scrna_evidence_tier,
      liana_support = df$liana_support,
      multinichenet_support = df$multinichenet_support,
      receiver_deg_support = df$receiver_deg_support,
      cellchat_hypothesis = df$cellchat_hypothesis,
      status = "candidate",
      reason = "",
      stringsAsFactors = FALSE
    )
    out[!duplicated(out[, c("comm_candidate_id"), drop = FALSE]), , drop = FALSE]
  })
  do.call(rbind, rows)
}

st08_collect_deconv_props <- function(cfg) {
  methods <- collect_deconv_method_manifests_st(cfg)
  ok_rows <- methods[methods$status == "ok" & nzchar(methods$proportion_tsv), , drop = FALSE]
  if (nrow(ok_rows) == 0) {
    return(data.frame(method = character(), section = character(), spot_id = character(), cell_type = character(), proportion = numeric(), stringsAsFactors = FALSE))
  }
  rows <- lapply(seq_len(nrow(ok_rows)), function(i) {
    row <- ok_rows[i, , drop = FALSE]
    df <- st07_read_tsv(st07_abs_path(row$proportion_tsv[[1]], cfg))
    if (nrow(df) == 0 || !all(c("spot_id", "cell_type", "proportion") %in% colnames(df))) {
      return(data.frame(method = character(), section = character(), spot_id = character(), cell_type = character(), proportion = numeric(), stringsAsFactors = FALSE))
    }
    df$proportion <- suppressWarnings(as.numeric(df$proportion))
    data.frame(
      method = row$method[[1]],
      section = row$section[[1]] %||% "all",
      spot_id = df$spot_id,
      cell_type = df$cell_type,
      proportion = df$proportion,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

st08_deconv_validation_status <- function(cfg) {
  manifest_tsv <- file.path(cfg$spatial_deconv_validation_table_dir, "validation_manifest.tsv")
  df <- st07_read_tsv(manifest_tsv)
  if (nrow(df) == 0 || !"status" %in% colnames(df)) return("missing")
  st08_scalar(df$status, "missing")
}

st08_deconv_support_level <- function(status) {
  status <- tolower(st08_scalar(status, "missing"))
  if (status %in% c("ok", "pass")) return("primary")
  if (status %in% c("warn")) return("supporting")
  if (status %in% c("ok_smoke")) return("smoke_only")
  "blocked"
}

st08_scrna_support_level <- function(df) {
  tier <- tolower(as.character(df$scrna_evidence_tier %||% "candidate"))
  bool_vec <- function(x) tolower(as.character(x %||% "no")) %in% c("yes", "true", "1", "ok", "support", "supported", "primary")
  liana <- bool_vec(df$liana_support)
  multinichenet <- bool_vec(df$multinichenet_support)
  receiver_deg <- bool_vec(df$receiver_deg_support)
  support_n <- as.integer(liana) + as.integer(multinichenet) + as.integer(receiver_deg)
  cellchat <- bool_vec(df$cellchat_hypothesis)
  out <- rep("blocked", length(tier))
  out[tier %in% c("primary") & support_n >= 2] <- "primary"
  out[tier %in% c("primary", "supporting", "exploratory") & support_n >= 1 & out == "blocked"] <- "supporting"
  out[(cellchat | tier %in% c("candidate", "hypothesis_only", "exploratory")) & out == "blocked"] <- "hypothesis_only"
  out
}

st08_gate_row <- function(question_id, module, gate_status, interpretation_allowed, reason, evidence_tier = "") {
  data.frame(
    question_id = question_id,
    module = module,
    gate_status = gate_status,
    interpretation_allowed = interpretation_allowed,
    evidence_tier = evidence_tier,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

st08_support_status <- function(value) {
  tolower(st08_scalar(value, "no")) %in% c("yes", "true", "1", "ok", "support", "supported")
}
