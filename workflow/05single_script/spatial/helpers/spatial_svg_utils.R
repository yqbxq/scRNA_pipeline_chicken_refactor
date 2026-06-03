st09_svg_method_columns <- function() {
  c(
    "section_id", "condition", "method", "gene_id", "gene_symbol", "score", "pvalue", "qvalue", "rank",
    "mean_expression", "n_spots", "n_genes_tested", "status", "reason",
    "is_ligand_candidate", "is_receptor_candidate", "is_receiver_target_candidate",
    "is_region_marker_candidate", "is_regulator_target_candidate", "linked_question_ids"
  )
}

st09_empty_method_df <- function() {
  out <- as.data.frame(setNames(rep(list(character()), length(st09_svg_method_columns())), st09_svg_method_columns()), stringsAsFactors = FALSE)
  out$score <- numeric()
  out$pvalue <- numeric()
  out$qvalue <- numeric()
  out$rank <- integer()
  out$mean_expression <- numeric()
  out$n_spots <- integer()
  out$n_genes_tested <- integer()
  out
}

st09_empty_method_row <- function(section_id = "all", condition = "all", method = "unknown", status = "skipped", reason = "") {
  data.frame(
    section_id = section_id,
    condition = condition,
    method = method,
    gene_id = "",
    gene_symbol = "",
    score = NA_real_,
    pvalue = NA_real_,
    qvalue = NA_real_,
    rank = NA_integer_,
    mean_expression = NA_real_,
    n_spots = 0L,
    n_genes_tested = 0L,
    status = status,
    reason = reason,
    is_ligand_candidate = "no",
    is_receptor_candidate = "no",
    is_receiver_target_candidate = "no",
    is_region_marker_candidate = "no",
    is_regulator_target_candidate = "no",
    linked_question_ids = "",
    stringsAsFactors = FALSE
  )
}

st09_as_yes_no <- function(value) {
  value <- tolower(trimws(as.character(value)))
  ifelse(value %in% c("yes", "true", "1", "on", "ok", "supported"), "yes", "no")
}

st09_read_tsv <- function(path) {
  if (exists("spatial_read_tsv", mode = "function")) {
    return(spatial_read_tsv(path))
  }
  if (!file.exists(path) || file.info(path)$size == 0) return(data.frame(stringsAsFactors = FALSE))
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
}

st09_write_tsv <- function(df, path) {
  if (exists("spatial_write_tsv", mode = "function")) {
    return(spatial_write_tsv(df, path))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

st09_pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) "" else hit
}

st09_scalar <- function(x, default = "") {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]])) return(default)
  value <- trimws(as.character(x[[1]]))
  if (!nzchar(value)) default else value
}

st09_normalize_gene <- function(value) {
  toupper(trimws(as.character(value)))
}

st09_normalize_method_df <- function(df, default_method = "unknown") {
  if (is.null(df) || nrow(df) == 0) return(st09_empty_method_df())
  out <- as.data.frame(df, stringsAsFactors = FALSE)
  for (col in st09_svg_method_columns()) {
    if (!col %in% colnames(out)) {
      out[[col]] <- if (col %in% c("score", "pvalue", "qvalue", "mean_expression")) NA_real_ else if (col %in% c("rank", "n_spots", "n_genes_tested")) NA_integer_ else ""
    }
  }
  out$method[!nzchar(trimws(out$method))] <- default_method
  for (col in c("score", "pvalue", "qvalue", "mean_expression")) out[[col]] <- suppressWarnings(as.numeric(out[[col]]))
  for (col in c("rank", "n_spots", "n_genes_tested")) out[[col]] <- suppressWarnings(as.integer(out[[col]]))
  for (col in c("is_ligand_candidate", "is_receptor_candidate", "is_receiver_target_candidate", "is_region_marker_candidate", "is_regulator_target_candidate")) {
    out[[col]] <- st09_as_yes_no(out[[col]])
  }
  out[, st09_svg_method_columns(), drop = FALSE]
}

st09_method_status_class <- function(status) {
  status <- tolower(trimws(as.character(status)))
  ifelse(status %in% c("ok"), "formal",
    ifelse(status %in% c("ok_proxy"), "proxy",
      ifelse(grepl("^skipped|^failed", status), "unusable", "unusable")
    )
  )
}

st09_collect_method_outputs <- function(svg_root) {
  paths <- c(
    Sys.getenv("SPATIALDE2_SVG_TSV", unset = ""),
    list.files(file.path(svg_root, "spatialde2"), pattern = "_svg\\.tsv$", full.names = TRUE),
    Sys.getenv("SPARKX_SVG_TSV", unset = ""),
    list.files(file.path(svg_root, "sparkx"), pattern = "_svg\\.tsv$", full.names = TRUE)
  )
  paths <- unique(paths[nzchar(paths) & file.exists(paths)])
  rows <- lapply(paths, function(path) st09_normalize_method_df(st09_read_tsv(path)))
  if (length(rows) == 0) return(st09_empty_method_df())
  do.call(rbind, rows)
}

st09_condition_key <- function(x) {
  value <- trimws(as.character(x))
  value[!nzchar(value)] <- "all"
  value
}

st09_build_consensus <- function(method_df) {
  if (nrow(method_df) == 0) {
    return(data.frame(
      section_id = character(), condition = character(), gene_id = character(), gene_symbol = character(),
      best_score = numeric(), best_rank = integer(), support_methods = character(), formal_method_n = integer(),
      proxy_method_n = integer(), condition_n = integer(), svg_tier = character(), evidence_tier = character(),
      status = character(), reason = character(), stringsAsFactors = FALSE
    ))
  }
  df <- st09_normalize_method_df(method_df)
  df <- df[nzchar(df$gene_symbol) | nzchar(df$gene_id), , drop = FALSE]
  if (nrow(df) == 0) return(st09_build_consensus(st09_empty_method_df()))
  df$gene_key <- st09_normalize_gene(ifelse(nzchar(df$gene_symbol), df$gene_symbol, df$gene_id))
  df$condition <- st09_condition_key(df$condition)
  df$status_class <- st09_method_status_class(df$status)
  keys <- unique(df[, c("section_id", "condition", "gene_key"), drop = FALSE])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    hit <- df[df$section_id == key$section_id & df$condition == key$condition & df$gene_key == key$gene_key, , drop = FALSE]
    formal <- hit[hit$status_class == "formal", , drop = FALSE]
    proxy <- hit[hit$status_class == "proxy", , drop = FALSE]
    formal_n <- length(unique(formal$method))
    proxy_n <- length(unique(proxy$method))
    tier <- if (formal_n >= 2) "high_confidence_svg" else if (formal_n == 1) "single_method_svg" else if (proxy_n > 0) "proxy_only_svg" else "unsupported"
    data.frame(
      section_id = key$section_id,
      condition = key$condition,
      gene_id = st09_scalar(hit$gene_id, ""),
      gene_symbol = st09_scalar(hit$gene_symbol, st09_scalar(hit$gene_id, "")),
      best_score = suppressWarnings(max(hit$score, na.rm = TRUE)),
      best_rank = suppressWarnings(min(hit$rank, na.rm = TRUE)),
      support_methods = paste(sort(unique(hit$method[nzchar(hit$method)])), collapse = ";"),
      formal_method_n = formal_n,
      proxy_method_n = proxy_n,
      condition_n = 1L,
      svg_tier = tier,
      evidence_tier = ifelse(tier == "proxy_only_svg", "exploratory_caution", ifelse(tier == "unsupported", "no", "exploratory_yes")),
      status = ifelse(tier == "unsupported", "unsupported", "ok"),
      reason = ifelse(tier == "proxy_only_svg", "proxy-only SVG evidence cannot satisfy I06/I07 PASS", ""),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$best_score[!is.finite(out$best_score)] <- NA_real_
  out$best_rank[!is.finite(out$best_rank)] <- NA_integer_

  by_gene <- split(out, st09_normalize_gene(out$gene_symbol))
  stable_genes <- names(by_gene)[vapply(by_gene, function(x) length(unique(x$condition[x$svg_tier %in% c("high_confidence_svg", "single_method_svg")])) >= 2, logical(1))]
  out$svg_tier[st09_normalize_gene(out$gene_symbol) %in% stable_genes & out$svg_tier %in% c("high_confidence_svg", "single_method_svg")] <- "condition_stable_svg"
  out
}

st09_support_level_from_tier <- function(tier) {
  tier <- as.character(tier)
  ifelse(tier %in% c("high_confidence_svg", "condition_specific_svg"), "strong",
    ifelse(tier %in% c("condition_stable_svg", "single_method_svg"), "moderate",
      ifelse(tier == "proxy_only_svg", "weak", "none")
    )
  )
}

st09_consensus_lookup <- function(consensus) {
  if (nrow(consensus) == 0) return(setNames(character(), character()))
  keys <- st09_normalize_gene(c(consensus$gene_symbol, consensus$gene_id))
  values <- c(consensus$svg_tier, consensus$svg_tier)
  values[!nzchar(keys)] <- ""
  stats::setNames(values[nzchar(keys)], keys[nzchar(keys)])
}

st09_rank_lookup <- function(consensus) {
  if (nrow(consensus) == 0) return(setNames(character(), character()))
  keys <- st09_normalize_gene(c(consensus$gene_symbol, consensus$gene_id))
  values <- as.character(c(consensus$best_rank, consensus$best_rank))
  stats::setNames(values[nzchar(keys)], keys[nzchar(keys)])
}

st09_tier_for_gene <- function(gene, lookup, default = "unsupported") {
  key <- st09_normalize_gene(gene)
  hit <- unname(lookup[key])
  ifelse(is.na(hit) | !nzchar(hit), default, hit)
}

st09_rank_for_gene <- function(gene, lookup) {
  key <- st09_normalize_gene(gene)
  hit <- unname(lookup[key])
  ifelse(is.na(hit), "", hit)
}

st09_svg_gene_sets <- function(consensus) {
  cols <- c("gene_set_id", "gene_id", "gene_symbol", "condition", "svg_tier", "support_methods", "source", "intended_use")
  if (nrow(consensus) == 0) return(as.data.frame(setNames(rep(list(character()), length(cols)), cols), stringsAsFactors = FALSE))
  make_set <- function(id, rows, use) {
    if (nrow(rows) == 0) return(NULL)
    data.frame(
      gene_set_id = id,
      gene_id = rows$gene_id,
      gene_symbol = rows$gene_symbol,
      condition = rows$condition,
      svg_tier = rows$svg_tier,
      support_methods = rows$support_methods,
      source = "ST09_SVG",
      intended_use = use,
      stringsAsFactors = FALSE
    )
  }
  rows <- list(
    make_set("SVG_condition_stable", consensus[consensus$svg_tier == "condition_stable_svg", , drop = FALSE], "condition-stable spatial gene evidence"),
    make_set("SVG_single_method_supplemental", consensus[consensus$svg_tier == "single_method_svg", , drop = FALSE], "supplemental single-method spatial gene evidence")
  )
  for (cond in unique(consensus$condition[nzchar(consensus$condition)])) {
    cond_df <- consensus[consensus$condition == cond, , drop = FALSE]
    rows[[length(rows) + 1L]] <- make_set(paste0("SVG_high_confidence_", cond), cond_df[cond_df$svg_tier == "high_confidence_svg", , drop = FALSE], "high-confidence condition-specific spatial genes")
    rows[[length(rows) + 1L]] <- make_set(paste0("SVG_", cond, "_specific"), cond_df[cond_df$svg_tier %in% c("high_confidence_svg", "single_method_svg"), , drop = FALSE], "condition-specific spatial genes")
  }
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(as.data.frame(setNames(rep(list(character()), length(cols)), cols), stringsAsFactors = FALSE))
  out <- do.call(rbind, rows)
  out[, cols, drop = FALSE]
}

st09_enrichment_handoff <- function(gene_sets) {
  cols <- c("gene_set_id", "condition", "n_genes", "recommended_enrichment_mode", "status", "reason")
  if (nrow(gene_sets) == 0) return(as.data.frame(setNames(rep(list(character()), length(cols)), cols), stringsAsFactors = FALSE))
  keys <- unique(gene_sets[, c("gene_set_id", "condition"), drop = FALSE])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    hit <- gene_sets[gene_sets$gene_set_id == keys$gene_set_id[[i]] & gene_sets$condition == keys$condition[[i]], , drop = FALSE]
    n <- length(unique(st09_normalize_gene(hit$gene_symbol)))
    data.frame(
      gene_set_id = keys$gene_set_id[[i]],
      condition = keys$condition[[i]],
      n_genes = n,
      recommended_enrichment_mode = ifelse(n >= 10, "over_representation", "descriptive_only"),
      status = ifelse(n > 0, "ok", "skipped_empty_gene_set"),
      reason = ifelse(n >= 10, "", "Too few SVG genes for formal enrichment; use as descriptive handoff only."),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)[, cols, drop = FALSE]
}

st09_lr_support <- function(consensus, candidates) {
  cols <- c(
    "comm_candidate_id", "lr_axis_id", "ligand", "receptor", "sender_cell_type", "receiver_cell_type", "condition",
    "ligand_svg_tier", "receptor_svg_tier", "ligand_svg_rank", "receptor_svg_rank",
    "receiver_target_svg_support_n", "receiver_target_svg_support_genes", "svg_spatial_support_level", "svg_support_reason"
  )
  if (nrow(candidates) == 0) return(as.data.frame(setNames(rep(list(character()), length(cols)), cols), stringsAsFactors = FALSE))
  lookup <- st09_consensus_lookup(consensus)
  ranks <- st09_rank_lookup(consensus)
  for (col in c("comm_candidate_id", "lr_axis_id", "ligand", "receptor", "sender_cell_type", "receiver_cell_type", "condition_value", "condition")) {
    if (!col %in% colnames(candidates)) candidates[[col]] <- ""
  }
  rows <- lapply(seq_len(nrow(candidates)), function(i) {
    row <- candidates[i, , drop = FALSE]
    ligand_tier <- st09_tier_for_gene(row$ligand, lookup)
    receptor_tier <- st09_tier_for_gene(row$receptor, lookup)
    level <- "none"
    if (ligand_tier %in% c("high_confidence_svg", "condition_specific_svg") || receptor_tier %in% c("high_confidence_svg", "condition_specific_svg")) {
      level <- "strong"
    } else if (ligand_tier %in% c("single_method_svg", "condition_stable_svg") || receptor_tier %in% c("single_method_svg", "condition_stable_svg")) {
      level <- "moderate"
    } else if (ligand_tier == "proxy_only_svg" || receptor_tier == "proxy_only_svg") {
      level <- "weak"
    }
    data.frame(
      comm_candidate_id = row$comm_candidate_id,
      lr_axis_id = row$lr_axis_id,
      ligand = row$ligand,
      receptor = row$receptor,
      sender_cell_type = row$sender_cell_type,
      receiver_cell_type = row$receiver_cell_type,
      condition = st09_scalar(row$condition_value, st09_scalar(row$condition, "")),
      ligand_svg_tier = ligand_tier,
      receptor_svg_tier = receptor_tier,
      ligand_svg_rank = st09_rank_for_gene(row$ligand, ranks),
      receptor_svg_rank = st09_rank_for_gene(row$receptor, ranks),
      receiver_target_svg_support_n = 0L,
      receiver_target_svg_support_genes = "",
      svg_spatial_support_level = level,
      svg_support_reason = ifelse(level == "none", "No SVG overlap for ligand/receptor or receiver targets.", paste("SVG gene-level support:", level)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)[, cols, drop = FALSE]
}

st09_regulation_support <- function(consensus, scenic_targets) {
  cols <- c("regulator", "target_gene", "condition", "target_svg_tier", "target_svg_rank", "target_svg_method_support", "regulation_svg_support_level")
  if (nrow(scenic_targets) == 0) return(as.data.frame(setNames(rep(list(character()), length(cols)), cols), stringsAsFactors = FALSE))
  lookup <- st09_consensus_lookup(consensus)
  ranks <- st09_rank_lookup(consensus)
  regulator_col <- st09_pick_col(scenic_targets, c("regulator", "tf", "target_id", "source_question_id"))
  target_col <- st09_pick_col(scenic_targets, c("target_gene", "gene_symbol", "gene", "receiver_target_gene"))
  condition_col <- st09_pick_col(scenic_targets, c("condition", "condition_value", "stage"))
  if (!nzchar(target_col)) return(as.data.frame(setNames(rep(list(character()), length(cols)), cols), stringsAsFactors = FALSE))
  out <- data.frame(
    regulator = if (nzchar(regulator_col)) scenic_targets[[regulator_col]] else "",
    target_gene = scenic_targets[[target_col]],
    condition = if (nzchar(condition_col)) scenic_targets[[condition_col]] else "all",
    stringsAsFactors = FALSE
  )
  out$target_svg_tier <- st09_tier_for_gene(out$target_gene, lookup)
  out$target_svg_rank <- st09_rank_for_gene(out$target_gene, ranks)
  out$target_svg_method_support <- out$target_svg_tier
  out$regulation_svg_support_level <- st09_support_level_from_tier(out$target_svg_tier)
  out[, cols, drop = FALSE]
}

st09_trajectory_support <- function(consensus, trajectory_pairs, gene_program_targets) {
  cols <- c("trajectory_id", "question_id", "gene_id", "gene_symbol", "trajectory_gene_role", "condition", "svg_tier", "svg_rank", "support_methods", "trajectory_svg_support_level")
  lookup <- st09_consensus_lookup(consensus)
  ranks <- st09_rank_lookup(consensus)
  rows <- list()
  if (nrow(gene_program_targets) > 0) {
    gene_col <- st09_pick_col(gene_program_targets, c("gene_symbol", "gene", "target_gene"))
    if (nzchar(gene_col)) {
      id_col <- st09_pick_col(gene_program_targets, c("trajectory_id", "comparison_id", "source_question_id"))
      q_col <- st09_pick_col(gene_program_targets, c("source_question_id", "question_id"))
      cond_col <- st09_pick_col(gene_program_targets, c("condition", "condition_value", "stage"))
      rows[[length(rows) + 1L]] <- data.frame(
        trajectory_id = if (nzchar(id_col)) gene_program_targets[[id_col]] else "",
        question_id = if (nzchar(q_col)) gene_program_targets[[q_col]] else "",
        gene_id = "",
        gene_symbol = gene_program_targets[[gene_col]],
        trajectory_gene_role = "gene_program_target",
        condition = if (nzchar(cond_col)) gene_program_targets[[cond_col]] else "all",
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(as.data.frame(setNames(rep(list(character()), length(cols)), cols), stringsAsFactors = FALSE))
  out <- do.call(rbind, rows)
  out$svg_tier <- st09_tier_for_gene(out$gene_symbol, lookup)
  out$svg_rank <- st09_rank_for_gene(out$gene_symbol, ranks)
  out$support_methods <- out$svg_tier
  out$trajectory_svg_support_level <- st09_support_level_from_tier(out$svg_tier)
  out[, cols, drop = FALSE]
}

st09_question_gate_status <- function(consensus) {
  formal_high <- any(consensus$svg_tier %in% c("high_confidence_svg", "condition_stable_svg"))
  any_single <- any(consensus$svg_tier == "single_method_svg")
  any_proxy <- any(consensus$svg_tier == "proxy_only_svg")
  conditions <- unique(consensus$condition[consensus$svg_tier %in% c("high_confidence_svg", "condition_stable_svg", "single_method_svg")])
  i06_status <- if (formal_high) "PASS" else if (any_single || any_proxy) "WARN" else "FAIL"
  i07_status <- if (length(conditions) >= 2 && formal_high) "PASS" else if (length(conditions) >= 1 || any_proxy) "WARN" else "FAIL"
  data.frame(
    question_id = c("I06_SVG", "I07_SVG_by_stage"),
    gate_status = c(i06_status, i07_status),
    gate_level = "exploratory",
    interpretation_allowed = c(
      ifelse(i06_status == "PASS", "exploratory_yes", ifelse(i06_status == "WARN", "exploratory_caution", "no")),
      ifelse(i07_status == "PASS", "exploratory_yes", ifelse(i07_status == "WARN", "exploratory_caution", "no"))
    ),
    evidence_tiers = c(paste(sort(unique(consensus$svg_tier)), collapse = ";"), paste(sort(unique(consensus$svg_tier)), collapse = ";")),
    source_module = "spatial_09c_svg_consensus",
    reason = c(
      ifelse(i06_status == "PASS", "Formal SVG evidence is available.", ifelse(i06_status == "WARN", "Only single-method or proxy SVG evidence is available.", "No usable SVG evidence.")),
      ifelse(i07_status == "PASS", "Both stage/condition comparison and formal SVG evidence are available.", ifelse(i07_status == "WARN", "Only one condition/method or proxy SVG evidence is available.", "No stage-level SVG comparison evidence."))
    ),
    stringsAsFactors = FALSE
  )
}
