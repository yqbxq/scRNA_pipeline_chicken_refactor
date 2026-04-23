#!/usr/bin/env Rscript

load_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      sprintf("缺少 R 包: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

load_required_packages(c("dplyr", "ggplot2", "jsonlite"))

reference_dir <- "/home/user_test/syf_f5/01shared_resources/reference/chicken"
clean_gtf <- file.path(reference_dir, "GRCg7b_genomic_clean.gtf")
reference_gtf <- file.path(reference_dir, "GRCg7b_genomic.gtf")
output_dir <- "/home/user_test/syf_f5/05projects/scrna_improve/ortholog_cache"
figure_dir <- file.path(output_dir, "figures")
manifest_path <- file.path(output_dir, "_manifest.json")
target_species <- c("human", "mouse")
module_name <- "00b_ortholog_report"
module_contract <- "00_ortholog_module"
module_version <- "1.0"
key_genes <- c("DRGX", "CREB3L2", "EMX2", "CEBPB", "FOSL2", "JUN", "PBX3", "HMGA1", "SREBF2")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

strip_ensembl_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

extract_gtf_attr <- function(attr_vec, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  matches <- regexec(pattern, attr_vec, perl = TRUE)
  values <- regmatches(attr_vec, matches)
  vapply(values, function(hit) {
    if (length(hit) >= 2) {
      hit[2]
    } else {
      ""
    }
  }, character(1))
}

select_gtf_path <- function(clean_gtf, reference_gtf) {
  for (path in c(clean_gtf, reference_gtf)) {
    if (nzchar(path) && file.exists(path)) {
      return(path)
    }
  }
  stop(
    sprintf("GTF 不存在: %s 或 %s", clean_gtf, reference_gtf),
    call. = FALSE
  )
}

read_reference_annotation_local <- function(clean_gtf, reference_gtf) {
  gtf_path <- select_gtf_path(clean_gtf, reference_gtf)
  con <- if (grepl("\\.gz$", gtf_path)) gzfile(gtf_path, open = "rt") else file(gtf_path, open = "rt")
  on.exit(close(con), add = TRUE)
  gtf_df <- tryCatch(
    read.delim(
      con,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE,
      comment.char = "#",
      quote = ""
    ),
    error = function(e) NULL
  )
  if (is.null(gtf_df) || nrow(gtf_df) == 0 || ncol(gtf_df) < 9) {
    stop(sprintf("无法解析 GTF: %s", gtf_path), call. = FALSE)
  }

  gene_rows <- gtf_df[gtf_df[[3]] == "gene", c(1, 9), drop = FALSE]
  if (nrow(gene_rows) == 0) {
    stop(sprintf("GTF 中没有 gene 记录: %s", gtf_path), call. = FALSE)
  }

  attr_field <- gene_rows[[2]]
  gene_id <- extract_gtf_attr(attr_field, "gene_id")
  gene_name <- extract_gtf_attr(attr_field, "gene_name")
  gene_name[!nzchar(gene_name)] <- extract_gtf_attr(attr_field, "Name")[!nzchar(gene_name)]
  gene_biotype <- extract_gtf_attr(attr_field, "gene_biotype")
  gene_biotype[!nzchar(gene_biotype)] <- extract_gtf_attr(attr_field, "gene_type")[!nzchar(gene_biotype)]
  gene_biotype[!nzchar(gene_biotype)] <- extract_gtf_attr(attr_field, "biotype")[!nzchar(gene_biotype)]

  annotation_df <- unique(data.frame(
    seqname = as.character(gene_rows[[1]]),
    gene_id = gene_id,
    gene_id_stripped = strip_ensembl_version(gene_id),
    gene_name = gene_name,
    gene_name_upper = toupper(gene_name),
    gene_biotype = gene_biotype,
    stringsAsFactors = FALSE
  ))
  annotation_df <- annotation_df[nzchar(annotation_df$gene_id) | nzchar(annotation_df$gene_name), , drop = FALSE]
  list(gtf_path = gtf_path, annotation_df = annotation_df)
}

normalize_key <- function(x) {
  toupper(trimws(as.character(x)))
}

orthology_bucket <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    grepl("one2one$", x) ~ "one2one",
    grepl("one2many$", x) ~ "one2many",
    grepl("many2many$", x) ~ "many2many",
    TRUE ~ "other"
  )
}

read_manifest_local <- function(manifest_path) {
  if (!file.exists(manifest_path)) {
    stop(sprintf("缺少 manifest: %s", manifest_path), call. = FALSE)
  }
  jsonlite::read_json(manifest_path, simplifyVector = FALSE)
}

resolve_output_local <- function(manifest, key) {
  entry <- manifest$outputs[[key]]
  if (is.null(entry) || is.null(entry$path)) {
    stop(sprintf("manifest 缺少输出键: %s", key), call. = FALSE)
  }
  raw_path <- as.character(entry$path)
  if (grepl("^/", raw_path)) {
    normalizePath(raw_path, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(as.character(manifest$base_dir), raw_path), winslash = "/", mustWork = FALSE)
  }
}

write_manifest_local <- function(manifest_path, new_outputs, module_name = NULL, base_dir = NULL, inputs = NULL) {
  if (file.exists(manifest_path)) {
    existing <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
    existing$outputs <- modifyList(existing$outputs %||% list(), new_outputs)
    existing$timestamp <- timestamp_now()
    manifest <- existing
  } else {
    if (is.null(module_name) || is.null(base_dir) || is.null(inputs)) {
      stop("manifest 不存在且缺少初始化参数", call. = FALSE)
    }
    manifest <- list(
      module = module_name,
      version = module_version,
      timestamp = timestamp_now(),
      base_dir = base_dir,
      inputs = inputs,
      outputs = new_outputs
    )
  }
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

relative_path_local <- function(path, base_dir) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  normalized_base <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
  prefix <- paste0(normalized_base, "/")
  if (startsWith(normalized_path, prefix)) {
    substring(normalized_path, nchar(prefix) + 1L)
  } else {
    normalized_path
  }
}

build_output_entry <- function(path, type, produced_by, row_semantics) {
  list(
    path = relative_path_local(path, output_dir),
    type = type,
    produced_by = produced_by,
    row_semantics = row_semantics
  )
}

write_markdown_local <- function(lines, path) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

escape_markdown_cell <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("\\|", "\\\\|", x)
  gsub("\n", "<br>", x, fixed = TRUE)
}

render_markdown_table_local <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return("No data available.")
  }
  headers <- paste(sprintf(" %s ", vapply(colnames(df), escape_markdown_cell, character(1))), collapse = "|")
  divider <- paste(rep(" --- ", ncol(df)), collapse = "|")
  rows <- apply(df, 1, function(row) {
    paste(sprintf(" %s ", vapply(row, escape_markdown_cell, character(1))), collapse = "|")
  })
  c(paste0("|", headers, "|"), paste0("|", divider, "|"), paste0("|", rows, "|"))
}

fmt_int <- function(x) {
  if (is.na(x)) {
    return("NA")
  }
  prettyNum(round(x), big.mark = ",", scientific = FALSE)
}

fmt_pct <- function(x) {
  if (is.na(x)) {
    return("NA")
  }
  sprintf("%.1f%%", 100 * x)
}

fmt_n_pct <- function(n, total) {
  if (is.na(total) || total <= 0) {
    return(fmt_int(n))
  }
  sprintf("%s (%s)", fmt_int(n), fmt_pct(n / total))
}

fmt_num <- function(x, digits = 1) {
  if (is.na(x)) {
    return("NA")
  }
  sprintf(paste0("%.", digits, "f"), x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE))
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  stats::median(x, na.rm = TRUE)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  min(x, na.rm = TRUE)
}

safe_rate <- function(numerator, denominator) {
  if (is.na(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  numerator / denominator
}

save_plot_local <- function(plot_obj, path, width, height) {
  ggplot2::ggsave(path, plot_obj, width = width, height = height, dpi = 300, bg = "white")
}

read_csv_required <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("缺少输入文件: %s", path), call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

prepare_best_table <- function(df, species_key) {
  df$target_species <- species_key
  df$type_bucket <- orthology_bucket(df$orthology_type)
  df
}

build_species_summary <- function(best_df, annotation_df, input_genes) {
  mapped_genes <- unique(best_df$external_gene_name)
  unmapped_genes <- setdiff(input_genes, mapped_genes)

  type_counts <- table(best_df$type_bucket)
  one2one_count <- if ("one2one" %in% names(type_counts)) unname(type_counts["one2one"]) else 0
  one2many_count <- if ("one2many" %in% names(type_counts)) unname(type_counts["one2many"]) else 0
  many2many_count <- if ("many2many" %in% names(type_counts)) unname(type_counts["many2many"]) else 0

  gold_count <- sum(best_df$pair_quality == "gold", na.rm = TRUE)
  silver_count <- sum(best_df$pair_quality == "silver", na.rm = TRUE)
  ambiguous_count <- sum(best_df$pair_quality == "ambiguous", na.rm = TRUE)

  one2one_df <- best_df[best_df$type_bucket == "one2one", , drop = FALSE]
  dnds <- one2one_df$dn / one2one_df$ds
  dnds[!is.finite(dnds) | one2one_df$ds <= 0] <- NA_real_

  goc_available_n <- sum(!is.na(best_df$goc_score))
  wga_available_n <- sum(!is.na(best_df$wga_coverage))

  unmapped_biotype <- annotation_df[annotation_df$gene_name %in% unmapped_genes, , drop = FALSE]
  if (nrow(unmapped_biotype) > 0) {
    unmapped_biotype <- unmapped_biotype %>%
      dplyr::mutate(gene_biotype = ifelse(is.na(gene_biotype) | gene_biotype == "", "unknown", gene_biotype)) %>%
      dplyr::count(gene_biotype, sort = TRUE, name = "count")
  } else {
    unmapped_biotype <- data.frame(gene_biotype = character(0), count = integer(0), stringsAsFactors = FALSE)
  }

  list(
    total_input = length(input_genes),
    mapped_count = length(mapped_genes),
    unmapped_count = length(unmapped_genes),
    coverage_fraction = safe_rate(length(mapped_genes), length(input_genes)),
    one2one_count = one2one_count,
    one2many_count = one2many_count,
    many2many_count = many2many_count,
    gold_count = gold_count,
    silver_count = silver_count,
    ambiguous_count = ambiguous_count,
    query_id_median = safe_median(one2one_df$perc_id),
    query_id_q25 = safe_quantile(one2one_df$perc_id, 0.25),
    query_id_q75 = safe_quantile(one2one_df$perc_id, 0.75),
    query_id_min = safe_min(one2one_df$perc_id),
    target_id_median = safe_median(one2one_df$perc_id_r1),
    target_id_q25 = safe_quantile(one2one_df$perc_id_r1, 0.25),
    target_id_q75 = safe_quantile(one2one_df$perc_id_r1, 0.75),
    target_id_min = safe_min(one2one_df$perc_id_r1),
    goc_available_rate = safe_rate(goc_available_n, nrow(best_df)),
    goc_support_rate = safe_rate(sum(!is.na(best_df$goc_score) & best_df$goc_score >= 75), goc_available_n),
    wga_available_rate = safe_rate(wga_available_n, nrow(best_df)),
    wga_support_rate = safe_rate(sum(!is.na(best_df$wga_coverage) & best_df$wga_coverage >= 50), wga_available_n),
    dnds_median = safe_median(dnds),
    positive_selection_count = sum(dnds > 1, na.rm = TRUE),
    unmapped_biotype = unmapped_biotype,
    unmapped_genes = unmapped_genes,
    best_df = best_df
  )
}

build_key_hits <- function(best_tables, key_genes) {
  combined <- dplyr::bind_rows(best_tables)
  if (nrow(combined) == 0) {
    return(combined)
  }
  combined %>%
    dplyr::filter(
      normalize_key(external_gene_name) %in% key_genes |
        normalize_key(target_gene_name) %in% key_genes
    ) %>%
    dplyr::select(
      target_species,
      external_gene_name,
      target_gene_name,
      ensembl_gene_id,
      target_ensembl_gene,
      orthology_type,
      orthology_confidence,
      perc_id,
      perc_id_r1,
      goc_score,
      wga_coverage,
      dn,
      ds,
      pair_quality
    )
}

build_species_comparison_long <- function(summaries) {
  dplyr::bind_rows(lapply(names(summaries), function(species_key) {
    summary_item <- summaries[[species_key]]
    data.frame(
      target_species = species_key,
      metric = c("one2one_fraction", "gold_fraction", "median_query_identity", "goc_support_rate"),
      value = c(
        100 * safe_rate(summary_item$one2one_count, summary_item$total_input),
        100 * safe_rate(summary_item$gold_count, summary_item$total_input),
        summary_item$query_id_median,
        100 * summary_item$goc_support_rate
      ),
      stringsAsFactors = FALSE
    )
  }))
}

recommend_species <- function(summaries) {
  species_keys <- names(summaries)
  if (length(species_keys) < 2) {
    return(list(
      recommended = species_keys[1],
      reasons = c("仅查询单一目标物种，未做跨物种推荐比较。")
    ))
  }

  score_matrix <- data.frame(
    target_species = species_keys,
    gold_count = vapply(summaries, function(x) x$gold_count, numeric(1)),
    one2one_count = vapply(summaries, function(x) x$one2one_count, numeric(1)),
    median_query_identity = vapply(summaries, function(x) ifelse(is.na(x$query_id_median), -Inf, x$query_id_median), numeric(1)),
    goc_support_rate = vapply(summaries, function(x) ifelse(is.na(x$goc_support_rate), -Inf, x$goc_support_rate), numeric(1)),
    stringsAsFactors = FALSE
  )

  ordered_metrics <- c("gold_count", "one2one_count", "median_query_identity", "goc_support_rate")
  recommended <- species_keys[1]
  for (metric_name in ordered_metrics) {
    metric_values <- score_matrix[[metric_name]]
    if (length(unique(metric_values)) > 1) {
      recommended <- score_matrix$target_species[which.max(metric_values)]
      break
    }
  }

  human <- summaries[["human"]]
  mouse <- summaries[["mouse"]]
  reasons <- c(
    sprintf(
      "Gold 映射数量: human=%s, mouse=%s",
      fmt_int(human$gold_count),
      fmt_int(mouse$gold_count)
    ),
    sprintf(
      "one2one 数量: human=%s, mouse=%s",
      fmt_int(human$one2one_count),
      fmt_int(mouse$one2one_count)
    ),
    sprintf(
      "中位 Query %%ID: human=%s, mouse=%s",
      fmt_num(human$query_id_median),
      fmt_num(mouse$query_id_median)
    ),
    sprintf(
      "GOC>=75 支持率: human=%s, mouse=%s",
      fmt_pct(human$goc_support_rate),
      fmt_pct(mouse$goc_support_rate)
    )
  )

  list(recommended = recommended, reasons = reasons)
}

ensure_dir(output_dir)
ensure_dir(figure_dir)

annotation_payload <- read_reference_annotation_local(clean_gtf, reference_gtf)
gtf_path <- annotation_payload$gtf_path
annotation_df <- annotation_payload$annotation_df
input_genes <- sort(unique(annotation_df$gene_name[nzchar(annotation_df$gene_name)]))

manifest <- read_manifest_local(manifest_path)
best_tables <- list()
all_tables <- list()
summaries <- list()

for (species_key in target_species) {
  best_path <- resolve_output_local(manifest, paste0(species_key, "_best"))
  all_path <- resolve_output_local(manifest, paste0(species_key, "_all"))
  best_df <- prepare_best_table(read_csv_required(best_path), species_key)
  all_df <- read_csv_required(all_path)
  best_tables[[species_key]] <- best_df
  all_tables[[species_key]] <- all_df
  summaries[[species_key]] <- build_species_summary(best_df, annotation_df, input_genes)
}

key_hits <- build_key_hits(all_tables, key_genes)
key_hits_path <- file.path(output_dir, "ortholog_key_gene_hits.csv")
write.csv(key_hits, key_hits_path, row.names = FALSE)

coverage_df <- data.frame(Metric = c(
  "Total input genes",
  "one2one",
  "one2many",
  "many2many",
  "Unmapped",
  "Coverage"
), stringsAsFactors = FALSE)
for (species_key in target_species) {
  summary_item <- summaries[[species_key]]
  coverage_df[[tools::toTitleCase(species_key)]] <- c(
    fmt_int(summary_item$total_input),
    fmt_n_pct(summary_item$one2one_count, summary_item$total_input),
    fmt_n_pct(summary_item$one2many_count, summary_item$total_input),
    fmt_n_pct(summary_item$many2many_count, summary_item$total_input),
    fmt_n_pct(summary_item$unmapped_count, summary_item$total_input),
    fmt_pct(summary_item$coverage_fraction)
  )
}

quality_df <- data.frame(
  Quality = c("Gold", "Silver", "Ambiguous"),
  Criteria = c(
    "one2one + orthology_confidence == 1",
    "one2one + orthology_confidence in {0, NA}",
    "one2many / many2many"
  ),
  stringsAsFactors = FALSE
)
for (species_key in target_species) {
  summary_item <- summaries[[species_key]]
  quality_df[[tools::toTitleCase(species_key)]] <- c(
    fmt_n_pct(summary_item$gold_count, summary_item$total_input),
    fmt_n_pct(summary_item$silver_count, summary_item$total_input),
    fmt_n_pct(summary_item$ambiguous_count, summary_item$total_input)
  )
}

identity_df <- dplyr::bind_rows(lapply(target_species, function(species_key) {
  summary_item <- summaries[[species_key]]
  data.frame(
    target_species = species_key,
    Direction = c("Chicken->Target", "Target->Chicken"),
    Median = c(summary_item$query_id_median, summary_item$target_id_median),
    Q25 = c(summary_item$query_id_q25, summary_item$target_id_q25),
    Q75 = c(summary_item$query_id_q75, summary_item$target_id_q75),
    Min = c(summary_item$query_id_min, summary_item$target_id_min),
    stringsAsFactors = FALSE
  )
}))
identity_df$Median <- vapply(identity_df$Median, fmt_num, character(1))
identity_df$Q25 <- vapply(identity_df$Q25, fmt_num, character(1))
identity_df$Q75 <- vapply(identity_df$Q75, fmt_num, character(1))
identity_df$Min <- vapply(identity_df$Min, fmt_num, character(1))

genomic_df <- data.frame(
  Metric = c("GOC score available", "GOC >= 75", "WGA coverage available", "WGA >= 50"),
  stringsAsFactors = FALSE
)
for (species_key in target_species) {
  summary_item <- summaries[[species_key]]
  genomic_df[[tools::toTitleCase(species_key)]] <- c(
    fmt_pct(summary_item$goc_available_rate),
    fmt_pct(summary_item$goc_support_rate),
    fmt_pct(summary_item$wga_available_rate),
    fmt_pct(summary_item$wga_support_rate)
  )
}

selection_df <- data.frame(
  Statistic = c("Median dN/dS", "dN/dS > 1"),
  stringsAsFactors = FALSE
)
for (species_key in target_species) {
  summary_item <- summaries[[species_key]]
  selection_df[[tools::toTitleCase(species_key)]] <- c(
    fmt_num(summary_item$dnds_median, digits = 2),
    fmt_int(summary_item$positive_selection_count)
  )
}

coverage_path <- file.path(output_dir, "ortholog_coverage_summary.csv")
write.csv(coverage_df, coverage_path, row.names = FALSE)

identity_plot_df <- dplyr::bind_rows(lapply(best_tables, function(best_df) {
  best_df[best_df$type_bucket == "one2one" & is.finite(best_df$perc_id) & is.finite(best_df$perc_id_r1), , drop = FALSE]
}))
identity_plot <- ggplot2::ggplot(identity_plot_df, ggplot2::aes(x = perc_id, y = perc_id_r1, color = pair_quality)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +
  ggplot2::geom_point(alpha = 0.7, size = 1.5) +
  ggplot2::facet_wrap(~target_species) +
  ggplot2::scale_color_manual(values = c(gold = "#2b8cbe", silver = "#969696", ambiguous = "#d95f0e")) +
  ggplot2::labs(
    title = "Identity Scatter (one2one)",
    x = "Query %ID",
    y = "Target %ID",
    color = "Pair quality"
  ) +
  ggplot2::theme_bw(base_size = 11)
identity_plot_path <- file.path(figure_dir, "ortholog_identity_scatter.png")
save_plot_local(identity_plot, identity_plot_path, width = 10, height = 5)

goc_plot_df <- dplyr::bind_rows(lapply(best_tables, function(best_df) {
  best_df[!is.na(best_df$goc_score), , drop = FALSE]
}))
goc_plot <- ggplot2::ggplot(goc_plot_df, ggplot2::aes(x = goc_score, fill = target_species)) +
  ggplot2::geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  ggplot2::geom_vline(xintercept = 75, linetype = "dashed", color = "firebrick") +
  ggplot2::facet_wrap(~target_species) +
  ggplot2::labs(
    title = "GOC Score Distribution",
    x = "GOC score",
    y = "Count",
    fill = "Target species"
  ) +
  ggplot2::theme_bw(base_size = 11)
goc_plot_path <- file.path(figure_dir, "ortholog_goc_distribution.png")
save_plot_local(goc_plot, goc_plot_path, width = 10, height = 5)

unmapped_plot_df <- dplyr::bind_rows(lapply(names(summaries), function(species_key) {
  df <- summaries[[species_key]]$unmapped_biotype
  if (nrow(df) == 0) {
    return(data.frame(target_species = species_key, gene_biotype = "none", count = 0, stringsAsFactors = FALSE))
  }
  df$target_species <- species_key
  df
}))
unmapped_plot <- ggplot2::ggplot(unmapped_plot_df, ggplot2::aes(
  x = stats::reorder(gene_biotype, count),
  y = count,
  fill = target_species
)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Unmapped Gene Biotype Distribution",
    x = "Gene biotype",
    y = "Count",
    fill = "Target species"
  ) +
  ggplot2::theme_bw(base_size = 11)
unmapped_plot_path <- file.path(figure_dir, "ortholog_unmapped_biotype.png")
save_plot_local(unmapped_plot, unmapped_plot_path, width = 10, height = 6)

comparison_plot_df <- build_species_comparison_long(summaries)
comparison_plot <- ggplot2::ggplot(comparison_plot_df, ggplot2::aes(x = target_species, y = value, fill = target_species)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::facet_wrap(~metric, scales = "free_y") +
  ggplot2::labs(
    title = "Species Comparison",
    x = "Target species",
    y = "Value",
    fill = "Target species"
  ) +
  ggplot2::theme_bw(base_size = 11)
comparison_plot_path <- file.path(figure_dir, "ortholog_species_comparison.png")
save_plot_local(comparison_plot, comparison_plot_path, width = 10, height = 6)

recommendation <- recommend_species(summaries)

report_lines <- c(
  "# Ortholog Mapping Quality Report",
  "",
  sprintf("Generated: %s | Source: `%s`", timestamp_now(), gtf_path),
  "",
  "## Coverage Summary",
  "",
  render_markdown_table_local(coverage_df),
  "",
  "## Per-pair Quality Grading",
  "",
  render_markdown_table_local(quality_df),
  "",
  "## Sequence Identity Distribution (one2one only)",
  "",
  render_markdown_table_local(identity_df),
  "",
  sprintf("Identity scatter plot: `%s`", identity_plot_path),
  "",
  "## Genomic Context Support",
  "",
  render_markdown_table_local(genomic_df),
  "",
  sprintf("GOC distribution plot: `%s`", goc_plot_path),
  "",
  "## Selection Pressure (dN/dS)",
  "",
  render_markdown_table_local(selection_df),
  "",
  "## Unmapped Gene Biotype Distribution",
  ""
)

for (species_key in target_species) {
  unmapped_df <- summaries[[species_key]]$unmapped_biotype
  if (nrow(unmapped_df) > 0) {
    unmapped_df$Fraction <- vapply(safe_rate(unmapped_df$count, sum(unmapped_df$count)), fmt_pct, character(1))
    names(unmapped_df) <- c("Biotype", "Count", "Fraction")
    report_lines <- c(
      report_lines,
      sprintf("### %s", tools::toTitleCase(species_key)),
      "",
      render_markdown_table_local(unmapped_df),
      ""
    )
  }
}

report_lines <- c(
  report_lines,
  sprintf("Unmapped biotype plot: `%s`", unmapped_plot_path),
  "",
  "## Key Gene Audit",
  ""
)

if (nrow(key_hits) > 0) {
  key_hits_md <- key_hits
  key_hits_md$orthology_confidence <- vapply(key_hits_md$orthology_confidence, fmt_num, character(1))
  key_hits_md$perc_id <- vapply(key_hits_md$perc_id, fmt_num, character(1))
  key_hits_md$perc_id_r1 <- vapply(key_hits_md$perc_id_r1, fmt_num, character(1))
  key_hits_md$goc_score <- vapply(key_hits_md$goc_score, fmt_num, character(1))
  key_hits_md$wga_coverage <- vapply(key_hits_md$wga_coverage, fmt_num, character(1))
  key_hits_md$dn <- vapply(key_hits_md$dn, fmt_num, character(1))
  key_hits_md$ds <- vapply(key_hits_md$ds, fmt_num, character(1))
  report_lines <- c(report_lines, render_markdown_table_local(key_hits_md), "")
} else {
  report_lines <- c(report_lines, "No key genes matched.", "")
}

report_lines <- c(
  report_lines,
  "## Recommendation",
  "",
  sprintf("**Recommended target species: %s**", recommendation$recommended),
  ""
)
report_lines <- c(report_lines, paste0("- ", recommendation$reasons), "")
report_lines <- c(
  report_lines,
  sprintf("Species comparison plot: `%s`", comparison_plot_path)
)

report_path <- file.path(output_dir, "ortholog_quality_report.md")
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = manifest_path,
  new_outputs = list(
    quality_report = build_output_entry(report_path, "md", module_name, "ortholog quality markdown report"),
    coverage_summary = build_output_entry(coverage_path, "csv", module_name, "cross-species coverage summary"),
    key_gene_hits = build_output_entry(key_hits_path, "csv", module_name, "key gene audit rows across target species"),
    identity_scatter = build_output_entry(identity_plot_path, "png", module_name, "identity scatter plot for one2one pairs"),
    goc_distribution = build_output_entry(goc_plot_path, "png", module_name, "GOC score distribution plot"),
    unmapped_biotype_plot = build_output_entry(unmapped_plot_path, "png", module_name, "unmapped biotype distribution plot"),
    species_comparison_plot = build_output_entry(comparison_plot_path, "png", module_name, "cross-species comparison bar chart")
  ),
  module_name = module_contract,
  base_dir = output_dir,
  inputs = list(gtf_path = gtf_path, target_species = target_species)
)

message("00b 完成。报告目录: ", output_dir)
