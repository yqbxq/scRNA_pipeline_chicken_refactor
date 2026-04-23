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

load_required_packages(c("dplyr", "biomaRt", "jsonlite"))

reference_dir <- "/home/user_test/syf_f5/01shared_resources/reference/chicken"
clean_gtf <- file.path(reference_dir, "GRCg7b_genomic_clean.gtf")
reference_gtf <- file.path(reference_dir, "GRCg7b_genomic.gtf")
output_dir <- "/home/user_test/syf_f5/05projects/scrna_improve/ortholog_cache"
figure_dir <- file.path(output_dir, "figures")
manifest_path <- file.path(output_dir, "_manifest.json")
target_species <- c("human", "mouse")
mirrors <- c("asia", "www", "useast")
chunk_size <- 500L
module_name <- "00a_build_ortholog_cache"
module_version <- "1.0"
module_contract <- "00_ortholog_module"

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

orthology_rank <- function(x) {
  dplyr::case_when(
    x == "ortholog_one2one" ~ 1L,
    x == "apparent_ortholog_one2one" ~ 2L,
    x == "ortholog_one2many" ~ 3L,
    x == "apparent_ortholog_one2many" ~ 4L,
    x == "ortholog_many2many" ~ 5L,
    TRUE ~ 99L
  )
}

species_prefix <- function(species_key) {
  switch(
    species_key,
    human = "hsapiens",
    mouse = "mmusculus",
    stop(sprintf("不支持的目标物种: %s", species_key), call. = FALSE)
  )
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

connect_chicken_mart <- function(mirror_name) {
  message("连接 Ensembl 鸡数据库镜像: ", mirror_name)
  biomaRt::useEnsembl(
    biomart = "genes",
    dataset = "ggallus_gene_ensembl",
    mirror = mirror_name
  )
}

query_one_chunk <- function(mart, gene_chunk, species_key) {
  prefix <- species_prefix(species_key)
  biomaRt::getBM(
    attributes = c(
      "external_gene_name",
      "ensembl_gene_id",
      paste0(prefix, "_homolog_associated_gene_name"),
      paste0(prefix, "_homolog_ensembl_gene"),
      paste0(prefix, "_homolog_orthology_type"),
      paste0(prefix, "_homolog_orthology_confidence"),
      paste0(prefix, "_homolog_perc_id"),
      paste0(prefix, "_homolog_perc_id_r1"),
      paste0(prefix, "_homolog_goc_score"),
      paste0(prefix, "_homolog_wga_coverage"),
      paste0(prefix, "_homolog_dn"),
      paste0(prefix, "_homolog_ds")
    ),
    filters = "external_gene_name",
    values = gene_chunk,
    mart = mart
  )
}

query_one_chunk_with_retry <- function(gene_chunk, chunk_idx, total_chunks, species_key, mirrors, max_attempts = 6L) {
  last_err <- NULL
  for (attempt in seq_len(max_attempts)) {
    mirror_name <- mirrors[((attempt - 1L) %% length(mirrors)) + 1L]
    message(sprintf(
      "  [%s] chunk %d/%d attempt %d/%d via %s",
      species_key,
      chunk_idx,
      total_chunks,
      attempt,
      max_attempts,
      mirror_name
    ))
    res <- tryCatch({
      mart <- connect_chicken_mart(mirror_name)
      query_one_chunk(mart, gene_chunk, species_key)
    }, error = function(e) {
      last_err <<- e
      NULL
    })
    if (!is.null(res)) {
      return(res)
    }
    if (!is.null(last_err)) {
      message("    failed: ", conditionMessage(last_err))
    }
    Sys.sleep(min(5L * attempt, 20L))
  }
  stop(
    sprintf(
      "[%s] chunk %d 在 %d 次尝试后仍失败。最后错误: %s",
      species_key,
      chunk_idx,
      max_attempts,
      if (is.null(last_err)) "unknown" else conditionMessage(last_err)
    ),
    call. = FALSE
  )
}

classify_pair_quality <- function(df) {
  bucket <- orthology_bucket(df$orthology_type)
  df$pair_quality <- dplyr::case_when(
    bucket == "one2one" & !is.na(df$orthology_confidence) & df$orthology_confidence == 1 ~ "gold",
    bucket == "one2one" ~ "silver",
    bucket %in% c("one2many", "many2many") ~ "ambiguous",
    TRUE ~ "ambiguous"
  )
  df
}

empty_candidate_table <- function() {
  data.frame(
    external_gene_name = character(0),
    ensembl_gene_id = character(0),
    target_species = character(0),
    target_gene_name = character(0),
    target_ensembl_gene = character(0),
    orthology_type = character(0),
    orthology_confidence = numeric(0),
    perc_id = numeric(0),
    perc_id_r1 = numeric(0),
    goc_score = numeric(0),
    wga_coverage = numeric(0),
    dn = numeric(0),
    ds = numeric(0),
    same_symbol = logical(0),
    orthology_rank = integer(0),
    pair_quality = character(0),
    stringsAsFactors = FALSE
  )
}

standardize_candidates <- function(raw_df, species_key) {
  if (is.null(raw_df) || nrow(raw_df) == 0) {
    return(empty_candidate_table())
  }

  prefix <- species_prefix(species_key)
  name_col <- paste0(prefix, "_homolog_associated_gene_name")
  ensembl_col <- paste0(prefix, "_homolog_ensembl_gene")
  type_col <- paste0(prefix, "_homolog_orthology_type")
  conf_col <- paste0(prefix, "_homolog_orthology_confidence")
  perc_id_col <- paste0(prefix, "_homolog_perc_id")
  perc_id_r1_col <- paste0(prefix, "_homolog_perc_id_r1")
  goc_col <- paste0(prefix, "_homolog_goc_score")
  wga_col <- paste0(prefix, "_homolog_wga_coverage")
  dn_col <- paste0(prefix, "_homolog_dn")
  ds_col <- paste0(prefix, "_homolog_ds")

  standardized <- raw_df %>%
    dplyr::transmute(
      external_gene_name = trimws(as.character(.data$external_gene_name)),
      ensembl_gene_id = trimws(as.character(.data$ensembl_gene_id)),
      target_species = species_key,
      target_gene_name = trimws(as.character(.data[[name_col]])),
      target_ensembl_gene = trimws(as.character(.data[[ensembl_col]])),
      orthology_type = trimws(as.character(.data[[type_col]])),
      orthology_confidence = safe_numeric(.data[[conf_col]]),
      perc_id = safe_numeric(.data[[perc_id_col]]),
      perc_id_r1 = safe_numeric(.data[[perc_id_r1_col]]),
      goc_score = safe_numeric(.data[[goc_col]]),
      wga_coverage = safe_numeric(.data[[wga_col]]),
      dn = safe_numeric(.data[[dn_col]]),
      ds = safe_numeric(.data[[ds_col]])
    ) %>%
    dplyr::filter(
      !is.na(external_gene_name),
      external_gene_name != "",
      (!is.na(target_gene_name) & target_gene_name != "") |
        (!is.na(target_ensembl_gene) & target_ensembl_gene != "")
    ) %>%
    dplyr::mutate(
      same_symbol = normalize_key(external_gene_name) == normalize_key(target_gene_name),
      orthology_rank = orthology_rank(orthology_type)
    ) %>%
    dplyr::distinct(
      external_gene_name,
      ensembl_gene_id,
      target_species,
      target_gene_name,
      target_ensembl_gene,
      orthology_type,
      orthology_confidence,
      perc_id,
      perc_id_r1,
      goc_score,
      wga_coverage,
      dn,
      ds,
      .keep_all = TRUE
    )

  classify_pair_quality(standardized)
}

select_best_per_gene <- function(df) {
  if (nrow(df) == 0) {
    return(df)
  }

  df %>%
    dplyr::arrange(
      external_gene_name,
      dplyr::desc(same_symbol),
      orthology_rank,
      dplyr::desc(dplyr::coalesce(orthology_confidence, -1)),
      dplyr::desc(dplyr::coalesce(perc_id, -1)),
      dplyr::desc(dplyr::coalesce(perc_id_r1, -1)),
      target_gene_name
    ) %>%
    dplyr::group_by(external_gene_name) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
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

build_output_entry <- function(path, type, produced_by, row_semantics, schema = NULL) {
  entry <- list(
    path = relative_path_local(path, output_dir),
    type = type,
    produced_by = produced_by,
    row_semantics = row_semantics
  )
  if (!is.null(schema)) {
    entry$schema <- schema
  }
  entry
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

schema_ortholog <- list(
  external_gene_name = "character",
  ensembl_gene_id = "character",
  target_species = "character",
  target_gene_name = "character",
  target_ensembl_gene = "character",
  orthology_type = "character",
  orthology_confidence = "numeric",
  perc_id = "numeric",
  perc_id_r1 = "numeric",
  goc_score = "numeric",
  wga_coverage = "numeric",
  dn = "numeric",
  ds = "numeric",
  same_symbol = "logical",
  orthology_rank = "integer",
  pair_quality = "character"
)

ensure_dir(output_dir)
ensure_dir(figure_dir)

annotation_payload <- read_reference_annotation_local(clean_gtf, reference_gtf)
gtf_path <- annotation_payload$gtf_path
annotation_df <- annotation_payload$annotation_df
genes <- sort(unique(annotation_df$gene_name[nzchar(annotation_df$gene_name)]))

if (length(genes) == 0) {
  stop("GTF 中没有可用的 gene_name。", call. = FALSE)
}

chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))
message("开始构建同源映射缓存。")
message("GTF: ", gtf_path)
message("需要查询的鸡基因数: ", length(genes))

manifest_outputs <- list()

for (species_key in target_species) {
  message("处理目标物种: ", species_key)
  res_list <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    res_list[[i]] <- query_one_chunk_with_retry(
      gene_chunk = chunks[[i]],
      chunk_idx = i,
      total_chunks = length(chunks),
      species_key = species_key,
      mirrors = mirrors
    )
  }

  ortholog_raw <- dplyr::bind_rows(res_list)
  ortholog_all <- standardize_candidates(ortholog_raw, species_key)
  ortholog_best <- select_best_per_gene(ortholog_all)

  all_path <- file.path(output_dir, sprintf("chicken_%s_orthologs_all_candidates.csv", species_key))
  best_path <- file.path(output_dir, sprintf("chicken_%s_orthologs.csv", species_key))

  write.csv(ortholog_all, all_path, row.names = FALSE)
  write.csv(ortholog_best, best_path, row.names = FALSE)

  manifest_outputs[[paste0(species_key, "_best")]] <- build_output_entry(
    path = best_path,
    type = "csv",
    produced_by = module_name,
    row_semantics = "one row per chicken gene (best hit)",
    schema = schema_ortholog
  )
  manifest_outputs[[paste0(species_key, "_all")]] <- build_output_entry(
    path = all_path,
    type = "csv",
    produced_by = module_name,
    row_semantics = "one row per chicken-target ortholog pair",
    schema = schema_ortholog
  )

  message(sprintf(
    "[%s] best hits: %d; candidate pairs: %d",
    species_key,
    nrow(ortholog_best),
    nrow(ortholog_all)
  ))
}

write_manifest_local(
  manifest_path = manifest_path,
  new_outputs = manifest_outputs,
  module_name = module_contract,
  base_dir = output_dir,
  inputs = list(
    gtf_path = gtf_path,
    target_species = target_species,
    mirrors = mirrors,
    chunk_size = chunk_size
  )
)

message("00a 完成。输出目录: ", output_dir)
