#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
)

source_utf8 <- function(path) {
  source(path, encoding = "UTF-8")
}

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "gtf_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ortholog_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "biomart_utils.R"))

load_required_packages(c("dplyr", "biomaRt", "jsonlite"))

cfg <- get_single_script_config()
module_name <- "00a_build_ortholog_cache"

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

standardize_candidates <- function(raw_df, species_key, annotation_df) {
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
  required_cols <- c(
    "external_gene_name",
    "ensembl_gene_id",
    name_col,
    ensembl_col,
    type_col,
    conf_col,
    perc_id_col,
    perc_id_r1_col,
    goc_col,
    wga_col,
    dn_col,
    ds_col
  )
  missing_cols <- setdiff(required_cols, names(raw_df))
  for (col_name in missing_cols) {
    raw_df[[col_name]] <- NA_character_
  }

  standardized <- raw_df %>%
    dplyr::transmute(
      external_gene_name = trimws(as.character(.data$external_gene_name)),
      ensembl_gene_id = strip_ensembl_version(trimws(as.character(.data$ensembl_gene_id))),
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

  annotation_lookup <- unique(annotation_df[, c("gene_id_stripped", "preferred_gene_label"), drop = FALSE])
  standardized <- standardized %>%
    dplyr::left_join(annotation_lookup, by = c("ensembl_gene_id" = "gene_id_stripped")) %>%
    dplyr::mutate(
      external_gene_name = dplyr::if_else(
        is.na(external_gene_name) | external_gene_name == "",
        preferred_gene_label,
        external_gene_name
      ),
      external_gene_name = dplyr::if_else(
        is.na(external_gene_name) | external_gene_name == "",
        ensembl_gene_id,
        external_gene_name
      )
    ) %>%
    dplyr::select(-preferred_gene_label)

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

ensure_dir(cfg$output_dir)
ensure_dir(cfg$figure_dir)

annotation_payload <- read_reference_annotation_local(cfg$clean_gtf, cfg$reference_gtf)
gtf_path <- annotation_payload$gtf_path
annotation_df <- annotation_payload$annotation_df
genes <- sort(unique(annotation_df$gene_id_stripped[nzchar(annotation_df$gene_id_stripped)]))

if (length(genes) == 0) {
  stop("GTF 中没有可用的 gene_id。", call. = FALSE)
}

chunks <- split(genes, ceiling(seq_along(genes) / cfg$chunk_size))
message("开始构建同源映射缓存。")
message("GTF: ", gtf_path)
message("需要查询的鸡基因数: ", length(genes))

manifest_outputs <- list()
schema_ortholog <- infer_schema_from_df(empty_candidate_table())

for (species_key in cfg$target_species) {
  message("处理目标物种: ", species_key)
  res_list <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    res_list[[i]] <- query_one_chunk_with_retry(
      gene_chunk = chunks[[i]],
      chunk_idx = i,
      total_chunks = length(chunks),
      species_key = species_key,
      mirrors = cfg$mirrors,
      filter_name = "ensembl_gene_id"
    )
  }

  ortholog_raw <- dplyr::bind_rows(res_list)
  ortholog_all <- standardize_candidates(ortholog_raw, species_key, annotation_df)
  ortholog_best <- select_best_per_gene(ortholog_all)

  all_path <- file.path(cfg$output_dir, sprintf("chicken_%s_orthologs_all_candidates.csv", species_key))
  best_path <- file.path(cfg$output_dir, sprintf("chicken_%s_orthologs.csv", species_key))

  write.csv(ortholog_all, all_path, row.names = FALSE)
  write.csv(ortholog_best, best_path, row.names = FALSE)

  manifest_outputs[[paste0(species_key, "_best")]] <- build_output_entry(
    path = best_path,
    type = "csv",
    produced_by = module_name,
    row_semantics = "one row per chicken gene (best hit)",
    base_dir = cfg$output_dir,
    schema = schema_ortholog
  )
  manifest_outputs[[paste0(species_key, "_all")]] <- build_output_entry(
    path = all_path,
    type = "csv",
    produced_by = module_name,
    row_semantics = "one row per chicken-target ortholog pair",
    base_dir = cfg$output_dir,
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
  manifest_path = cfg$manifest_path,
  new_outputs = manifest_outputs,
  module_name = cfg$module_contract,
  base_dir = cfg$output_dir,
  inputs = list(
    gtf_path = gtf_path,
    target_species = cfg$target_species,
    mirrors = cfg$mirrors,
    chunk_size = cfg$chunk_size
  ),
  version = cfg$module_version,
  depends_on = list()
)

message("00a 完成。输出目录: ", cfg$output_dir)
