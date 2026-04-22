source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(biomaRt)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

input_rds <- file.path(cfg$checkpoint_dir, "03_after_annotation.rds")
if (!file.exists(input_rds)) {
  stop("缺少输入检查点: 03_after_annotation.rds", call. = FALSE)
}

ortholog_mode <- tolower(trimws(Sys.getenv("SCENIC_ORTHOLOG_MODE", "strict")))
if (!ortholog_mode %in% c("strict", "relaxed")) {
  stop("SCENIC_ORTHOLOG_MODE 仅支持 strict 或 relaxed", call. = FALSE)
}

output_dir <- Sys.getenv("ORTHOLOG_CACHE_DIR", file.path(cfg$results_dir, "ortholog_cache"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
obj <- maybe_join_layers(obj)
counts <- get_assay_matrix(obj, type = "counts")
genes <- unique(rownames(counts))
genes <- genes[!is.na(genes) & genes != ""]

query_one_chunk <- function(mart, gene_chunk) {
  getBM(
    attributes = c(
      "external_gene_name",
      "ensembl_gene_id",
      "hsapiens_homolog_associated_gene_name",
      "hsapiens_homolog_ensembl_gene",
      "hsapiens_homolog_orthology_type",
      "hsapiens_homolog_orthology_confidence"
    ),
    filters = "external_gene_name",
    values = gene_chunk,
    mart = mart
  )
}

connect_chicken_mart <- function(mirror_name) {
  message("连接 Ensembl 鸡数据库镜像: ", mirror_name)
  useEnsembl(
    biomart = "genes",
    dataset = "ggallus_gene_ensembl",
    mirror = mirror_name
  )
}

query_one_chunk_with_retry <- function(gene_chunk, chunk_idx, total_chunks, mirrors, max_attempts = 6L) {
  last_err <- NULL
  for (attempt in seq_len(max_attempts)) {
    mirror_name <- mirrors[((attempt - 1L) %% length(mirrors)) + 1L]
    message(sprintf(
      "  chunk %d/%d attempt %d/%d via %s",
      chunk_idx, total_chunks, attempt, max_attempts, mirror_name
    ))
    res <- tryCatch({
      mart <- connect_chicken_mart(mirror_name)
      query_one_chunk(mart, gene_chunk)
    }, error = function(e) {
      last_err <<- e
      NULL
    })
    if (!is.null(res)) {
      return(res)
    }
    message("    failed: ", conditionMessage(last_err))
    Sys.sleep(min(5L * attempt, 20L))
  }
  stop(
    sprintf(
      "chunk %d 在 %d 次尝试后仍失败。最后错误: %s",
      chunk_idx, max_attempts, conditionMessage(last_err)
    ),
    call. = FALSE
  )
}

normalize_key <- function(x) {
  toupper(trimws(as.character(x)))
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

chunk_size <- as.integer(Sys.getenv("ORTHOLOG_QUERY_CHUNK_SIZE", "500"))
chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))
mirrors <- unique(c(cfg$ensembl_mirror, "www", "useast"))

message(sprintf("开始构建同源映射缓存，模式: %s", ortholog_mode))
message("需要查询的鸡基因数: ", length(genes))
res_list <- vector("list", length(chunks))
for (i in seq_along(chunks)) {
  res_list[[i]] <- query_one_chunk_with_retry(
    gene_chunk = chunks[[i]],
    chunk_idx = i,
    total_chunks = length(chunks),
    mirrors = mirrors
  )
}

ortholog_raw <- bind_rows(res_list)

ortholog_all <- ortholog_raw %>%
  transmute(
    external_gene_name = trimws(external_gene_name),
    ensembl_gene_id = trimws(ensembl_gene_id),
    hsapiens_homolog_associated_gene_name = trimws(hsapiens_homolog_associated_gene_name),
    hsapiens_homolog_ensembl_gene = trimws(hsapiens_homolog_ensembl_gene),
    orthology_type = trimws(hsapiens_homolog_orthology_type),
    orthology_confidence = suppressWarnings(as.numeric(trimws(as.character(hsapiens_homolog_orthology_confidence))))
  ) %>%
  filter(
    external_gene_name != "",
    !is.na(external_gene_name),
    hsapiens_homolog_associated_gene_name != "",
    !is.na(hsapiens_homolog_associated_gene_name)
  ) %>%
  mutate(
    chicken_symbol_upper = normalize_key(external_gene_name),
    human_symbol_upper = normalize_key(hsapiens_homolog_associated_gene_name),
    same_symbol = chicken_symbol_upper == human_symbol_upper,
    orthology_rank = orthology_rank(orthology_type)
  ) %>%
  distinct(
    external_gene_name,
    hsapiens_homolog_associated_gene_name,
    ensembl_gene_id,
    hsapiens_homolog_ensembl_gene,
    orthology_type,
    orthology_confidence,
    .keep_all = TRUE
  )

selected_tbl <- if (ortholog_mode == "strict") {
  ortholog_all %>%
    filter(
      orthology_type == "ortholog_one2one",
      ensembl_gene_id != "",
      !is.na(ensembl_gene_id),
      hsapiens_homolog_ensembl_gene != "",
      !is.na(hsapiens_homolog_ensembl_gene)
    ) %>%
    distinct(external_gene_name, hsapiens_homolog_associated_gene_name, .keep_all = TRUE) %>%
    distinct(external_gene_name, .keep_all = TRUE) %>%
    distinct(hsapiens_homolog_associated_gene_name, .keep_all = TRUE)
} else {
  ortholog_all %>%
    arrange(
      external_gene_name,
      desc(same_symbol),
      orthology_rank,
      desc(ifelse(is.na(orthology_confidence), -1, orthology_confidence)),
      hsapiens_homolog_associated_gene_name
    ) %>%
    group_by(external_gene_name) %>%
    slice_head(n = 1) %>%
    ungroup()
}

coverage_tbl <- data.frame(
  ortholog_mode = ortholog_mode,
  total_input_genes = length(unique(genes)),
  mapped_genes = nrow(selected_tbl),
  coverage_fraction = nrow(selected_tbl) / length(unique(genes)),
  stringsAsFactors = FALSE
)

type_summary <- selected_tbl %>%
  count(orthology_type, sort = TRUE)

pipeline_map <- selected_tbl %>%
  select(external_gene_name, hsapiens_homolog_associated_gene_name)

write.csv(
  selected_tbl %>%
    select(
      external_gene_name,
      hsapiens_homolog_associated_gene_name,
      ensembl_gene_id,
      hsapiens_homolog_ensembl_gene,
      orthology_type,
      orthology_confidence
    ),
  file.path(output_dir, "chicken_human_orthologs.csv"),
  row.names = FALSE
)
saveRDS(
  selected_tbl %>%
    select(
      external_gene_name,
      hsapiens_homolog_associated_gene_name,
      ensembl_gene_id,
      hsapiens_homolog_ensembl_gene,
      orthology_type,
      orthology_confidence
    ),
  file.path(output_dir, "chicken_human_orthologs.rds")
)
write.csv(
  pipeline_map,
  file.path(output_dir, "chicken_human_orthologs_for_pipeline.csv"),
  row.names = FALSE
)
write.csv(coverage_tbl, file.path(output_dir, "ortholog_coverage.csv"), row.names = FALSE)
write.csv(type_summary, file.path(output_dir, "ortholog_type_summary.csv"), row.names = FALSE)

if (ortholog_mode == "relaxed") {
  write.csv(
    ortholog_all %>%
      select(
        external_gene_name,
        hsapiens_homolog_associated_gene_name,
        ensembl_gene_id,
        hsapiens_homolog_ensembl_gene,
        orthology_type,
        orthology_confidence
      ),
    file.path(output_dir, "chicken_human_orthologs_all_candidates.csv"),
    row.names = FALSE
  )
}

coverage_txt <- c(
  sprintf("ortholog_mode=%s", ortholog_mode),
  sprintf("total_input_genes=%d", coverage_tbl$total_input_genes),
  sprintf("mapped_genes=%d", coverage_tbl$mapped_genes),
  sprintf("coverage_fraction=%.6f", coverage_tbl$coverage_fraction)
)
writeLines(coverage_txt, file.path(output_dir, "ortholog_coverage.txt"))

key_genes <- c("DRGX", "CREB3L2", "EMX2", "CEBPB", "FOSL2", "JUN", "PBX3", "HMGA1", "SREBF2")
key_hits <- selected_tbl %>%
  filter(
    normalize_key(external_gene_name) %in% key_genes |
      normalize_key(hsapiens_homolog_associated_gene_name) %in% key_genes
  ) %>%
  select(
    external_gene_name,
    hsapiens_homolog_associated_gene_name,
    ensembl_gene_id,
    hsapiens_homolog_ensembl_gene,
    orthology_type,
    orthology_confidence
  )
write.csv(key_hits, file.path(output_dir, "ortholog_key_gene_hits.csv"), row.names = FALSE)

message("同源映射缓存构建完成。")
message(sprintf(
  "Mode=%s, Coverage=%d / %d = %.4f",
  ortholog_mode,
  coverage_tbl$mapped_genes,
  coverage_tbl$total_input_genes,
  coverage_tbl$coverage_fraction
))
