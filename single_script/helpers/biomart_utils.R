species_prefix <- function(species_key) {
  switch(
    species_key,
    human = "hsapiens",
    mouse = "mmusculus",
    stop(sprintf("不支持的目标物种: %s", species_key), call. = FALSE)
  )
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
