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
  gene_name[!nzchar(gene_name)] <- extract_gtf_attr(attr_field, "gene")[!nzchar(gene_name)]
  gene_name[!nzchar(gene_name)] <- extract_gtf_attr(attr_field, "Name")[!nzchar(gene_name)]
  gene_biotype <- extract_gtf_attr(attr_field, "gene_biotype")
  gene_biotype[!nzchar(gene_biotype)] <- extract_gtf_attr(attr_field, "gene_type")[!nzchar(gene_biotype)]
  gene_biotype[!nzchar(gene_biotype)] <- extract_gtf_attr(attr_field, "biotype")[!nzchar(gene_biotype)]
  gene_id_stripped <- strip_ensembl_version(gene_id)

  n_total <- nrow(gene_rows)
  n_with_id <- sum(nzchar(gene_id))
  n_with_name <- sum(nzchar(gene_name))
  n_with_biotype <- sum(nzchar(gene_biotype))
  message(sprintf("GTF gene rows: %d", n_total))
  message(sprintf(
    "GTF gene_id populated: %d / %d (%.1f%%)",
    n_with_id,
    n_total,
    100 * n_with_id / n_total
  ))
  message(sprintf(
    "GTF gene_name populated: %d / %d (%.1f%%)",
    n_with_name,
    n_total,
    100 * n_with_name / n_total
  ))
  message(sprintf(
    "GTF gene_biotype populated: %d / %d (%.1f%%)",
    n_with_biotype,
    n_total,
    100 * n_with_biotype / n_total
  ))
  biotype_table <- sort(table(ifelse(nzchar(gene_biotype), gene_biotype, "unknown")), decreasing = TRUE)
  top_biotypes <- head(biotype_table, 10)
  message(sprintf(
    "GTF top biotypes: %s",
    paste(sprintf("%s=%s", names(top_biotypes), as.integer(top_biotypes)), collapse = ", ")
  ))
  if (n_with_id < n_total * 0.95) {
    stop("GTF gene_id 填充率低于 95%，请检查 GTF 格式。", call. = FALSE)
  }
  if (n_with_name < n_total * 0.5) {
    message("GTF gene_name 填充率低于 50%，将使用 Ensembl gene_id 作为查询主键，并用 gene_id 回填缺失 gene_name。")
  }

  preferred_gene_label <- gene_name
  preferred_gene_label[!nzchar(preferred_gene_label)] <- gene_id[!nzchar(preferred_gene_label)]
  preferred_gene_label[!nzchar(preferred_gene_label)] <- gene_id_stripped[!nzchar(preferred_gene_label)]

  annotation_df <- unique(data.frame(
    seqname = as.character(gene_rows[[1]]),
    gene_id = gene_id,
    gene_id_stripped = gene_id_stripped,
    gene_name = gene_name,
    gene_name_upper = toupper(gene_name),
    preferred_gene_label = preferred_gene_label,
    preferred_gene_label_upper = toupper(preferred_gene_label),
    gene_biotype = gene_biotype,
    stringsAsFactors = FALSE
  ))
  annotation_df <- annotation_df[nzchar(annotation_df$gene_id) | nzchar(annotation_df$gene_name), , drop = FALSE]

  list(gtf_path = gtf_path, annotation_df = annotation_df)
}
