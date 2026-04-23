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
