read_ortholog_lut_07 <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(data.frame(chicken_symbol = character(0), human_symbol = character(0), stringsAsFactors = FALSE))
  }
  lut <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")
  chicken_col <- intersect(c("chicken_symbol", "chicken_gene", "source_symbol", "external_gene_name"), colnames(lut))[1]
  human_col <- intersect(c("human_symbol", "human_gene", "target_symbol", "hsapiens_homolog_associated_gene_name"), colnames(lut))[1]
  if (is.na(chicken_col) || is.na(human_col)) {
    stop(sprintf("Ortholog LUT missing chicken/human columns: %s", paste(colnames(lut), collapse = ", ")), call. = FALSE)
  }
  out <- data.frame(
    chicken_symbol = trimws(as.character(lut[[chicken_col]])),
    human_symbol = trimws(as.character(lut[[human_col]])),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$chicken_symbol) & nzchar(out$human_symbol) & out$human_symbol != "-", , drop = FALSE]
  out <- out[!duplicated(toupper(out$chicken_symbol)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

ortholog_chicken_to_human_vec <- function(symbols, lut) {
  values <- trimws(as.character(symbols))
  out <- values
  if (is.character(lut)) {
    lut <- read_ortholog_lut_07(lut)
  }
  if (is.data.frame(lut) && nrow(lut) > 0) {
    lookup <- stats::setNames(lut$human_symbol, toupper(lut$chicken_symbol))
    hit <- unname(lookup[toupper(values)])
    out[!is.na(hit) & nzchar(hit)] <- hit[!is.na(hit) & nzchar(hit)]
  }
  out[is.na(out) | out == "-"] <- ""
  out
}

ortholog_chicken_to_human_complex_vec <- function(values, lut) {
  vapply(values, function(value) {
    value <- trimws(as.character(value))
    if (!nzchar(value) || is.na(value) || value == "-") {
      return("")
    }
    parts <- unlist(strsplit(value, "[_+/]", perl = TRUE), use.names = FALSE)
    parts <- trimws(parts[nzchar(trimws(parts))])
    paste(ortholog_chicken_to_human_vec(parts, lut), collapse = "_")
  }, character(1), USE.NAMES = FALSE)
}
