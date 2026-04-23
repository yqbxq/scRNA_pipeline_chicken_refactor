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

fmt_num <- function(x, digits = 1) {
  if (is.na(x)) {
    return("NA")
  }
  sprintf(paste0("%.", digits, "f"), x)
}

safe_rate <- function(numerator, denominator) {
  if (length(numerator) > 1) {
    return(ifelse(is.na(denominator) | denominator <= 0, NA_real_, numerator / denominator))
  }
  if (is.na(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  numerator / denominator
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  stats::median(x, na.rm = TRUE)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE))
}

save_plot_local <- function(plot_obj, path, width, height) {
  ggplot2::ggsave(path, plot_obj, width = width, height = height, dpi = 300, bg = "white")
}
