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

format_report_bullets_v04 <- function(items) {
  items <- as.character(items)
  items <- items[nzchar(items)]
  if (length(items) == 0) {
    return(character())
  }
  ifelse(grepl("^\\s*-\\s+", items), items, sprintf("- %s", items))
}

format_key_files_v04 <- function(key_files) {
  if (is.null(key_files) || length(key_files) == 0) {
    return("- No key files recorded.")
  }
  key_names <- names(key_files)
  if (is.null(key_names)) {
    key_names <- rep("", length(key_files))
  }
  lines <- character()
  for (i in seq_along(key_files)) {
    label <- key_names[[i]]
    path <- key_files[[i]]
    path <- paste(as.character(path), collapse = ", ")
    if (nzchar(label)) {
      lines <- c(lines, sprintf("- `%s`: `%s`", label, path))
    } else {
      lines <- c(lines, sprintf("- `%s`", path))
    }
  }
  lines
}

render_triage_section_v04 <- function(triage_df) {
  if (is.null(triage_df)) {
    return(character())
  }
  if (exists("render_triage_markdown", mode = "function")) {
    return(render_triage_markdown(triage_df))
  }
  render_markdown_table_local(triage_df)
}

build_report_lines_v04 <- function(
  title,
  header_bullets = character(),
  key_files = list(),
  review_focus = character(),
  extra_sections = list(),
  triage_df = NULL
) {
  title_line <- as.character(title)[[1]]
  if (!grepl("^#\\s+", title_line)) {
    title_line <- sprintf("# %s", title_line)
  }

  lines <- c(title_line)
  header_lines <- format_report_bullets_v04(header_bullets)
  if (length(header_lines) > 0) {
    lines <- c(lines, "", header_lines)
  }

  lines <- c(lines, "", "## Key Files", format_key_files_v04(key_files))

  focus_lines <- format_report_bullets_v04(review_focus)
  if (length(focus_lines) == 0) {
    focus_lines <- "- No specific review focus recorded."
  }
  lines <- c(lines, "", "## Review Focus", focus_lines)

  if (!is.null(extra_sections) && length(extra_sections) > 0) {
    section_names <- names(extra_sections)
    if (is.null(section_names)) {
      section_names <- rep("", length(extra_sections))
    }
    for (i in seq_along(extra_sections)) {
      section_title <- section_names[[i]]
      section_lines <- as.character(extra_sections[[i]])
      if (!nzchar(section_title)) {
        next
      }
      lines <- c(lines, "", sprintf("## %s", section_title), section_lines)
    }
  }

  triage_lines <- render_triage_section_v04(triage_df)
  if (length(triage_lines) > 0) {
    lines <- c(lines, "", "## Triage Summary", triage_lines)
  }

  lines
}
