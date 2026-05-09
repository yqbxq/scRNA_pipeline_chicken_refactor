project_feature_palette <- function() {
  c(
    "#2C7BB6", "#00A6CA", "#00CCBC", "#90EB9D",
    "#FFFF8C", "#F9D057", "#F29E2E", "#E76818", "#D7191C"
  )
}

paper_feature_palette <- project_feature_palette

save_plot_dual <- function(plot_obj, png_path, width, height, dpi = 300, pdf_path = NULL) {
  if (is.null(pdf_path)) {
    pdf_path <- sub("\\.png$", ".pdf", png_path)
  }
  ensure_dir(dirname(png_path))
  ggplot2::ggsave(png_path, plot_obj, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(pdf_path, plot_obj, width = width, height = height, bg = "white")
}
