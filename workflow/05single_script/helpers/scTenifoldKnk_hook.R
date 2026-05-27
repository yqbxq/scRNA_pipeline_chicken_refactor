scTenifoldKnk_hook_path <- function(cfg) {
  cfg$scTenifoldKnk_hook_results_tsv %||%
    file.path(cfg$communication_table_dir %||% file.path(cfg$table_dir, "communication"), "scTenifoldKnk_hook_results.tsv")
}

render_scTenifoldKnk_section_with_caveats <- function(path) {
  stop(
    sprintf(
      "scTenifoldKnk hook results were found at %s, but scTenifoldKnk interpretation is not in scope for this pipeline. Review docs/scTenifoldKnk_hook_disclaimer.md before enabling this section.",
      path
    ),
    call. = FALSE
  )
}

render_panel_scTenifoldKnk_hook <- function(consensus_df, cfg) {
  path <- scTenifoldKnk_hook_path(cfg)
  if (file.exists(path)) {
    render_scTenifoldKnk_section_with_caveats(path)
  }
  tsv <- data.frame(
    status = "not_enabled",
    path = path,
    reason = "scTenifoldKnk is not part of the current communication evidence chain.",
    stringsAsFactors = FALSE
  )
  list(
    name = "scTenifoldKnk_hook",
    status = "not_enabled",
    md_lines = c(
      "## Panel 4: scTenifoldKnk Hook",
      "",
      "scTenifoldKnk is not enabled. If hook results appear, 07e fails loudly until the method caveat is reviewed.",
      "",
      communication_report_md_table(tsv, max_rows = 5L)
    ),
    plots = list(),
    tsv = tsv
  )
}
