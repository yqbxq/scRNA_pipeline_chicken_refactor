read_reference_inventory <- function(path = Sys.getenv("SPATIAL_REFERENCE_INVENTORY_FILE", unset = "")) {
  if (!nzchar(path) || !file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  read_tsv_optional(path)
}

resolve_local_reference_path <- function(inv) {
  if (is.null(inv) || nrow(inv) == 0 || !"selected_reference" %in% colnames(inv)) {
    stop("spatial_reference_inventory.tsv is empty or missing selected_reference", call. = FALSE)
  }
  selected <- trimws(as.character(inv$selected_reference[[1]]))
  if (!nzchar(selected)) {
    stop("spatial_reference_inventory.tsv has empty selected_reference", call. = FALSE)
  }
  if (!grepl("^/", selected)) {
    project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
    selected <- file.path(project_root, selected)
  }
  if (!file.exists(selected)) {
    stop(
      "scRNA reference not found; run scRNA pipeline first or update spatial_reference_inventory.tsv",
      call. = FALSE
    )
  }
  normalizePath(selected, winslash = "/", mustWork = FALSE)
}
