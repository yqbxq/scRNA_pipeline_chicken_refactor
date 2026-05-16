source_utf8 <- function(path) {
  source(path, encoding = "UTF-8")
}

pipeline_root_from_env <- function() {
  root <- Sys.getenv("PIPELINE_ROOT", unset = "")
  if (nzchar(root)) {
    return(normalizePath(root, winslash = "/", mustWork = FALSE))
  }
  if (exists("PIPELINE_ROOT", inherits = TRUE)) {
    root <- get("PIPELINE_ROOT", inherits = TRUE)
    if (length(root) > 0 && nzchar(root[[1]])) {
      return(normalizePath(root[[1]], winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(file.path(getwd()), winslash = "/", mustWork = FALSE)
}

PIPELINE_ROOT <- pipeline_root_from_env()
single_script_helper_dir <- file.path(PIPELINE_ROOT, "workflow", "05single_script", "helpers")

for (helper in c("runtime_utils.R", "manifest_utils.R", "report_utils.R", "metadata_io.R")) {
  helper_path <- file.path(single_script_helper_dir, helper)
  if (file.exists(helper_path)) {
    source_utf8(helper_path)
  }
}
