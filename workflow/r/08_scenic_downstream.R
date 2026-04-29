message("[DEPRECATED] workflow/r/08_scenic_downstream.R now delegates to shell/05single_script/08d_scenic_downstream.R")

pipeline_root <- Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT", unset = ""))
if (!nzchar(pipeline_root)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    pipeline_root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..", ".."), winslash = "/", mustWork = FALSE)
  } else {
    pipeline_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
}

source(file.path(pipeline_root, "shell", "05single_script", "08d_scenic_downstream.R"), encoding = "UTF-8")
