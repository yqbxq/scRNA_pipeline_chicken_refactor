options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

ensure_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

install_local_git_repo <- function(repo) {
  clone_dir <- tempfile(pattern = "gitpkg_")
  on.exit(unlink(clone_dir, recursive = TRUE, force = TRUE), add = TRUE)

  clone_url <- sprintf("https://github.com/%s.git", repo)
  cloned <- FALSE
  last_status <- NULL

  for (attempt in seq_len(3)) {
    unlink(clone_dir, recursive = TRUE, force = TRUE)
    last_status <- system2(
      "git",
      c("-c", "http.version=HTTP/1.1", "clone", "--depth", "1", clone_url, clone_dir)
    )
    if (identical(last_status, 0L) && dir.exists(clone_dir)) {
      cloned <- TRUE
      break
    }
    Sys.sleep(attempt)
  }

  if (!cloned) {
    stop(sprintf("Failed to clone %s after repeated attempts", clone_url))
  }

  remotes::install_local(
    clone_dir,
    upgrade = "never",
    dependencies = c("Depends", "Imports", "LinkingTo"),
    build_vignettes = FALSE
  )
}

assert_pkgs <- function(pkgs, label) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      sprintf("缺少%s依赖包: %s", label, paste(missing_pkgs, collapse = ", ")),
      call. = FALSE
    )
  }
}

ensure_cran("remotes")
ensure_cran("BiocManager")

for (pkg in c(
  "Seurat", "SeuratObject", "dplyr", "tibble", "tidyr", "readr",
  "ggplot2", "patchwork", "cowplot", "data.table", "SingleCellExperiment"
)) {
  if (pkg == "SingleCellExperiment") {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  } else {
    ensure_cran(pkg)
  }
}

if (!requireNamespace("scDesign3", quietly = TRUE)) {
  install_local_git_repo("SONGDONGYUAN1994/scDesign3")
}

required_pkgs <- c(
  "Seurat", "SeuratObject", "dplyr", "tibble", "tidyr", "readr",
  "ggplot2", "patchwork", "cowplot", "data.table", "SingleCellExperiment",
  "scDesign3"
)
assert_pkgs(required_pkgs, "验证分析")

installed <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
installed <- installed[order(installed$Package), ]
write.csv(installed, file.path(Sys.getenv("LOG_DIR"), "R_validation_installed_packages.csv"), row.names = FALSE)
writeLines(utils::capture.output(sessionInfo()), file.path(Sys.getenv("LOG_DIR"), "R_validation_sessionInfo.txt"))
