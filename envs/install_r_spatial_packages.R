options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

ensure_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

ensure_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
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

ensure_github <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install_local_git_repo(repo)
  }
}

ensure_decoupleR <- function() {
  if (requireNamespace("decoupleR", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  tryCatch(
    {
      install_local_git_repo("saezlab/decoupleR")
    },
    error = function(github_err) {
      message("Local git install for decoupleR failed, retrying via BiocManager...")
      BiocManager::install("decoupleR", ask = FALSE, update = FALSE)
    }
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
  "Seurat", "SeuratObject", "hdf5r", "reticulate", "yaml", "dplyr",
  "tibble", "tidyr", "readr", "ggplot2", "patchwork", "cowplot",
  "ggrepel", "data.table", "broom", "dbscan", "sctransform"
)) {
  ensure_cran(pkg)
}

for (pkg in c(
  "SingleCellExperiment", "SpatialExperiment", "SummarizedExperiment",
  "BiocParallel", "BiocNeighbors", "BayesSpace", "spatialDE"
)) {
  ensure_bioc(pkg)
}

ensure_decoupleR()

ensure_github("spacexr", "dmcable/spacexr")
ensure_github("CARD", "YingMa0107/CARD")
ensure_github("SPARK", "xzhoulab/SPARK")

required_pkgs <- c(
  "Seurat", "SeuratObject", "hdf5r", "reticulate", "yaml", "dplyr",
  "tibble", "tidyr", "readr", "ggplot2", "patchwork", "cowplot",
  "ggrepel", "data.table", "dbscan", "sctransform", "SingleCellExperiment", "SpatialExperiment",
  "SummarizedExperiment", "BiocParallel", "BayesSpace", "decoupleR",
  "spatialDE", "spacexr", "CARD", "SPARK"
)
assert_pkgs(required_pkgs, "空间分析")

installed <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
installed <- installed[order(installed$Package), ]
write.csv(installed, file.path(Sys.getenv("LOG_DIR"), "R_spatial_installed_packages.csv"), row.names = FALSE)
writeLines(utils::capture.output(sessionInfo()), file.path(Sys.getenv("LOG_DIR"), "R_spatial_sessionInfo.txt"))
