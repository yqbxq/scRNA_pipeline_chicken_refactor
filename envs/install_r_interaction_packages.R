options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

ensure_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

ensure_cran_version <- function(pkg, min_version) {
  if (!requireNamespace(pkg, quietly = TRUE) || packageVersion(pkg) < min_version) {
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
ensure_cran("gprofiler2")
ensure_cran_version("NMF", "0.23.0")

for (pkg in c(
  "Seurat", "SeuratObject", "hdf5r", "yaml", "dplyr", "tibble", "tidyr",
  "readr", "ggplot2", "patchwork", "cowplot", "ggrepel", "data.table",
  "future", "future.apply", "igraph", "circlize", "RColorBrewer", "NMF",
  "broom", "collapse", "ggalluvial", "ggnetwork", "ggpubr", "sna",
  "svglite", "tidyverse", "fdrtool", "Hmisc", "caret", "randomForest",
  "DiagrammeR", "mlrMBO", "parallelMap", "emoa", "DiceKriging", "e1071",
  "shadowtext"
)) {
  ensure_cran(pkg)
}

for (pkg in c(
  "SingleCellExperiment", "clusterProfiler", "org.Gg.eg.db",
  "biomaRt", "ComplexHeatmap", "BiocNeighbors"
)) {
  ensure_bioc(pkg)
}

ensure_decoupleR()

ensure_github("presto", "immunogenomics/presto")
ensure_github("CellChat", "jinworks/CellChat")
ensure_github("nichenetr", "saeyslab/nichenetr")

required_pkgs <- c(
  "Seurat", "SeuratObject", "hdf5r", "yaml", "dplyr", "tibble", "tidyr",
  "readr", "ggplot2", "patchwork", "cowplot", "ggrepel", "data.table",
  "future", "future.apply", "igraph", "circlize", "RColorBrewer", "NMF",
  "SingleCellExperiment", "clusterProfiler", "org.Gg.eg.db", "biomaRt",
  "ComplexHeatmap", "decoupleR", "gprofiler2", "presto", "CellChat", "nichenetr"
)
assert_pkgs(required_pkgs, "通讯分析")

installed <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
installed <- installed[order(installed$Package), ]
write.csv(installed, file.path(Sys.getenv("LOG_DIR"), "R_interaction_installed_packages.csv"), row.names = FALSE)
writeLines(utils::capture.output(sessionInfo()), file.path(Sys.getenv("LOG_DIR"), "R_interaction_sessionInfo.txt"))
