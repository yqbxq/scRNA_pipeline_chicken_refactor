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

cran_pkgs <- c(
  "Seurat",
  "hdf5r",
  "dplyr",
  "tibble",
  "tidyr",
  "readr",
  "ggplot2",
  "patchwork",
  "pheatmap",
  "yaml",
  "igraph",
  "VGAM"
)
for (pkg in cran_pkgs) {
  ensure_cran(pkg)
}

ensure_bioc("Biobase")

if (!requireNamespace("DDRTree", quietly = TRUE)) {
  remotes::install_github(
    "cole-trapnell-lab/DDRTree",
    upgrade = "never",
    dependencies = FALSE,
    build_vignettes = FALSE
  )
}

if (!requireNamespace("monocle", quietly = TRUE)) {
  BiocManager::install("monocle", ask = FALSE, update = FALSE)
}

if (!requireNamespace("topGO", quietly = TRUE)) {
  BiocManager::install("topGO", ask = FALSE, update = FALSE)
}

if (!requireNamespace("phateR", quietly = TRUE)) {
  remotes::install_github(
    "KlugerLab/phateR",
    upgrade = "never",
    dependencies = FALSE,
    build_vignettes = FALSE
  )
}

required_pkgs <- c(
  "Seurat",
  "hdf5r",
  "dplyr",
  "tibble",
  "tidyr",
  "readr",
  "ggplot2",
  "patchwork",
  "pheatmap",
  "yaml",
  "igraph",
  "VGAM",
  "Biobase",
  "DDRTree",
  "monocle",
  "topGO",
  "phateR"
)
assert_pkgs(required_pkgs, "legacy 轨迹")

installed <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
installed <- installed[order(installed$Package), ]
write.csv(installed, file.path(Sys.getenv("LOG_DIR"), "R_legacy_installed_packages.csv"), row.names = FALSE)
writeLines(utils::capture.output(sessionInfo()), file.path(Sys.getenv("LOG_DIR"), "R_legacy_sessionInfo.txt"))
