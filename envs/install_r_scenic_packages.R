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

ensure_optional_github <- function(pkg, repo, ref = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }

  tryCatch(
    {
      remotes::install_github(
        repo,
        ref = ref,
        upgrade = "never",
        dependencies = FALSE,
        build_vignettes = FALSE
      )
    },
    error = function(e) {
      warning(sprintf("Failed to install optional package %s: %s", pkg, conditionMessage(e)))
      log_dir <- Sys.getenv("LOG_DIR", "")
      if (nzchar(log_dir)) {
        dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
        writeLines(
          c(
            sprintf("Optional package installation failed: %s", pkg),
            conditionMessage(e)
          ),
          file.path(log_dir, sprintf("optional_%s_install_warning.txt", pkg))
        )
      }
    }
  )
}

assert_pkgs <- function(pkgs, label) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      sprintf(
        "Missing %s packages after installation: %s",
        label,
        paste(missing_pkgs, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

ensure_cran("remotes")
ensure_cran("BiocManager")

cran_pkgs <- c(
  "Seurat",
  "dplyr",
  "tibble",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "circlize",
  "viridis",
  "data.table",
  "getopt",
  "philentropy"
)
for (pkg in cran_pkgs) {
  ensure_cran(pkg)
}

bioc_pkgs <- c(
  "AUCell",
  "RcisTarget",
  "GENIE3",
  "ComplexHeatmap",
  "fgsea",
  "Biobase"
)
for (pkg in bioc_pkgs) {
  ensure_bioc(pkg)
}

if (!requireNamespace("SCENIC", quietly = TRUE)) {
  remotes::install_github(
    "aertslab/SCENIC",
    ref = "v1.3.0",
    upgrade = "never",
    dependencies = FALSE,
    build_vignettes = FALSE
  )
}

ensure_optional_github(
  "scFunctions",
  "FloWuenne/scFunctions",
  ref = "07c34da4db10d63d8123b21685d8974f98740948"
)

core_pkgs <- c(
  "Seurat",
  "dplyr",
  "tibble",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "circlize",
  "viridis",
  "data.table",
  "philentropy",
  "AUCell",
  "RcisTarget",
  "GENIE3",
  "ComplexHeatmap",
  "fgsea",
  "Biobase",
  "SCENIC"
)
assert_pkgs(core_pkgs, "core SCENIC")

installed <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
installed <- installed[order(installed$Package), ]
write.csv(installed, file.path(Sys.getenv("LOG_DIR"), "R_scenic_installed_packages.csv"), row.names = FALSE)
writeLines(utils::capture.output(sessionInfo()), file.path(Sys.getenv("LOG_DIR"), "R_scenic_sessionInfo.txt"))
