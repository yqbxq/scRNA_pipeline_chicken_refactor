options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

ensure_pkg <- function(pkg, installer) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    installer()
  }
}

require_pkgs <- function(pkgs, label) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      sprintf(
        "缺少%s依赖包: %s。请先重新运行 workflow/01_install_envs.sh 更新 conda 环境。",
        label,
        paste(missing_pkgs, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

ensure_pkg("remotes", function() install.packages("remotes"))
ensure_pkg("gprofiler2", function() install.packages("gprofiler2"))
ensure_pkg(
  "presto",
  function() remotes::install_github(
    "immunogenomics/presto",
    upgrade = "never",
    dependencies = FALSE,
    build_vignettes = FALSE
  )
)

# 主流程核心依赖优先由 conda 提供，避免在 R 里追最新版造成代系漂移。
required_main_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "harmony",
  "dplyr",
  "tibble",
  "tidyr",
  "readr",
  "patchwork",
  "cowplot",
  "ggrepel",
  "readxl",
  "hdf5r",
  "reticulate",
  "yaml",
  "presto",
  "RColorBrewer",
  "circlize",
  "pheatmap",
  "future",
  "future.apply",
  "BiocParallel",
  "data.table",
  "SingleCellExperiment",
  "scater",
  "scDblFinder",
  "celda",
  "edgeR",
  "limma",
  "muscat",
  "speckle",
  "slingshot",
  "tradeSeq",
  "monocle3",
  "clusterProfiler",
  "gprofiler2",
  "org.Gg.eg.db",
  "biomaRt",
  "fields",
  "ROCR"
)
require_pkgs(required_main_pkgs, "主流程")

if (!requireNamespace("SoupX", quietly = TRUE)) {
  remotes::install_github(
    "constantAmateur/SoupX",
    upgrade = "never",
    dependencies = FALSE,
    build_vignettes = FALSE
  )
}
require_pkgs("SoupX", "SoupX")

# DoubletFinder 在 conda 中没有现成包，固定到当前可解析的 GitHub 提交。
if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  remotes::install_github(
    "chris-mcginnis-ucsf/DoubletFinder",
    ref = "1b244d8f0d54b4b1cb4365639931bbb16f01e1cd",
    upgrade = "never",
    dependencies = FALSE,
    build_vignettes = FALSE
  )
}
require_pkgs("DoubletFinder", "DoubletFinder")

installed <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
installed <- installed[order(installed$Package), ]
write.csv(installed, file.path(Sys.getenv("LOG_DIR"), "R_main_installed_packages.csv"), row.names = FALSE)
writeLines(utils::capture.output(sessionInfo()), file.path(Sys.getenv("LOG_DIR"), "R_main_sessionInfo.txt"))
