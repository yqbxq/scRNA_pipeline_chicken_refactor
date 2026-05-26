env_or_default_07c <- function(name, default = "") {
  if (exists("env_or_default_03", mode = "function")) {
    return(env_or_default_03(name, default))
  }
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

multinichenet_available_07c <- function() {
  requireNamespace("multinichenetr", quietly = TRUE) &&
    requireNamespace("SingleCellExperiment", quietly = TRUE)
}

determine_07c_mode <- function(cfg, pairs = NULL) {
  mode <- tolower(trimws(env_or_default_07c("NICHENET_MODE", "auto")))
  enabled <- tolower(trimws(env_or_default_07c("MULTINICHENET_ENABLED", "auto")))
  if (identical(mode, "skip") || identical(enabled, "no")) {
    return("skip")
  }
  if (identical(mode, "nichenet_legacy")) {
    return("nichenet_legacy")
  }
  if (identical(mode, "multinichenet")) {
    if (!multinichenet_available_07c()) {
      stop("NICHENET_MODE=multinichenet but multinichenetr/SingleCellExperiment is unavailable.", call. = FALSE)
    }
    return("multinichenet")
  }
  if (multinichenet_available_07c()) "multinichenet" else "nichenet_legacy"
}

run_multinichenet_for_pair <- function(pair_row, cfg) {
  stop("run_multinichenet_for_pair requires the full P-R03-C runtime fixture and is not invoked by the legacy fallback path yet.", call. = FALSE)
}

compute_min_samples_per_group_07c <- function(sample_sheet) {
  if (!file.exists(sample_sheet)) {
    return(NA_integer_)
  }
  samples <- read.delim(sample_sheet, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  group_col <- intersect(c("condition", "group_id", "stage"), colnames(samples))[1]
  if (is.na(group_col)) {
    return(NA_integer_)
  }
  min(tabulate(match(samples[[group_col]], unique(samples[[group_col]]))))
}
