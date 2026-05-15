#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else sys.frame(1)$ofile
script_dir <- dirname(normalizePath(script_path))
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
helper_dir <- file.path(repo_root, "workflow", "05single_script", "helpers")

source(file.path(helper_dir, "runtime_utils.R"), encoding = "UTF-8")
source(file.path(helper_dir, "scdesign3_engine_helpers.R"), encoding = "UTF-8")

tmp <- tempfile("scdesign3_mock_")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

counts <- matrix(
  c(
    8, 9, 7, 8, 1, 1, 1, 2, 0, 0, 1, 0,
    7, 8, 8, 9, 1, 0, 1, 1, 0, 1, 0, 0,
    0, 1, 0, 0, 8, 8, 9, 7, 1, 1, 0, 1,
    1, 0, 1, 0, 7, 9, 8, 8, 0, 0, 1, 1,
    0, 0, 1, 0, 1, 0, 1, 0, 8, 9, 8, 7,
    1, 0, 0, 1, 0, 1, 0, 1, 9, 8, 7, 8
  ),
  nrow = 6,
  byrow = TRUE
)
rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
truth <- rep(c("A", "B", "C"), each = 4)
input_path <- file.path(tmp, "mock_input.rds")
saveRDS(list(counts = counts, meta = data.frame(cell_type = truth, row.names = colnames(counts))), input_path)

scd_missing_engine_packages <- function() character(0)
scd_prepare_target_data <- function(path, truth_col, max_cells_per_label, n_hvg, seed) {
  stopifnot(identical(path, input_path))
  stopifnot(max_cells_per_label == 12L)
  stopifnot(n_hvg == 6L)
  list(
    counts = counts,
    truth = truth,
    truth_col = truth_col,
    n_real_cells = ncol(counts),
    n_truth_labels = length(unique(truth))
  )
}
fit_scdesign3 <- function(prepared, truth_col, family_use = "nb", n_cores = 1L) {
  list(prepared = prepared, truth_col = truth_col)
}
simulate_synthetic_counts <- function(fit, sim_i, seed) {
  list(counts = fit$prepared$counts, truth = fit$prepared$truth)
}
recluster_synthetic <- function(counts, resolution, expected_clusters, seed, n_pcs = 30L) {
  stopifnot(n_pcs == 4L)
  truth
}

stopifnot(identical(scd_gate_decision(c(ARI = 0.90, NMI = 0.80, min_Jaccard = 0.40), "ARI>=0.80", "ARI>=0.60;NMI>=0.60;min_Jaccard>=0.45", "ARI<0.60"), "PASS"))
stopifnot(identical(scd_gate_decision(c(ARI = 0.70, NMI = 0.65, min_Jaccard = 0.50), "ARI>=0.80", "ARI>=0.60;NMI>=0.60;min_Jaccard>=0.45", "ARI<0.60"), "WARN"))
stopifnot(identical(scd_gate_decision(c(ARI = 0.55, NMI = 0.90, min_Jaccard = 0.90), "ARI>=0.80", "ARI>=0.60;NMI>=0.60;min_Jaccard>=0.45", "ARI<0.60"), "FAIL"))

target <- data.frame(
  target_id = "SCD_CLUSTER_MOCK",
  target_type = "cluster_robustness",
  layer_id = "panorama",
  input_object = "mock",
  truth_col = "cell_type",
  questions_covered = "MOCK",
  n_simulations = "2",
  resolution_grid = "0.4",
  mixture_design = "observed_balanced",
  primary_metric = "ARI",
  pass_threshold = "ARI>=0.80",
  warn_threshold = "ARI>=0.60;NMI>=0.60;min_Jaccard>=0.45",
  fail_threshold = "ARI<0.60",
  max_cells_per_label = "12",
  n_hvg = "6",
  n_pcs = "4",
  output_dir = "",
  status = "active",
  resolved_input_path = input_path,
  stringsAsFactors = FALSE
)
engine_cfg <- list(
  seed = 42L,
  n_simulations_default = 1L,
  n_cores = 1L,
  family_use = "nb",
  max_cells_per_label = 2000L,
  n_hvg = 2000L,
  n_pcs = 30L,
  resolution_default = c(0.6),
  resolution_grid_override = "",
  checkpoint_dir = tmp,
  figure_root = file.path(tmp, "figures")
)

result <- run_cluster_target(target, engine_cfg)
stopifnot(identical(result$target_metrics$status[[1]], "ok"))
stopifnot(identical(result$target_metrics$gate_status[[1]], "PASS"))
stopifnot(identical(result$target_metrics$n_simulations_requested[[1]], "2"))
stopifnot(identical(result$target_metrics$n_simulations_done[[1]], "2"))
stopifnot(nrow(result$per_simulation) == 2L)
stopifnot(nrow(result$per_label) > 0)
stopifnot(file.exists(file.path(engine_cfg$figure_root, "SCD_CLUSTER_MOCK", "ARI_distribution.png")))

skip_target <- target
skip_target$target_type <- "composition_robustness"
skip_result <- run_cluster_target(skip_target, engine_cfg)
stopifnot(identical(skip_result$engine_status$status[[1]], "unsupported_m1"))
stopifnot(length(skip_result$figures) == 0L)

cat("smoke_scdesign3_engine_mock_ok\n")
