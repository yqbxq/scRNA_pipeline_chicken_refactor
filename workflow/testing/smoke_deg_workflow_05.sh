#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [[ -z "${R_BIN:-}" ]]; then
  if [[ -n "${R_MAIN_ENV_PREFIX:-}" && -x "${R_MAIN_ENV_PREFIX}/bin/Rscript" ]]; then
    R_BIN="${R_MAIN_ENV_PREFIX}/bin/Rscript"
  elif command -v Rscript >/dev/null 2>&1; then
    R_BIN="$(command -v Rscript)"
  else
    R_BIN=""
  fi
fi

[[ -x "${R_BIN}" ]] || {
  echo "Missing Rscript for smoke_deg_workflow_05: ${R_BIN}" >&2
  exit 1
}

TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/scrna_deg05_smoke.XXXXXX")"
trap 'rm -rf "${TMP_ROOT}"' EXIT

PROJECT_ROOT="${TMP_ROOT}/project"
RESULTS_DIR="${PROJECT_ROOT}/results"
CHECKPOINT_DIR="${RESULTS_DIR}/checkpoints"
TABLE_DIR="${RESULTS_DIR}/tables"
FIGURE_DIR="${RESULTS_DIR}/figures"
REPORT_DIR="${PROJECT_ROOT}/reports"
EDA_REPORT_DIR="${REPORT_DIR}/eda"
MANIFEST_DIR="${RESULTS_DIR}/manifests"
METADATA_DIR="${PROJECT_ROOT}/metadata"
PROJECT_CONFIG_DIR="${PROJECT_ROOT}/config"
LOG_DIR="${PROJECT_ROOT}/logs"

mkdir -p \
  "${CHECKPOINT_DIR}/layers/panorama" \
  "${TABLE_DIR}" \
  "${FIGURE_DIR}" \
  "${EDA_REPORT_DIR}" \
  "${MANIFEST_DIR}" \
  "${METADATA_DIR}" \
  "${PROJECT_CONFIG_DIR}" \
  "${LOG_DIR}"

export PIPELINE_ROOT
export PROJECT_ROOT
export RESULTS_DIR
export CHECKPOINT_DIR
export TABLE_DIR
export FIGURE_DIR
export REPORT_DIR
export EDA_REPORT_DIR
export MANIFEST_DIR
export METADATA_DIR
export PROJECT_CONFIG_DIR
export LOG_DIR
export SAMPLE_SHEET="${METADATA_DIR}/samples.tsv"
export CANONICAL_SAMPLE_SHEET="${METADATA_DIR}/samples.canonical.tsv"
export COMPARISON_SHEET="${METADATA_DIR}/comparisons.tsv"
export OBJECT_LAYER_CONFIG_FILE="${PROJECT_CONFIG_DIR}/object_layers.tsv"
export LAYER_STATUS_FILE="${TABLE_DIR}/layer_status.tsv"
export SAMPLE_NAMES="S1,S2"
export DEG_IDENT_1="Ctrl"
export DEG_IDENT_2="Treat"
export MIN_BIOLOGICAL_REPLICATES="2"
export DEG_MIN_CELLS_PER_GROUP="3"
export DEG_LOGFC_THRESHOLD="0"
export RANDOM_SEED="42"

cat > "${SAMPLE_SHEET}" <<'EOF'
sample_id	condition	group_id	biological_replicate	technical_replicate	batch	input_mode	input_source	source_path	platform	gene_id_type	reference_version	timepoint	tissue	chemistry	run_main	run_velocity	run_scenic
S1	Ctrl	Ctrl	Ctrl_rep1	1	batch1	matrix	smoke	/tmp	10x_cellranger	symbol	GRCg7b	T1	tissue_x	v3	yes	no	no
S2	Treat	Treat	Treat_rep1	1	batch2	matrix	smoke	/tmp	10x_cellranger	symbol	GRCg7b	T1	tissue_x	v3	yes	no	no
EOF

cat > "${COMPARISON_SHEET}" <<'EOF'
comparison_id	ident_1	ident_2	enabled	group_var	batch_var	layer_scope	min_biological_replicates	subset_column	subset_value	force_exploratory	min_cells_per_group	logfc_threshold
Ctrl_vs_Treat	Ctrl	Treat	yes	group_id	batch	*	2			no	3	0
GC_forced	Ctrl	Treat	yes	group_id	batch	*	2	cell_type	GC	yes	3	0
EOF

"${R_BIN}" - <<'EOF'
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

project_root <- Sys.getenv("PROJECT_ROOT")
checkpoint_dir <- Sys.getenv("CHECKPOINT_DIR")
table_dir <- Sys.getenv("TABLE_DIR")

genes <- paste0("gene", seq_len(12))
cells <- paste0("cell", seq_len(16))
counts <- matrix(1L, nrow = length(genes), ncol = length(cells), dimnames = list(genes, cells))
counts[1:3, 1:8] <- counts[1:3, 1:8] + 5L
counts[4:6, 9:16] <- counts[4:6, 9:16] + 5L

meta <- data.frame(
  sample_id = rep(c("S1", "S2"), each = 8),
  orig.ident = rep(c("S1", "S2"), each = 8),
  group_id = rep(c("Ctrl", "Treat"), each = 8),
  condition = rep(c("Ctrl", "Treat"), each = 8),
  biological_replicate = rep(c("Ctrl_rep1", "Treat_rep1"), each = 8),
  batch = rep(c("batch1", "batch2"), each = 8),
  panorama_cluster = rep(c("C0", "C1", "C0", "C1"), each = 4),
  cluster_id = rep(c("C0", "C1", "C0", "C1"), each = 4),
  cell_type = rep(c("GC", "Other", "GC", "Other"), each = 4),
  annotation_label = rep(c("GC", "Other", "GC", "Other"), each = 4),
  row.names = cells,
  stringsAsFactors = FALSE
)

obj <- CreateSeuratObject(counts = Matrix(counts, sparse = TRUE), meta.data = meta, min.features = 0)
obj <- NormalizeData(obj, verbose = FALSE)

layer_dir <- file.path(checkpoint_dir, "layers", "panorama")
dir.create(layer_dir, recursive = TRUE, showWarnings = FALSE)
annotated_rds <- file.path(layer_dir, "panorama_after_annotation.rds")
saveRDS(obj, annotated_rds)
saveRDS(obj, file.path(checkpoint_dir, "03_after_annotation.rds"))

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  data.frame(
    layer_id = "panorama",
    status = "annotated",
    layer_role = "panorama",
    cluster_column = "panorama_cluster",
    clustered_rds = annotated_rds,
    annotated_rds = annotated_rds,
    stringsAsFactors = FALSE
  ),
  file = file.path(table_dir, "layer_status.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
EOF

"${R_BIN}" "${PIPELINE_ROOT}/workflow/05single_script/05a_marker_discovery.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/05single_script/05b_pseudobulk_de.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/05single_script/05c_composition.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/05single_script/05d_deg_eda.R"

"${R_BIN}" - <<'EOF'
table_dir <- Sys.getenv("TABLE_DIR")
report_dir <- file.path(Sys.getenv("REPORT_DIR"), "eda")

marker_manifest <- read.delim(file.path(table_dir, "marker_discovery", "marker_discovery_manifest.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
pb_manifest <- read.delim(file.path(table_dir, "pseudobulk_ds", "pseudobulk_manifest.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
comp_manifest <- read.delim(file.path(table_dir, "composition", "composition_manifest.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
status_matrix <- read.delim(file.path(table_dir, "deg", "deg_status_matrix.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
deg_summary <- read.delim(file.path(table_dir, "deg", "deg_summary.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(any(pb_manifest$comparison_id == "Ctrl_vs_Treat" & pb_manifest$inference_status == "exploratory_only"))
stopifnot(any(pb_manifest$comparison_id == "GC_forced" & pb_manifest$inference_status == "exploratory_forced"))
stopifnot(any(comp_manifest$comparison_id == "Ctrl_vs_Treat" & comp_manifest$inference_status == "exploratory_only"))
stopifnot(any(comp_manifest$comparison_id == "GC_forced" & comp_manifest$inference_status == "exploratory_forced"))
stopifnot(file.exists(file.path(report_dir, "deg", "report.md")))
stopifnot(any(status_matrix$comparison_id == "GC_forced" & status_matrix$pseudobulk_status == "exploratory_forced"))
stopifnot(all(c("exploratory_marker_total_n", "exploratory_marker_significant_n") %in% colnames(deg_summary)))
stopifnot(any(deg_summary$comparison_id == "GC_forced" & deg_summary$exploratory_marker_total_n > 0))

gc_summary_path <- marker_manifest$exploratory_summary_tsv[marker_manifest$comparison_id == "GC_forced"][[1]]
gc_summary <- read.delim(gc_summary_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(gc_summary$cluster_id %in% "C0"))

cat("smoke_deg_workflow_05: PASS\n")
EOF
