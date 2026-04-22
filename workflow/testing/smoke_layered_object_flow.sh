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
  echo "Missing Rscript for smoke_layered_object_flow: ${R_BIN}" >&2
  exit 1
}

TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/scrna_layered_smoke.XXXXXX")"
trap 'rm -rf "${TMP_ROOT}"' EXIT

PROJECT_ROOT="${TMP_ROOT}/project"
PROJECT_CONFIG_DIR="${PROJECT_ROOT}/config"
REPORT_DIR="${PROJECT_ROOT}/reports"
RESULTS_DIR="${PROJECT_ROOT}/results"
CHECKPOINT_DIR="${RESULTS_DIR}/checkpoints"
FIGURE_DIR="${RESULTS_DIR}/figures"
TABLE_DIR="${RESULTS_DIR}/tables"
LOG_DIR="${PROJECT_ROOT}/logs"
CONFIG_FILE="${TMP_ROOT}/config.sh"
OBJECT_LAYER_CONFIG_FILE="${PROJECT_CONFIG_DIR}/object_layers.tsv"

mkdir -p \
  "${PROJECT_CONFIG_DIR}" \
  "${PROJECT_ROOT}/metadata" \
  "${REPORT_DIR}/eda/pre_qc" \
  "${REPORT_DIR}/eda/post_qc" \
  "${REPORT_DIR}/eda/integration" \
  "${REPORT_DIR}/eda/annotation" \
  "${REPORT_DIR}/intake" \
  "${RESULTS_DIR}" \
  "${CHECKPOINT_DIR}" \
  "${FIGURE_DIR}" \
  "${TABLE_DIR}" \
  "${LOG_DIR}" \
  "${PROJECT_ROOT}/status"

cat > "${CONFIG_FILE}" <<EOF
#!/usr/bin/env bash
export PIPELINE_ROOT="${PIPELINE_ROOT}"
export PROJECT_ROOT="${PROJECT_ROOT}"
export DATA_DIR="${PROJECT_ROOT}/data"
export RESULTS_DIR="${RESULTS_DIR}"
export CHECKPOINT_DIR="${CHECKPOINT_DIR}"
export FIGURE_DIR="${FIGURE_DIR}"
export TABLE_DIR="${TABLE_DIR}"
export LOG_DIR="${LOG_DIR}"
export REPORT_DIR="${REPORT_DIR}"
export EDA_REPORT_DIR="${REPORT_DIR}/eda"
export PRE_QC_REPORT_DIR="${REPORT_DIR}/eda/pre_qc"
export POST_QC_REPORT_DIR="${REPORT_DIR}/eda/post_qc"
export INTEGRATION_REPORT_DIR="${REPORT_DIR}/eda/integration"
export ANNOTATION_REPORT_DIR="${REPORT_DIR}/eda/annotation"
export INTAKE_REPORT_DIR="${REPORT_DIR}/intake"
export PROJECT_CONFIG_DIR="${PROJECT_CONFIG_DIR}"
export OBJECT_LAYER_CONFIG_FILE="${OBJECT_LAYER_CONFIG_FILE}"
export SAMPLE_SHEET="${PROJECT_ROOT}/metadata/samples.tsv"
export CANONICAL_SAMPLE_SHEET="${PROJECT_ROOT}/metadata/samples.canonical.tsv"
export COMPARISON_SHEET="${PROJECT_ROOT}/metadata/comparisons.tsv"
export INPUT_INVENTORY_FILE="${REPORT_DIR}/intake/input_inventory.tsv"
export BRANCH_READINESS_FILE="${REPORT_DIR}/intake/branch_readiness.tsv"
export QC_THRESHOLD_FILE="${PROJECT_CONFIG_DIR}/qc_thresholds.tsv"
export EDA_GATE_FILE="${PROJECT_CONFIG_DIR}/eda_gates.tsv"
export SAMPLE_NAMES="S1,S2"
export ANALYSIS_GROUP_1_NAME="Ctrl"
export ANALYSIS_GROUP_1_SAMPLES="S1"
export ANALYSIS_GROUP_2_NAME="Treat"
export ANALYSIS_GROUP_2_SAMPLES="S2"
export RAW_GROUP_1_NAME="Ctrl"
export RAW_GROUP_1_SAMPLES="S1"
export RAW_GROUP_2_NAME="Treat"
export RAW_GROUP_2_SAMPLES="S2"
export RAW_SAMPLES="S1,S2"
export DEG_IDENT_1="Ctrl"
export DEG_IDENT_2="Treat"
export RANDOM_SEED="42"
export HVG_NFEATURES="8"
export PCA_DIMS="1:6"
export TARGET_CLUSTERS="2"
export RES_RANGE="0.10,0.20,0.30,0.40"
export RES_FINE_STEP="0.05"
export INTEGRATION_MODE="none"
export QC_MIN_NFEATURE="0"
export QC_MIN_NCOUNT="0"
export QC_MIN_LOG10UMI="0"
export QC_MAX_MITO_PCT="100"
EOF

{
  printf 'layer_id\tlayer_role\tenabled\tparent_layer\tsample_include\tsample_exclude\tselection_column\tselection_values\trebuild_normalization\thvg_nfeatures\tpca_dims\ttarget_clusters\tres_range\tres_fine_step\tintegration_mode\tdescription\n'
  printf 'panorama\tpanorama\tyes\t\t\t\t\t\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "8" \
    "1:6" \
    "2" \
    "0.10,0.20,0.30,0.40" \
    "0.05" \
    "none" \
    "Panorama test layer"
  printf 'subcluster_1\tsubcluster\tyes\tpanorama\t\t\tlineage_seed\tlineage_a\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "6" \
    "1:5" \
    "2" \
    "0.20,0.40,0.60,0.80,1.00" \
    "0.05" \
    "none" \
    "First test subcluster"
  printf 'subcluster_2\tsubcluster\tyes\tpanorama\t\t\tlineage_seed\tlineage_b\tyes\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "6" \
    "1:5" \
    "2" \
    "0.20,0.40,0.60,0.80,1.00" \
    "0.05" \
    "none" \
    "Second test subcluster"
} > "${OBJECT_LAYER_CONFIG_FILE}"

cat > "${PROJECT_ROOT}/metadata/samples.tsv" <<'EOF'
sample_id	condition	biological_replicate	technical_replicate	batch	input_mode	input_source	source_path	platform	gene_id_type	reference_version	timepoint	tissue	chemistry	run_main	run_velocity	run_scenic
S1	Ctrl	Ctrl_1	1	batch1	matrix	smoke	/tmp	10x_cellranger	symbol	GRCg7b	E12	ovary	v3	yes	no	no
S2	Treat	Treat_1	1	batch2	matrix	smoke	/tmp	10x_cellranger	symbol	GRCg7b	E12	ovary	v3	yes	no	no
EOF

cat > "${PROJECT_ROOT}/metadata/comparisons.tsv" <<'EOF'
comparison_id	ident_1	ident_2	enabled
Ctrl_vs_Treat	Ctrl	Treat	yes
EOF

cat > "${PROJECT_CONFIG_DIR}/qc_thresholds.tsv" <<'EOF'
sample_id	qc_min_nfeature	qc_min_ncount	qc_min_log10umi	qc_max_mito_pct
EOF

cat > "${PROJECT_CONFIG_DIR}/eda_gates.tsv" <<'EOF'
gate_id	status	approved_by	notes
pre_qc	approved		
post_qc	approved		
integration	approved		
annotation	approved		
EOF

set -a
# shellcheck disable=SC1090
source "${CONFIG_FILE}"
set +a

"${R_BIN}" - <<'EOF'
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

set.seed(42)

genes <- c(
  "LINEAGE_A_CORE1", "LINEAGE_A_CORE2", "LINEAGE_B_CORE1", "LINEAGE_B_CORE2",
  "STATE_A1", "STATE_A2", "STATE_B1", "STATE_B2",
  "HK1", "HK2", "NOISE1", "NOISE2"
)

make_group_counts <- function(sample_id, group_name, lineage_seed, sublineage_seed, n_cells) {
  counts <- matrix(rpois(length(genes) * n_cells, lambda = 1), nrow = length(genes))
  rownames(counts) <- genes
  colnames(counts) <- sprintf("%s_%s_%02d", sample_id, group_name, seq_len(n_cells))

  if (identical(lineage_seed, "lineage_a")) {
    counts[genes %in% c("LINEAGE_A_CORE1", "LINEAGE_A_CORE2"), ] <- counts[genes %in% c("LINEAGE_A_CORE1", "LINEAGE_A_CORE2"), ] + 18L
  } else {
    counts[genes %in% c("LINEAGE_B_CORE1", "LINEAGE_B_CORE2"), ] <- counts[genes %in% c("LINEAGE_B_CORE1", "LINEAGE_B_CORE2"), ] + 18L
  }

  if (identical(sublineage_seed, "state_a1")) {
    counts[genes == "STATE_A1", ] <- counts[genes == "STATE_A1", ] + 28L
  } else if (identical(sublineage_seed, "state_a2")) {
    counts[genes == "STATE_A2", ] <- counts[genes == "STATE_A2", ] + 28L
  } else if (identical(sublineage_seed, "state_b1")) {
    counts[genes == "STATE_B1", ] <- counts[genes == "STATE_B1", ] + 28L
  } else if (identical(sublineage_seed, "state_b2")) {
    counts[genes == "STATE_B2", ] <- counts[genes == "STATE_B2", ] + 28L
  }

  counts[genes %in% c("HK1", "HK2"), ] <- counts[genes %in% c("HK1", "HK2"), ] + 6L

  meta <- data.frame(
    row.names = colnames(counts),
    orig.ident = rep(sample_id, n_cells),
    analysis_group = rep(ifelse(sample_id == "S1", "Ctrl", "Treat"), n_cells),
    lineage_seed = rep(lineage_seed, n_cells),
    sublineage_seed = rep(sublineage_seed, n_cells),
    stringsAsFactors = FALSE
  )

  list(counts = counts, meta = meta)
}

group_defs <- list(
  c("S1", "a1", "lineage_a", "state_a1", 12L),
  c("S1", "a2", "lineage_a", "state_a2", 12L),
  c("S1", "b1", "lineage_b", "state_b1", 12L),
  c("S1", "b2", "lineage_b", "state_b2", 12L),
  c("S2", "a1", "lineage_a", "state_a1", 12L),
  c("S2", "a2", "lineage_a", "state_a2", 12L),
  c("S2", "b1", "lineage_b", "state_b1", 12L),
  c("S2", "b2", "lineage_b", "state_b2", 12L)
)

payloads <- lapply(group_defs, function(def) {
  make_group_counts(
    sample_id = def[[1]],
    group_name = def[[2]],
    lineage_seed = def[[3]],
    sublineage_seed = def[[4]],
    n_cells = as.integer(def[[5]])
  )
})

counts <- do.call(cbind, lapply(payloads, `[[`, "counts"))
meta <- do.call(rbind, lapply(payloads, `[[`, "meta"))

obj <- CreateSeuratObject(counts = Matrix(counts, sparse = TRUE), meta.data = meta, min.features = 0)
obj$orig.ident <- meta$orig.ident
obj$analysis_group <- meta$analysis_group
obj$lineage_seed <- meta$lineage_seed
obj$sublineage_seed <- meta$sublineage_seed

saveRDS(obj, file.path(Sys.getenv("CHECKPOINT_DIR"), "01_after_qc_doublet.rds"))
EOF

"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/02_build_reductions.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/02a_integration_eda.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/02b_finalize_clustering.R"

PANORAMA_CLUSTERED_RDS="${CHECKPOINT_DIR}/layers/panorama/panorama_after_clustering.rds"
PANORAMA_SIGNATURE_TSV="${TABLE_DIR}/layers/panorama/build_signature.tsv"
PANORAMA_MTIME_BEFORE="$(stat -c %Y "${PANORAMA_CLUSTERED_RDS}")"

sleep 1
"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/02b_finalize_clustering.R"

[[ -f "${PANORAMA_SIGNATURE_TSV}" ]] || {
  echo "Missing panorama cache signature: ${PANORAMA_SIGNATURE_TSV}" >&2
  exit 1
}

PANORAMA_MTIME_AFTER="$(stat -c %Y "${PANORAMA_CLUSTERED_RDS}")"
[[ "${PANORAMA_MTIME_BEFORE}" == "${PANORAMA_MTIME_AFTER}" ]] || {
  echo "Panorama checkpoint was unexpectedly rebuilt on second 02b run" >&2
  exit 1
}

"${R_BIN}" - <<'EOF'
suppressPackageStartupMessages({
  library(Seurat)
})

project_root <- Sys.getenv("PROJECT_ROOT")
checkpoint_dir <- file.path(project_root, "results", "checkpoints")
table_dir <- file.path(project_root, "results", "tables")
eda_dir <- file.path(project_root, "reports", "eda", "integration")

required_paths <- c(
  file.path(checkpoint_dir, "02_reduction_candidates.rds"),
  file.path(checkpoint_dir, "02_after_clustering.rds"),
  file.path(checkpoint_dir, "layers", "panorama", "panorama_reduction_candidates.rds"),
  file.path(checkpoint_dir, "layers", "panorama", "panorama_after_clustering.rds"),
  file.path(table_dir, "layers", "panorama", "build_signature.tsv"),
  file.path(checkpoint_dir, "layers", "subcluster_1", "subcluster_1_reduction_candidates.rds"),
  file.path(checkpoint_dir, "layers", "subcluster_1", "subcluster_1_after_clustering.rds"),
  file.path(checkpoint_dir, "layers", "subcluster_2", "subcluster_2_reduction_candidates.rds"),
  file.path(checkpoint_dir, "layers", "subcluster_2", "subcluster_2_after_clustering.rds"),
  file.path(table_dir, "layer_status.tsv"),
  file.path(eda_dir, "panorama", "report.md"),
  file.path(eda_dir, "subcluster_1", "report.md"),
  file.path(eda_dir, "subcluster_2", "report.md")
)

missing <- required_paths[!file.exists(required_paths)]
stopifnot(length(missing) == 0)

status_df <- read.delim(file.path(table_dir, "layer_status.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
status_map <- setNames(status_df$status, status_df$layer_id)
stopifnot(all(unname(status_map[c("panorama", "subcluster_1", "subcluster_2")]) == c("reused", "built", "built")))

panorama <- readRDS(file.path(checkpoint_dir, "layers", "panorama", "panorama_after_clustering.rds"))
subcluster_1_obj <- readRDS(file.path(checkpoint_dir, "layers", "subcluster_1", "subcluster_1_after_clustering.rds"))
subcluster_2_obj <- readRDS(file.path(checkpoint_dir, "layers", "subcluster_2", "subcluster_2_after_clustering.rds"))

for (obj in list(panorama, subcluster_1_obj, subcluster_2_obj)) {
  stopifnot("pca" %in% Reductions(obj))
  stopifnot("umap" %in% Reductions(obj))
  stopifnot(length(VariableFeatures(obj)) > 0)
  stopifnot(length(obj@graphs) > 0)
}

stopifnot(length(levels(panorama$seurat_clusters)) == 2)
stopifnot(length(levels(subcluster_1_obj$seurat_clusters)) == 2)
stopifnot(length(levels(subcluster_2_obj$seurat_clusters)) == 2)
stopifnot("panorama_cluster" %in% colnames(subcluster_1_obj@meta.data))
stopifnot("panorama_cluster" %in% colnames(subcluster_2_obj@meta.data))
stopifnot("subcluster_1_cluster" %in% colnames(subcluster_1_obj@meta.data))
stopifnot("subcluster_2_cluster" %in% colnames(subcluster_2_obj@meta.data))

cat("smoke_layered_object_flow: PASS\n")
EOF
