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
  echo "Missing Rscript for smoke_annotation_stats_split: ${R_BIN}" >&2
  exit 1
}

TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/scrna_p03_smoke.XXXXXX")"
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
MARKER_PANEL_DIR="${PROJECT_CONFIG_DIR}/marker_panels"

mkdir -p \
  "${PROJECT_CONFIG_DIR}" \
  "${MARKER_PANEL_DIR}" \
  "${PROJECT_ROOT}/metadata" \
  "${REPORT_DIR}/eda/pre_qc" \
  "${REPORT_DIR}/eda/post_qc" \
  "${REPORT_DIR}/eda/integration" \
  "${REPORT_DIR}/eda/annotation" \
  "${REPORT_DIR}/eda/marker_discovery" \
  "${REPORT_DIR}/eda/pseudobulk_ds" \
  "${REPORT_DIR}/eda/composition" \
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
export MARKER_PANEL_DIR="${MARKER_PANEL_DIR}"
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
export HVG_NFEATURES="10"
export PCA_DIMS="1:6"
export TARGET_CLUSTERS="2"
export RES_RANGE="0.10,0.20,0.30,0.40"
export RES_FINE_STEP="0.05"
export INTEGRATION_MODE="none"
export QC_MIN_NFEATURE="0"
export QC_MIN_NCOUNT="0"
export QC_MIN_LOG10UMI="0"
export QC_MAX_MITO_PCT="100"
export MIN_BIOLOGICAL_REPLICATES="2"
EOF

{
  printf 'layer_id\tlayer_role\tenabled\tparent_layer\tsample_include\tsample_exclude\tselection_column\tselection_values\trebuild_normalization\thvg_nfeatures\tpca_dims\ttarget_clusters\tres_range\tres_fine_step\tintegration_mode\tdescription\n'
  printf 'panorama\tpanorama\tyes\t\t\t\t\t\tyes\t10\t1:6\t2\t0.10,0.20,0.30,0.40\t0.05\tnone\tNeutral panorama layer\n'
  printf 'subcluster_1\tsubcluster\tyes\tpanorama\t\t\tlineage_seed\tlineage_a\tyes\t8\t1:5\t2\t0.20,0.40,0.60,0.80\t0.05\tnone\tNeutral subcluster A\n'
  printf 'subcluster_2\tsubcluster\tyes\tpanorama\t\t\tlineage_seed\tlineage_b\tyes\t8\t1:5\t2\t0.20,0.40,0.60,0.80\t0.05\tnone\tNeutral subcluster B\n'
} > "${OBJECT_LAYER_CONFIG_FILE}"

cat > "${PROJECT_ROOT}/metadata/samples.tsv" <<'EOF'
sample_id	condition	group_id	biological_replicate	technical_replicate	batch	input_mode	input_source	source_path	platform	gene_id_type	reference_version	timepoint	tissue	chemistry	run_main	run_velocity	run_scenic
S1	Ctrl	Ctrl	Ctrl_rep1	1	batch1	matrix	smoke	/tmp	10x_cellranger	symbol	GRCg7b	T1	tissue_x	v3	yes	no	no
S2	Treat	Treat	Treat_rep1	1	batch2	matrix	smoke	/tmp	10x_cellranger	symbol	GRCg7b	T1	tissue_x	v3	yes	no	no
EOF

cat > "${PROJECT_ROOT}/metadata/comparisons.tsv" <<'EOF'
comparison_id	ident_1	ident_2	enabled	group_var	batch_var	layer_scope	min_biological_replicates
Ctrl_vs_Treat	Ctrl	Treat	yes	group_id	batch	*	2
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
  "LINEAGE-A-CORE1", "LINEAGE-A-CORE2",
  "LINEAGE-B-CORE1", "LINEAGE-B-CORE2",
  "A1-MARK1", "A1-MARK2",
  "A2-MARK1", "A2-MARK2",
  "B1-MARK1", "B1-MARK2",
  "B2-MARK1", "B2-MARK2",
  "HK1", "HK2", "NOISE1", "NOISE2"
)

make_group_counts <- function(sample_id, group_id, lineage_seed, sublineage_seed, n_cells) {
  counts <- matrix(rpois(length(genes) * n_cells, lambda = 1), nrow = length(genes))
  rownames(counts) <- genes
  colnames(counts) <- sprintf("%s_%s_%02d", sample_id, sublineage_seed, seq_len(n_cells))

  if (identical(lineage_seed, "lineage_a")) {
    counts[genes %in% c("LINEAGE-A-CORE1", "LINEAGE-A-CORE2"), ] <- counts[genes %in% c("LINEAGE-A-CORE1", "LINEAGE-A-CORE2"), ] + 20L
  } else {
    counts[genes %in% c("LINEAGE-B-CORE1", "LINEAGE-B-CORE2"), ] <- counts[genes %in% c("LINEAGE-B-CORE1", "LINEAGE-B-CORE2"), ] + 20L
  }

  if (identical(sublineage_seed, "a1")) {
    counts[genes %in% c("A1-MARK1", "A1-MARK2"), ] <- counts[genes %in% c("A1-MARK1", "A1-MARK2"), ] + 30L
  } else if (identical(sublineage_seed, "a2")) {
    counts[genes %in% c("A2-MARK1", "A2-MARK2"), ] <- counts[genes %in% c("A2-MARK1", "A2-MARK2"), ] + 30L
  } else if (identical(sublineage_seed, "b1")) {
    counts[genes %in% c("B1-MARK1", "B1-MARK2"), ] <- counts[genes %in% c("B1-MARK1", "B1-MARK2"), ] + 30L
  } else if (identical(sublineage_seed, "b2")) {
    counts[genes %in% c("B2-MARK1", "B2-MARK2"), ] <- counts[genes %in% c("B2-MARK1", "B2-MARK2"), ] + 30L
  }

  counts[genes %in% c("HK1", "HK2"), ] <- counts[genes %in% c("HK1", "HK2"), ] + 6L

  meta <- data.frame(
    row.names = colnames(counts),
    orig.ident = rep(sample_id, n_cells),
    sample_id = rep(sample_id, n_cells),
    group_id = rep(group_id, n_cells),
    condition = rep(group_id, n_cells),
    analysis_group = rep(group_id, n_cells),
    biological_replicate = rep(paste0(group_id, "_rep1"), n_cells),
    batch = rep(ifelse(sample_id == "S1", "batch1", "batch2"), n_cells),
    tissue = rep("tissue_x", n_cells),
    lineage_seed = rep(lineage_seed, n_cells),
    sublineage_seed = rep(sublineage_seed, n_cells),
    stringsAsFactors = FALSE
  )

  list(counts = counts, meta = meta)
}

payloads <- list(
  make_group_counts("S1", "Ctrl", "lineage_a", "a1", 12L),
  make_group_counts("S1", "Ctrl", "lineage_a", "a2", 12L),
  make_group_counts("S1", "Ctrl", "lineage_b", "b1", 12L),
  make_group_counts("S1", "Ctrl", "lineage_b", "b2", 12L),
  make_group_counts("S2", "Treat", "lineage_a", "a1", 12L),
  make_group_counts("S2", "Treat", "lineage_a", "a2", 12L),
  make_group_counts("S2", "Treat", "lineage_b", "b1", 12L),
  make_group_counts("S2", "Treat", "lineage_b", "b2", 12L)
)

counts <- do.call(cbind, lapply(payloads, `[[`, "counts"))
meta <- do.call(rbind, lapply(payloads, `[[`, "meta"))

obj <- CreateSeuratObject(counts = Matrix(counts, sparse = TRUE), meta.data = meta, min.features = 0)
saveRDS(obj, file.path(Sys.getenv("CHECKPOINT_DIR"), "01_after_qc_doublet.rds"))
EOF

"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/02_build_reductions.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/02a_integration_eda.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/02b_finalize_clustering.R"

# First run without any active marker panel TSV.
"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/03_annotation.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/03a_annotation_eda.R"

"${R_BIN}" - <<'EOF'
summary_df <- read.delim(file.path(Sys.getenv("TABLE_DIR"), "annotation", "layer_annotation_summary.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(summary_df) == 3)
stopifnot(all(summary_df$undetermined_n >= 1))

for (layer_id in c("panorama", "subcluster_1", "subcluster_2")) {
  ann_path <- file.path(Sys.getenv("TABLE_DIR"), "annotation", "layers", layer_id, "annotation_table.tsv")
  stopifnot(file.exists(ann_path))
  ann_df <- read.delim(ann_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("data_evidence", "literature_evidence", "module_score_check", "confidence") %in% colnames(ann_df)))
  stopifnot(all(ann_df$confidence == "未定"))
}
EOF

cat > "${MARKER_PANEL_DIR}/neutral_panel.tsv" <<'EOF'
layer_id	celltype	gene	evidence_source	tissue	panel_name	evidence_note	confidence_ceiling
panorama	Type_A	LINEAGE-A-CORE1	literature	tissue_x	neutral_panel	smoke direct evidence	确定
panorama	Type_A	LINEAGE-A-CORE2	literature	tissue_x	neutral_panel	smoke direct evidence	确定
panorama	Type_B	LINEAGE-B-CORE1	literature	tissue_x	neutral_panel	smoke direct evidence	确定
panorama	Type_B	LINEAGE-B-CORE2	literature	tissue_x	neutral_panel	smoke direct evidence	确定
subcluster_1	Subtype_A1	A1-MARK1	literature	tissue_x	neutral_panel	smoke direct evidence	确定
subcluster_1	Subtype_A1	A1-MARK2	literature	tissue_x	neutral_panel	smoke direct evidence	确定
subcluster_1	Subtype_A2	A2-MARK1	literature	tissue_x	neutral_panel	smoke direct evidence	确定
subcluster_1	Subtype_A2	A2-MARK2	literature	tissue_x	neutral_panel	smoke direct evidence	确定
subcluster_2	State_B1_candidate	B1-MARK1	ortholog	tissue_x	neutral_panel	smoke indirect evidence	暂定
subcluster_2	State_B2_candidate	B2-MARK1	ortholog	tissue_x	neutral_panel	smoke indirect evidence	暂定
EOF

"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/03_annotation.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/03a_annotation_eda.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/04_marker_discovery.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/04a_pseudobulk_ds.R"
"${R_BIN}" "${PIPELINE_ROOT}/workflow/r/04b_composition.R"

"${R_BIN}" - <<'EOF'
table_dir <- Sys.getenv("TABLE_DIR")
report_dir <- file.path(Sys.getenv("REPORT_DIR"), "eda")

summary_df <- read.delim(file.path(table_dir, "annotation", "layer_annotation_summary.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(any(summary_df$determined_n > 0))
stopifnot(any(summary_df$tentative_n > 0))

agg_meta_path <- file.path(table_dir, "pseudobulk_ds", "layers", "panorama", "pseudobulk_column_metadata.tsv")
stopifnot(file.exists(agg_meta_path))
agg_meta <- read.delim(agg_meta_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("cluster_id", "sample_id", "group_id", "biological_replicate") %in% colnames(agg_meta)))
stopifnot(nrow(agg_meta) == length(unique(paste(agg_meta$cluster_id, agg_meta$sample_id, sep = "__"))))

pb_status <- read.delim(file.path(table_dir, "pseudobulk_ds", "comparisons", "Ctrl_vs_Treat", "panorama", "ds_status.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(pb_status$inference_status == "exploratory_only"))
stopifnot(!file.exists(file.path(table_dir, "pseudobulk_ds", "comparisons", "Ctrl_vs_Treat", "panorama", "formal_pseudobulk_ds.tsv")))

comp_props <- file.path(table_dir, "composition", "layers", "panorama", "sample_level_proportions.tsv")
stopifnot(file.exists(comp_props))
comp_status <- read.delim(file.path(table_dir, "composition", "comparisons", "Ctrl_vs_Treat", "panorama", "composition_status.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(comp_status$inference_status == "exploratory_only"))
stopifnot(!file.exists(file.path(table_dir, "composition", "comparisons", "Ctrl_vs_Treat", "panorama", "formal_propeller.tsv")))

pb_report <- readLines(file.path(report_dir, "pseudobulk_ds", "report.md"), warn = FALSE, encoding = "UTF-8")
stopifnot(any(grepl("每组 N=1", pb_report, fixed = TRUE)))
comp_report <- readLines(file.path(report_dir, "composition", "report.md"), warn = FALSE, encoding = "UTF-8")
stopifnot(any(grepl("每组 N=1", comp_report, fixed = TRUE)))

md_manifest <- file.path(table_dir, "marker_discovery", "marker_discovery_manifest.tsv")
stopifnot(file.exists(md_manifest))

cat("smoke_annotation_stats_split: PASS\n")
EOF
