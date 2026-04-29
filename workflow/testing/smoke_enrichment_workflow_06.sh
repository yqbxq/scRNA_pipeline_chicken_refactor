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
  echo "Missing Rscript for smoke_enrichment_workflow_06: ${R_BIN}" >&2
  exit 1
}

if [[ -z "${R_LIBS_USER:-}" ]]; then
  if [[ -n "${R_LIBS_MAIN:-}" ]]; then
    export R_LIBS_USER="${R_LIBS_MAIN}"
  elif [[ -n "${R_MAIN_ENV_PREFIX:-}" ]]; then
    SHARED_ENV_ROOT="$(cd "${R_MAIN_ENV_PREFIX}/../.." && pwd)"
    if [[ -d "${SHARED_ENV_ROOT}/R_libs_main" ]]; then
      export R_LIBS_USER="${SHARED_ENV_ROOT}/R_libs_main"
    fi
  fi
fi

TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/scrna_enrichment06_smoke.XXXXXX")"
trap 'rm -rf "${TMP_ROOT}"' EXIT

PROJECT_ROOT="${TMP_ROOT}/project"
RESULTS_DIR="${PROJECT_ROOT}/results"
TABLE_DIR="${RESULTS_DIR}/tables"
FIGURE_DIR="${RESULTS_DIR}/figures"
REPORT_DIR="${PROJECT_ROOT}/reports"
EDA_REPORT_DIR="${REPORT_DIR}/eda"
MANIFEST_DIR="${RESULTS_DIR}/manifests"
METADATA_DIR="${PROJECT_ROOT}/metadata"
PROJECT_CONFIG_DIR="${PROJECT_ROOT}/config"

mkdir -p \
  "${TABLE_DIR}/enrichment/go/panorama/Ctrl_vs_Treat/C0" \
  "${TABLE_DIR}/enrichment/go/panorama/Ctrl_vs_Treat/C1" \
  "${TABLE_DIR}/enrichment/kegg/panorama/Ctrl_vs_Treat/C0" \
  "${TABLE_DIR}/deg" \
  "${FIGURE_DIR}" \
  "${EDA_REPORT_DIR}" \
  "${MANIFEST_DIR}" \
  "${METADATA_DIR}" \
  "${PROJECT_CONFIG_DIR}"

export PIPELINE_ROOT
export PROJECT_ROOT
export RESULTS_DIR
export TABLE_DIR
export FIGURE_DIR
export REPORT_DIR
export EDA_REPORT_DIR
export MANIFEST_DIR
export METADATA_DIR
export PROJECT_CONFIG_DIR
export ENRICHMENT_PVALUE_CUTOFF="0.05"
export ENRICHMENT_KEGG_TIMEOUT_SEC="5"

if "${R_BIN}" - <<'EOF'
required <- c("dplyr", "tibble", "jsonlite", "ggplot2", "AnnotationDbi", "org.Gg.eg.db")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
quit(status = if (length(missing) == 0) 0 else 1)
EOF
then
  "${R_BIN}" "${PIPELINE_ROOT}/workflow/05single_script/06a_go_enrichment.R"
  "${R_BIN}" "${PIPELINE_ROOT}/workflow/05single_script/06b_kegg_enrichment.R"

  "${R_BIN}" - <<'EOF'
project_root <- Sys.getenv("PROJECT_ROOT")
results_dir <- Sys.getenv("RESULTS_DIR")
go_manifest <- file.path(results_dir, "tables", "enrichment", "go", "go_enrichment_manifest.tsv")
kegg_manifest <- file.path(results_dir, "tables", "enrichment", "kegg", "kegg_enrichment_manifest.tsv")
go_report <- file.path(project_root, "reports", "eda", "enrichment", "06a_go_report.md")
kegg_report <- file.path(project_root, "reports", "eda", "enrichment", "06b_kegg_report.md")
stopifnot(file.exists(go_manifest), file.exists(kegg_manifest), file.exists(go_report), file.exists(kegg_report))
stopifnot(nrow(read.delim(go_manifest, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)) == 0)
stopifnot(nrow(read.delim(kegg_manifest, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)) == 0)
EOF

  SYN_DEG="${TABLE_DIR}/deg/smoke_deg.tsv"
  cat > "${SYN_DEG}" <<EOF
gene	cluster_id	avg_log2FC	p_val_adj
SMOKE_GENE_01	C0	2.0	0.001
SMOKE_GENE_02	C0	1.5	0.002
SMOKE_GENE_03	C0	1.0	0.003
SMOKE_GENE_04	C0	0.5	0.004
SMOKE_GENE_05	C0	0.25	0.005
SMOKE_GENE_06	C0	-0.4	0.006
SMOKE_GENE_07	C0	-0.8	0.007
SMOKE_GENE_08	C0	-1.2	0.008
EOF

  cat > "${TABLE_DIR}/deg/deg_status_matrix.tsv" <<EOF
layer_id	comparison_id	marker_status	pseudobulk_status	composition_status	marker_results_tsv	pseudobulk_results_tsv	composition_results_tsv
panorama	Ctrl_vs_Treat	ok	exploratory_only	exploratory_only	${SYN_DEG}		
EOF

  ENRICHMENT_MIN_INPUT_GENES=9999 "${R_BIN}" "${PIPELINE_ROOT}/workflow/05single_script/06a_go_enrichment.R"
  ENRICHMENT_MIN_INPUT_GENES=9999 "${R_BIN}" "${PIPELINE_ROOT}/workflow/05single_script/06b_kegg_enrichment.R"

  "${R_BIN}" - <<'EOF'
results_dir <- Sys.getenv("RESULTS_DIR")
go_manifest <- file.path(results_dir, "tables", "enrichment", "go", "go_enrichment_manifest.tsv")
kegg_manifest <- file.path(results_dir, "tables", "enrichment", "kegg", "kegg_enrichment_manifest.tsv")
go_df <- read.delim(go_manifest, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
kegg_df <- read.delim(kegg_manifest, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(go_df) == 9)
stopifnot(nrow(kegg_df) == 3)
stopifnot(all(go_df$status == "too_few_genes"))
stopifnot(all(kegg_df$status == "too_few_genes"))
EOF

  if "${R_BIN}" - <<'EOF'
required <- c("clusterProfiler", "AnnotationDbi", "org.Gg.eg.db")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
quit(status = if (length(missing) == 0) 0 else 1)
EOF
  then
    "${R_BIN}" - <<'EOF'
script_dir <- file.path(Sys.getenv("PIPELINE_ROOT"), "shell", "05single_script")
source_utf8 <- function(path) source(path, encoding = "UTF-8")
source_utf8(file.path(script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(script_dir, "helpers", "config.R"))
source_utf8(file.path(script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(script_dir, "helpers", "project_paths_05.R"))
source_utf8(file.path(script_dir, "helpers", "project_paths_06.R"))
source_utf8(file.path(script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(script_dir, "helpers", "comparison_utils.R"))
source_utf8(file.path(script_dir, "helpers", "deg_utils.R"))
source_utf8(file.path(script_dir, "helpers", "enrichment_utils.R"))
org_db <- get_org_db_06("org.Gg.eg.db")
entrez <- AnnotationDbi::keys(org_db, keytype = "ENTREZID")
entrez <- unique(entrez[!is.na(entrez) & nzchar(entrez)])
stopifnot(length(entrez) >= 30)
entrez <- entrez[seq_len(30)]
go_res <- run_go_enrichment_single(entrez, org_db, ont = "BP", pvalue_cutoff = 1, qvalue_cutoff = 1, min_gs_size = 1, max_gs_size = 5000)
stopifnot(go_res$status %in% c("ok", "no_terms"))
kegg_res <- run_kegg_enrichment_single(entrez, organism = "gga", pvalue_cutoff = 1, qvalue_cutoff = 1, min_gs_size = 1, max_gs_size = 5000, timeout_sec = 5)
stopifnot(kegg_res$status %in% c("ok", "no_terms", "timeout", "error"))
EOF
  else
    echo "smoke_enrichment_workflow_06: skipping direct clusterProfiler runner check; clusterProfiler is not available" >&2
  fi
else
  echo "smoke_enrichment_workflow_06: skipping 06a/06b script checks; enrichment R packages are not available" >&2
fi

GO_C0="${TABLE_DIR}/enrichment/go/panorama/Ctrl_vs_Treat/C0/go_BP_all_chicken.tsv"
GO_C1="${TABLE_DIR}/enrichment/go/panorama/Ctrl_vs_Treat/C1/go_BP_all_chicken.tsv"
KEGG_C0="${TABLE_DIR}/enrichment/kegg/panorama/Ctrl_vs_Treat/C0/kegg_all_chicken.tsv"

cat > "${GO_C0}" <<EOF
layer_id	comparison_id	cluster_id	gene_direction	analysis_type	ontology	source_species	enrichment_source	deg_source	deg_inference_status	ID	Description	GeneRatio	BgRatio	pvalue	p.adjust	qvalue	geneID	Count
panorama	Ctrl_vs_Treat	C0	all	go	BP	chicken	org.Gg.eg.db	exploratory_marker	exploratory_only	GO:0001	shared developmental process	4/20	40/1000	0.001	0.010	0.020	G1/G2/G3/G4	4
EOF

cat > "${GO_C1}" <<EOF
layer_id	comparison_id	cluster_id	gene_direction	analysis_type	ontology	source_species	enrichment_source	deg_source	deg_inference_status	ID	Description	GeneRatio	BgRatio	pvalue	p.adjust	qvalue	geneID	Count
panorama	Ctrl_vs_Treat	C1	all	go	BP	chicken	org.Gg.eg.db	exploratory_marker	exploratory_only	GO:0001	shared developmental process	3/18	40/1000	0.002	0.020	0.030	G2/G5/G6	3
EOF

cat > "${KEGG_C0}" <<EOF
layer_id	comparison_id	cluster_id	gene_direction	analysis_type	ontology	source_species	enrichment_source	deg_source	deg_inference_status	ID	Description	GeneRatio	BgRatio	pvalue	p.adjust	qvalue	geneID	Count
panorama	Ctrl_vs_Treat	C0	all	kegg	KEGG	chicken	org.Gg.eg.db_to_KEGG_gga	exploratory_marker	exploratory_only	gga00010	Glycolysis / Gluconeogenesis	2/20	30/1000	0.003	0.030	0.040	G7/G8	2
EOF

cat > "${TABLE_DIR}/enrichment/go/go_enrichment_manifest.tsv" <<EOF
layer_id	comparison_id	cluster_id	gene_direction	analysis_type	ontology	source_species	deg_source	deg_inference_status	status	reason	input_gene_n	mapped_gene_n	mapping_rate	significant_term_n	enrichment_tsv	enrichment_significant_tsv	dotplot_png	barplot_png
panorama	Ctrl_vs_Treat	C0	all	go	BP	chicken	exploratory_marker	exploratory_only	ok		20	18	0.90	1	${GO_C0}			
panorama	Ctrl_vs_Treat	C1	all	go	BP	chicken	exploratory_marker	exploratory_only	ok		18	16	0.89	1	${GO_C1}			
EOF

cat > "${TABLE_DIR}/enrichment/kegg/kegg_enrichment_manifest.tsv" <<EOF
layer_id	comparison_id	cluster_id	gene_direction	analysis_type	ontology	source_species	deg_source	deg_inference_status	status	reason	input_gene_n	mapped_gene_n	mapping_rate	significant_term_n	enrichment_tsv	enrichment_significant_tsv	dotplot_png	barplot_png
panorama	Ctrl_vs_Treat	C0	all	kegg	KEGG	chicken	exploratory_marker	exploratory_only	ok		20	18	0.90	1	${KEGG_C0}			
EOF

cat > "${TABLE_DIR}/deg/deg_status_matrix.tsv" <<EOF
layer_id	comparison_id	marker_status	pseudobulk_status	composition_status	marker_results_tsv	pseudobulk_results_tsv	composition_results_tsv
panorama	Ctrl_vs_Treat	ok	exploratory_only	exploratory_only	${GO_C0}		
EOF

"${R_BIN}" "${PIPELINE_ROOT}/workflow/05single_script/06c_enrichment_eda.R"

"${R_BIN}" - <<'EOF'
project_root <- Sys.getenv("PROJECT_ROOT")
results_dir <- Sys.getenv("RESULTS_DIR")
report_path <- file.path(project_root, "reports", "eda", "enrichment", "report.md")
summary_path <- file.path(results_dir, "tables", "enrichment", "enrichment_summary.tsv")
shared_path <- file.path(project_root, "reports", "eda", "enrichment", "shared_pathways.tsv")
manifest_path <- file.path(results_dir, "manifests", "06c_enrichment_eda", "_manifest.json")
go_heatmap <- file.path(results_dir, "figures", "enrichment", "cross_cluster_go_bp_heatmap.png")
kegg_heatmap <- file.path(results_dir, "figures", "enrichment", "cross_cluster_kegg_heatmap.png")

stopifnot(file.exists(report_path))
stopifnot(file.exists(summary_path))
stopifnot(file.exists(shared_path))
stopifnot(file.exists(manifest_path))
stopifnot(file.exists(go_heatmap))
stopifnot(file.exists(kegg_heatmap))

summary_df <- read.delim(summary_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
shared_df <- read.delim(shared_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(any(summary_df$cluster_id == "C0" & summary_df$total_significant_n >= 2))
stopifnot(any(shared_df$ID == "GO:0001" & shared_df$cluster_n == 2))

cat("smoke_enrichment_workflow_06: PASS\n")
EOF
