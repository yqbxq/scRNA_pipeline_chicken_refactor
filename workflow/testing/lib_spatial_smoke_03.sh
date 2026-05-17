#!/usr/bin/env bash

spatial_smoke03_require_r() {
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "skip: Rscript is not available"
    exit 0
  fi
  if ! Rscript -e 'quit(status = if (requireNamespace("Seurat", quietly = TRUE) && requireNamespace("Matrix", quietly = TRUE)) 0 else 42)' >/dev/null 2>&1; then
    echo "skip: Seurat and Matrix are required for spatial module 03 smoke"
    exit 0
  fi
}

spatial_smoke03_setup() {
  local tmp_root="$1"
  local pipeline_root="$2"
  mkdir -p "${tmp_root}/refs" "${tmp_root}/matrix" "${tmp_root}/project/config/marker_panels" \
    "${tmp_root}/project/metadata" "${tmp_root}/project/reports/intake"

  cat > "${tmp_root}/refs/genes.gtf" <<'EOF'
chrM	source	gene	1	4	.	+	.	gene_id "MT-CO1"; gene_name "MT-CO1";
EOF

  Rscript - <<RSCRIPT
if (!requireNamespace("Matrix", quietly = TRUE)) quit(status = 42)
set.seed(31)
root <- "${tmp_root}/matrix"
features <- paste0("gene", seq_len(80))
features[1:8] <- c("MT-CO1", "AMH", "FOXL2", "STAR", "CYP17A1", "COL1A1", "DCN", "PECAM1")
for (sec in c("sec_1", "sec_2")) {
  dir.create(file.path(root, sec), recursive = TRUE, showWarnings = FALSE)
  counts <- matrix(rpois(80 * 60, lambda = 2), nrow = 80)
  counts[1, ] <- rpois(60, lambda = 1)
  counts[2:3, 1:30] <- counts[2:3, 1:30] + 8
  counts[4:5, 31:60] <- counts[4:5, 31:60] + 8
  counts[6:7, 10:20] <- counts[6:7, 10:20] + 5
  counts[8, 45:55] <- counts[8, 45:55] + 5
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  Matrix::writeMM(counts, file.path(root, sec, "matrix.mtx"))
  write.table(
    data.frame(gene_id = features, gene_name = features, type = "Gene Expression"),
    file.path(root, sec, "features.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  barcodes <- paste0(sec, "_spot", seq_len(60))
  writeLines(barcodes, file.path(root, sec, "barcodes.tsv"))
  write.csv(
    data.frame(
      barcode = barcodes,
      row = c(rep(seq_len(6), each = 5), rep(seq_len(6), each = 5)),
      col = c(rep(seq_len(5), times = 6), rep(seq_len(5) + 10, times = 6)),
      in_tissue = TRUE
    ),
    file.path(root, sec, "coords.csv"),
    row.names = FALSE
  )
}
RSCRIPT

  cat > "${tmp_root}/project/config/project_config.sh" <<EOF
export PROJECT_ROOT="${tmp_root}/project"
export PIPELINE_ROOT="${pipeline_root}"
export DATA_DIR="${tmp_root}/project/data"
export RESULTS_DIR="${tmp_root}/project/results"
export REPORT_DIR="${tmp_root}/project/reports"
export EDA_REPORT_DIR="${tmp_root}/project/reports/eda"
export INTAKE_REPORT_DIR="${tmp_root}/project/reports/intake"
export MANIFEST_DIR="${tmp_root}/project/results/manifests"
export METADATA_DIR="${tmp_root}/project/metadata"
export PROJECT_CONFIG_DIR="${tmp_root}/project/config"
export SAMPLE_SHEET="${tmp_root}/project/metadata/samples.tsv"
export CANONICAL_SAMPLE_SHEET="${tmp_root}/project/metadata/samples.canonical.tsv"
export SECTION_SHEET="${tmp_root}/project/metadata/sections.tsv"
export SPATIAL_INPUT_INVENTORY_FILE="${tmp_root}/project/reports/intake/spatial_input_inventory.tsv"
export SPATIAL_INTAKE_CONTRACT_FILE="${tmp_root}/project/config/spatial_intake_contract.tsv"
export SPATIAL_QC_THRESHOLD_FILE="${tmp_root}/project/config/spatial_qc_thresholds.tsv"
export SPATIAL_OBJECT_LAYER_FILE="${tmp_root}/project/config/spatial_object_layers.tsv"
export EDA_GATE_FILE="${tmp_root}/project/config/eda_gates.tsv"
export MARKER_PANEL_DIR="${tmp_root}/project/config/marker_panels"
export SPATIAL_RESULTS_DIR="${tmp_root}/project/results/spatial"
export SPATIAL_CHECKPOINT_DIR="${tmp_root}/project/results/spatial/checkpoints"
export SPATIAL_TABLE_DIR="${tmp_root}/project/results/spatial/tables"
export SPATIAL_FIGURE_DIR="${tmp_root}/project/results/spatial/figures"
export SPATIAL_PRE_QC_REPORT_DIR="${tmp_root}/project/reports/eda/spatial_pre_qc"
export SPATIAL_POST_QC_REPORT_DIR="${tmp_root}/project/reports/eda/spatial_post_qc"
export SPATIAL_NORMALIZATION_VARIANTS_DIR="${tmp_root}/project/results/spatial/checkpoints/normalization_variants"
export SPATIAL_NORMALIZATION_COMPARE_DIR="${tmp_root}/project/reports/spatial/normalization_compare"
export SPATIAL_NORMALIZATION_OVERRIDE_FILE="${tmp_root}/project/config/spatial_normalization_override.tsv"
export SPATIAL_INTEGRATION_COMPARE_DIR="${tmp_root}/project/reports/eda/spatial_integration"
export SPATIAL_CLUSTERING_VARIANTS_DIR="${tmp_root}/project/results/spatial/checkpoints/clustering_variants"
export SPATIAL_CLUSTERING_COMPARE_DIR="${tmp_root}/project/reports/spatial/clustering_compare"
export SPATIAL_CLUSTERING_OVERRIDE_FILE="${tmp_root}/project/config/spatial_clustering_override.tsv"
export SPATIAL_REGION_ANNOTATION_DIR="${tmp_root}/project/reports/eda/spatial_region_annotation"
export SPATIAL_REGION_ANNOTATION_TABLE_DIR="${tmp_root}/project/results/spatial/tables/spatial_region_annotation"
export SPATIAL_REGION_ANNOTATION_OVERRIDE_FILE="${tmp_root}/project/config/spatial_region_annotation_override.tsv"
export SPATIAL_DEFAULT_NORMALIZATION_METHOD="m1_lognormalize"
export SPATIAL_DEFAULT_INTEGRATION_MODE="none"
export SPATIAL_INTEGRATION_RUN_MODES="none,harmony,cca"
export SPATIAL_CLUSTERING_BACKENDS="b1_seurat_snn,b2_bayesspace,b3_spagcn,b4_stagate"
export SPATIAL_DEFAULT_CLUSTERING_BACKEND="b1_seurat_snn"
export SPATIAL_CLUSTER_TARGET="4"
export PY_SPATIAL_BIN="${tmp_root}/missing_py_spatial/bin/python"
export SPAGCN_PY_BIN="${tmp_root}/missing_py_spatial_legacy/bin/python"
export STAGATE_PY_BIN="${tmp_root}/missing_py_spatial/bin/python"
export CLEAN_GTF="${tmp_root}/refs/genes.gtf"
EOF

  {
    printf 'sample_id\tcondition\tbiological_replicate\ttechnical_replicate\tbatch\tinput_mode\tinput_source\tsource_path\tplatform\tgene_id_type\treference_version\tgroup_id\ttimepoint\ttissue\tchemistry\trun_main\trun_velocity\trun_scenic\tmodality\tsection_id\tchip_id\tbundle_layout\timage_path\trun_spatial\trun_deconv\trun_joint\tnotes\n'
    for sec in sec_1 sec_2; do
      printf 'st_%s\tsmoke\t%s\t1\tst_batch\tspatial_matrix\tsmoke\t%s/matrix/%s\tgeneric\tauto\tcustom_reference\tsmoke\tD0\tchicken_ovary\tauto\tno\tno\tno\tspatial\t%s\tchip_%s\tgeneric_spatial_matrix\t\tyes\tauto\tauto\tsmoke\n' "${sec}" "${sec}" "${tmp_root}" "${sec}" "${sec}" "${sec}"
    done
  } > "${tmp_root}/project/metadata/samples.canonical.tsv"
  cp "${tmp_root}/project/metadata/samples.canonical.tsv" "${tmp_root}/project/metadata/samples.tsv"

  {
    printf 'section_id\tchip_id\ttissue_block\tcondition\tplatform\tbundle_root\timage_lowres\timage_hires\ttissue_positions\tscalefactors\tenabled\tnotes\n'
    for sec in sec_1 sec_2; do
      printf '%s\tchip_%s\tblock\tsmoke\tgeneric\t%s/matrix/%s\t\t\t\t\tyes\tsmoke\n' "${sec}" "${sec}" "${tmp_root}" "${sec}"
    done
  } > "${tmp_root}/project/metadata/sections.tsv"

  {
    printf 'sample_id\tsection_id\tmodality\tcondition\tplatform\tbundle_layout\tsource_path\tsection_bundle_root\tstandardized_outs\tready\tready_status\trequired_relpaths\tmissing_required_relpaths\toptional_present_relpaths\toptional_missing_relpaths\tnotes\n'
    for sec in sec_1 sec_2; do
      printf 'st_%s\t%s\tspatial\tsmoke\tgeneric\tgeneric_spatial_matrix\t%s/matrix/%s\t%s/matrix/%s\t%s/matrix/%s\ttrue\ttrue\t\t\t\t\tsmoke\n' "${sec}" "${sec}" "${tmp_root}" "${sec}" "${tmp_root}" "${sec}" "${tmp_root}" "${sec}"
    done
  } > "${tmp_root}/project/reports/intake/spatial_input_inventory.tsv"

  cat > "${tmp_root}/project/config/spatial_intake_contract.tsv" <<'EOF'
platform	bundle_layout	required_relpaths	optional_relpaths	standardized_outs_relpath	ready_status_when_present	notes
generic	generic_spatial_matrix	matrix.mtx,features.tsv,barcodes.tsv	coords.csv	.	true	Smoke generic spatial matrix.
EOF

  cat > "${tmp_root}/project/config/spatial_qc_thresholds.tsv" <<'EOF'
section_id	qc_min_nfeature	qc_max_nfeature	qc_min_ncount	qc_max_ncount	qc_max_mito_pct	mito_set_override	spatial_aware_filter	excessive_drop_threshold
__DEFAULT__	5		20		95		true	0.5
sec_1	5		20		95		true	0.5
sec_2	5		20		95		true	0.5
EOF
  cp "${pipeline_root}/config/spatial_object_layers.tsv.template" "${tmp_root}/project/config/spatial_object_layers.tsv"
  cp "${pipeline_root}/config/spatial_normalization_override.tsv.template" "${tmp_root}/project/config/spatial_normalization_override.tsv"
  cp "${pipeline_root}/config/spatial_clustering_override.tsv.template" "${tmp_root}/project/config/spatial_clustering_override.tsv"
  cp "${pipeline_root}/config/spatial_region_annotation_override.tsv.template" "${tmp_root}/project/config/spatial_region_annotation_override.tsv"
  cp "${pipeline_root}/config/eda_gates.tsv.template" "${tmp_root}/project/config/eda_gates.tsv"
  cp "${pipeline_root}/config/marker_panels/spatial_region_panel.tsv.template" "${tmp_root}/project/config/marker_panels/spatial_region_panel.tsv"

  export SCRNA_PIPELINE_CONFIG="${tmp_root}/project/config/project_config.sh"
  set -a
  # shellcheck disable=SC1090
  source "${SCRNA_PIPELINE_CONFIG}"
  set +a
}

spatial_smoke03_run_to_normalized() {
  local pipeline_root="$1"
  Rscript "${pipeline_root}/workflow/05single_script/spatial/01_build_spatial_objects.R"
  Rscript "${pipeline_root}/workflow/05single_script/spatial/01a_pre_spot_qc_eda.R"
  Rscript "${pipeline_root}/workflow/05single_script/spatial/01b_spot_qc_filter.R"
  Rscript "${pipeline_root}/workflow/05single_script/spatial/01c_post_spot_qc_eda.R"
  Rscript "${pipeline_root}/workflow/05single_script/spatial/02_normalize_spatial.R" --methods m0,m1
}
