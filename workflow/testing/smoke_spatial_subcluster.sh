#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/spatial_subcluster_smoke.XXXXXX")"
trap 'rm -rf "${TMP_ROOT}"' EXIT

source "${SCRIPT_DIR}/lib_spatial_smoke_03.sh"
spatial_smoke03_require_r
spatial_smoke03_setup "${TMP_ROOT}" "${PIPELINE_ROOT}"
spatial_smoke03_run_to_normalized "${PIPELINE_ROOT}"

export SPATIAL_SUBCLUSTER_MIN_SPOTS="20"

Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/02a_spatial_integration_eda.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/02b_finalize_spatial_clustering.R" --backends b1_seurat_snn
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/03_region_annotation.R" "${MARKER_PANEL_DIR}"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/03a_region_annotation_eda.R"

Rscript - <<'RSCRIPT'
source(file.path(Sys.getenv("PIPELINE_ROOT"), "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(Sys.getenv("PIPELINE_ROOT"), "workflow/05single_script/spatial/helpers/project_paths_spatial.R"), encoding = "UTF-8")
cfg <- get_spatial_script_config()
panorama <- readRDS(cfg$spatial_panorama_annotated_rds)
spot_n <- suppressWarnings(as.integer(sub(".*_spot", "", colnames(panorama))))
spot_n[is.na(spot_n)] <- seq_along(spot_n[is.na(spot_n)])
region <- ifelse(spot_n <= 30, "GC_rich", "TC_rich")
panorama$region <- factor(region, levels = c("GC_rich", "TC_rich", "stroma", "vasculature", "mixed_or_uncertain"))
saveRDS(panorama, cfg$spatial_panorama_annotated_rds)
RSCRIPT

Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/04a_subcluster_build.R" "${SPATIAL_OBJECT_LAYER_FILE}"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/04b_subcluster_annotate.R" "${MARKER_PANEL_DIR}"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/04c_subcluster_eda.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/testing/smoke_spatial_subcluster.R" default_disabled

Rscript - <<'RSCRIPT'
path <- Sys.getenv("SPATIAL_OBJECT_LAYER_FILE")
layers <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
enable <- layers$layer_id %in% c("region_gc")
layers$enabled[enable] <- "yes"
layers$normalization_methods[enable] <- "log"
layers$hvg_nfeatures[enable] <- "80"
layers$pca_dims[enable] <- "1:10"
layers$target_clusters[enable] <- "2"
layers$res_range[enable] <- "0.2,0.4"
write.table(layers, path, sep = "\t", quote = FALSE, row.names = FALSE)
RSCRIPT

cat > "${MARKER_PANEL_DIR}/region_gc.tsv" <<'EOF'
layer_id	celltype	gene	evidence_source	tissue	panel_name	evidence_note	confidence_ceiling
region_gc	inner_ring	COL1A1	smoke	chicken_ovary	region_gc_smoke	Subregion smoke marker	confirmed
region_gc	inner_ring	DCN	smoke	chicken_ovary	region_gc_smoke	Subregion smoke marker	confirmed
region_gc	outer_ring	AMH	smoke	chicken_ovary	region_gc_smoke	Subregion smoke marker	confirmed
region_gc	outer_ring	FOXL2	smoke	chicken_ovary	region_gc_smoke	Subregion smoke marker	confirmed
EOF

Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/04a_subcluster_build.R" "${SPATIAL_OBJECT_LAYER_FILE}"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/04b_subcluster_annotate.R" "${MARKER_PANEL_DIR}"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/04c_subcluster_eda.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/testing/smoke_spatial_subcluster.R" one_enabled

rm -f "${MARKER_PANEL_DIR}/region_gc.tsv"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/04b_subcluster_annotate.R" "${MARKER_PANEL_DIR}"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/04c_subcluster_eda.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/testing/smoke_spatial_subcluster.R" one_enabled_no_panel

cat > "${MARKER_PANEL_DIR}/region_gc.tsv" <<'EOF'
layer_id	celltype	gene	evidence_source	tissue	panel_name	evidence_note	confidence_ceiling
region_gc	inner_ring	COL1A1	smoke	chicken_ovary	region_gc_smoke	Subregion smoke marker	confirmed
region_gc	inner_ring	DCN	smoke	chicken_ovary	region_gc_smoke	Subregion smoke marker	confirmed
region_gc	outer_ring	AMH	smoke	chicken_ovary	region_gc_smoke	Subregion smoke marker	confirmed
region_gc	outer_ring	FOXL2	smoke	chicken_ovary	region_gc_smoke	Subregion smoke marker	confirmed
EOF
cat > "${MARKER_PANEL_DIR}/region_tc.tsv" <<'EOF'
layer_id	celltype	gene	evidence_source	tissue	panel_name	evidence_note	confidence_ceiling
region_tc	inner_cuff	PECAM1	smoke	chicken_ovary	region_tc_smoke	Subregion smoke marker	confirmed
region_tc	outer_cuff	STAR	smoke	chicken_ovary	region_tc_smoke	Subregion smoke marker	confirmed
region_tc	outer_cuff	CYP17A1	smoke	chicken_ovary	region_tc_smoke	Subregion smoke marker	confirmed
EOF
Rscript - <<'RSCRIPT'
path <- Sys.getenv("SPATIAL_OBJECT_LAYER_FILE")
layers <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
enable <- layers$layer_id %in% c("region_gc", "region_tc")
layers$enabled[enable] <- "yes"
layers$normalization_methods[enable] <- "log"
layers$hvg_nfeatures[enable] <- "80"
layers$pca_dims[enable] <- "1:10"
layers$target_clusters[enable] <- "2"
layers$res_range[enable] <- "0.2,0.4"
write.table(layers, path, sep = "\t", quote = FALSE, row.names = FALSE)
RSCRIPT

Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/04a_subcluster_build.R" "${SPATIAL_OBJECT_LAYER_FILE}"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/04b_subcluster_annotate.R" "${MARKER_PANEL_DIR}"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/04c_subcluster_eda.R"
Rscript "${PIPELINE_ROOT}/workflow/05single_script/spatial/testing/smoke_spatial_subcluster.R" two_enabled

echo "smoke_spatial_subcluster_ok"
