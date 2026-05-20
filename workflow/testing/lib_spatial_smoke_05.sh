#!/usr/bin/env bash

spatial_smoke05_prepare_region_fixture() {
  local tmp_root="$1"
  local pipeline_root="$2"
  source "${pipeline_root}/workflow/testing/lib_spatial_smoke_03.sh"
  spatial_smoke03_require_r
  spatial_smoke03_setup "${tmp_root}" "${pipeline_root}"
  spatial_smoke03_run_to_normalized "${pipeline_root}"
  export COMPARISON_SHEET="${METADATA_DIR}/comparisons.tsv"

  cat > "${COMPARISON_SHEET}" <<'EOF'
comparison_id	source_question_id	display_question_id	output_alias	report_title	layer_scope	contrast_axis	analysis_mode	analysis_modality	analysis_unit	stat_level	group_var	ident_1	ident_2	subset_column	subset_value	aggregation_group_var	composition_group_var	batch_var	enabled	min_biological_replicates	force_exploratory	min_cells_per_group	logfc_threshold	produces_gene_program	gene_program_role	notes
syf_vs_f5_region	syf_vs_f5_region	syf_vs_f5_region	syf_vs_f5_region	Spatial smoke region comparison	panorama_st	condition_split	condition_pairwise	spatial	region	section_level	condition	syf	f5			section_id	condition	batch	yes	3	no	3	0.1	no		smoke
EOF

  Rscript "${pipeline_root}/workflow/05single_script/spatial/02a_spatial_integration_eda.R"
  Rscript "${pipeline_root}/workflow/05single_script/spatial/02b_finalize_spatial_clustering.R" --backends b1_seurat_snn
  Rscript "${pipeline_root}/workflow/05single_script/spatial/03_region_annotation.R" "${MARKER_PANEL_DIR}"
  Rscript "${pipeline_root}/workflow/05single_script/spatial/03a_region_annotation_eda.R"

  Rscript - <<'RSCRIPT'
source(file.path(Sys.getenv("PIPELINE_ROOT"), "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(Sys.getenv("PIPELINE_ROOT"), "workflow/05single_script/spatial/helpers/project_paths_spatial.R"), encoding = "UTF-8")
cfg <- get_spatial_script_config()
panorama <- readRDS(cfg$spatial_panorama_annotated_rds)
meta <- panorama@meta.data
sec <- as.character(meta$section_id)
panorama$condition <- ifelse(sec == "sec_1", "syf", "f5")
spot_n <- suppressWarnings(as.integer(sub(".*_spot", "", colnames(panorama))))
spot_n[is.na(spot_n)] <- seq_along(spot_n[is.na(spot_n)])
panorama$region <- factor(ifelse(spot_n <= 30, "GC_rich", "TC_rich"), levels = c("GC_rich", "TC_rich", "mixed_or_uncertain"))
saveRDS(panorama, cfg$spatial_panorama_annotated_rds)
RSCRIPT
}

spatial_smoke05_prepare_passgate_fixture() {
  local tmp_root="$1"
  local pipeline_root="$2"
  spatial_smoke05_prepare_region_fixture "${tmp_root}" "${pipeline_root}"
  Rscript - <<'RSCRIPT'
source(file.path(Sys.getenv("PIPELINE_ROOT"), "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(Sys.getenv("PIPELINE_ROOT"), "workflow/05single_script/spatial/helpers/project_paths_spatial.R"), encoding = "UTF-8")
cfg <- get_spatial_script_config()
panorama <- readRDS(cfg$spatial_panorama_annotated_rds)
meta <- panorama@meta.data
splits <- ave(seq_len(ncol(panorama)), as.character(meta$condition), FUN = seq_along)
rep_id <- ((splits - 1) %% 3) + 1
panorama$sample_id <- paste0(as.character(meta$condition), "_rep", rep_id)
panorama$batch <- ifelse(rep_id == 2, "B2", "B1")
saveRDS(panorama, cfg$spatial_panorama_annotated_rds)
RSCRIPT
  sed -i 's/	3	no	3	0.1/	2	no	3	0.1/' "${COMPARISON_SHEET}"
}

spatial_smoke05_assert_manifest_status() {
  local manifest_tsv="$1"
  local expected="$2"
  local strict="${3:-no}"
  Rscript - "${manifest_tsv}" "${expected}" "${strict}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
expected <- args[[2]]
strict <- identical(args[[3]], "yes")
if (!file.exists(path)) stop(sprintf("missing manifest tsv: %s", path), call. = FALSE)
df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(df) == 0) stop("manifest has zero rows", call. = FALSE)
if (!expected %in% df$status) stop(sprintf("expected status %s missing; got: %s", expected, paste(unique(df$status), collapse = ",")), call. = FALSE)
if (strict && !all(df$status == expected)) stop(sprintf("strict mode: not all rows are %s; got: %s", expected, paste(unique(df$status), collapse = ",")), call. = FALSE)
bad <- grep("^failed_", df$status, value = TRUE)
if (length(bad) > 0) stop(sprintf("smoke saw failure rows: %s", paste(unique(bad), collapse = ",")), call. = FALSE)
RSCRIPT
}
