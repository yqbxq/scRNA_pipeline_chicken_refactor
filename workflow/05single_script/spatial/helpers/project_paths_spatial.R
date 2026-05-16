spatial_env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

spatial_env_numeric <- function(name, default) {
  value <- suppressWarnings(as.numeric(spatial_env_or_default(name, as.character(default))))
  if (is.na(value)) default else value
}

spatial_env_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(spatial_env_or_default(name, as.character(default))))
  if (is.na(value)) default else value
}

get_spatial_script_config <- function() {
  project_root <- spatial_env_or_default("PROJECT_ROOT", normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  pipeline_root <- spatial_env_or_default("PIPELINE_ROOT", normalizePath(file.path(getwd()), winslash = "/", mustWork = FALSE))
  metadata_dir <- spatial_env_or_default("METADATA_DIR", file.path(project_root, "metadata"))
  config_dir <- spatial_env_or_default("PROJECT_CONFIG_DIR", file.path(project_root, "config"))
  results_dir <- spatial_env_or_default("RESULTS_DIR", file.path(project_root, "results"))
  report_dir <- spatial_env_or_default("REPORT_DIR", file.path(project_root, "reports"))
  eda_report_dir <- spatial_env_or_default("EDA_REPORT_DIR", file.path(report_dir, "eda"))
  intake_report_dir <- spatial_env_or_default("INTAKE_REPORT_DIR", file.path(report_dir, "intake"))
  manifest_dir <- spatial_env_or_default("MANIFEST_DIR", file.path(results_dir, "manifests"))
  spatial_results_dir <- spatial_env_or_default("SPATIAL_RESULTS_DIR", file.path(results_dir, "spatial"))
  spatial_table_dir <- spatial_env_or_default("SPATIAL_TABLE_DIR", file.path(spatial_results_dir, "tables"))
  spatial_figure_dir <- spatial_env_or_default("SPATIAL_FIGURE_DIR", file.path(spatial_results_dir, "figures"))
  spatial_checkpoint_dir <- spatial_env_or_default("SPATIAL_CHECKPOINT_DIR", file.path(spatial_results_dir, "checkpoints"))
  spatial_pre_qc_report_dir <- spatial_env_or_default("SPATIAL_PRE_QC_REPORT_DIR", file.path(eda_report_dir, "spatial_pre_qc"))
  spatial_post_qc_report_dir <- spatial_env_or_default("SPATIAL_POST_QC_REPORT_DIR", file.path(eda_report_dir, "spatial_post_qc"))
  normalization_compare_dir <- spatial_env_or_default("SPATIAL_NORMALIZATION_COMPARE_DIR", file.path(report_dir, "spatial", "normalization_compare"))
  py_spatial_prefix <- spatial_env_or_default("PY_SPATIAL_ENV_PREFIX", "")
  py_spatial_bin_default <- if (nzchar(py_spatial_prefix)) file.path(py_spatial_prefix, "bin", "python") else "python"

  list(
    project_root = project_root,
    pipeline_root = pipeline_root,
    data_dir = spatial_env_or_default("DATA_DIR", file.path(project_root, "data")),
    metadata_dir = metadata_dir,
    config_dir = config_dir,
    results_dir = results_dir,
    report_dir = report_dir,
    eda_report_dir = eda_report_dir,
    intake_report_dir = intake_report_dir,
    manifest_dir = manifest_dir,
    spatial_results_dir = spatial_results_dir,
    spatial_checkpoint_dir = spatial_checkpoint_dir,
    spatial_table_dir = spatial_table_dir,
    spatial_figure_dir = spatial_figure_dir,
    spatial_pre_qc_report_dir = spatial_pre_qc_report_dir,
    spatial_post_qc_report_dir = spatial_post_qc_report_dir,
    spatial_pre_qc_table_dir = file.path(spatial_table_dir, "spatial_pre_qc"),
    spatial_pre_qc_figure_dir = file.path(spatial_pre_qc_report_dir, "figures"),
    spatial_post_qc_table_dir = file.path(spatial_table_dir, "spatial_post_qc"),
    spatial_post_qc_figure_dir = file.path(spatial_post_qc_report_dir, "figures"),
    spatial_post_qc_mask_dir = file.path(spatial_post_qc_report_dir, "figures", "qc_filter_mask"),
    spatial_post_qc_overlay_dir = file.path(spatial_post_qc_report_dir, "figures", "raw_vs_post_qc"),
    normalization_variants_dir = spatial_env_or_default("SPATIAL_NORMALIZATION_VARIANTS_DIR", file.path(spatial_checkpoint_dir, "normalization_variants")),
    normalization_compare_dir = normalization_compare_dir,
    sample_sheet = spatial_env_or_default("SAMPLE_SHEET", file.path(metadata_dir, "samples.tsv")),
    canonical_sample_sheet = spatial_env_or_default("CANONICAL_SAMPLE_SHEET", file.path(metadata_dir, "samples.canonical.tsv")),
    section_sheet = spatial_env_or_default("SECTION_SHEET", file.path(metadata_dir, "sections.tsv")),
    spatial_reference_inventory_file = spatial_env_or_default("SPATIAL_REFERENCE_INVENTORY_FILE", file.path(metadata_dir, "spatial_reference_inventory.tsv")),
    spatial_input_inventory_file = spatial_env_or_default("SPATIAL_INPUT_INVENTORY_FILE", file.path(intake_report_dir, "spatial_input_inventory.tsv")),
    spatial_intake_contract_file = spatial_env_or_default("SPATIAL_INTAKE_CONTRACT_FILE", file.path(config_dir, "spatial_intake_contract.tsv")),
    spatial_qc_threshold_file = spatial_env_or_default("SPATIAL_QC_THRESHOLD_FILE", file.path(config_dir, "spatial_qc_thresholds.tsv")),
    spatial_object_layer_file = spatial_env_or_default("SPATIAL_OBJECT_LAYER_FILE", file.path(config_dir, "spatial_object_layers.tsv")),
    eda_gate_file = spatial_env_or_default("EDA_GATE_FILE", file.path(config_dir, "eda_gates.tsv")),
    clean_gtf = spatial_env_or_default("CLEAN_GTF", file.path(project_root, "reference", "genes.clean.gtf")),
    module_01_manifest_path = file.path(manifest_dir, "spatial_01_build_objects", "_manifest.json"),
    module_01a_manifest_path = file.path(manifest_dir, "spatial_01a_pre_qc_eda", "_manifest.json"),
    module_02_filter_manifest_path = file.path(manifest_dir, "spatial_01b_qc_filter", "_manifest.json"),
    module_02_eda_manifest_path = file.path(manifest_dir, "spatial_01c_post_qc_eda", "_manifest.json"),
    module_02_norm_manifest_path = file.path(manifest_dir, "spatial_02_normalize", "_manifest.json"),
    load_summary_csv = file.path(spatial_table_dir, "load_summary.csv"),
    load_diagnostics_tsv = file.path(spatial_table_dir, "load_diagnostics.tsv"),
    pre_qc_report_md = file.path(spatial_pre_qc_report_dir, "report.md"),
    pre_qc_triage_tsv = file.path(spatial_table_dir, "spatial_pre_qc", "triage_signals.tsv"),
    qc_filter_summary_tsv = file.path(spatial_table_dir, "spatial_post_qc", "qc_filter_summary.tsv"),
    post_qc_report_md = file.path(spatial_post_qc_report_dir, "report.md"),
    post_filter_triage_tsv = file.path(spatial_table_dir, "spatial_post_qc", "post_filter_triage.tsv"),
    selected_method_tsv = file.path(normalization_compare_dir, "selected_method.tsv"),
    selected_method_override_file = spatial_env_or_default("SPATIAL_NORMALIZATION_OVERRIDE_FILE", file.path(config_dir, "spatial_normalization_override.tsv")),
    pearson_residuals_py = spatial_env_or_default("SPATIAL_PEARSON_RESIDUALS_PY", file.path(pipeline_root, "workflow", "04python", "normalize_pearson_residuals.py")),
    py_spatial_bin = spatial_env_or_default("PY_SPATIAL_BIN", py_spatial_bin_default),
    random_seed = spatial_env_integer("RANDOM_SEED", 42L),
    qc_min_nfeature_default = spatial_env_numeric("SPATIAL_QC_MIN_NFEATURE", 200),
    qc_min_ncount_default = spatial_env_numeric("SPATIAL_QC_MIN_NCOUNT", 500),
    qc_max_mito_pct_default = spatial_env_numeric("SPATIAL_QC_MAX_MITO_PCT", 20),
    qc_max_nfeature_default = spatial_env_numeric("SPATIAL_QC_MAX_NFEATURE", Inf),
    qc_max_ncount_default = spatial_env_numeric("SPATIAL_QC_MAX_NCOUNT", Inf),
    excessive_drop_threshold = spatial_env_numeric("SPATIAL_EXCESSIVE_DROP_THRESHOLD", 0.5),
    spatial_dbscan_eps_factor = spatial_env_numeric("SPATIAL_DBSCAN_EPS_FACTOR", 1.5),
    default_normalization_method = spatial_env_or_default("SPATIAL_DEFAULT_NORMALIZATION_METHOD", "m3_sct_v2"),
    default_hvg_n = spatial_env_integer("SPATIAL_HVG_NFEATURES", spatial_env_integer("HVG_NFEATURES", 2000L)),
    module_version = spatial_env_or_default("MODULE_SPATIAL_01_VERSION", "1.0"),
    module_02_version = spatial_env_or_default("MODULE_SPATIAL_02_VERSION", "1.0")
  )
}

prepare_dirs_spatial <- function(cfg) {
  dirs <- c(
    cfg$spatial_results_dir,
    cfg$spatial_checkpoint_dir,
    cfg$spatial_table_dir,
    cfg$spatial_figure_dir,
    cfg$spatial_pre_qc_report_dir,
    cfg$spatial_post_qc_report_dir,
    cfg$spatial_pre_qc_table_dir,
    cfg$spatial_pre_qc_figure_dir,
    cfg$spatial_post_qc_table_dir,
    cfg$spatial_post_qc_figure_dir,
    cfg$spatial_post_qc_mask_dir,
    cfg$spatial_post_qc_overlay_dir,
    cfg$normalization_variants_dir,
    cfg$normalization_compare_dir,
    dirname(cfg$module_01_manifest_path),
    dirname(cfg$module_01a_manifest_path),
    dirname(cfg$module_02_filter_manifest_path),
    dirname(cfg$module_02_eda_manifest_path),
    dirname(cfg$module_02_norm_manifest_path),
    dirname(cfg$load_summary_csv),
    dirname(cfg$load_diagnostics_tsv),
    dirname(cfg$pre_qc_report_md),
    dirname(cfg$pre_qc_triage_tsv),
    dirname(cfg$qc_filter_summary_tsv),
    dirname(cfg$post_qc_report_md),
    dirname(cfg$post_filter_triage_tsv),
    dirname(cfg$selected_method_tsv)
  )
  invisible(lapply(dirs[nzchar(dirs)], ensure_dir))
}

spatial_safe_id <- function(value) {
  value <- trimws(as.character(value))
  value <- gsub("[^A-Za-z0-9._-]+", "_", value, perl = TRUE)
  value <- gsub("^_+|_+$", "", value, perl = TRUE)
  if (!nzchar(value)) "value" else value
}
