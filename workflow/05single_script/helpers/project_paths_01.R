env_or_default_01 <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

env_numeric_01 <- function(name, default) {
  value <- suppressWarnings(as.numeric(env_or_default_01(name, as.character(default))))
  if (is.na(value)) default else value
}

env_integer_01 <- function(name, default) {
  value <- suppressWarnings(as.integer(env_or_default_01(name, as.character(default))))
  if (is.na(value)) default else value
}

get_single_script_config_01 <- function() {
  project_root <- env_or_default_01("PROJECT_ROOT", normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  metadata_dir <- env_or_default_01("METADATA_DIR", file.path(project_root, "metadata"))
  config_dir <- env_or_default_01("PROJECT_CONFIG_DIR", file.path(project_root, "config"))
  results_dir <- env_or_default_01("RESULTS_DIR", file.path(project_root, "results"))
  report_dir <- env_or_default_01("REPORT_DIR", file.path(project_root, "reports"))
  intake_report_dir <- env_or_default_01("INTAKE_REPORT_DIR", file.path(report_dir, "intake"))
  eda_report_dir <- env_or_default_01("EDA_REPORT_DIR", file.path(report_dir, "eda"))
  manifest_dir <- env_or_default_01("MANIFEST_DIR", file.path(results_dir, "manifests"))
  reference_root <- env_or_default_01("REFERENCE_ROOT", file.path(project_root, "reference", "chicken"))
  reference_dir <- env_or_default_01("REFERENCE_DIR", file.path(reference_root, "ensembl_release112"))
  ortholog_cache_dir <- env_or_default_01("ORTHOLOG_CACHE_DIR", file.path(results_dir, "ortholog_cache"))

  list(
    project_root = project_root,
    data_dir = env_or_default_01("DATA_DIR", file.path(project_root, "data")),
    metadata_dir = metadata_dir,
    config_dir = config_dir,
    results_dir = results_dir,
    checkpoint_dir = env_or_default_01("CHECKPOINT_DIR", file.path(results_dir, "checkpoints")),
    table_dir = env_or_default_01("TABLE_DIR", file.path(results_dir, "tables")),
    figure_dir = env_or_default_01("FIGURE_DIR", file.path(results_dir, "figures")),
    manifest_dir = manifest_dir,
    module_01a_manifest_dir = file.path(manifest_dir, "01a_build_raw_objects"),
    module_01b_manifest_dir = file.path(manifest_dir, "01b_pre_qc_eda"),
    module_01a_manifest_path = file.path(manifest_dir, "01a_build_raw_objects", "_manifest.json"),
    module_01b_manifest_path = file.path(manifest_dir, "01b_pre_qc_eda", "_manifest.json"),
    report_dir = report_dir,
    intake_report_dir = intake_report_dir,
    eda_report_dir = eda_report_dir,
    pre_qc_report_dir = env_or_default_01("PRE_QC_REPORT_DIR", file.path(eda_report_dir, "pre_qc")),
    sample_sheet = env_or_default_01("SAMPLE_SHEET", file.path(metadata_dir, "samples.tsv")),
    canonical_sample_sheet = env_or_default_01("CANONICAL_SAMPLE_SHEET", file.path(metadata_dir, "samples.canonical.tsv")),
    input_inventory_file = env_or_default_01("INPUT_INVENTORY_FILE", file.path(intake_report_dir, "input_inventory.tsv")),
    branch_readiness_file = env_or_default_01("BRANCH_READINESS_FILE", file.path(intake_report_dir, "branch_readiness.tsv")),
    qc_threshold_file = env_or_default_01("QC_THRESHOLD_FILE", file.path(config_dir, "qc_thresholds.tsv")),
    eda_gate_file = env_or_default_01("EDA_GATE_FILE", file.path(config_dir, "eda_gates.tsv")),
    mito_gene_list_file = env_or_default_01("MITO_GENE_LIST_FILE", file.path(config_dir, "mito_gene_list.txt")),
    rbc_gene_list_file = env_or_default_01("RBC_GENE_LIST_FILE", file.path(config_dir, "rbc_gene_list.txt")),
    reference_root = reference_root,
    reference_dir = reference_dir,
    reference_version = env_or_default_01("REFERENCE_VERSION", "ensembl_release112"),
    clean_gtf = env_or_default_01("CLEAN_GTF", file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112_clean.gtf")),
    reference_gtf = env_or_default_01("REFERENCE_GTF", file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf")),
    ortholog_cache_dir = ortholog_cache_dir,
    ortholog_manifest_path = env_or_default_01("ORTHOLOG_MANIFEST", file.path(ortholog_cache_dir, "_manifest.json")),
    cellranger_out_dir = env_or_default_01("CELLRANGER_OUT_DIR", ""),
    random_seed = env_integer_01("RANDOM_SEED", 42L),
    qc_min_nfeature = env_numeric_01("QC_MIN_NFEATURE", 200),
    qc_min_ncount = env_numeric_01("QC_MIN_NCOUNT", 1000),
    qc_min_log10umi = env_numeric_01("QC_MIN_LOG10UMI", 0.7),
    qc_max_mito_pct = env_numeric_01("QC_MAX_MITO_PCT", 20),
    triage_frac_below_cutoff = env_numeric_01("TRIAGE_FRAC_BELOW_CUTOFF", 0.35),
    triage_frac_above_mito = env_numeric_01("TRIAGE_FRAC_ABOVE_MITO", 0.25),
    triage_density_peaks = env_integer_01("TRIAGE_DENSITY_PEAKS", 2L),
    module_version = env_or_default_01("MODULE_01_VERSION", "1.0")
  )
}

prepare_dirs_01 <- function(cfg) {
  dirs <- c(
    cfg$results_dir,
    cfg$checkpoint_dir,
    cfg$table_dir,
    cfg$figure_dir,
    cfg$manifest_dir,
    cfg$module_01a_manifest_dir,
    cfg$module_01b_manifest_dir,
    cfg$report_dir,
    cfg$intake_report_dir,
    cfg$eda_report_dir,
    cfg$pre_qc_report_dir,
    cfg$config_dir
  )
  invisible(lapply(dirs[nzchar(dirs)], ensure_dir))
}
