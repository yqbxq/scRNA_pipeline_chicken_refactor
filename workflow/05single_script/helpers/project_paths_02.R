env_or_default_02 <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

env_numeric_02 <- function(name, default) {
  value <- suppressWarnings(as.numeric(env_or_default_02(name, as.character(default))))
  if (is.na(value)) default else value
}

env_integer_02 <- function(name, default) {
  value <- suppressWarnings(as.integer(env_or_default_02(name, as.character(default))))
  if (is.na(value)) default else value
}

env_bool_02 <- function(name, default = TRUE) {
  value <- tolower(env_or_default_02(name, if (isTRUE(default)) "yes" else "no"))
  value %in% c("yes", "true", "1", "on")
}

split_csv_02 <- function(value) {
  if (length(value) == 0) {
    return(character(0))
  }
  value <- trimws(as.character(value[[1]]))
  if (is.na(value) || !nzchar(value)) {
    return(character(0))
  }
  parts <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  parts[nzchar(parts)]
}

parse_index_spec_02 <- function(value) {
  value <- trimws(as.character(value[[1]]))
  if (length(value) == 0 || is.na(value) || !nzchar(value)) {
    return(integer(0))
  }
  if (grepl("^[0-9]+:[0-9]+$", value)) {
    parts <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
    return(seq.int(parts[1], parts[2]))
  }
  as.integer(split_csv_02(value))
}

get_single_script_config_02 <- function() {
  project_root <- env_or_default_02("PROJECT_ROOT", normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  metadata_dir <- env_or_default_02("METADATA_DIR", file.path(project_root, "metadata"))
  config_dir <- env_or_default_02("PROJECT_CONFIG_DIR", file.path(project_root, "config"))
  results_dir <- env_or_default_02("RESULTS_DIR", file.path(project_root, "results"))
  report_dir <- env_or_default_02("REPORT_DIR", file.path(project_root, "reports"))
  intake_report_dir <- env_or_default_02("INTAKE_REPORT_DIR", file.path(report_dir, "intake"))
  eda_report_dir <- env_or_default_02("EDA_REPORT_DIR", file.path(report_dir, "eda"))
  manifest_dir <- env_or_default_02("MANIFEST_DIR", file.path(results_dir, "manifests"))
  reference_root <- env_or_default_02("REFERENCE_ROOT", file.path(project_root, "reference", "chicken"))
  reference_dir <- env_or_default_02("REFERENCE_DIR", file.path(reference_root, "ensembl_release112"))
  ortholog_cache_dir <- env_or_default_02("ORTHOLOG_CACHE_DIR", file.path(results_dir, "ortholog_cache"))

  list(
    project_root = project_root,
    data_dir = env_or_default_02("DATA_DIR", file.path(project_root, "data")),
    metadata_dir = metadata_dir,
    config_dir = config_dir,
    results_dir = results_dir,
    checkpoint_dir = env_or_default_02("CHECKPOINT_DIR", file.path(results_dir, "checkpoints")),
    table_dir = env_or_default_02("TABLE_DIR", file.path(results_dir, "tables")),
    figure_dir = env_or_default_02("FIGURE_DIR", file.path(results_dir, "figures")),
    manifest_dir = manifest_dir,
    module_01a_manifest_path = file.path(manifest_dir, "01a_build_raw_objects", "_manifest.json"),
    module_01b_manifest_path = file.path(manifest_dir, "01b_pre_qc_eda", "_manifest.json"),
    module_02a_manifest_dir = file.path(manifest_dir, "02a_ambient_branch"),
    module_02b1_manifest_dir = file.path(manifest_dir, "02b1_qc_filter"),
    module_02b2_manifest_dir = file.path(manifest_dir, "02b2_doublet"),
    module_02c_manifest_dir = file.path(manifest_dir, "02c_post_qc_eda"),
    module_02a_manifest_path = file.path(manifest_dir, "02a_ambient_branch", "_manifest.json"),
    module_02b1_manifest_path = file.path(manifest_dir, "02b1_qc_filter", "_manifest.json"),
    module_02b2_manifest_path = file.path(manifest_dir, "02b2_doublet", "_manifest.json"),
    module_02c_manifest_path = file.path(manifest_dir, "02c_post_qc_eda", "_manifest.json"),
    report_dir = report_dir,
    intake_report_dir = intake_report_dir,
    eda_report_dir = eda_report_dir,
    ambient_report_dir = env_or_default_02("AMBIENT_REPORT_DIR", file.path(eda_report_dir, "ambient")),
    pre_qc_report_dir = env_or_default_02("PRE_QC_REPORT_DIR", file.path(eda_report_dir, "pre_qc")),
    post_qc_report_dir = env_or_default_02("POST_QC_REPORT_DIR", file.path(eda_report_dir, "post_qc")),
    sample_sheet = env_or_default_02("SAMPLE_SHEET", file.path(metadata_dir, "samples.tsv")),
    canonical_sample_sheet = env_or_default_02("CANONICAL_SAMPLE_SHEET", file.path(metadata_dir, "samples.canonical.tsv")),
    input_inventory_file = env_or_default_02("INPUT_INVENTORY_FILE", file.path(intake_report_dir, "input_inventory.tsv")),
    branch_readiness_file = env_or_default_02("BRANCH_READINESS_FILE", file.path(intake_report_dir, "branch_readiness.tsv")),
    qc_threshold_file = env_or_default_02("QC_THRESHOLD_FILE", file.path(config_dir, "qc_thresholds.tsv")),
    eda_gate_file = env_or_default_02("EDA_GATE_FILE", file.path(config_dir, "eda_gates.tsv")),
    mito_gene_list_file = env_or_default_02("MITO_GENE_LIST_FILE", file.path(config_dir, "mito_gene_list.txt")),
    rbc_gene_list_file = env_or_default_02("RBC_GENE_LIST_FILE", file.path(config_dir, "rbc_gene_list.txt")),
    reference_root = reference_root,
    reference_dir = reference_dir,
    reference_version = env_or_default_02("REFERENCE_VERSION", "ensembl_release112"),
    clean_gtf = env_or_default_02("CLEAN_GTF", file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112_clean.gtf")),
    reference_gtf = env_or_default_02("REFERENCE_GTF", file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf")),
    ortholog_cache_dir = ortholog_cache_dir,
    ortholog_manifest_path = env_or_default_02("ORTHOLOG_MANIFEST", file.path(ortholog_cache_dir, "_manifest.json")),
    cellranger_out_dir = env_or_default_02("CELLRANGER_OUT_DIR", ""),
    random_seed = env_integer_02("RANDOM_SEED", 42L),
    qc_min_nfeature = env_numeric_02("QC_MIN_NFEATURE", 200),
    qc_min_ncount = env_numeric_02("QC_MIN_NCOUNT", 1000),
    qc_min_log10umi = env_numeric_02("QC_MIN_LOG10UMI", 0.7),
    qc_max_mito_pct = env_numeric_02("QC_MAX_MITO_PCT", 20),
    triage_frac_below_cutoff = env_numeric_02("TRIAGE_FRAC_BELOW_CUTOFF", 0.35),
    triage_frac_above_mito = env_numeric_02("TRIAGE_FRAC_ABOVE_MITO", 0.25),
    triage_density_peaks = env_integer_02("TRIAGE_DENSITY_PEAKS", 2L),
    ambient_primary_method = tolower(env_or_default_02("AMBIENT_PRIMARY_METHOD", "soupx")),
    ambient_fallback_method = tolower(env_or_default_02("AMBIENT_FALLBACK_METHOD", "decontx")),
    ambient_apply_policy = tolower(env_or_default_02("AMBIENT_APPLY_POLICY", "manual")),
    ambient_min_cells = env_integer_02("AMBIENT_MIN_CELLS", 50L),
    ambient_cluster_dims = parse_index_spec_02(env_or_default_02("AMBIENT_CLUSTER_DIMS", "1:20")),
    ambient_cluster_resolution = env_numeric_02("AMBIENT_CLUSTER_RESOLUTION", 0.4),
    ambient_marker_top_n = env_integer_02("AMBIENT_MARKER_TOP_N", 3L),
    ambient_recommend_min_contamination = env_numeric_02("AMBIENT_RECOMMEND_MIN_CONTAMINATION", 0.05),
    cellbender_mode = tolower(env_or_default_02("CELLBENDER_MODE", "stub")),
    cellbender_fpr = env_numeric_02("CELLBENDER_FPR", 0.01),
    cellbender_cuda = tolower(env_or_default_02("CELLBENDER_CUDA", "yes")),
    cellbender_extra_args = env_or_default_02("CELLBENDER_EXTRA_ARGS", ""),
    doublet_rate = env_numeric_02("DOUBLET_RATE", 0.008),
    doublet_rate_per_1k = env_numeric_02("DOUBLET_RATE_PER_1K", env_numeric_02("DOUBLET_RATE", 0.008)),
    doublet_primary_caller = tolower(env_or_default_02("DOUBLET_PRIMARY_CALLER", "scDblFinder")),
    doublet_secondary_caller = env_or_default_02("DOUBLET_SECONDARY_CALLER", "DoubletFinder"),
    doublet_secondary_enabled = env_bool_02("DOUBLET_SECONDARY_ENABLED", TRUE),
    doublet_min_cells = env_integer_02("DOUBLET_MIN_CELLS", 50L),
    doublet_dims = parse_index_spec_02(env_or_default_02("DOUBLET_DIMS", "1:20")),
    hvg_nfeatures = env_integer_02("HVG_NFEATURES", 2000L),
    pca_dims = parse_index_spec_02(env_or_default_02("PCA_DIMS", "1:30")),
    target_clusters = env_integer_02("TARGET_CLUSTERS", 15L),
    res_range = suppressWarnings(as.numeric(split_csv_02(env_or_default_02("RES_RANGE", "0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60")))),
    module_version = env_or_default_02("MODULE_02_VERSION", "1.0")
  )
}

prepare_dirs_02 <- function(cfg) {
  dirs <- c(
    cfg$results_dir,
    cfg$checkpoint_dir,
    cfg$table_dir,
    cfg$figure_dir,
    cfg$manifest_dir,
    cfg$module_02a_manifest_dir,
    cfg$module_02b1_manifest_dir,
    cfg$module_02b2_manifest_dir,
    cfg$module_02c_manifest_dir,
    cfg$report_dir,
    cfg$intake_report_dir,
    cfg$eda_report_dir,
    cfg$ambient_report_dir,
    cfg$pre_qc_report_dir,
    cfg$post_qc_report_dir,
    cfg$config_dir
  )
  invisible(lapply(dirs[nzchar(dirs)], ensure_dir))
}
