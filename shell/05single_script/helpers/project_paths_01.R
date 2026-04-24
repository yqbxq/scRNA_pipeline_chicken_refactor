get_single_script_config_01 <- function() {
  project_root <- "/home/user_test/syf_f5/05projects/scrna_improve"
  reference_root <- "/home/user_test/syf_f5/01shared_resources/reference/chicken"
  reference_dir <- file.path(reference_root, "ensembl_release112")
  metadata_dir <- file.path(project_root, "metadata")
  config_dir <- file.path(project_root, "config")
  results_dir <- file.path(project_root, "results")
  report_dir <- file.path(project_root, "reports")
  intake_report_dir <- file.path(report_dir, "intake")
  eda_report_dir <- file.path(report_dir, "eda")
  manifest_dir <- file.path(results_dir, "manifests")
  ortholog_cache_dir <- file.path(project_root, "ortholog_cache")

  list(
    project_root = project_root,
    data_dir = file.path(project_root, "data"),
    metadata_dir = metadata_dir,
    config_dir = config_dir,
    results_dir = results_dir,
    checkpoint_dir = file.path(results_dir, "checkpoints"),
    table_dir = file.path(results_dir, "tables"),
    figure_dir = file.path(results_dir, "figures"),
    manifest_dir = manifest_dir,
    module_01a_manifest_dir = file.path(manifest_dir, "01a_build_raw_objects"),
    module_01b_manifest_dir = file.path(manifest_dir, "01b_pre_qc_eda"),
    module_01a_manifest_path = file.path(manifest_dir, "01a_build_raw_objects", "_manifest.json"),
    module_01b_manifest_path = file.path(manifest_dir, "01b_pre_qc_eda", "_manifest.json"),
    report_dir = report_dir,
    intake_report_dir = intake_report_dir,
    eda_report_dir = eda_report_dir,
    pre_qc_report_dir = file.path(eda_report_dir, "pre_qc"),
    sample_sheet = file.path(metadata_dir, "samples.tsv"),
    canonical_sample_sheet = file.path(metadata_dir, "samples.canonical.tsv"),
    input_inventory_file = file.path(intake_report_dir, "input_inventory.tsv"),
    branch_readiness_file = file.path(intake_report_dir, "branch_readiness.tsv"),
    qc_threshold_file = file.path(config_dir, "qc_thresholds.tsv"),
    eda_gate_file = file.path(config_dir, "eda_gates.tsv"),
    mito_gene_list_file = file.path(config_dir, "mito_gene_list.txt"),
    reference_root = reference_root,
    reference_dir = reference_dir,
    reference_version = "ensembl_release112",
    clean_gtf = file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112_clean.gtf"),
    reference_gtf = file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf"),
    ortholog_cache_dir = ortholog_cache_dir,
    ortholog_manifest_path = file.path(ortholog_cache_dir, "_manifest.json"),
    cellranger_out_dir = "",
    random_seed = 42L,
    qc_min_nfeature = 200,
    qc_min_ncount = 1000,
    qc_min_log10umi = 0.7,
    qc_max_mito_pct = 20,
    triage_frac_below_cutoff = 0.35,
    triage_frac_above_mito = 0.25,
    triage_density_peaks = 2L,
    module_version = "1.0"
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
