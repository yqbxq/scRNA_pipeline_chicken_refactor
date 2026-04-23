get_single_script_config <- function() {
  reference_dir <- "/home/user_test/syf_f5/01shared_resources/reference/chicken"
  output_dir <- "/home/user_test/syf_f5/05projects/scrna_improve/ortholog_cache"

  list(
    reference_dir = reference_dir,
    clean_gtf = file.path(reference_dir, "GRCg7b_genomic_clean.gtf"),
    reference_gtf = file.path(reference_dir, "GRCg7b_genomic.gtf"),
    output_dir = output_dir,
    figure_dir = file.path(output_dir, "figures"),
    manifest_path = file.path(output_dir, "_manifest.json"),
    target_species = c("human", "mouse"),
    mirrors = c("asia", "www", "useast"),
    chunk_size = 500L,
    key_genes = c("DRGX", "CREB3L2", "EMX2", "CEBPB", "FOSL2", "JUN", "PBX3", "HMGA1", "SREBF2"),
    module_contract = "00_ortholog_module",
    module_version = "1.0"
  )
}
