get_single_script_config <- function() {
  reference_root <- "/home/user_test/syf_f5/01shared_resources/reference/chicken"
  reference_dir <- file.path(reference_root, "ensembl_release112")
  output_dir <- "/home/user_test/syf_f5/05projects/scrna_improve/ortholog_cache"

  list(
    reference_root = reference_root,
    reference_dir = reference_dir,
    clean_gtf = file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112_clean.gtf"),
    reference_gtf = file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf"),
    ensembl_gtf_gz = file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf.gz"),
    ensembl_fasta_gz = file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz"),
    ensembl_fasta = file.path(reference_dir, "genome.fa"),
    cellranger_ref_ensembl = file.path(reference_dir, "cellranger_ref_ensembl_r112"),
    output_dir = output_dir,
    figure_dir = file.path(output_dir, "figures"),
    manifest_path = file.path(output_dir, "_manifest.json"),
    target_species = c("human", "mouse"),
    mirrors = c("www", "useast", "asia"),
    chunk_size = 500L,
    key_genes = c("DRGX", "CREB3L2", "EMX2", "CEBPB", "FOSL2", "JUN", "PBX3", "HMGA1", "SREBF2"),
    module_contract = "00_ortholog_module",
    module_version = "1.0"
  )
}
