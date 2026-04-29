get_single_script_config <- function() {
  env_or_default <- function(name, default) {
    value <- Sys.getenv(name, unset = "")
    if (nzchar(value)) value else default
  }

  split_env_csv <- function(name, default) {
    value <- env_or_default(name, default)
    parts <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
    parts[nzchar(parts)]
  }

  project_root <- env_or_default("PROJECT_ROOT", normalizePath(file.path(getwd()), mustWork = FALSE))
  results_dir <- env_or_default("RESULTS_DIR", file.path(project_root, "results"))
  reference_root <- env_or_default("REFERENCE_ROOT", file.path(project_root, "reference", "chicken"))
  reference_dir <- env_or_default("REFERENCE_DIR", file.path(reference_root, "ensembl_release112"))
  output_dir <- env_or_default("ORTHOLOG_CACHE_DIR", file.path(results_dir, "ortholog_cache"))

  list(
    reference_root = reference_root,
    reference_dir = reference_dir,
    clean_gtf = env_or_default("CLEAN_GTF", file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112_clean.gtf")),
    reference_gtf = env_or_default("REFERENCE_GTF", file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf")),
    ensembl_gtf_gz = env_or_default("ENSEMBL_GTF_GZ", file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf.gz")),
    ensembl_fasta_gz = env_or_default("GENOME_FASTA_GZ", file.path(reference_dir, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz")),
    ensembl_fasta = env_or_default("GENOME_FASTA", file.path(reference_dir, "genome.fa")),
    star_index_dir = env_or_default("STAR_INDEX_DIR", file.path(reference_dir, "star_index")),
    output_dir = output_dir,
    figure_dir = file.path(output_dir, "figures"),
    manifest_path = file.path(output_dir, "_manifest.json"),
    target_species = split_env_csv("ORTHOLOG_TARGET_SPECIES", "human,mouse"),
    mirrors = split_env_csv("ENSEMBL_MIRRORS", "www,useast,asia"),
    chunk_size = as.integer(env_or_default("ORTHOLOG_QUERY_CHUNK_SIZE", "500")),
    key_genes = split_env_csv("ORTHOLOG_KEY_GENES", "DRGX,CREB3L2,EMX2,CEBPB,FOSL2,JUN,PBX3,HMGA1,SREBF2"),
    module_contract = "00_ortholog_module",
    module_version = "1.0"
  )
}
