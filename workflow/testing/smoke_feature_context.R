args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript smoke_feature_context.R /path/to/pipeline_root", call. = FALSE)
}

pipeline_root <- normalizePath(args[[1]], mustWork = TRUE)
helper_dir <- file.path(pipeline_root, "workflow", "05single_script", "helpers")
source(file.path(helper_dir, "metadata_io.R"))
source(file.path(helper_dir, "gtf_utils.R"))
source(file.path(helper_dir, "qc_utils.R"))

tmp_dir <- tempfile("scrna_feature_contract_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

load_cell_cycle_genes_local <- function(cfg) {
  list(
    s.genes = "FOO1",
    g2m.genes = "BAR1",
    source = "smoke_stub"
  )
}

make_cfg <- function(reference_gtf = "", clean_gtf = "", mito_gene_list_file = "", rbc_gene_list_file = "", reference_version = "") {
  list(
    reference_gtf = reference_gtf,
    clean_gtf = clean_gtf,
    mito_gene_list_file = mito_gene_list_file,
    rbc_gene_list_file = rbc_gene_list_file,
    reference_version = reference_version
  )
}

gtf_path <- file.path(tmp_dir, "mini.gtf")
writeLines(
  c(
    'MT\tmock\tgene\t1\t1000\t.\t+\t.\tgene_id "ENSGALG00000000001"; gene_name "MT-CO1"; gene_biotype "protein_coding";',
    '1\tmock\tgene\t1\t1000\t.\t+\t.\tgene_id "ENSGALG00000000002"; gene_name "RPLP0"; gene_biotype "protein_coding";',
    '1\tmock\tgene\t1\t1000\t.\t+\t.\tgene_id "ENSGALG00000000003"; gene_name "FOO1"; gene_biotype "protein_coding";',
    '1\tmock\tgene\t1\t1000\t.\t+\t.\tgene_id "ENSGALG00000000004"; gene_name "HBAA"; gene_biotype "protein_coding";'
  ),
  con = gtf_path
)

gene_names <- c(
  "ENSGALG00000000001",
  "ENSGALG00000000002",
  "ENSGALG00000000003",
  "ENSGALG00000000004"
)

payload <- resolve_feature_context(
  gene_names,
  cfg = make_cfg(reference_gtf = gtf_path),
  declared_gene_id_type = "ensembl"
)

stopifnot(identical(payload$feature_name_profile, "ensembl"))
stopifnot(identical(payload$mito_detection_method, "gtf"))
stopifnot(payload$mito_feature_count >= 1)
stopifnot(payload$ribo_feature_count >= 1)
stopifnot(identical(payload$rbc_detection_method, "not_configured"))
stopifnot(payload$rbc_feature_count == 0)

rbc_empty_example <- read_identifier_list(file.path(pipeline_root, "config", "rbc_gene_list.txt"))
stopifnot(length(rbc_empty_example) == 0)

rbc_list_path <- file.path(tmp_dir, "rbc_gene_list.txt")
writeLines(c("# one gene per line", "ENSGALG00000000004"), con = rbc_list_path)
payload_rbc <- resolve_feature_context(
  gene_names,
  cfg = make_cfg(reference_gtf = gtf_path, rbc_gene_list_file = rbc_list_path),
  declared_gene_id_type = "ensembl"
)

stopifnot(identical(payload_rbc$rbc_detection_method, "user_list"))
stopifnot(payload_rbc$rbc_feature_count == 1)
stopifnot(identical(payload_rbc$rbc_features, "ENSGALG00000000004"))
stopifnot(identical(payload_rbc$rbc_detected_gene_names_text, "HBAA"))

rbc_no_match_path <- file.path(tmp_dir, "rbc_gene_list_no_match.txt")
writeLines(c("# one gene per line", "ENSGALG00000099999"), con = rbc_no_match_path)
payload_rbc_no_match <- resolve_feature_context(
  gene_names,
  cfg = make_cfg(reference_gtf = gtf_path, rbc_gene_list_file = rbc_no_match_path),
  declared_gene_id_type = "ensembl"
)

stopifnot(identical(payload_rbc_no_match$rbc_detection_method, "configured_no_match"))
stopifnot(payload_rbc_no_match$rbc_feature_count == 0)

mito_list_path <- file.path(tmp_dir, "mito_gene_list.txt")
writeLines(c("# one gene per line", "CUSTOM_MT"), con = mito_list_path)
payload_user_mito <- resolve_feature_context(
  c("CUSTOM-MT", "RPLP0", "FOO1"),
  cfg = make_cfg(mito_gene_list_file = mito_list_path),
  declared_gene_id_type = "symbol"
)

stopifnot(identical(payload_user_mito$mito_detection_method, "user_list"))
stopifnot(payload_user_mito$mito_feature_count == 1)
stopifnot(payload_user_mito$ribo_feature_count >= 1)

payload_prefix <- resolve_feature_context(
  c("ND1", "COX3", "ATP6AP1"),
  cfg = make_cfg(reference_version = "GRCg7b"),
  declared_gene_id_type = "symbol"
)

stopifnot(identical(payload_prefix$mito_detection_method, "prefix_fallback"))
stopifnot(payload_prefix$mito_feature_count == 2)
stopifnot(!"ATP6AP1" %in% payload_prefix$mito_features)
stopifnot(grepl("mito_prefix_fallback", payload_prefix$mito_warning_codes_text, fixed = TRUE))

cat("smoke_feature_context: PASS\n")
