.pipeline_cache_01 <- new.env(parent = emptyenv())

safe_get_field <- function(inventory_row, sample_row, inv_field, sample_field = inv_field, default = "") {
  if (!is.null(inventory_row) && inv_field %in% colnames(inventory_row)) {
    return(normalize_scalar_value(inventory_row[[inv_field]][1], default))
  }
  if (!is.null(sample_row) && sample_field %in% colnames(sample_row)) {
    return(normalize_scalar_value(sample_row[[sample_field]][1], default))
  }
  default
}

read_reference_annotation_cached <- function(cfg) {
  candidates <- c(cfg$clean_gtf, cfg$reference_gtf)
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  cache_key <- paste0("annotation::", candidates[1])
  if (exists(cache_key, envir = .pipeline_cache_01, inherits = FALSE)) {
    return(get(cache_key, envir = .pipeline_cache_01, inherits = FALSE))
  }

  annotation_df <- tryCatch(
    read_reference_annotation_local(cfg$clean_gtf, cfg$reference_gtf)$annotation_df,
    error = function(e) {
      warning(sprintf("无法读取参考注释，feature contract 将降级: %s", conditionMessage(e)), call. = FALSE)
      data.frame(stringsAsFactors = FALSE)
    }
  )
  assign(cache_key, annotation_df, envir = .pipeline_cache_01)
  annotation_df
}

guess_feature_name_profile <- function(gene_names) {
  gene_names <- as.character(gene_names)
  gene_names <- gene_names[nzchar(gene_names)]
  if (length(gene_names) == 0) {
    return("unknown")
  }
  gene_names <- head(gene_names, 2000)
  stripped <- strip_ensembl_version(gene_names)
  ensembl_like <- grepl("^ENS[A-Z0-9]*G[0-9]+$", stripped)
  symbol_like <- !ensembl_like & grepl("^[A-Za-z][A-Za-z0-9._-]*$", gene_names)
  ensembl_fraction <- mean(ensembl_like)
  symbol_fraction <- mean(symbol_like)
  if (symbol_fraction >= 0.80 && ensembl_fraction <= 0.10) {
    return("symbol")
  }
  if (ensembl_fraction >= 0.80 && symbol_fraction <= 0.10) {
    return("ensembl")
  }
  if (any(ensembl_like) && any(symbol_like)) {
    return("mixed")
  }
  "unknown"
}

normalize_feature_token <- function(x) {
  x <- toupper(strip_ensembl_version(trimws(as.character(x))))
  gsub("[-_]+", "-", x)
}

match_features_from_annotation <- function(gene_names, annotation_df, target_rows) {
  if (nrow(annotation_df) == 0 || nrow(target_rows) == 0) {
    return(character(0))
  }

  gene_names <- as.character(gene_names)
  gene_names_upper <- toupper(gene_names)
  gene_ids_stripped <- strip_ensembl_version(gene_names)
  keep <- gene_ids_stripped %in% target_rows$gene_id_stripped
  if (any(nzchar(target_rows$gene_name_upper))) {
    keep <- keep | gene_names_upper %in% target_rows$gene_name_upper
  }
  if ("preferred_gene_label_upper" %in% colnames(target_rows) && any(nzchar(target_rows$preferred_gene_label_upper))) {
    keep <- keep | gene_names_upper %in% target_rows$preferred_gene_label_upper
  }
  unique(gene_names[keep])
}

matched_gene_names_from_annotation <- function(features, annotation_df) {
  if (length(features) == 0 || nrow(annotation_df) == 0) {
    return(character(0))
  }
  feature_upper <- toupper(as.character(features))
  feature_stripped <- toupper(strip_ensembl_version(as.character(features)))
  matched <- annotation_df[
    toupper(annotation_df$gene_id_stripped) %in% feature_stripped |
      toupper(annotation_df$gene_name) %in% feature_upper,
    ,
    drop = FALSE
  ]
  if ("preferred_gene_label" %in% colnames(matched)) {
    labels <- matched$preferred_gene_label[nzchar(matched$preferred_gene_label)]
    if (length(labels) > 0) {
      return(unique(labels))
    }
  }
  unique(matched$gene_name[nzchar(matched$gene_name)])
}

match_features_from_identifier_list <- function(gene_names, identifiers, annotation_df) {
  identifiers <- unique(trimws(as.character(identifiers)))
  identifiers <- identifiers[nzchar(identifiers)]
  if (length(identifiers) == 0) {
    return(character(0))
  }

  gene_names <- as.character(gene_names)
  gene_names_upper <- normalize_feature_token(gene_names)
  gene_ids_stripped <- normalize_feature_token(gene_names)
  id_upper <- normalize_feature_token(identifiers)
  id_stripped <- normalize_feature_token(identifiers)

  keep <- gene_names_upper %in% id_upper | gene_ids_stripped %in% id_stripped
  if (nrow(annotation_df) > 0) {
    target_rows <- annotation_df[
      toupper(annotation_df$gene_id_stripped) %in% id_stripped |
        toupper(annotation_df$gene_name) %in% id_upper,
      ,
      drop = FALSE
    ]
    keep <- keep | gene_names %in% match_features_from_annotation(gene_names, annotation_df, target_rows)
  }

  unique(gene_names[keep])
}

canonical_vertebrate_mito_genes <- function() {
  c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB")
}

normalize_mito_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  sub("^MT-", "", x)
}

infer_mito_seqnames_from_annotation <- function(annotation_df) {
  if (nrow(annotation_df) == 0) {
    return(character(0))
  }

  seqname_upper <- toupper(annotation_df$seqname)
  direct_hits <- unique(annotation_df$seqname[
    seqname_upper %in% c("M", "MT", "CHRM", "CHRMT", "MITOCHONDRION_GENOME") |
      grepl("MITO", seqname_upper)
  ])
  if (length(direct_hits) > 0) {
    return(direct_hits)
  }

  canonical_hits <- annotation_df[normalize_mito_symbol(annotation_df$gene_name) %in% canonical_vertebrate_mito_genes(), , drop = FALSE]
  if (nrow(canonical_hits) == 0 && "preferred_gene_label" %in% colnames(annotation_df)) {
    canonical_hits <- annotation_df[
      normalize_mito_symbol(annotation_df$preferred_gene_label) %in% canonical_vertebrate_mito_genes(),
      ,
      drop = FALSE
    ]
  }
  if (nrow(canonical_hits) == 0) {
    return(character(0))
  }

  seq_counts <- sort(table(canonical_hits$seqname), decreasing = TRUE)
  best_count <- as.integer(seq_counts[[1]])
  if (is.na(best_count) || best_count < 10) {
    return(character(0))
  }
  names(seq_counts)[seq_counts == best_count]
}

infer_reference_species <- function(cfg, annotation_df, gene_names) {
  candidates <- tolower(c(
    normalize_scalar_value(cfg$reference_version),
    normalize_scalar_value(cfg$reference_gtf),
    normalize_scalar_value(cfg$clean_gtf),
    head(strip_ensembl_version(gene_names), 200)
  ))
  if (nrow(annotation_df) > 0) {
    candidates <- c(
      candidates,
      tolower(head(annotation_df$gene_id, 200)),
      tolower(head(annotation_df$gene_name, 200))
    )
  }
  candidates <- candidates[nzchar(candidates)]
  joined <- paste(candidates, collapse = " ")

  if (grepl("ensgalg|grcg|galgal|gallus|chicken|bgalgal", joined)) {
    return("chicken")
  }
  if (grepl("ensmusg|grcm|mm[0-9]+|mus musculus|mouse", joined)) {
    return("mouse")
  }
  if (grepl("ensg|grch|hg[0-9]+|homo sapiens|human", joined)) {
    return("human")
  }
  "unknown"
}

resolve_mito_feature_set <- function(gene_names, cfg, annotation_df, species) {
  method <- "failed"
  detail <- ""
  features <- character(0)
  detected_gene_names <- character(0)
  seqnames <- character(0)

  if (nrow(annotation_df) > 0) {
    mito_seqnames <- infer_mito_seqnames_from_annotation(annotation_df)
    if (length(mito_seqnames) > 0) {
      target_rows <- annotation_df[annotation_df$seqname %in% mito_seqnames, , drop = FALSE]
      features <- match_features_from_annotation(gene_names, annotation_df, target_rows)
      detected_gene_names <- matched_gene_names_from_annotation(features, annotation_df)
      if (length(features) > 0) {
        method <- "gtf"
        detail <- sprintf("seqname=%s", paste(mito_seqnames, collapse = ","))
        seqnames <- mito_seqnames
      }
    }
  }

  if (length(features) == 0) {
    user_genes <- read_identifier_list(cfg$mito_gene_list_file)
    features <- match_features_from_identifier_list(gene_names, user_genes, annotation_df)
    detected_gene_names <- matched_gene_names_from_annotation(features, annotation_df)
    if (length(features) > 0) {
      method <- "user_list"
      detail <- normalize_scalar_value(cfg$mito_gene_list_file)
    }
  }

  if (length(features) == 0) {
    gene_names_upper <- toupper(as.character(gene_names))
    fallback_patterns <- switch(
      species,
      human = c("^MT-"),
      mouse = c("^MT-"),
      chicken = c("^J6367-", "^(?:MT-)?(?:ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|ATP6|ATP8|CYTB)$"),
      c("^MT-", "^J6367-", "^(?:MT-)?(?:ND1|ND2|ND3|ND4|ND4L|ND5|ND6|COX1|COX2|COX3|ATP6|ATP8|CYTB)$")
    )
    for (pattern in fallback_patterns) {
      hits <- unique(as.character(gene_names)[grepl(pattern, gene_names_upper, perl = TRUE)])
      if (length(hits) > 0) {
        features <- unique(c(features, hits))
      }
    }
    detected_gene_names <- if (nrow(annotation_df) > 0) matched_gene_names_from_annotation(features, annotation_df) else features
    if (length(features) > 0) {
      method <- "prefix_fallback"
      detail <- sprintf("%s:%s", species, paste(fallback_patterns, collapse = ";"))
    }
  }

  warning_codes <- character(0)
  warning_messages <- character(0)
  vertebrate_like <- species %in% c("human", "mouse", "chicken")

  if (identical(method, "prefix_fallback")) {
    warning_codes <- c(warning_codes, "mito_prefix_fallback")
    warning_messages <- c(warning_messages, "使用 prefix fallback 识别线粒体基因，准确性不保证。")
  }
  if (length(features) == 0) {
    method <- "failed"
    detail <- if (nzchar(detail)) detail else "no_matching_mito_features"
    warning_codes <- c(warning_codes, "mito_detection_failed")
    warning_messages <- c(warning_messages, "未检测到线粒体基因，percent.mito 不可信。")
  } else if (vertebrate_like && length(features) < 10) {
    warning_codes <- c(warning_codes, "mito_feature_count_low")
    warning_messages <- c(
      warning_messages,
      sprintf("脊椎动物参考通常应命中约 13 个线粒体蛋白编码基因；当前仅检测到 %d 个。", length(features))
    )
  }

  list(
    mito_features = unique(features),
    mito_feature_count = length(unique(features)),
    mito_detection_method = method,
    mito_detection_detail = detail,
    mito_detected_gene_names = unique(detected_gene_names),
    mito_reference_seqnames = unique(seqnames),
    mito_warning_codes = unique(warning_codes),
    mito_warning_messages = unique(warning_messages),
    species_guess = species
  )
}

load_cell_cycle_genes_local <- function(cfg) {
  cc_path <- require_ortholog_cc_genes_local(cfg)
  mapped <- tryCatch(readRDS(cc_path), error = function(e) {
    stop(sprintf("无法读取 cc_genes.rds: %s。请先运行 workflow/03stages/00_ortholog.sh", cc_path), call. = FALSE)
  })
  if (is.list(mapped) && all(c("s.genes", "g2m.genes") %in% names(mapped))) {
    mapped$source <- "ortholog_manifest_cc_genes"
    return(mapped)
  }

  stop(
    sprintf(
      "01a 需要 00c 通过 ortholog manifest 暴露 cc_genes 输出。请先完成 00a -> 00b -> 00c: %s",
      cfg$ortholog_manifest_path
    ),
    call. = FALSE
  )
}

require_ortholog_cc_genes_local <- function(cfg) {
  if (!file.exists(cfg$ortholog_manifest_path)) {
    stop(sprintf("缺少 cc_genes.rds，请先运行 workflow/03stages/00_ortholog.sh。缺少 manifest: %s", cfg$ortholog_manifest_path), call. = FALSE)
  }
  manifest <- read_manifest_local(cfg$ortholog_manifest_path)
  cc_path <- tryCatch(resolve_output_local(manifest, "cc_genes"), error = function(e) {
    stop(sprintf("缺少 cc_genes.rds，请先运行 workflow/03stages/00_ortholog.sh。ortholog manifest 未暴露 cc_genes: %s", cfg$ortholog_manifest_path), call. = FALSE)
  })
  if (!file.exists(cc_path)) {
    stop(sprintf("缺少 cc_genes.rds，请先运行 workflow/03stages/00_ortholog.sh。manifest 指向的文件不存在: %s", cc_path), call. = FALSE)
  }
  cc_path
}

resolve_feature_context <- function(gene_names, cfg, declared_gene_id_type = "auto") {
  annotation_df <- read_reference_annotation_cached(cfg)
  feature_name_profile <- guess_feature_name_profile(gene_names)
  gene_id_type <- tolower(normalize_scalar_value(declared_gene_id_type, "auto"))
  if (!gene_id_type %in% c("symbol", "ensembl", "mixed", "unknown")) {
    gene_id_type <- feature_name_profile
  }
  if (!gene_id_type %in% c("symbol", "ensembl", "mixed", "unknown")) {
    gene_id_type <- "unknown"
  }

  ribo_features <- character(0)
  species <- infer_reference_species(cfg, annotation_df, gene_names)
  mito_payload <- resolve_mito_feature_set(gene_names, cfg, annotation_df, species)

  if (nrow(annotation_df) > 0) {
    biotype_upper <- toupper(annotation_df$gene_biotype)
    ribo_rows <- annotation_df[
      grepl("^(RPL|RPS|MRPL|MRPS)", annotation_df$gene_name_upper) |
        grepl("RIBOSOM", biotype_upper),
      ,
      drop = FALSE
    ]
    ribo_features <- match_features_from_annotation(gene_names, annotation_df, ribo_rows)
  }
  if (length(ribo_features) == 0 && gene_id_type %in% c("symbol", "mixed", "unknown")) {
    ribo_features <- unique(gene_names[grepl("^(RPL|RPS|Mrpl|Mrps|MRPL|MRPS)", gene_names)])
  }

  rbc_identifiers <- read_identifier_list(cfg$rbc_gene_list_file)
  rbc_features <- match_features_from_identifier_list(
    gene_names,
    rbc_identifiers,
    annotation_df
  )
  rbc_detected_gene_names <- if (nrow(annotation_df) > 0) {
    matched_gene_names_from_annotation(rbc_features, annotation_df)
  } else {
    rbc_features
  }
  rbc_detection_method <- if (length(rbc_identifiers) == 0) {
    "not_configured"
  } else if (length(rbc_features) == 0) {
    "configured_no_match"
  } else {
    "user_list"
  }

  cc_genes <- load_cell_cycle_genes_local(cfg)
  s_features <- match_features_from_identifier_list(gene_names, cc_genes$s.genes, annotation_df)
  g2m_features <- match_features_from_identifier_list(gene_names, cc_genes$g2m.genes, annotation_df)

  list(
    feature_name_profile = feature_name_profile,
    gene_id_type = gene_id_type,
    mito_features = unique(mito_payload$mito_features),
    ribo_features = unique(ribo_features),
    mito_feature_count = mito_payload$mito_feature_count,
    mito_detection_method = mito_payload$mito_detection_method,
    mito_detection_detail = mito_payload$mito_detection_detail,
    mito_detected_gene_names = unique(mito_payload$mito_detected_gene_names),
    mito_detected_gene_names_text = collapse_unique_values(mito_payload$mito_detected_gene_names),
    mito_detected_feature_names_text = collapse_unique_values(mito_payload$mito_features),
    mito_reference_seqnames = unique(mito_payload$mito_reference_seqnames),
    mito_reference_seqnames_text = collapse_unique_values(mito_payload$mito_reference_seqnames),
    mito_warning_codes = unique(mito_payload$mito_warning_codes),
    mito_warning_messages = unique(mito_payload$mito_warning_messages),
    mito_warning_codes_text = collapse_unique_values(mito_payload$mito_warning_codes),
    mito_warning_messages_text = collapse_unique_values(mito_payload$mito_warning_messages),
    species_guess = mito_payload$species_guess,
    ribo_feature_count = length(unique(ribo_features)),
    rbc_features = unique(rbc_features),
    rbc_feature_count = length(unique(rbc_features)),
    rbc_detected_gene_names = unique(rbc_detected_gene_names),
    rbc_detected_gene_names_text = collapse_unique_values(rbc_detected_gene_names),
    rbc_detection_method = rbc_detection_method,
    rbc_gene_list_file = normalize_scalar_value(cfg$rbc_gene_list_file),
    cell_cycle_s_features = unique(s_features),
    cell_cycle_g2m_features = unique(g2m_features),
    cell_cycle_s_feature_count = length(unique(s_features)),
    cell_cycle_g2m_feature_count = length(unique(g2m_features)),
    cell_cycle_gene_source = cc_genes$source,
    annotation_available = nrow(annotation_df) > 0
  )
}

add_basic_qc_metrics <- function(seu, cfg, declared_gene_id_type = "auto") {
  feature_context <- resolve_feature_context(rownames(seu), cfg, declared_gene_id_type = declared_gene_id_type)

  seu$percent.mito <- if (length(feature_context$mito_features) > 0) {
    Seurat::PercentageFeatureSet(seu, features = feature_context$mito_features)
  } else {
    0
  }
  seu$percent.ribo <- if (length(feature_context$ribo_features) > 0) {
    Seurat::PercentageFeatureSet(seu, features = feature_context$ribo_features)
  } else {
    0
  }
  seu$percent.rbc <- if (length(feature_context$rbc_features) > 0) {
    Seurat::PercentageFeatureSet(seu, features = feature_context$rbc_features)
  } else {
    0
  }
  seu$log10GenesPerUMI <- log10(seu$nFeature_RNA + 1) / log10(seu$nCount_RNA + 1)
  seu$log10GenesPerUMI[!is.finite(seu$log10GenesPerUMI)] <- 0
  seu@misc$feature_contract <- list(
    feature_name_profile = feature_context$feature_name_profile,
    gene_id_type = feature_context$gene_id_type,
    mito_feature_count = feature_context$mito_feature_count,
    mito_detection_method = feature_context$mito_detection_method,
    mito_detection_detail = feature_context$mito_detection_detail,
    mito_detected_gene_names = feature_context$mito_detected_gene_names,
    mito_detected_gene_names_text = feature_context$mito_detected_gene_names_text,
    mito_detected_feature_names_text = feature_context$mito_detected_feature_names_text,
    mito_reference_seqnames_text = feature_context$mito_reference_seqnames_text,
    mito_warning_codes = feature_context$mito_warning_codes,
    mito_warning_messages = feature_context$mito_warning_messages,
    mito_warning_codes_text = feature_context$mito_warning_codes_text,
    mito_warning_messages_text = feature_context$mito_warning_messages_text,
    species_guess = feature_context$species_guess,
    ribo_feature_count = feature_context$ribo_feature_count,
    rbc_feature_count = feature_context$rbc_feature_count,
    rbc_detected_gene_names = feature_context$rbc_detected_gene_names,
    rbc_detected_gene_names_text = feature_context$rbc_detected_gene_names_text,
    rbc_detection_method = feature_context$rbc_detection_method,
    rbc_gene_list_file = feature_context$rbc_gene_list_file,
    cell_cycle_s_feature_count = feature_context$cell_cycle_s_feature_count,
    cell_cycle_g2m_feature_count = feature_context$cell_cycle_g2m_feature_count,
    cell_cycle_gene_source = feature_context$cell_cycle_gene_source,
    annotation_available = feature_context$annotation_available
  )
  list(object = seu, feature_context = feature_context)
}

maybe_join_layers <- function(seu) {
  if (exists("JoinLayers", mode = "function")) {
    seu <- JoinLayers(seu)
  }
  seu
}

density_peaks <- function(x) {
  x <- stats::na.omit(as.numeric(x))
  if (length(x) < 100 || length(unique(x)) < 20) {
    return(1L)
  }
  dens <- stats::density(x, na.rm = TRUE)
  y <- dens$y
  max(1L, sum(diff(sign(diff(y))) < 0))
}
