suppressPackageStartupMessages({
  library(Seurat)
})

must_getenv <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || identical(value, "")) {
    stop(sprintf("缺少环境变量: %s", name), call. = FALSE)
  }
  value
}

split_csv <- function(x) {
  if (length(x) == 0) {
    return(character(0))
  }
  x <- as.character(x[[1]])
  if (is.na(x)) {
    return(character(0))
  }
  x <- trimws(x)
  if (!nzchar(x)) {
    return(character(0))
  }
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_index_spec <- function(x) {
  x <- trimws(x)
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    parts <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    return(seq.int(parts[1], parts[2]))
  }
  as.integer(split_csv(x))
}

read_cfg <- function() {
  list(
    project_root = must_getenv("PROJECT_ROOT"),
    data_dir = must_getenv("DATA_DIR"),
    results_dir = must_getenv("RESULTS_DIR"),
    checkpoint_dir = must_getenv("CHECKPOINT_DIR"),
    figure_dir = must_getenv("FIGURE_DIR"),
    table_dir = must_getenv("TABLE_DIR"),
    log_dir = must_getenv("LOG_DIR"),
    report_dir = Sys.getenv("REPORT_DIR", file.path(must_getenv("PROJECT_ROOT"), "reports")),
    intake_report_dir = Sys.getenv("INTAKE_REPORT_DIR", file.path(Sys.getenv("REPORT_DIR", file.path(must_getenv("PROJECT_ROOT"), "reports")), "intake")),
    eda_report_dir = Sys.getenv("EDA_REPORT_DIR", file.path(Sys.getenv("REPORT_DIR", file.path(must_getenv("PROJECT_ROOT"), "reports")), "eda")),
    cellranger_out_dir = Sys.getenv("CELLRANGER_OUT_DIR", ""),
    velocity_input_dir = Sys.getenv("VELOCITY_INPUT_DIR", ""),
    velocity_loom_dir = Sys.getenv("VELOCITY_LOOM_DIR", ""),
    velocity_output_dir = Sys.getenv("VELOCITY_OUTPUT_DIR", ""),
    scenic_input_dir = Sys.getenv("SCENIC_INPUT_DIR", ""),
    scenic_output_dir = Sys.getenv("SCENIC_OUTPUT_DIR", ""),
    sample_sheet = Sys.getenv("SAMPLE_SHEET", ""),
    canonical_sample_sheet = Sys.getenv("CANONICAL_SAMPLE_SHEET", ""),
    comparison_sheet = Sys.getenv("COMPARISON_SHEET", ""),
    input_inventory_file = Sys.getenv("INPUT_INVENTORY_FILE", ""),
    branch_readiness_file = Sys.getenv("BRANCH_READINESS_FILE", ""),
    reference_dir = Sys.getenv("REFERENCE_DIR", ""),
    reference_gtf = Sys.getenv("REFERENCE_GTF", ""),
    clean_gtf = Sys.getenv("CLEAN_GTF", ""),
    qc_threshold_file = Sys.getenv("QC_THRESHOLD_FILE", ""),
    eda_gate_file = Sys.getenv("EDA_GATE_FILE", ""),
    object_layer_config_file = Sys.getenv("OBJECT_LAYER_CONFIG_FILE", file.path(must_getenv("PROJECT_ROOT"), "config", "object_layers.tsv")),
    marker_panel_dir = Sys.getenv("MARKER_PANEL_DIR", file.path(must_getenv("PROJECT_ROOT"), "config", "marker_panels")),
    mito_gene_list_file = Sys.getenv("MITO_GENE_LIST_FILE", file.path(must_getenv("PROJECT_ROOT"), "config", "mito_gene_list.txt")),
    min_biological_replicates = as.integer(Sys.getenv("MIN_BIOLOGICAL_REPLICATES", "2")),
    sample_names = split_csv(must_getenv("SAMPLE_NAMES")),
    analysis_group_1_name = Sys.getenv("ANALYSIS_GROUP_1_NAME", ""),
    analysis_group_1_samples = split_csv(Sys.getenv("ANALYSIS_GROUP_1_SAMPLES", "")),
    analysis_group_2_name = Sys.getenv("ANALYSIS_GROUP_2_NAME", ""),
    analysis_group_2_samples = split_csv(Sys.getenv("ANALYSIS_GROUP_2_SAMPLES", "")),
    raw_group_1_name = Sys.getenv("RAW_GROUP_1_NAME", "Group1"),
    raw_group_1_samples = split_csv(Sys.getenv("RAW_GROUP_1_SAMPLES", "")),
    raw_group_2_name = Sys.getenv("RAW_GROUP_2_NAME", "Group2"),
    raw_group_2_samples = split_csv(Sys.getenv("RAW_GROUP_2_SAMPLES", "")),
    raw_samples = split_csv(Sys.getenv("RAW_SAMPLES", "")),
    deg_ident_1 = Sys.getenv("DEG_IDENT_1", ""),
    deg_ident_2 = Sys.getenv("DEG_IDENT_2", ""),
    random_seed = as.integer(Sys.getenv("RANDOM_SEED", "42")),
    qc_min_nfeature = as.numeric(Sys.getenv("QC_MIN_NFEATURE", "200")),
    qc_min_ncount = as.numeric(Sys.getenv("QC_MIN_NCOUNT", "1000")),
    qc_min_log10umi = as.numeric(Sys.getenv("QC_MIN_LOG10UMI", "0.7")),
    qc_max_mito_pct = as.numeric(Sys.getenv("QC_MAX_MITO_PCT", "20")),
    ambient_primary_method = tolower(Sys.getenv("AMBIENT_PRIMARY_METHOD", "soupx")),
    ambient_fallback_method = tolower(Sys.getenv("AMBIENT_FALLBACK_METHOD", "decontx")),
    ambient_apply_policy = tolower(Sys.getenv("AMBIENT_APPLY_POLICY", "manual")),
    ambient_min_cells = as.integer(Sys.getenv("AMBIENT_MIN_CELLS", "50")),
    ambient_cluster_dims = parse_index_spec(Sys.getenv("AMBIENT_CLUSTER_DIMS", "1:20")),
    ambient_cluster_resolution = as.numeric(Sys.getenv("AMBIENT_CLUSTER_RESOLUTION", "0.4")),
    ambient_marker_top_n = as.integer(Sys.getenv("AMBIENT_MARKER_TOP_N", "3")),
    ambient_recommend_min_contamination = as.numeric(Sys.getenv("AMBIENT_RECOMMEND_MIN_CONTAMINATION", "0.05")),
    cellbender_mode = tolower(Sys.getenv("CELLBENDER_MODE", "stub")),
    cellbender_fpr = as.numeric(Sys.getenv("CELLBENDER_FPR", "0.01")),
    cellbender_cuda = tolower(Sys.getenv("CELLBENDER_CUDA", "yes")),
    cellbender_extra_args = Sys.getenv("CELLBENDER_EXTRA_ARGS", ""),
    doublet_rate = as.numeric(Sys.getenv("DOUBLET_RATE", "0.008")),
    doublet_rate_per_1k = as.numeric(Sys.getenv("DOUBLET_RATE_PER_1K", Sys.getenv("DOUBLET_RATE", "0.008"))),
    doublet_primary_caller = tolower(Sys.getenv("DOUBLET_PRIMARY_CALLER", "scDblFinder")),
    doublet_secondary_caller = Sys.getenv("DOUBLET_SECONDARY_CALLER", "DoubletFinder"),
    doublet_secondary_enabled = tolower(Sys.getenv("DOUBLET_SECONDARY_ENABLED", "yes")) %in% c("yes", "true", "1", "on"),
    doublet_min_cells = as.integer(Sys.getenv("DOUBLET_MIN_CELLS", "50")),
    doublet_dims = parse_index_spec(Sys.getenv("DOUBLET_DIMS", "1:20")),
    hvg_nfeatures = as.integer(Sys.getenv("HVG_NFEATURES", "2000")),
    pca_dims = parse_index_spec(Sys.getenv("PCA_DIMS", "1:30")),
    target_clusters = as.integer(Sys.getenv("TARGET_CLUSTERS", "15")),
    res_range = as.numeric(split_csv(Sys.getenv("RES_RANGE", "0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60"))),
    res_fine_step = as.numeric(Sys.getenv("RES_FINE_STEP", "0.005")),
    triage_frac_below_cutoff = as.numeric(Sys.getenv("TRIAGE_FRAC_BELOW_CUTOFF", "0.35")),
    triage_frac_above_mito = as.numeric(Sys.getenv("TRIAGE_FRAC_ABOVE_MITO", "0.25")),
    triage_density_peaks = as.integer(Sys.getenv("TRIAGE_DENSITY_PEAKS", "2")),
    integration_mode = tolower(Sys.getenv("INTEGRATION_MODE", "harmony")),
    trajectory_start = Sys.getenv("TRAJECTORY_START", ""),
    trajectory_coarse_label = Sys.getenv("TRAJECTORY_COARSE_LABEL", "cell_type"),
    trajectory_fine_label = Sys.getenv("TRAJECTORY_FINE_LABEL", "seurat_clusters"),
    trajectory_fine_start_cluster = Sys.getenv("TRAJECTORY_FINE_START_CLUSTER", ""),
    trajectory_fine_top_n = as.integer(Sys.getenv("TRAJECTORY_FINE_TOP_N", "4")),
    trajectory_tradeseq_label = Sys.getenv("TRAJECTORY_TRADESEQ_LABEL", "cell_type"),
    tradeseq_knots = as.integer(Sys.getenv("TRADESEQ_KNOTS", "6")),
    ensembl_mirror = Sys.getenv("ENSEMBL_MIRROR", "asia")
  )
}

prepare_dirs <- function(cfg) {
  dirs <- c(
    cfg$results_dir,
    cfg$checkpoint_dir,
    cfg$figure_dir,
    cfg$table_dir,
    cfg$log_dir,
    cfg$report_dir,
    cfg$intake_report_dir,
    cfg$eda_report_dir,
    file.path(cfg$eda_report_dir, "pre_qc"),
    file.path(cfg$eda_report_dir, "ambient"),
    file.path(cfg$eda_report_dir, "post_qc"),
    file.path(cfg$eda_report_dir, "integration"),
    file.path(cfg$eda_report_dir, "annotation"),
    cfg$velocity_input_dir,
    cfg$velocity_output_dir,
    cfg$scenic_input_dir,
    cfg$scenic_output_dir
  )
  dirs <- dirs[nzchar(dirs)]
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

normalize_flag <- function(x, default = "yes") {
  if (length(x) == 0) {
    return(character(0))
  }
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x[x == ""] <- default
  tolower(x)
}

active_sample_sheet_path <- function(cfg) {
  if (nzchar(cfg$canonical_sample_sheet) && file.exists(cfg$canonical_sample_sheet)) {
    return(cfg$canonical_sample_sheet)
  }
  cfg$sample_sheet
}

read_tsv_optional <- function(path) {
  if (!nzchar(path) || !file.exists(path) || isTRUE(file.info(path)$size == 0)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  # Filter comment/blank lines before parsing so project templates with
  # explanatory headers remain readable and load consistently across locales.
  lines <- sub("^\ufeff", "", lines)
  trimmed <- trimws(lines)
  keep <- nzchar(trimmed) & !startsWith(trimmed, "#")
  if (!any(keep)) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  con <- textConnection(lines[keep])
  on.exit(close(con), add = TRUE)
  read.delim(
    con,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = "",
    fill = TRUE
  )
}

read_sample_sheet <- function(cfg) {
  path <- active_sample_sheet_path(cfg)
  df <- read_tsv_optional(path)
  if (nrow(df) == 0) {
    df <- data.frame(
      sample_id = cfg$sample_names,
      condition = cfg$sample_names,
      biological_replicate = cfg$sample_names,
      technical_replicate = "1",
      batch = "default",
      input_mode = "",
      input_source = "",
      source_path = "",
      platform = "auto",
      gene_id_type = "auto",
      reference_version = "",
      group_id = "",
      timepoint = "",
      tissue = "",
      chemistry = "",
      run_main = "yes",
      run_velocity = "yes",
      run_scenic = "yes",
      stringsAsFactors = FALSE
    )
  }

  expected_cols <- c(
    "sample_id", "condition", "biological_replicate", "technical_replicate",
    "batch", "input_mode", "input_source", "source_path",
    "platform", "gene_id_type", "reference_version", "group_id", "timepoint", "tissue", "chemistry",
    "run_main", "run_velocity", "run_scenic"
  )
  for (col in expected_cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  df
}

main_sample_sheet <- function(cfg) {
  df <- read_sample_sheet(cfg)
  df[normalize_flag(df$run_main, "yes") != "no", , drop = FALSE]
}

main_sample_names <- function(cfg) {
  df <- main_sample_sheet(cfg)
  sample_ids <- unique(df$sample_id)
  sample_ids[nzchar(sample_ids)]
}

read_input_inventory <- function(cfg) {
  read_tsv_optional(cfg$input_inventory_file)
}

read_branch_readiness <- function(cfg) {
  read_tsv_optional(cfg$branch_readiness_file)
}

read_qc_threshold_overrides <- function(cfg) {
  df <- read_tsv_optional(cfg$qc_threshold_file)
  if (nrow(df) == 0) {
    return(df)
  }

  required <- c("sample_id", "qc_min_nfeature", "qc_min_ncount", "qc_min_log10umi", "qc_max_mito_pct")
  for (col in required) {
    if (!col %in% colnames(df)) {
      df[[col]] <- NA_character_
    }
  }

  numeric_cols <- setdiff(required, "sample_id")
  for (col in numeric_cols) {
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }
  df
}

get_sample_qc_thresholds <- function(cfg, sample_id, overrides = NULL) {
  thresholds <- list(
    qc_min_nfeature = cfg$qc_min_nfeature,
    qc_min_ncount = cfg$qc_min_ncount,
    qc_min_log10umi = cfg$qc_min_log10umi,
    qc_max_mito_pct = cfg$qc_max_mito_pct
  )

  if (is.null(overrides)) {
    overrides <- read_qc_threshold_overrides(cfg)
  }
  if (nrow(overrides) == 0 || !"sample_id" %in% colnames(overrides)) {
    return(thresholds)
  }

  hit <- overrides[overrides$sample_id == sample_id, , drop = FALSE]
  if (nrow(hit) == 0) {
    return(thresholds)
  }

  for (name in names(thresholds)) {
    override_value <- suppressWarnings(as.numeric(hit[[name]][1]))
    if (!is.na(override_value)) {
      thresholds[[name]] <- override_value
    }
  }
  thresholds
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

write_markdown <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, con = path, useBytes = TRUE)
}

ensure_eda_stage_dir <- function(cfg, stage_name) {
  stage_dir <- file.path(cfg$eda_report_dir, stage_name)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  stage_dir
}

resolve_inventory_row <- function(cfg, sample_id, inventory_df = NULL) {
  if (is.null(inventory_df)) {
    inventory_df <- read_input_inventory(cfg)
  }
  if (nrow(inventory_df) == 0 || !"sample_id" %in% colnames(inventory_df)) {
    return(NULL)
  }
  hit <- inventory_df[inventory_df$sample_id == sample_id, , drop = FALSE]
  if (nrow(hit) == 0) {
    return(NULL)
  }
  hit[1, , drop = FALSE]
}

.pipeline_cache <- new.env(parent = emptyenv())

first_existing_path <- function(paths) {
  paths <- paths[nzchar(paths)]
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return("")
  }
  existing[1]
}

normalize_scalar_value <- function(x, default = "") {
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(default)
  }
  trimws(as.character(x))
}

display_scalar_value <- function(x, default = "NA") {
  value <- normalize_scalar_value(x)
  if (nzchar(value)) {
    return(value)
  }
  default
}

safe_get_field <- function(inventory_row, sample_row, inv_field, sample_field = inv_field, default = "") {
  if (!is.null(inventory_row) && inv_field %in% colnames(inventory_row)) {
    return(normalize_scalar_value(inventory_row[[inv_field]][1], default))
  }
  if (!is.null(sample_row) && sample_field %in% colnames(sample_row)) {
    return(normalize_scalar_value(sample_row[[sample_field]][1], default))
  }
  default
}

resolve_doubletfinder_fns <- function(pkg_env = NULL) {
  if (is.null(pkg_env)) {
    pkg_env <- asNamespace("DoubletFinder")
  }

  sweep_name <- if (exists("paramSweep_v3", envir = pkg_env, mode = "function")) {
    "paramSweep_v3"
  } else {
    "paramSweep"
  }
  doublet_name <- if (exists("doubletFinder_v3", envir = pkg_env, mode = "function")) {
    "doubletFinder_v3"
  } else {
    "doubletFinder"
  }

  list(
    paramSweep = get(sweep_name, envir = pkg_env),
    summarizeSweep = get("summarizeSweep", envir = pkg_env),
    find.pK = get("find.pK", envir = pkg_env),
    doubletFinder = get(doublet_name, envir = pkg_env)
  )
}

safe_merge_meta <- function(meta_df, new_df, by = "cell_id", cols = NULL) {
  if (is.null(new_df) || nrow(new_df) == 0) {
    return(meta_df)
  }
  if (!is.null(cols)) {
    keep_cols <- unique(c(by, cols))
    keep_cols <- keep_cols[keep_cols %in% colnames(new_df)]
    new_df <- new_df[, keep_cols, drop = FALSE]
  }

  result <- dplyr::left_join(meta_df, new_df, by = by, suffix = c("", ".new"))
  new_cols <- grep("\\.new$", colnames(result), value = TRUE)
  for (new_col in new_cols) {
    base_col <- sub("\\.new$", "", new_col)
    if (base_col %in% colnames(meta_df)) {
      use_new <- !is.na(result[[new_col]])
      result[[base_col]][use_new] <- result[[new_col]][use_new]
    } else {
      result[[base_col]] <- result[[new_col]]
    }
  }

  result[, !grepl("\\.new$", colnames(result)), drop = FALSE]
}

safe_process_sample <- function(sample_id, fn, on_error = NULL) {
  tryCatch(
    fn(),
    error = function(e) {
      warning(
        sprintf("样本 %s 处理失败: %s", sample_id, conditionMessage(e)),
        call. = FALSE
      )
      if (is.function(on_error)) {
        return(on_error(e))
      }
      NULL
    }
  )
}

TRIAGE_SIGNAL_REGISTRY <- list(
  low_complexity_burden = list(
    display_name = "大量细胞未达标",
    suspected_issue = "超过当前 QC 最低阈值的细胞比例偏高，说明样本整体质量可能较差，或当前阈值对这个样本来说过严。",
    recommended_action = "先打开 violin 图和 cutoff burden 图查看分布。如果低质量细胞形成独立峰，沿用当前阈值；如果分布连续，考虑为该样本单独调整阈值。"
  ),
  mito_detection_failed = list(
    display_name = "线粒体基因识别失败",
    suspected_issue = "未能稳定识别线粒体基因，当前 percent.mito 不可信。",
    recommended_action = "优先使用 reference GTF 染色体位点或 config/mito_gene_list.txt 明确线粒体基因列表，再继续解释或依赖 percent.mito。"
  ),
  mito_prefix_fallback = list(
    display_name = "线粒体识别仅依赖前缀匹配",
    suspected_issue = "当前样本没有直接使用 reference/GTF 命中线粒体基因，而是退回到前缀匹配。",
    recommended_action = "结果可用于初步 QC，但进入正式解释前应尽快用 reference GTF 或用户显式列表复核。"
  ),
  mito_feature_count_low = list(
    display_name = "线粒体基因命中数偏少",
    suspected_issue = "虽然识别到了线粒体基因，但命中数偏少，percent.mito 的稳定性仍然可疑。",
    recommended_action = "先检查 GTF 染色体命名、features 列选择和 mito_gene_list.txt 是否需要显式补充。"
  ),
  high_mito_burden = list(
    display_name = "高线粒体负担",
    suspected_issue = "高 mt 细胞比例偏高，样本可能存在明显的低质量负担，也可能混入特定高 mt 群体。",
    recommended_action = "先判断这是全样本质量问题还是特定细胞群体的生物信号，再决定是否收紧 mt cutoff。"
  ),
  possible_multimodal_qc = list(
    display_name = "QC 分布疑似多峰",
    suspected_issue = "QC 指标分布不太像单一总体，说明同一样本里可能混有不同质量层或不同输入状态。",
    recommended_action = "先标记为人工复核，避免直接套单一固定 cutoff 清洗。"
  ),
  low_genes_low_saturation = list(
    display_name = "低基因数且饱和度仍偏低",
    suspected_issue = "每细胞基因数偏低，同时测序饱和度也偏低，说明追加测序可能仍有收益。",
    recommended_action = "继续核查文库复杂度，同时评估是否值得追加测序。"
  ),
  low_genes_high_saturation = list(
    display_name = "低基因数但饱和度已高",
    suspected_issue = "每细胞基因数偏低，但测序饱和度已接近平稳，问题更像是样本或文库复杂度本身受限。",
    recommended_action = "优先回头检查样本与文库复杂度，不要默认把问题归因于测序深度不足。"
  ),
  ambient_soupx_ready = list(
    display_name = "具备 SoupX 主路径条件",
    suspected_issue = "该样本已经具备 raw droplets 等前提，ambient 分支可以按主路径执行。",
    recommended_action = "进入 ambient branch 时优先走 SoupX，并输出 contamination summary 与前后 marker leakage 对比。"
  ),
  ambient_decontx_fallback = list(
    display_name = "仅具备 DecontX fallback 条件",
    suspected_issue = "该样本不满足 SoupX 主路径要求，但仍可走对象级的 DecontX fallback。",
    recommended_action = "在 ambient 报告中明确标记为 DecontX fallback，不要与 SoupX 等价解释。"
  ),
  gene_id_type_unrecognized = list(
    display_name = "gene_id_type 未识别",
    suspected_issue = "feature naming profile 无法稳定识别 gene_id_type，后续 percent.mito、marker 命中和同源映射解释都存在风险。",
    recommended_action = "先确认 features.tsv 的命名列和 reference GTF 是否匹配，再继续解释 percent.mito 或 marker 命中。"
  ),
  raw_matrix_missing_for_soupx = list(
    display_name = "缺少 SoupX 所需原始矩阵",
    suspected_issue = "该样本没有可用于 SoupX 的 raw droplets 输入，ambient 主路径不可用。",
    recommended_action = "在 ambient 报告里确认是否降级到 DecontX fallback，并把“可做 ambient”与“已完成 correction”明确区分。"
  ),
  input_contract_limits = list(
    display_name = "输入契约存在限制",
    suspected_issue = "当前样本的输入契约仍有边界条件，后续解释可能受到约束。",
    recommended_action = "结合 sample_qc_summary.tsv 一起确认这些限制是否影响继续推进。"
  ),
  project_readiness_limits = list(
    display_name = "项目级 readiness 仍有限制",
    suspected_issue = "项目级 intake/readiness 还没有完全放行，继续推进前需要确认这些限制是否可接受。",
    recommended_action = "在进入后续整合和分支模块前，先显式记录这些限制及接受理由。"
  ),
  heavy_qc_loss = list(
    display_name = "QC 保留率过低",
    suspected_issue = "该样本在 QC 后保留率偏低，说明阈值或样本质量都需要重新审视。",
    recommended_action = "优先回看 pre-QC 报告中的 cutoff burden 和样本级质量分布。"
  ),
  high_primary_doublet_rate = list(
    display_name = "primary doublet rate 偏高",
    suspected_issue = "默认主调用器识别出的 doublet 比例偏高，后续聚类和整合容易受到影响。",
    recommended_action = "重点审阅 doublet UMAP、score 分布与高风险 cluster。"
  ),
  primary_secondary_discordance = list(
    display_name = "主次 doublet 调用器不一致",
    suspected_issue = "primary 与 secondary 调用器对同一批细胞给出了冲突结论。",
    recommended_action = "在保留 primary 为默认主调用器的前提下，重点复核冲突区域。"
  ),
  sample_collapse = list(
    display_name = "样本在 QC/doublet 后接近塌缩",
    suspected_issue = "经过 QC 和 doublet 过滤后，该样本剩余细胞过少，后续结构解释会明显变得不稳定。",
    recommended_action = "在继续整合和注释前先评估这个样本是否仍适合保留。"
  ),
  ambient_recommended_not_applied = list(
    display_name = "ambient 建议替换但未应用",
    suspected_issue = "ambient 分支认为应替换 counts，但主对象当前仍保留原始 counts。",
    recommended_action = "进入后续整合和注释前，应把这个差异显式记录在报告中。"
  ),
  doublet_high_risk_cluster = list(
    display_name = "存在 doublet 高风险 cluster",
    suspected_issue = "某些 provisional cluster 中的 doublet 富集明显高于样本背景。",
    recommended_action = "结合 UMAP 空间和 marker 混合模式复核这些 cluster。"
  ),
  batch_dominant_pca_axes = list(
    display_name = "PCA 主轴被样本来源主导",
    suspected_issue = "未整合空间中的主导变异更像 sample/batch，而不是目标生物学结构。",
    recommended_action = "优先检查是否需要 batch correction 或 integration，再决定是否沿用未整合空间。"
  ),
  possible_overcorrection = list(
    display_name = "可能发生过度校正",
    suspected_issue = "整合后条件结构被明显抹平，可能已经开始损失你真正关心的生物学差异。",
    recommended_action = "如果整合后条件结构被明显抹平，优先改用更保守策略或退回未整合分析。"
  ),
  annotation_conflict = list(
    display_name = "注释证据存在冲突",
    suspected_issue = "同一层中存在 marker 证据与最终命名不一致的 cluster，需要单独审阅。",
    recommended_action = "回看 annotation evidence 表，检查 marker panel 是否重叠或 cluster 是否需要进一步拆分。"
  ),
  high_undetermined_ratio = list(
    display_name = "未定注释比例偏高",
    suspected_issue = "该层里未定 cluster 占比偏高，说明当前 panel 或证据还不足以稳定命名。",
    recommended_action = "优先补充 marker panel，或把这一层保留为待定结果而不是强行命名。"
  ),
  no_marker_panel = list(
    display_name = "缺少 marker panel",
    suspected_issue = "该层当前没有 marker panel，注释主要依赖数据驱动证据。",
    recommended_action = "在 config/marker_panels/ 中补充组织特异性 marker panel，再复核该层注释。"
  )
)

TRIAGE_SIGNAL_PREFIX_REGISTRY <- list(
  cell_cycle_axis_ = list(
    display_name = "细胞周期主导轴",
    suspected_issue = "当前主导降维轴中包含较多细胞周期高载荷基因，cell cycle 可能正在主导结构变化。",
    recommended_action = "只有在 cell cycle 明显遮蔽了你的核心科学问题时，再考虑回归，而不是默认直接回归。"
  ),
  stress_axis_ = list(
    display_name = "应激/即时早期基因主导轴",
    suspected_issue = "当前主导降维轴中富集 stress/IEG 高载荷基因，更像 dissociation stress 或技术扰动。",
    recommended_action = "优先怀疑应激来源，不要立刻把它解释成新的稳定细胞状态。"
  )
)

resolve_triage_signal_spec <- function(signal_id) {
  signal_id <- normalize_scalar_value(signal_id)
  if (!nzchar(signal_id)) {
    return(list(
      display_name = "未命名信号",
      suspected_issue = "",
      recommended_action = ""
    ))
  }

  if (signal_id %in% names(TRIAGE_SIGNAL_REGISTRY)) {
    return(TRIAGE_SIGNAL_REGISTRY[[signal_id]])
  }

  matched_prefixes <- names(TRIAGE_SIGNAL_PREFIX_REGISTRY)[
    startsWith(signal_id, names(TRIAGE_SIGNAL_PREFIX_REGISTRY))
  ]
  if (length(matched_prefixes) > 0) {
    prefix <- matched_prefixes[[which.max(nchar(matched_prefixes))]]
    return(TRIAGE_SIGNAL_PREFIX_REGISTRY[[prefix]])
  }

  list(
    display_name = signal_id,
    suspected_issue = "",
    recommended_action = ""
  )
}

empty_triage_df <- function(include_sample = TRUE) {
  cols <- list(
    severity = character(0),
    signal_id = character(0),
    display_name = character(0),
    suspected_issue = character(0),
    evidence = character(0),
    recommended_action = character(0),
    manual_review_required = character(0)
  )
  if (isTRUE(include_sample)) {
    cols <- c(list(sample_id = character(0)), cols)
  }
  do.call(
    data.frame,
    c(cols, list(stringsAsFactors = FALSE, check.names = FALSE))
  )
}

make_triage_row <- function(sample_id = NULL, severity, signal_id, evidence = "", suspected_issue = "", recommended_action = "", manual_review_required = "yes") {
  spec <- resolve_triage_signal_spec(signal_id)
  row <- data.frame(
    severity = severity,
    signal_id = signal_id,
    display_name = display_scalar_value(spec$display_name, signal_id),
    suspected_issue = display_scalar_value(suspected_issue, display_scalar_value(spec$suspected_issue, "待补充")),
    evidence = display_scalar_value(evidence, "待补充"),
    recommended_action = display_scalar_value(recommended_action, display_scalar_value(spec$recommended_action, "待补充")),
    manual_review_required = display_scalar_value(manual_review_required, "yes"),
    stringsAsFactors = FALSE
  )
  if (!is.null(sample_id)) {
    row <- cbind(data.frame(sample_id = display_scalar_value(sample_id, "NA"), stringsAsFactors = FALSE), row)
  }
  row
}

triage_subject_label <- function(sample_id) {
  sample_id <- normalize_scalar_value(sample_id)
  if (!nzchar(sample_id)) {
    return("")
  }
  if (identical(sample_id, "__PROJECT__")) {
    return("项目级")
  }
  sprintf("`%s`", sample_id)
}

render_triage_markdown <- function(triage_df, include_sample = TRUE) {
  if (is.null(triage_df) || nrow(triage_df) == 0) {
    return("- 当前没有 triage 信号。")
  }

  lines <- character(0)
  for (i in seq_len(nrow(triage_df))) {
    row <- triage_df[i, , drop = FALSE]
    header_parts <- c()
    if (isTRUE(include_sample) && "sample_id" %in% colnames(row)) {
      subject_label <- triage_subject_label(row$sample_id[[1]])
      if (nzchar(subject_label)) {
        header_parts <- c(header_parts, subject_label)
      }
    }
    header_parts <- c(
      header_parts,
      sprintf("[%s]", display_scalar_value(row$severity[[1]], "info")),
      display_scalar_value(row$display_name[[1]], display_scalar_value(row$signal_id[[1]], "未命名信号"))
    )
    lines <- c(
      lines,
      sprintf("- %s", paste(header_parts, collapse = " ")),
      sprintf("  问题：%s", display_scalar_value(row$suspected_issue[[1]], "待补充")),
      sprintf("  证据：%s", display_scalar_value(row$evidence[[1]], "待补充")),
      sprintf("  建议：%s", display_scalar_value(row$recommended_action[[1]], "待补充"))
    )
  }
  lines
}

escape_markdown_cell <- function(x) {
  value <- display_scalar_value(x, "NA")
  value <- gsub("\\|", "\\\\|", value)
  value <- gsub("\n", "<br>", value, fixed = TRUE)
  value
}

render_markdown_table <- function(df) {
  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) {
    return("- 无")
  }

  header <- paste(sprintf(" %s ", colnames(df)), collapse = "|")
  separator <- paste(rep(" --- ", ncol(df)), collapse = "|")
  body <- apply(df, 1, function(row) {
    paste(sprintf(" %s ", vapply(row, escape_markdown_cell, character(1))), collapse = "|")
  })

  c(
    paste0("|", header, "|"),
    paste0("|", separator, "|"),
    paste0("|", body, "|")
  )
}

strip_ensembl_version <- function(x) {
  sub("\\.[0-9]+$", "", x)
}

collapse_unique_values <- function(x, sep = ", ") {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return("")
  }
  paste(unique(x), collapse = sep)
}

read_identifier_list <- function(path) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path) || !file.exists(path) || isTRUE(file.info(path)$size == 0)) {
    return(character(0))
  }

  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) == 0) {
    return(character(0))
  }

  lines <- sub("^\ufeff", "", lines)
  trimmed <- trimws(lines)
  trimmed <- trimmed[nzchar(trimmed) & !startsWith(trimmed, "#")]
  if (length(trimmed) == 0) {
    return(character(0))
  }

  tokens <- unlist(strsplit(trimmed, "[,\t ]+", perl = TRUE), use.names = FALSE)
  tokens <- trimws(tokens)
  tokens <- tokens[nzchar(tokens)]
  tokens <- tokens[!tolower(tokens) %in% c("gene_id", "gene_name", "feature", "feature_name")]
  unique(tokens)
}

resolve_metrics_summary_path <- function(cfg, sample_id, inventory_df = NULL) {
  inventory_row <- resolve_inventory_row(cfg, sample_id, inventory_df = inventory_df)
  candidates <- character(0)

  if (!is.null(inventory_row) && "metrics_path" %in% colnames(inventory_row)) {
    metrics_path <- normalize_scalar_value(inventory_row$metrics_path[1])
    if (nzchar(metrics_path) && file.exists(metrics_path)) {
      return(metrics_path)
    }
  }

  if (!is.null(inventory_row) && "cellranger_sample_dir" %in% colnames(inventory_row)) {
    sample_dir <- normalize_scalar_value(inventory_row$cellranger_sample_dir[1])
    if (nzchar(sample_dir)) {
      candidates <- c(
        candidates,
        file.path(sample_dir, "outs", "metrics_summary.csv"),
        file.path(sample_dir, "metrics_summary.csv"),
        file.path(sample_dir, "metrics_summary.xls")
      )
    }
  }

  if (nzchar(cfg$cellranger_out_dir)) {
    candidates <- c(
      candidates,
      file.path(cfg$cellranger_out_dir, sample_id, "outs", "metrics_summary.csv"),
      file.path(cfg$cellranger_out_dir, sample_id, "metrics_summary.xls")
    )
  }

  first_existing_path(candidates)
}

canonicalize_metric_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

extract_metric_value <- function(metrics_df, candidates) {
  if (is.null(metrics_df) || nrow(metrics_df) == 0) {
    return(NA_character_)
  }
  metric_map <- setNames(colnames(metrics_df), canonicalize_metric_name(colnames(metrics_df)))
  for (candidate in candidates) {
    key <- canonicalize_metric_name(candidate)
    if (key %in% names(metric_map)) {
      value <- metrics_df[[metric_map[[key]]]][1]
      if (length(value) > 0 && !is.na(value)) {
        return(as.character(value))
      }
    }
  }
  NA_character_
}

as_numeric_metric <- function(value) {
  if (length(value) == 0 || is.na(value) || !nzchar(value)) {
    return(NA_real_)
  }
  as.numeric(gsub("[^0-9.]+", "", value))
}

read_metrics_table <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    return(NULL)
  }

  readers <- list(
    function() read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    function() read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    function() utils::read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  )

  df <- NULL
  for (reader in readers) {
    df <- tryCatch(reader(), error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) {
      break
    }
  }
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  df
}

read_cellranger_metrics <- function(path) {
  df <- read_metrics_table(path)
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  list(
    estimated_cells = as_numeric_metric(extract_metric_value(df, c("Estimated Number of Cells"))),
    mean_reads_per_cell = as_numeric_metric(extract_metric_value(df, c("Mean Reads per Cell"))),
    median_genes_per_cell = as_numeric_metric(extract_metric_value(df, c("Median Genes per Cell", "Median Genes per Cell Associated with Cell"))),
    sequencing_saturation = as_numeric_metric(extract_metric_value(df, c("Sequencing Saturation"))),
    fraction_reads_in_cells = as_numeric_metric(extract_metric_value(df, c("Fraction Reads in Cells", "Reads Mapped Confidently to Transcriptome")))
  )
}

extract_gtf_attr <- function(attr_vec, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  matches <- regexec(pattern, attr_vec, perl = TRUE)
  values <- regmatches(attr_vec, matches)
  vapply(values, function(hit) {
    if (length(hit) >= 2) {
      hit[2]
    } else {
      ""
    }
  }, character(1))
}

read_reference_annotation <- function(cfg) {
  candidates <- c(cfg$clean_gtf, cfg$reference_gtf)
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  cache_key <- paste0("annotation::", candidates[1])
  if (exists(cache_key, envir = .pipeline_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .pipeline_cache, inherits = FALSE))
  }

  gtf_path <- candidates[1]
  con <- if (grepl("\\.gz$", gtf_path)) gzfile(gtf_path, open = "rt") else file(gtf_path, open = "rt")
  on.exit(close(con), add = TRUE)
  gtf_df <- tryCatch(
    read.delim(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char = "#", quote = ""),
    error = function(e) NULL
  )
  if (is.null(gtf_df) || nrow(gtf_df) == 0 || ncol(gtf_df) < 9) {
    empty_df <- data.frame(stringsAsFactors = FALSE)
    assign(cache_key, empty_df, envir = .pipeline_cache)
    return(empty_df)
  }

  gene_rows <- gtf_df[gtf_df[[3]] == "gene", c(1, 9), drop = FALSE]
  if (nrow(gene_rows) == 0) {
    empty_df <- data.frame(stringsAsFactors = FALSE)
    assign(cache_key, empty_df, envir = .pipeline_cache)
    return(empty_df)
  }

  attr_field <- gene_rows[[2]]
  gene_id <- extract_gtf_attr(attr_field, "gene_id")
  gene_name <- extract_gtf_attr(attr_field, "gene_name")
  gene_name[!nzchar(gene_name)] <- extract_gtf_attr(attr_field, "Name")[!nzchar(gene_name)]
  gene_biotype <- extract_gtf_attr(attr_field, "gene_biotype")
  gene_biotype[!nzchar(gene_biotype)] <- extract_gtf_attr(attr_field, "gene_type")[!nzchar(gene_biotype)]
  gene_biotype[!nzchar(gene_biotype)] <- extract_gtf_attr(attr_field, "biotype")[!nzchar(gene_biotype)]

  annotation_df <- unique(data.frame(
    seqname = as.character(gene_rows[[1]]),
    gene_id = gene_id,
    gene_id_stripped = strip_ensembl_version(gene_id),
    gene_name = gene_name,
    gene_name_upper = toupper(gene_name),
    gene_biotype = gene_biotype,
    stringsAsFactors = FALSE
  ))
  annotation_df <- annotation_df[nzchar(annotation_df$gene_id) | nzchar(annotation_df$gene_name), , drop = FALSE]
  assign(cache_key, annotation_df, envir = .pipeline_cache)
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
  unique(gene_names[keep])
}

canonical_vertebrate_mito_genes <- function() {
  c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB")
}

normalize_mito_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^MT-", "", x)
  x
}

normalize_feature_token <- function(x) {
  x <- toupper(strip_ensembl_version(trimws(as.character(x))))
  gsub("[-_]+", "-", x)
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

  if (grepl("ensgalg|grcg|galgal|gallus|chicken", joined)) {
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

resolve_feature_context <- function(gene_names, cfg, declared_gene_id_type = "auto") {
  annotation_df <- read_reference_annotation(cfg)
  feature_name_profile <- guess_feature_name_profile(gene_names)
  gene_id_type <- tolower(normalize_scalar_value(declared_gene_id_type, "auto"))
  if (!gene_id_type %in% c("symbol", "ensembl", "mixed", "unknown")) {
    gene_id_type <- feature_name_profile
  }
  if (!gene_id_type %in% c("symbol", "ensembl", "mixed", "unknown")) {
    gene_id_type <- "unknown"
  }

  gene_names_upper <- toupper(as.character(gene_names))
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

  cc_genes <- tryCatch(Seurat::cc.genes.updated.2019, error = function(e) NULL)
  s_features <- character(0)
  g2m_features <- character(0)
  if (!is.null(cc_genes)) {
    target_s <- unique(toupper(cc_genes$s.genes))
    target_g2m <- unique(toupper(cc_genes$g2m.genes))
    if (nrow(annotation_df) > 0) {
      s_rows <- annotation_df[annotation_df$gene_name_upper %in% target_s, , drop = FALSE]
      g2m_rows <- annotation_df[annotation_df$gene_name_upper %in% target_g2m, , drop = FALSE]
      s_features <- match_features_from_annotation(gene_names, annotation_df, s_rows)
      g2m_features <- match_features_from_annotation(gene_names, annotation_df, g2m_rows)
    } else {
      s_features <- gene_names[gene_names_upper %in% target_s]
      g2m_features <- gene_names[gene_names_upper %in% target_g2m]
    }
  }

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
    cell_cycle_s_features = unique(s_features),
    cell_cycle_g2m_features = unique(g2m_features),
    cell_cycle_s_feature_count = length(unique(s_features)),
    cell_cycle_g2m_feature_count = length(unique(g2m_features)),
    annotation_available = nrow(annotation_df) > 0
  )
}

add_basic_qc_metrics <- function(seu, cfg, declared_gene_id_type = "auto") {
  feature_context <- resolve_feature_context(rownames(seu), cfg, declared_gene_id_type = declared_gene_id_type)

  seu$percent.mito <- if (length(feature_context$mito_features) > 0) {
    PercentageFeatureSet(seu, features = feature_context$mito_features)
  } else {
    0
  }
  seu$percent.ribo <- if (length(feature_context$ribo_features) > 0) {
    PercentageFeatureSet(seu, features = feature_context$ribo_features)
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
    cell_cycle_s_feature_count = feature_context$cell_cycle_s_feature_count,
    cell_cycle_g2m_feature_count = feature_context$cell_cycle_g2m_feature_count,
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

get_assay_matrix <- function(seu, assay = "RNA", type = c("counts", "data")) {
  type <- match.arg(type)
  tryCatch(
    GetAssayData(seu, assay = assay, layer = type),
    error = function(e) GetAssayData(seu, assay = assay, slot = type)
  )
}

pick_expression_matrix <- function(x) {
  if (!is.list(x)) {
    return(x)
  }
  if ("Gene Expression" %in% names(x)) {
    return(x[["Gene Expression"]])
  }
  x[[1]]
}

read_matrix_from_path <- function(path) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path) || !file.exists(path)) {
    stop(sprintf("矩阵路径不存在: %s", path), call. = FALSE)
  }
  if (dir.exists(path)) {
    return(pick_expression_matrix(Read10X(data.dir = path)))
  }
  if (grepl("\\.h5$", path, ignore.case = TRUE)) {
    return(pick_expression_matrix(Read10X_h5(path)))
  }
  stop(sprintf("暂不支持的矩阵路径类型: %s", path), call. = FALSE)
}

resolve_sample_split_var <- function(seu) {
  if ("sample_id" %in% colnames(seu@meta.data)) {
    return("sample_id")
  }
  "orig.ident"
}

split_object_by_sample <- function(seu) {
  split.by <- resolve_sample_split_var(seu)
  split_objs <- SplitObject(seu, split.by = split.by)
  split_objs[order(names(split_objs))]
}

build_sample_branch_model <- function(sample_obj,
                                      dims,
                                      resolution,
                                      hvg_nfeatures,
                                      min_cells = 50L) {
  sample_obj <- maybe_join_layers(sample_obj)
  if (ncol(sample_obj) == 0) {
    return(list(
      status = "skipped_empty",
      object = sample_obj,
      usable_dims = integer(0),
      cluster_col = "provisional_cluster"
    ))
  }

  sample_obj$provisional_cluster <- "cluster_1"
  sample_obj$provisional_umap_1 <- NA_real_
  sample_obj$provisional_umap_2 <- NA_real_

  if (ncol(sample_obj) < min_cells) {
    return(list(
      status = "skipped_low_cells",
      object = sample_obj,
      usable_dims = integer(0),
      cluster_col = "provisional_cluster"
    ))
  }

  model_obj <- sample_obj
  model_obj <- NormalizeData(model_obj, verbose = FALSE)
  model_obj <- FindVariableFeatures(
    model_obj,
    selection.method = "vst",
    nfeatures = hvg_nfeatures,
    verbose = FALSE
  )
  model_obj <- ScaleData(model_obj, verbose = FALSE)
  model_obj <- RunPCA(model_obj, verbose = FALSE)

  usable_dims <- dims[dims <= ncol(Embeddings(model_obj, "pca"))]
  usable_dims <- usable_dims[is.finite(usable_dims)]
  if (length(usable_dims) == 0) {
    return(list(
      status = "skipped_no_pcs",
      object = sample_obj,
      usable_dims = integer(0),
      cluster_col = "provisional_cluster"
    ))
  }

  model_obj <- FindNeighbors(model_obj, dims = usable_dims, verbose = FALSE)
  model_obj <- FindClusters(model_obj, resolution = resolution, verbose = FALSE)
  model_obj$provisional_cluster <- as.character(Idents(model_obj))

  if (length(usable_dims) >= 2) {
    model_obj <- RunUMAP(model_obj, dims = usable_dims, verbose = FALSE)
    umap_mat <- Embeddings(model_obj, reduction = "umap")
    if (ncol(umap_mat) >= 2) {
      model_obj$provisional_umap_1 <- umap_mat[, 1]
      model_obj$provisional_umap_2 <- umap_mat[, 2]
    }
  }

  list(
    status = "built",
    object = model_obj,
    usable_dims = usable_dims,
    cluster_col = "provisional_cluster"
  )
}

estimate_expected_doublet_rate <- function(n_cells, rate_per_1k) {
  n_cells <- suppressWarnings(as.numeric(n_cells))
  rate_per_1k <- suppressWarnings(as.numeric(rate_per_1k))
  if (!is.finite(n_cells) || n_cells <= 0 || !is.finite(rate_per_1k) || rate_per_1k <= 0) {
    return(NA_real_)
  }
  rate <- rate_per_1k * (n_cells / 1000)
  rate <- min(max(rate, 0.001), 0.30)
  rate
}

estimate_expected_doublet_count <- function(n_cells, rate_per_1k) {
  rate <- estimate_expected_doublet_rate(n_cells, rate_per_1k)
  if (!is.finite(rate)) {
    return(NA_integer_)
  }
  as.integer(max(1, round(rate * n_cells)))
}

choose_ambient_method <- function(readiness_row, cfg) {
  if (is.null(readiness_row) || nrow(readiness_row) == 0) {
    return(list(
      preferred_method = "none",
      fallback_method = "none",
      apply_policy = cfg$ambient_apply_policy,
      soupx_ready = FALSE,
      decontx_ready = FALSE,
      cellbender_ready = FALSE
    ))
  }

  as_bool <- function(col_name) {
    if (!col_name %in% colnames(readiness_row)) {
      return(FALSE)
    }
    identical(tolower(normalize_scalar_value(readiness_row[[col_name]][1], "false")), "true")
  }

  preferred_method <- normalize_scalar_value(readiness_row$ambient_preferred_method[1], "none")
  fallback_method <- normalize_scalar_value(readiness_row$ambient_fallback_method[1], "none")
  list(
    preferred_method = preferred_method,
    fallback_method = fallback_method,
    apply_policy = normalize_scalar_value(readiness_row$ambient_apply_default[1], cfg$ambient_apply_policy),
    soupx_ready = as_bool("ambient_soupx_ready"),
    decontx_ready = as_bool("ambient_decontx_ready"),
    cellbender_ready = as_bool("ambient_cellbender_ready"),
    raw_matrix_available = as_bool("raw_matrix_available"),
    raw_matrix_kind = normalize_scalar_value(readiness_row$raw_matrix_kind[1]),
    ambient_notes = normalize_scalar_value(readiness_row$ambient_notes[1])
  )
}

ambient_recommend_apply <- function(contamination_fraction, cfg) {
  contamination_fraction <- suppressWarnings(as.numeric(contamination_fraction))
  if (!is.finite(contamination_fraction)) {
    return(FALSE)
  }
  contamination_fraction >= cfg$ambient_recommend_min_contamination
}

merge_named_objects <- function(object_list) {
  nonempty <- object_list[vapply(object_list, function(obj) as.integer(ncol(obj)), integer(1)) > 0]
  if (length(nonempty) == 0) {
    stop("对象列表为空，无法合并。", call. = FALSE)
  }
  merged <- nonempty[[1]]
  if (length(nonempty) > 1) {
    merged <- merge(x = nonempty[[1]], y = nonempty[2:length(nonempty)])
  }
  maybe_join_layers(merged)
}

safe_find_all_markers <- function(seu, cluster_col = "provisional_cluster", top_n = 3L) {
  if (!cluster_col %in% colnames(seu@meta.data) || length(unique(seu[[cluster_col, drop = TRUE]])) < 2) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  marker_df <- tryCatch(
    FindAllMarkers(
      seu,
      only.pos = TRUE,
      group.by = cluster_col,
      assay = DefaultAssay(seu),
      verbose = FALSE
    ),
    error = function(e) data.frame(stringsAsFactors = FALSE)
  )
  if (nrow(marker_df) == 0) {
    return(marker_df)
  }
  marker_df <- dplyr::group_by(marker_df, cluster)
  marker_df <- dplyr::arrange(marker_df, dplyr::desc(avg_log2FC), .by_group = TRUE)
  marker_df <- dplyr::slice_head(marker_df, n = top_n)
  dplyr::ungroup(marker_df)
}

calculate_marker_leakage <- function(before_counts, after_counts, cluster_labels, marker_df) {
  if (nrow(marker_df) == 0 || length(cluster_labels) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  cluster_labels <- as.character(cluster_labels)
  out_rows <- list()
  for (i in seq_len(nrow(marker_df))) {
    row <- marker_df[i, , drop = FALSE]
    gene <- as.character(row$gene[1])
    cluster_id <- as.character(row$cluster[1])
    if (!gene %in% rownames(before_counts) || !gene %in% rownames(after_counts)) {
      next
    }
    in_cluster <- cluster_labels == cluster_id
    out_cluster <- !in_cluster
    before_gene <- as.numeric(before_counts[gene, , drop = TRUE])
    after_gene <- as.numeric(after_counts[gene, , drop = TRUE])
    before_in <- if (any(in_cluster)) mean(before_gene[in_cluster]) else NA_real_
    before_out <- if (any(out_cluster)) mean(before_gene[out_cluster]) else NA_real_
    after_in <- if (any(in_cluster)) mean(after_gene[in_cluster]) else NA_real_
    after_out <- if (any(out_cluster)) mean(after_gene[out_cluster]) else NA_real_
    out_rows[[length(out_rows) + 1]] <- data.frame(
      gene = gene,
      source_cluster = cluster_id,
      mean_in_source_before = before_in,
      mean_outside_before = before_out,
      leakage_before = ifelse(is.finite(before_in) && before_in > 0, before_out / before_in, NA_real_),
      mean_in_source_after = after_in,
      mean_outside_after = after_out,
      leakage_after = ifelse(is.finite(after_in) && after_in > 0, after_out / after_in, NA_real_),
      stringsAsFactors = FALSE
    )
  }
  if (length(out_rows) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  dplyr::bind_rows(out_rows)
}

build_cellbender_stub <- function(sample_id, sample_obj, inventory_row, cfg, raw_counts = NULL) {
  raw_path <- if (!is.null(inventory_row) && "raw_matrix_dir" %in% colnames(inventory_row)) {
    normalize_scalar_value(inventory_row$raw_matrix_dir[1])
  } else {
    ""
  }
  output_path <- file.path(cfg$checkpoint_dir, "ambient", sample_id, paste0(sample_id, "_cellbender_filtered.h5"))
  expected_cells <- ncol(sample_obj)
  total_droplets <- if (!is.null(raw_counts)) ncol(raw_counts) else NA_integer_
  args <- c(
    "cellbender", "remove-background",
    "--input", raw_path,
    "--output", output_path,
    "--expected-cells", as.character(expected_cells)
  )
  if (is.finite(total_droplets) && !is.na(total_droplets)) {
    args <- c(args, "--total-droplets-included", as.character(total_droplets))
  }
  if (isTRUE(tolower(cfg$cellbender_cuda) %in% c("yes", "true", "1", "on"))) {
    args <- c(args, "--cuda")
  }
  if (is.finite(cfg$cellbender_fpr) && !is.na(cfg$cellbender_fpr)) {
    args <- c(args, "--fpr", as.character(cfg$cellbender_fpr))
  }
  extra_args <- split_csv(cfg$cellbender_extra_args)
  if (length(extra_args) > 0) {
    args <- c(args, extra_args)
  }
  data.frame(
    sample_id = sample_id,
    cellbender_mode = cfg$cellbender_mode,
    raw_matrix_path = raw_path,
    expected_cells = expected_cells,
    total_droplets_included = total_droplets,
    fpr = cfg$cellbender_fpr,
    cuda_requested = cfg$cellbender_cuda,
    output_path = output_path,
    command_stub = paste(shQuote(args), collapse = " "),
    stringsAsFactors = FALSE
  )
}

detect_mito_features <- function(gene_names) {
  resolve_feature_context(gene_names, read_cfg())$mito_features
}

detect_ribo_features <- function(gene_names) {
  resolve_feature_context(gene_names, read_cfg())$ribo_features
}

readiness_note_labels <- function() {
  c(
    main_inputs_incomplete = "主流程输入不完整",
    velocity_inputs_incomplete = "velocity 上游输入不完整",
    ambient_inputs_incomplete = "ambient 分支上游输入不完整",
    scenic_inputs_incomplete = "SCENIC 上游输入不完整",
    replicates_insufficient = "replicate 信息不足",
    reference_version_mixed = "reference_version 混用",
    gene_id_type_unrecognized = "gene_id_type 未识别",
    missing_filtered_matrix = "filtered matrix 缺失",
    missing_raw_matrix = "raw matrix 缺失",
    ambient_branch_unavailable = "ambient 分支不可用",
    soupx_requires_raw_matrix = "SoupX 需要 raw droplets",
    decontx_fallback_only = "无 raw droplets，ambient 默认走 DecontX fallback",
    cellbender_stub_requires_10x_raw = "CellBender stub 目前仅标记为 10X raw 直接 ready",
    cellbender_stub_only = "CellBender 本轮仅提供 readiness/config stub",
    missing_bam = "BAM 缺失",
    missing_fastq = "FASTQ 缺失",
    missing_cellranger_out = "cellranger_out 缺失",
    replicate_fields_missing = "replicate/batch 字段缺失",
    velocity_requires_cellranger = "当前 velocity 仍要求 Cell Ranger 风格输入"
  )
}

describe_readiness_notes <- function(raw_notes) {
  labels <- readiness_note_labels()
  notes <- trimws(unlist(strsplit(normalize_scalar_value(raw_notes), ";", fixed = TRUE)))
  notes <- notes[nzchar(notes)]
  if (length(notes) == 0) {
    return(character(0))
  }
  vapply(notes, function(note) {
    if (note %in% names(labels)) labels[[note]] else note
  }, character(1))
}

celltype_marker_list <- function(marker_panel_dir = NULL) {
  if (is.null(marker_panel_dir) || !nzchar(marker_panel_dir)) {
    marker_panel_dir <- Sys.getenv("MARKER_PANEL_DIR", file.path(must_getenv("PROJECT_ROOT"), "config", "marker_panels"))
  }

  if (!dir.exists(marker_panel_dir)) {
    message("No marker panel found, skipping literature validation")
    return(list())
  }

  panel_files <- list.files(marker_panel_dir, pattern = "\\.tsv$", full.names = TRUE, ignore.case = TRUE)
  if (length(panel_files) == 0) {
    message("No marker panel found, skipping literature validation")
    return(list())
  }

  panel_rows <- lapply(panel_files, function(panel_path) {
    panel_df <- read_tsv_optional(panel_path)
    if (nrow(panel_df) == 0) {
      return(NULL)
    }
    if (!"layer_id" %in% colnames(panel_df)) {
      panel_df$layer_id <- "*"
    }
    if (!"celltype" %in% colnames(panel_df) && "cell_type" %in% colnames(panel_df)) {
      panel_df$celltype <- panel_df$cell_type
    }
    if (!"evidence_source" %in% colnames(panel_df)) {
      panel_df$evidence_source <- "unspecified"
    }
    required_cols <- c("celltype", "gene")
    if (!all(required_cols %in% colnames(panel_df))) {
      warning(sprintf("跳过 marker panel（缺少 layer_id/celltype/gene/evidence_source 所需字段）: %s", panel_path), call. = FALSE)
      return(NULL)
    }
    panel_df <- panel_df[, c("layer_id", "celltype", "gene", "evidence_source"), drop = FALSE]
    panel_df$celltype <- trimws(as.character(panel_df$celltype))
    panel_df$gene <- trimws(as.character(panel_df$gene))
    panel_df <- panel_df[nzchar(panel_df$celltype) & nzchar(panel_df$gene), , drop = FALSE]
    if (nrow(panel_df) == 0) {
      return(NULL)
    }
    panel_df
  })

  panel_rows <- Filter(Negate(is.null), panel_rows)
  if (length(panel_rows) == 0) {
    message("No marker panel found, skipping literature validation")
    return(list())
  }

  panel_df <- unique(do.call(rbind, panel_rows))
  split(panel_df$gene, panel_df$celltype)
}
