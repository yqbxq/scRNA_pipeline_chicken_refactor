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
    suspected_issue = "该样本在基础 QC 后细胞保留率偏低，可能是样本质量、文库复杂度或阈值设置问题。",
    recommended_action = "回看 pre-QC 分布和样本级阈值；若分布连续，优先人工复核后再决定是否放宽阈值。"
  ),
  high_primary_doublet_rate = list(
    display_name = "Primary doublet 检出率偏高",
    suspected_issue = "scDblFinder 检出的 doublet rate 高于预期，可能存在上样浓度过高或局部混合群。",
    recommended_action = "重点审阅 doublet score 分布、UMAP 空间位置和高风险 cluster，再决定是否调整过滤策略。"
  ),
  primary_secondary_discordance = list(
    display_name = "Primary/secondary doublet 调用不一致",
    suspected_issue = "scDblFinder 与 DoubletFinder 对部分细胞给出不同判断，说明 doublet 边界存在不确定性。",
    recommended_action = "保留 scDblFinder 作为默认主调用器，同时人工复核冲突区域和 score 分布。"
  ),
  sample_collapse = list(
    display_name = "QC/doublet 后样本细胞数过少",
    suspected_issue = "该样本过滤后剩余细胞太少，后续聚类、整合或差异分析稳定性较差。",
    recommended_action = "确认该样本是否继续纳入主线；必要时记录为低细胞数风险或从后续分析中排除。"
  ),
  ambient_recommended_not_applied = list(
    display_name = "Ambient 建议替换但未应用",
    suspected_issue = "ambient 分支估计污染比例达到建议阈值，但主对象仍未替换为校正 counts。",
    recommended_action = "进入整合前明确是否接受原始 counts；如不接受，调整 AMBIENT_APPLY_POLICY 后重跑 02。"
  ),
  doublet_high_risk_cluster = list(
    display_name = "Cluster doublet 风险偏高",
    suspected_issue = "某个 provisional cluster 中 primary doublet fraction 偏高，可能代表混合细胞或边界区域。",
    recommended_action = "结合 marker、UMAP 和 doublet score 复核该 cluster，必要时在后续整合前剔除或标注。"
  ),
  integration_overcorrection_suspected = list(
    display_name = "整合疑似过度校正",
    suspected_issue = "整合后样本/批次效应被显著抑制，但 group/condition 间的真实生物差异也可能被模糊。",
    recommended_action = "在 03b 报告中对比未整合 baseline 的 sample/group 关联指标与 marker 表现，必要时选用更弱的整合策略。"
  ),
  integration_undercorrection_suspected = list(
    display_name = "整合疑似不足",
    suspected_issue = "整合后样本仍呈样本驱动结构，未达到合并分析的预期。",
    recommended_action = "检查 group_var 和样本构成，必要时调整 integration_modes 或增加更强的整合候选后重跑 03a2/03b。"
  ),
  integration_method_unavailable = list(
    display_name = "整合方法环境缺失",
    suspected_issue = "请求的整合方法在当前 R/Python 环境不可用或运行失败，模块已记录降级原因。",
    recommended_action = "在 03a2 manifest 和候选索引中确认降级原因；补齐环境后重跑，或从 integration_modes 中移除该方法。"
  ),
  resolution_fallback_used = list(
    display_name = "Resolution 回退到默认值",
    suspected_issue = "粗+细 resolution 搜索没有精确命中 target_clusters，03c 已选择最接近目标的结果。",
    recommended_action = "在 03e 报告中确认 cluster 数是否合理；若不合理，调整 target_clusters 或 res_range 后重跑 03c。"
  ),
  annotation_low_confidence = list(
    display_name = "注释置信度偏低",
    suspected_issue = "确定级别 cluster 占比 < 50%，多数 cluster 仍为暂定或未定。",
    recommended_action = "补充 marker panel 或检查 module score 阈值；必要时回到 03c 调整 resolution。"
  ),
  annotation_panel_conflict = list(
    display_name = "Panel 之间存在冲突",
    suspected_issue = "至少 1 个 cluster 同时显著命中两个候选 panel，注释证据不唯一。",
    recommended_action = "审阅 annotation_evidence.tsv，确认冲突 panel 是否需要合并/拆分，或 cluster 是否需要进一步 subcluster。"
  ),
  subcluster_skipped_no_filter = list(
    display_name = "子层缺少过滤条件被跳过",
    suspected_issue = "object_layers.tsv 中该 subcluster 行未提供 sample_include/sample_exclude，也未提供 selection_column/selection_values。",
    recommended_action = "根据 panorama 注释结果回填 selection_column+selection_values，或用 sample_include 显式列出样本，再重跑 04a。"
  ),
  subcluster_candidate_capped = list(
    display_name = "子层整合候选被截断",
    suspected_issue = "该 subcluster 的 normalization_methods × integration_mode 候选数超过上限，04a 已按出现顺序截断。",
    recommended_action = "如需评估被截断的方法，请提高 MAX_INTEGRATION_CANDIDATES_PER_LAYER 或缩小该层候选列表后重跑。"
  ),
  subcluster_low_confidence = list(
    display_name = "子层注释置信度偏低",
    suspected_issue = "确定级别 subcluster 占比 < 50%，多数 cluster 仍为暂定或未定。",
    recommended_action = "补充该层 marker panel，或回到 04a 调整 target_clusters/res_range 后重跑。"
  ),
  subcluster_resolution_fallback_used = list(
    display_name = "子层 resolution 回退",
    suspected_issue = "子层粗+细 resolution 搜索没有精确命中 target_clusters，04a 已选择最接近目标的结果。",
    recommended_action = "在 04c 报告中确认 cluster 数是否合理；若不合理，调整该层 target_clusters 或 res_range 后重跑 04a。"
  )
)

resolve_triage_signal_spec <- function(signal_id) {
  signal_id <- normalize_scalar_value(signal_id)
  if (!nzchar(signal_id)) {
    return(list(display_name = "未命名信号", suspected_issue = "", recommended_action = ""))
  }
  if (signal_id %in% names(TRIAGE_SIGNAL_REGISTRY)) {
    return(TRIAGE_SIGNAL_REGISTRY[[signal_id]])
  }
  list(display_name = signal_id, suspected_issue = "", recommended_action = "")
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
  do.call(data.frame, c(cols, list(stringsAsFactors = FALSE, check.names = FALSE)))
}

make_triage_row <- function(sample_id = NULL,
                            severity,
                            signal_id,
                            evidence = "",
                            suspected_issue = "",
                            recommended_action = "",
                            manual_review_required = "yes") {
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
    row <- cbind(
      data.frame(sample_id = display_scalar_value(sample_id, "NA"), stringsAsFactors = FALSE),
      row
    )
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
    velocity_requires_bam = "velocity 需要 BAM 输入"
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

sample_contract_rows_pre_qc <- function(row) {
  data.frame(
    主题 = c(
      "输入契约", "输入契约", "输入契约", "输入契约",
      "标准路径", "标准路径", "标准路径", "标准路径",
      "Ambient", "Ambient", "Ambient", "Ambient", "Ambient", "Ambient",
      "Feature/QC", "Feature/QC", "Feature/QC", "Feature/QC", "Feature/QC",
      "Feature/QC", "Feature/QC", "Feature/QC", "Feature/QC", "Feature/QC"
    ),
    字段 = c(
      "platform", "feature_name_profile", "gene_id_type", "reference_version",
      "filtered_matrix_dir", "raw_matrix_dir", "metrics_path", "bam_path",
      "preferred_method", "fallback_method", "soupx_ready", "decontx_ready", "cellbender_ready", "ambient_notes",
      "mito_detection_method", "mito_detection_detail", "species_guess", "mito_reference_seqnames", "mito_feature_count",
      "ribo_feature_count", "cell_cycle_s_features", "cell_cycle_g2m_features", "mito_gene_list", "mito_warnings"
    ),
    值 = c(
      display_scalar_value(row$platform),
      display_scalar_value(row$feature_name_profile),
      display_scalar_value(row$gene_id_type_resolved),
      display_scalar_value(row$reference_version),
      display_scalar_value(row$filtered_matrix_dir),
      display_scalar_value(row$raw_matrix_dir),
      display_scalar_value(row$metrics_summary_path),
      display_scalar_value(row$bam_path),
      display_scalar_value(row$ambient_preferred_method),
      display_scalar_value(row$ambient_fallback_method),
      display_scalar_value(row$ambient_soupx_ready),
      display_scalar_value(row$ambient_decontx_ready),
      display_scalar_value(row$ambient_cellbender_ready),
      display_scalar_value(row$ambient_notes),
      display_scalar_value(row$mito_detection_method, "failed"),
      display_scalar_value(row$mito_detection_detail),
      display_scalar_value(row$species_guess, "unknown"),
      display_scalar_value(row$mito_reference_seqnames),
      display_scalar_value(row$mito_feature_count),
      display_scalar_value(row$ribo_feature_count),
      display_scalar_value(row$cell_cycle_s_feature_count),
      display_scalar_value(row$cell_cycle_g2m_feature_count),
      display_scalar_value(row$mito_detected_gene_names),
      display_scalar_value(row$mito_warning_messages)
    ),
    stringsAsFactors = FALSE
  )
}
