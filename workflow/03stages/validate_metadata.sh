#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")/../02lib" && pwd)/common.sh"

prepare_project_state_dirs

[[ -f "${SAMPLE_SHEET}" ]] || die "缺少样本表: ${SAMPLE_SHEET}"
[[ -f "${COMPARISON_SHEET}" ]] || die "缺少比较设计表: ${COMPARISON_SHEET}"

python_bin="$(detect_python)"

env \
  SAMPLE_SHEET="${SAMPLE_SHEET}" \
  CANONICAL_SAMPLE_SHEET="${CANONICAL_SAMPLE_SHEET}" \
  COMPARISON_SHEET="${COMPARISON_SHEET}" \
  SECTION_SHEET="${SECTION_SHEET}" \
  SPATIAL_REFERENCE_INVENTORY_FILE="${SPATIAL_REFERENCE_INVENTORY_FILE}" \
  WAIVER_FILE="${WAIVER_FILE}" \
  "${python_bin}" - <<'PY'
import csv
import os
import re
from pathlib import Path

sample_path = Path(os.environ["SAMPLE_SHEET"])
canonical_path = Path(os.environ["CANONICAL_SAMPLE_SHEET"])
comparison_path = Path(os.environ["COMPARISON_SHEET"])
section_path = Path(os.environ["SECTION_SHEET"])
spatial_reference_inventory_path = Path(os.environ["SPATIAL_REFERENCE_INVENTORY_FILE"])
waiver_path = Path(os.environ["WAIVER_FILE"])


def cell(row, key, default=""):
    value = row.get(key, default)
    if value is None:
        return default
    return str(value)


required_sample_cols = [
    "sample_id",
    "condition",
    "biological_replicate",
    "technical_replicate",
    "batch",
    "input_mode",
    "input_source",
    "source_path",
]
optional_sample_defaults = {
    "platform": "auto",
    "gene_id_type": "auto",
    "reference_version": "",
    "group_id": "",
    "timepoint": "",
    "tissue": "",
    "chemistry": "",
    "run_main": "yes",
    "run_velocity": "auto",
    "run_scenic": "yes",
    "modality": "scrna",
    "section_id": "",
    "chip_id": "",
    "bundle_layout": "",
    "image_path": "",
    "run_spatial": "auto",
    "run_deconv": "auto",
    "run_joint": "auto",
    "notes": "",
}
canonical_columns = required_sample_cols + list(optional_sample_defaults.keys())
allowed_input_modes = {"fastq", "cellranger_out", "matrix", "visium_bundle", "saw_bundle", "stomics_bundle", "spatial_matrix"}
allowed_run_values = {"yes", "no", "auto"}
allowed_modalities = {"scrna", "spatial"}
allowed_platforms = {"10x_cellranger", "dnbelab_c", "generic_mex", "visium", "saw", "stomics", "generic", "auto", ""}
allowed_gene_id_types = {"symbol", "ensembl", "mixed", "unknown", "auto", ""}


with sample_path.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"样本表缺少表头: {sample_path}")
    fieldnames = [field[2:] if field.startswith("# ") else field for field in reader.fieldnames]
    missing = [col for col in required_sample_cols if col not in fieldnames]
    if missing:
        raise SystemExit(f"样本表缺少必需列: {', '.join(missing)}")
    sample_rows = [
        {fieldnames[i]: cell(row, reader.fieldnames[i]) for i in range(len(fieldnames))}
        for row in reader
    ]

if not sample_rows:
    raise SystemExit("样本表为空，至少需要一行数据。")

seen_samples = set()
conditions = set()
canonical_rows = []
spatial_section_ids = set()
for row in sample_rows:
    canonical = {}
    for col in required_sample_cols:
        canonical[col] = cell(row, col).strip()

    reference_version = cell(row, "reference_version").strip()
    if not reference_version:
        reference_version = cell(row, "reference").strip()
    for col, default in optional_sample_defaults.items():
        if col == "reference_version":
            canonical[col] = reference_version or default
        else:
            canonical[col] = cell(row, col, default).strip() or default

    sample_id = canonical["sample_id"]
    if not sample_id:
        raise SystemExit("样本表存在空 sample_id。")
    if sample_id in seen_samples:
        raise SystemExit(f"sample_id 重复: {sample_id}")
    seen_samples.add(sample_id)

    if not canonical["condition"]:
        raise SystemExit(f"{sample_id} 缺少 condition。")
    conditions.add(canonical["condition"])

    if canonical["input_mode"] not in allowed_input_modes:
        raise SystemExit(
            f"{sample_id} 的 input_mode 非法: {canonical['input_mode']} "
            f"(允许值: {', '.join(sorted(allowed_input_modes))})"
        )
    canonical["platform"] = canonical["platform"].lower()
    if canonical["platform"] not in allowed_platforms:
        raise SystemExit(
            f"{sample_id} 的 platform 非法: {canonical['platform']} "
            f"(允许值: {', '.join(sorted(v for v in allowed_platforms if v))})"
        )
    canonical["gene_id_type"] = canonical["gene_id_type"].lower()
    if canonical["gene_id_type"] not in allowed_gene_id_types:
        raise SystemExit(
            f"{sample_id} 的 gene_id_type 非法: {canonical['gene_id_type']} "
            f"(允许值: {', '.join(sorted(v for v in allowed_gene_id_types if v))})"
        )
    if not canonical["source_path"]:
        raise SystemExit(f"{sample_id} 缺少 source_path。")

    for run_col in ("run_main", "run_velocity", "run_scenic"):
        canonical[run_col] = canonical[run_col].lower()
        if canonical[run_col] not in allowed_run_values:
            raise SystemExit(f"{sample_id} 的 {run_col} 非法: {canonical[run_col]}")

    canonical["modality"] = canonical["modality"].lower()
    if canonical["modality"] not in allowed_modalities:
        raise SystemExit(
            f"{sample_id} 的 modality 非法: {canonical['modality']} "
            f"(允许值: {', '.join(sorted(allowed_modalities))})"
        )
    for run_col in ("run_spatial", "run_deconv", "run_joint"):
        canonical[run_col] = canonical[run_col].lower()
        if canonical[run_col] not in allowed_run_values:
            raise SystemExit(f"{sample_id} 的 {run_col} 非法: {canonical[run_col]}")
    if canonical["modality"] == "scrna":
        if canonical["run_spatial"] == "auto":
            canonical["run_spatial"] = "no"
        if canonical["run_deconv"] == "auto":
            canonical["run_deconv"] = "no"
        if canonical["run_joint"] == "auto":
            canonical["run_joint"] = "no"
    else:
        canonical["run_main"] = "no"
        canonical["run_velocity"] = "no"
        canonical["run_scenic"] = "no"
        if not canonical["section_id"]:
            raise SystemExit(f"{sample_id} 是 spatial 样本但缺少 section_id。")
        if canonical["run_spatial"] == "auto":
            canonical["run_spatial"] = "yes"
        spatial_section_ids.add(canonical["section_id"])

    canonical_rows.append(canonical)

if spatial_section_ids:
    if not section_path.exists():
        raise SystemExit(f"存在 spatial 样本，但缺少切片登记表: {section_path}")
    with section_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"切片登记表缺少表头: {section_path}")
        section_fieldnames = [field[2:] if field.startswith("# ") else field for field in reader.fieldnames]
        required_section_cols = [
            "section_id", "chip_id", "tissue_block", "condition", "platform",
            "bundle_root", "image_lowres", "image_hires", "tissue_positions",
            "scalefactors", "enabled",
        ]
        missing = [col for col in required_section_cols if col not in section_fieldnames]
        if missing:
            raise SystemExit(f"切片登记表缺少必需列: {', '.join(missing)}")
        section_rows = [
            {section_fieldnames[i]: cell(row, reader.fieldnames[i]) for i in range(len(section_fieldnames))}
            for row in reader
        ]
    registered_sections = {cell(row, "section_id").strip() for row in section_rows}
    missing_sections = sorted(section for section in spatial_section_ids if section not in registered_sections)
    if missing_sections:
        raise SystemExit(f"samples.tsv.section_id 未登记到 sections.tsv: {', '.join(missing_sections)}")

    if not spatial_reference_inventory_path.exists():
        spatial_reference_inventory_path.parent.mkdir(parents=True, exist_ok=True)
        spatial_reference_inventory_path.write_text(
            "inspected_file\tobject_class\tkey_columns_or_items\tselected_reference\tselected_annotation_col\tpseudotime_evidence\tcommunication_evidence\tfrozen_reference_id\treference_sha256\treference_freeze_timestamp\n"
            "results/checkpoints/03_after_annotation.rds\tSeurat\tcell_type,cell_subtype\tresults/checkpoints/03_after_annotation.rds\tcell_subtype\tnone\tnone\t\t\t\n",
            encoding="utf-8",
        )

canonical_path.parent.mkdir(parents=True, exist_ok=True)
with canonical_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=canonical_columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(canonical_rows)

required_comparison_cols = ["comparison_id", "ident_1", "ident_2", "enabled"]
safe_subset_column_re = re.compile(r"^[A-Za-z0-9_.:-]+$")
with comparison_path.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"比较设计表缺少表头: {comparison_path}")
    missing = [col for col in required_comparison_cols if col not in reader.fieldnames]
    if missing:
        raise SystemExit(f"比较设计表缺少必需列: {', '.join(missing)}")
    comparison_rows = list(reader)

for row in comparison_rows:
    comparison_id = cell(row, "comparison_id").strip()
    ident_1 = cell(row, "ident_1").strip()
    ident_2 = cell(row, "ident_2").strip()
    enabled = cell(row, "enabled").strip().lower()
    subset_column = cell(row, "subset_column").strip()
    subset_value = cell(row, "subset_value").strip()
    force_exploratory = cell(row, "force_exploratory").strip().lower()
    min_cells_per_group = cell(row, "min_cells_per_group").strip()
    logfc_threshold = cell(row, "logfc_threshold").strip()
    analysis_modality = cell(row, "analysis_modality", "scrna").strip().lower() or "scrna"
    group_var = cell(row, "group_var").strip()
    contrast_axis = cell(row, "contrast_axis").strip()

    if not comparison_id:
        raise SystemExit("比较设计表存在空 comparison_id。")
    if analysis_modality not in {"scrna", "spatial", "joint"}:
        raise SystemExit(f"{comparison_id} 的 analysis_modality 非法: {analysis_modality}")
    condition_like_comparison = group_var == "condition" or "condition" in contrast_axis
    if condition_like_comparison and ident_1 not in conditions:
        raise SystemExit(f"{comparison_id} 引用了不存在的 ident_1: {ident_1}")
    if condition_like_comparison and ident_2 not in conditions:
        raise SystemExit(f"{comparison_id} 引用了不存在的 ident_2: {ident_2}")
    if enabled not in {"yes", "no", "true", "false"}:
        raise SystemExit(f"{comparison_id} 的 enabled 非法: {enabled}")
    if bool(subset_column) != bool(subset_value):
        raise SystemExit(f"{comparison_id} 的 subset_column/subset_value 必须同时填写或同时留空。")
    if subset_column and not safe_subset_column_re.match(subset_column):
        raise SystemExit(
            f"{comparison_id} 的 subset_column 含非法字符: {subset_column} "
            "(仅允许字母、数字、下划线、点、冒号和短横线)"
        )
    if force_exploratory and force_exploratory not in {"yes", "no", "true", "false", "1", "0", "on", "off"}:
        raise SystemExit(f"{comparison_id} 的 force_exploratory 非法: {force_exploratory}")
    if min_cells_per_group:
        try:
            min_cells_value = int(min_cells_per_group)
        except ValueError:
            raise SystemExit(f"{comparison_id} 的 min_cells_per_group 必须是整数。")
        if min_cells_value < 3:
            raise SystemExit(f"{comparison_id} 的 min_cells_per_group 必须 >= 3（Wilcoxon 在 < 3 时 p 值不可信）。")
    if logfc_threshold:
        try:
            float(logfc_threshold)
        except ValueError:
            raise SystemExit(f"{comparison_id} 的 logfc_threshold 必须是数字。")

if not waiver_path.exists():
    waiver_path.parent.mkdir(parents=True, exist_ok=True)
    waiver_path.write_text("check_id\tscope\treason\tapproved_by\n", encoding="utf-8")

print(f"validated_samples={len(canonical_rows)}")
print(f"validated_conditions={len(conditions)}")
print(f"validated_comparisons={len(comparison_rows)}")
print(f"validated_spatial_sections={len(spatial_section_ids)}")
PY

update_workflow_status \
  "metadata_validated" \
  "bash ${PIPELINE_ROOT}/workflow/03stages/audit_inputs.sh" \
  "status.metadata_valid=true"

echo "metadata 校验完成"
echo "canonical sample sheet: ${CANONICAL_SAMPLE_SHEET}"
