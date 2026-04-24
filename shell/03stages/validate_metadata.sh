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
  WAIVER_FILE="${WAIVER_FILE}" \
  "${python_bin}" - <<'PY'
import csv
import os
from pathlib import Path

sample_path = Path(os.environ["SAMPLE_SHEET"])
canonical_path = Path(os.environ["CANONICAL_SAMPLE_SHEET"])
comparison_path = Path(os.environ["COMPARISON_SHEET"])
waiver_path = Path(os.environ["WAIVER_FILE"])

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
    "timepoint": "",
    "tissue": "",
    "chemistry": "",
    "run_main": "yes",
    "run_velocity": "auto",
    "run_scenic": "yes",
}
canonical_columns = required_sample_cols + list(optional_sample_defaults.keys())
allowed_input_modes = {"fastq", "cellranger_out", "matrix"}
allowed_run_values = {"yes", "no", "auto"}
allowed_platforms = {"10x_cellranger", "dnbelab_c", "generic_mex", "auto", ""}
allowed_gene_id_types = {"symbol", "ensembl", "mixed", "unknown", "auto", ""}

with sample_path.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"样本表缺少表头: {sample_path}")
    missing = [col for col in required_sample_cols if col not in reader.fieldnames]
    if missing:
        raise SystemExit(f"样本表缺少必需列: {', '.join(missing)}")
    sample_rows = list(reader)

if not sample_rows:
    raise SystemExit("样本表为空，至少需要一行数据。")

seen_samples = set()
conditions = set()
canonical_rows = []
for row in sample_rows:
    canonical = {}
    for col in required_sample_cols:
        canonical[col] = row.get(col, "").strip()

    reference_version = row.get("reference_version", "").strip()
    if not reference_version:
        reference_version = row.get("reference", "").strip()
    for col, default in optional_sample_defaults.items():
        if col == "reference_version":
            canonical[col] = reference_version or default
        else:
            canonical[col] = row.get(col, default).strip() or default

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

    canonical_rows.append(canonical)

canonical_path.parent.mkdir(parents=True, exist_ok=True)
with canonical_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=canonical_columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(canonical_rows)

required_comparison_cols = ["comparison_id", "ident_1", "ident_2", "enabled"]
with comparison_path.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"比较设计表缺少表头: {comparison_path}")
    missing = [col for col in required_comparison_cols if col not in reader.fieldnames]
    if missing:
        raise SystemExit(f"比较设计表缺少必需列: {', '.join(missing)}")
    comparison_rows = list(reader)

for row in comparison_rows:
    comparison_id = row.get("comparison_id", "").strip()
    ident_1 = row.get("ident_1", "").strip()
    ident_2 = row.get("ident_2", "").strip()
    enabled = row.get("enabled", "").strip().lower()

    if not comparison_id:
        raise SystemExit("比较设计表存在空 comparison_id。")
    if ident_1 not in conditions:
        raise SystemExit(f"{comparison_id} 引用了不存在的 ident_1: {ident_1}")
    if ident_2 not in conditions:
        raise SystemExit(f"{comparison_id} 引用了不存在的 ident_2: {ident_2}")
    if enabled not in {"yes", "no", "true", "false"}:
        raise SystemExit(f"{comparison_id} 的 enabled 非法: {enabled}")

if not waiver_path.exists():
    waiver_path.parent.mkdir(parents=True, exist_ok=True)
    waiver_path.write_text("check_id\tscope\treason\tapproved_by\n", encoding="utf-8")

print(f"validated_samples={len(canonical_rows)}")
print(f"validated_conditions={len(conditions)}")
print(f"validated_comparisons={len(comparison_rows)}")
PY

update_workflow_status \
  "metadata_validated" \
  "bash ${PIPELINE_ROOT}/shell/03stages/audit_inputs.sh" \
  "status.metadata_valid=true"

echo "metadata 校验完成"
echo "canonical sample sheet: ${CANONICAL_SAMPLE_SHEET}"
