#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

prepare_project_state_dirs

INPUT_SHEET="${CANONICAL_SAMPLE_SHEET}"
if [[ ! -f "${INPUT_SHEET}" ]]; then
  INPUT_SHEET="${SAMPLE_SHEET}"
fi
[[ -f "${INPUT_SHEET}" ]] || die "缺少样本表，请先运行 02_validate_metadata.sh"
[[ -f "${COMPARISON_SHEET}" ]] || die "缺少比较设计表: ${COMPARISON_SHEET}"

python_bin="$(detect_python)"

env \
  INPUT_SHEET="${INPUT_SHEET}" \
  COMPARISON_SHEET="${COMPARISON_SHEET}" \
  FASTQ_DIR="${FASTQ_DIR}" \
  CELLRANGER_OUT_DIR="${CELLRANGER_OUT_DIR}" \
  DATA_DIR="${DATA_DIR}" \
  INPUT_INVENTORY_FILE="${INPUT_INVENTORY_FILE}" \
  BRANCH_READINESS_FILE="${BRANCH_READINESS_FILE}" \
  PIPELINE_ROOT="${PIPELINE_ROOT}" \
  "${python_bin}" - <<'PY'
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ["PIPELINE_ROOT"]) / "workflow" / "python"))

from intake_contract import normalize_flag, normalize_value, resolve_sample_contract

input_sheet = Path(os.environ["INPUT_SHEET"])
comparison_sheet = Path(os.environ["COMPARISON_SHEET"])
fastq_dir = Path(os.environ["FASTQ_DIR"])
cellranger_out_dir = Path(os.environ["CELLRANGER_OUT_DIR"])
data_dir = Path(os.environ["DATA_DIR"])
inventory_path = Path(os.environ["INPUT_INVENTORY_FILE"])
readiness_path = Path(os.environ["BRANCH_READINESS_FILE"])

with input_sheet.open("r", encoding="utf-8", newline="") as handle:
    sample_rows = list(csv.DictReader(handle, delimiter="\t"))

with comparison_sheet.open("r", encoding="utf-8", newline="") as handle:
    comparison_rows = list(csv.DictReader(handle, delimiter="\t"))

conditions_to_reps: dict[str, set[str]] = {}
replicate_fields_complete = True
project_reference_versions: set[str] = set()
for row in sample_rows:
    if normalize_flag(row.get("run_main", "yes")) == "no":
        continue
    condition = normalize_value(row.get("condition"))
    biological_replicate = normalize_value(row.get("biological_replicate"))
    technical_replicate = normalize_value(row.get("technical_replicate"))
    batch = normalize_value(row.get("batch"))
    reference_version = normalize_value(row.get("reference_version"))
    if biological_replicate:
        conditions_to_reps.setdefault(condition, set()).add(biological_replicate)
    if not (biological_replicate and technical_replicate and batch):
        replicate_fields_complete = False
    if reference_version:
        project_reference_versions.add(reference_version)

inventory_rows = []
readiness_rows = []

main_project_ready = True
velocity_project_ready = True
ambient_project_ready = True
scenic_project_ready = True

main_requested = False
velocity_requested = False
scenic_requested = False
any_unknown_gene_id_type = False
any_raw_matrix_available = False

for row in sample_rows:
    sample_id = normalize_value(row.get("sample_id"))
    run_main = normalize_flag(row.get("run_main", "yes")) == "yes"
    run_velocity = normalize_flag(row.get("run_velocity", "auto")) in {"yes", "auto"}
    run_scenic = normalize_flag(row.get("run_scenic", "yes")) in {"yes", "auto"}
    contract = resolve_sample_contract(
        row=row,
        fastq_dir=fastq_dir,
        cellranger_out_dir=cellranger_out_dir,
        data_dir=data_dir,
    )

    notes = []
    if run_main and contract["has_filtered_matrix"] != "true":
        notes.append("missing_filtered_matrix")
    if run_velocity and contract["has_bam"] != "true":
        notes.append("missing_bam")
    if run_velocity and contract["platform_resolved"] != "10x_cellranger":
        notes.append("velocity_requires_cellranger")
    if run_velocity and contract["platform_resolved"] == "10x_cellranger" and contract["has_cellranger_out"] != "true":
        notes.append("missing_cellranger_out")
    if contract["input_mode"] == "fastq" and contract["has_fastq"] != "true":
        notes.append("missing_fastq")
    sample_decontx_ready = contract["has_filtered_matrix"] == "true"
    sample_soupx_ready = sample_decontx_ready and contract["has_raw_matrix"] == "true"
    sample_cellbender_ready = (
        sample_soupx_ready
        and contract["platform_resolved"] == "10x_cellranger"
        and contract.get("raw_matrix_kind", "") in {"mex_dir", "h5"}
    )
    sample_ambient_ready = sample_decontx_ready
    ambient_notes = []
    if contract["has_raw_matrix"] != "true":
        notes.append("missing_raw_matrix")
        ambient_notes.append("soupx_requires_raw_matrix")
    if not sample_ambient_ready:
        notes.append("ambient_branch_unavailable")
    if sample_ambient_ready and not sample_soupx_ready:
        ambient_notes.append("decontx_fallback_only")
    if contract["has_raw_matrix"] == "true" and not sample_cellbender_ready:
        ambient_notes.append("cellbender_stub_requires_10x_raw")
    if sample_cellbender_ready:
        ambient_notes.append("cellbender_stub_only")
    ambient_preferred_method = "soupx" if sample_soupx_ready else ("decontx" if sample_decontx_ready else "none")
    ambient_fallback_method = "decontx" if ambient_preferred_method == "soupx" and sample_decontx_ready else "none"
    if contract["gene_id_type_resolved"] == "unknown":
        notes.append("gene_id_type_unrecognized")
        any_unknown_gene_id_type = True
    biological_replicate = normalize_value(row.get("biological_replicate"))
    technical_replicate = normalize_value(row.get("technical_replicate"))
    batch = normalize_value(row.get("batch"))
    sample_replicate_fields_complete = bool(biological_replicate and technical_replicate and batch)
    if not sample_replicate_fields_complete:
        notes.append("replicate_fields_missing")

    inventory_rows.append(contract)

    sample_main_ready = contract["has_filtered_matrix"] == "true"
    sample_velocity_ready = (
        contract["platform_resolved"] == "10x_cellranger"
        and contract["has_cellranger_out"] == "true"
        and contract["has_bam"] == "true"
    )
    sample_scenic_ready = contract["has_filtered_matrix"] == "true"

    readiness_rows.append({
        "sample_id": sample_id,
        "platform_resolved": contract["platform_resolved"],
        "gene_id_type_resolved": contract["gene_id_type_resolved"],
        "reference_version": contract["reference_version"],
        "main_upstream_ready": str(sample_main_ready).lower(),
        "velocity_upstream_ready": str(sample_velocity_ready).lower(),
        "ambient_upstream_ready": str(sample_ambient_ready).lower(),
        "scenic_upstream_ready": str(sample_scenic_ready).lower(),
        "de_replicate_ready": "false",
        "raw_matrix_available": contract["has_raw_matrix"],
        "raw_matrix_kind": contract.get("raw_matrix_kind", ""),
        "ambient_any_ready": str(sample_ambient_ready).lower(),
        "ambient_soupx_ready": str(sample_soupx_ready).lower(),
        "ambient_decontx_ready": str(sample_decontx_ready).lower(),
        "ambient_cellbender_ready": str(sample_cellbender_ready).lower(),
        "ambient_preferred_method": ambient_preferred_method,
        "ambient_fallback_method": ambient_fallback_method,
        "ambient_apply_default": "manual",
        "ambient_notes": ";".join(ambient_notes),
        "gene_id_type_recognized": str(contract["gene_id_type_resolved"] != "unknown").lower(),
        "replicate_info_sufficient": str(sample_replicate_fields_complete).lower(),
        "notes": ";".join(notes),
    })

    if run_main:
        main_requested = True
        scenic_requested = scenic_requested or run_scenic
        main_project_ready = main_project_ready and sample_main_ready
        scenic_project_ready = scenic_project_ready and sample_scenic_ready
        ambient_project_ready = ambient_project_ready and sample_ambient_ready
        any_raw_matrix_available = any_raw_matrix_available or contract["has_raw_matrix"] == "true"
    if run_velocity:
        velocity_requested = True
        velocity_project_ready = velocity_project_ready and sample_velocity_ready

enabled_comparisons = [
    row for row in comparison_rows
    if row.get("enabled", "").strip().lower() in {"yes", "true"}
]
de_replicate_ready = True
if enabled_comparisons:
    for row in enabled_comparisons:
        ident_1 = row["ident_1"].strip()
        ident_2 = row["ident_2"].strip()
        if len(conditions_to_reps.get(ident_1, set())) < 2 or len(conditions_to_reps.get(ident_2, set())) < 2:
            de_replicate_ready = False
            break
else:
    de_replicate_ready = False

reference_version_consistent = len(project_reference_versions) <= 1
replicate_info_sufficient = replicate_fields_complete and de_replicate_ready

project_notes = []
if main_requested and not main_project_ready:
    project_notes.append("main_inputs_incomplete")
if velocity_requested and not velocity_project_ready:
    project_notes.append("velocity_inputs_incomplete")
if scenic_requested and not scenic_project_ready:
    project_notes.append("scenic_inputs_incomplete")
if main_requested and not ambient_project_ready:
    project_notes.append("ambient_inputs_incomplete")
if any_unknown_gene_id_type:
    project_notes.append("gene_id_type_unrecognized")
if not reference_version_consistent:
    project_notes.append("reference_version_mixed")
if not replicate_info_sufficient:
    project_notes.append("replicates_insufficient")

main_inventory_rows = [
    row for row, sample_row in zip(inventory_rows, sample_rows)
    if normalize_flag(sample_row.get("run_main", "yes")) != "no"
]
project_platforms = {row["platform_resolved"] for row in main_inventory_rows if row["platform_resolved"]}
project_gene_id_types = {row["gene_id_type_resolved"] for row in main_inventory_rows if row["gene_id_type_resolved"]}
project_ambient_methods = {row["ambient_preferred_method"] for row in readiness_rows if row["sample_id"] != "__PROJECT__" and row["ambient_preferred_method"]}
project_ambient_fallbacks = {row["ambient_fallback_method"] for row in readiness_rows if row["sample_id"] != "__PROJECT__" and row["ambient_fallback_method"] and row["ambient_fallback_method"] != "none"}
project_ambient_notes = sorted({note for row in readiness_rows for note in row.get("ambient_notes", "").split(";") if row["sample_id"] != "__PROJECT__" for note in [note.strip()] if note})

readiness_rows.append({
    "sample_id": "__PROJECT__",
    "platform_resolved": next(iter(project_platforms)) if len(project_platforms) == 1 else "mixed",
    "gene_id_type_resolved": next(iter(project_gene_id_types)) if len(project_gene_id_types) == 1 else "mixed",
    "reference_version": next(iter(project_reference_versions)) if len(project_reference_versions) == 1 else "",
    "main_upstream_ready": str(main_requested and main_project_ready).lower(),
    "velocity_upstream_ready": str(velocity_requested and velocity_project_ready).lower(),
    "ambient_upstream_ready": str(main_requested and ambient_project_ready).lower(),
    "scenic_upstream_ready": str(main_requested and scenic_project_ready).lower(),
    "de_replicate_ready": str(de_replicate_ready).lower(),
    "raw_matrix_available": str(any_raw_matrix_available).lower(),
    "raw_matrix_kind": "mixed",
    "ambient_any_ready": str(main_requested and ambient_project_ready).lower(),
    "ambient_soupx_ready": "mixed" if any_raw_matrix_available and any(row["ambient_soupx_ready"] == "false" for row in readiness_rows if row["sample_id"] != "__PROJECT__") else str(any_raw_matrix_available).lower(),
    "ambient_decontx_ready": str(main_requested and ambient_project_ready).lower(),
    "ambient_cellbender_ready": "mixed" if any(row["ambient_cellbender_ready"] == "true" for row in readiness_rows if row["sample_id"] != "__PROJECT__") and any(row["ambient_cellbender_ready"] == "false" for row in readiness_rows if row["sample_id"] != "__PROJECT__") else str(any(row["ambient_cellbender_ready"] == "true" for row in readiness_rows if row["sample_id"] != "__PROJECT__")).lower(),
    "ambient_preferred_method": next(iter(project_ambient_methods)) if len(project_ambient_methods) == 1 else "mixed",
    "ambient_fallback_method": next(iter(project_ambient_fallbacks)) if len(project_ambient_fallbacks) == 1 else ("mixed" if len(project_ambient_fallbacks) > 1 else "none"),
    "ambient_apply_default": "manual",
    "ambient_notes": ";".join(project_ambient_notes),
    "gene_id_type_recognized": str(not any_unknown_gene_id_type).lower(),
    "replicate_info_sufficient": str(replicate_info_sufficient).lower(),
    "notes": ";".join(project_notes),
})

inventory_path.parent.mkdir(parents=True, exist_ok=True)
with inventory_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sample_id",
            "input_mode",
            "source_path",
            "platform",
            "platform_resolved",
            "gene_id_type",
            "gene_id_type_resolved",
            "reference_version",
            "resolved_input_root",
            "cellranger_sample_dir",
            "filtered_matrix_dir",
            "raw_matrix_dir",
            "raw_matrix_kind",
            "metrics_path",
            "bam_path",
            "web_summary_path",
            "has_fastq",
            "has_cellranger_out",
            "has_filtered_matrix",
            "has_raw_matrix",
            "has_bam",
            "has_web_summary",
            "has_metrics_summary",
            "feature_name_profile",
            "feature_name_source",
            "n_features_checked",
            "n_symbol_like",
            "n_ensembl_like",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(inventory_rows)

with readiness_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sample_id",
            "platform_resolved",
            "gene_id_type_resolved",
            "reference_version",
            "main_upstream_ready",
            "velocity_upstream_ready",
            "ambient_upstream_ready",
            "scenic_upstream_ready",
            "de_replicate_ready",
            "raw_matrix_available",
            "raw_matrix_kind",
            "ambient_any_ready",
            "ambient_soupx_ready",
            "ambient_decontx_ready",
            "ambient_cellbender_ready",
            "ambient_preferred_method",
            "ambient_fallback_method",
            "ambient_apply_default",
            "ambient_notes",
            "gene_id_type_recognized",
            "replicate_info_sufficient",
            "notes",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(readiness_rows)
PY

update_workflow_status \
  "input_audited" \
  "bash ${PIPELINE_ROOT}/workflow/04_standardize_inputs.sh"

echo "输入审计完成"
echo "inventory: ${INPUT_INVENTORY_FILE}"
echo "readiness: ${BRANCH_READINESS_FILE}"
