#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


SPATIAL_PLATFORMS = {"visium", "saw", "stomics", "generic"}
SPATIAL_INPUT_MODES = {"visium_bundle", "saw_bundle", "stomics_bundle", "spatial_matrix"}


def norm(raw: str | None) -> str:
    return (raw or "").strip()


def norm_lower(raw: str | None) -> str:
    return norm(raw).lower()


def read_tsv(path: Path, *, allow_missing: bool = False) -> list[dict[str, str]]:
    if allow_missing and (not path.exists() or path.stat().st_size == 0):
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        fieldnames = [field[2:] if field.startswith("# ") else field for field in reader.fieldnames]
        rows = []
        for row in reader:
            rows.append({fieldnames[i]: row.get(reader.fieldnames[i], "") for i in range(len(fieldnames))})
        return rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def ordered_union(*field_groups: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for fields in field_groups:
        for field in fields:
            if field and field not in seen:
                seen.add(field)
                out.append(field)
    return out


def is_spatial_sample(row: dict[str, str]) -> bool:
    modality = norm_lower(row.get("modality")) or "scrna"
    platform = norm_lower(row.get("platform"))
    input_mode = norm_lower(row.get("input_mode"))
    run_spatial = norm_lower(row.get("run_spatial")) or "auto"
    return (
        modality == "spatial"
        or platform in SPATIAL_PLATFORMS
        or input_mode in SPATIAL_INPUT_MODES
        or run_spatial == "yes"
    )


def split_paths(raw: str) -> list[str]:
    return [item.strip() for item in (raw or "").split(",") if item.strip()]


def fill_tokens(raw: str, sample: dict[str, str]) -> str:
    replacements = {
        "<sample_id>": norm(sample.get("sample_id")),
        "<section_id>": norm(sample.get("section_id")),
        "<chip_id>": norm(sample.get("chip_id")) or norm(sample.get("section_id")),
    }
    out = raw
    for key, value in replacements.items():
        out = out.replace(key, value)
    return out


def candidate_relpath(root: Path, relpath: str) -> Path:
    candidate = root / relpath
    if candidate.exists():
        return candidate
    if relpath.startswith("outs/") and root.name == "outs":
        return root / relpath.removeprefix("outs/")
    return candidate


def find_contract(
    contracts: list[dict[str, str]],
    platform: str,
    bundle_layout: str,
) -> dict[str, str] | None:
    platform = platform or "generic"
    for row in contracts:
        if norm_lower(row.get("platform")) == platform and norm_lower(row.get("bundle_layout")) == bundle_layout:
            return row
    for row in contracts:
        if norm_lower(row.get("platform")) == platform:
            return row
    return None


def spatial_rows(samples: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in samples if is_spatial_sample(row)]


def audit_spatial(args: argparse.Namespace) -> int:
    samples = read_tsv(args.samples)
    sections = {norm(row.get("section_id")): row for row in read_tsv(args.sections, allow_missing=True)}
    contracts = read_tsv(args.contract)
    data_dir = args.data_dir

    inventory_rows: list[dict[str, str]] = []
    for sample in spatial_rows(samples):
        sample_id = norm(sample.get("sample_id"))
        section_id = norm(sample.get("section_id"))
        platform = norm_lower(sample.get("platform")) or "generic"
        if platform == "auto":
            platform = "generic"
        bundle_layout = norm_lower(sample.get("bundle_layout"))
        if not bundle_layout:
            bundle_layout = {
                "visium": "outs_visium",
                "saw": "saw_bin50",
                "stomics": "stomics_native",
            }.get(platform, "spatial_matrix")

        source_path = Path(norm(sample.get("source_path")) or ".")
        contract = find_contract(contracts, platform, bundle_layout)
        section_row = sections.get(section_id, {})
        notes: list[str] = []
        if not section_id:
            notes.append("missing_section_id")
        elif section_id not in sections:
            notes.append("section_not_registered")
        if contract is None:
            notes.append("missing_platform_contract")

        required_relpaths = []
        optional_relpaths = []
        standardized_outs_relpath = "outs"
        ready_status_when_present = "true"
        if contract is not None:
            required_relpaths = [fill_tokens(item, sample) for item in split_paths(contract.get("required_relpaths", ""))]
            optional_relpaths = [fill_tokens(item, sample) for item in split_paths(contract.get("optional_relpaths", ""))]
            standardized_outs_relpath = norm(contract.get("standardized_outs_relpath")) or "outs"
            ready_status_when_present = norm(contract.get("ready_status_when_present")) or "true"

        missing_required = []
        present_optional = []
        missing_optional = []
        if not source_path.exists():
            notes.append("missing_source_path")
            missing_required = required_relpaths[:]
            missing_optional = optional_relpaths[:]
        else:
            for relpath in required_relpaths:
                if not candidate_relpath(source_path, relpath).exists():
                    missing_required.append(relpath)
            for relpath in optional_relpaths:
                if candidate_relpath(source_path, relpath).exists():
                    present_optional.append(relpath)
                else:
                    missing_optional.append(relpath)

        ready_status = ready_status_when_present if not missing_required and contract is not None else "false"
        ready = "true" if ready_status.startswith("true") else "false"
        standardized_outs = data_dir / sample_id / standardized_outs_relpath

        inventory_rows.append(
            {
                "sample_id": sample_id,
                "section_id": section_id,
                "modality": "spatial",
                "condition": norm(sample.get("condition")),
                "platform": platform,
                "bundle_layout": bundle_layout,
                "source_path": str(source_path),
                "section_bundle_root": norm(section_row.get("bundle_root")),
                "standardized_outs": str(standardized_outs),
                "ready": ready,
                "ready_status": ready_status,
                "required_relpaths": ",".join(required_relpaths),
                "missing_required_relpaths": ",".join(missing_required),
                "optional_present_relpaths": ",".join(present_optional),
                "optional_missing_relpaths": ",".join(missing_optional),
                "notes": ";".join(notes),
            }
        )

    inventory_fields = [
        "sample_id",
        "section_id",
        "modality",
        "condition",
        "platform",
        "bundle_layout",
        "source_path",
        "section_bundle_root",
        "standardized_outs",
        "ready",
        "ready_status",
        "required_relpaths",
        "missing_required_relpaths",
        "optional_present_relpaths",
        "optional_missing_relpaths",
        "notes",
    ]
    write_tsv(args.output, inventory_rows, inventory_fields)
    if args.readiness is not None:
        merge_spatial_readiness(args.readiness, inventory_rows)
    print(f"spatial_samples={len(inventory_rows)}")
    print(f"spatial_ready={sum(1 for row in inventory_rows if row['ready'] == 'true')}")
    return 0


def merge_spatial_readiness(path: Path, inventory_rows: list[dict[str, str]]) -> None:
    existing_rows = read_tsv(path, allow_missing=True)
    existing_fields = []
    if path.exists() and path.stat().st_size > 0:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            existing_fields = reader.fieldnames or []

    spatial_fields = [
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
        "modality",
        "section_id",
        "spatial_upstream_ready",
        "spatial_ready_status",
        "spatial_standardized_outs",
        "notes",
    ]
    rows = [row for row in existing_rows if row.get("modality") != "spatial"]
    spatial_project_ready = bool(inventory_rows) and all(row["ready"] == "true" for row in inventory_rows)
    for inv in inventory_rows:
        notes = [note for note in inv.get("notes", "").split(";") if note]
        if inv["ready"] != "true":
            notes.append("spatial_inputs_incomplete")
        rows.append(
            {
                "sample_id": inv["sample_id"],
                "platform_resolved": inv["platform"],
                "gene_id_type_resolved": "auto",
                "reference_version": "",
                "main_upstream_ready": "false",
                "velocity_upstream_ready": "false",
                "ambient_upstream_ready": "false",
                "scenic_upstream_ready": "false",
                "de_replicate_ready": "false",
                "raw_matrix_available": "false",
                "raw_matrix_kind": "spatial_bundle",
                "ambient_any_ready": "false",
                "ambient_soupx_ready": "false",
                "ambient_decontx_ready": "false",
                "ambient_cellbender_ready": "false",
                "ambient_preferred_method": "none",
                "ambient_fallback_method": "none",
                "ambient_apply_default": "none",
                "ambient_notes": "",
                "gene_id_type_recognized": "true",
                "replicate_info_sufficient": "true",
                "modality": "spatial",
                "section_id": inv["section_id"],
                "spatial_upstream_ready": inv["ready"],
                "spatial_ready_status": inv["ready_status"],
                "spatial_standardized_outs": inv["standardized_outs"],
                "notes": ";".join(notes),
            }
        )

    project_row = next((row for row in rows if row.get("sample_id") == "__PROJECT__"), None)
    if project_row is not None:
        project_row["spatial_upstream_ready"] = str(spatial_project_ready).lower()
        existing_notes = [note for note in project_row.get("notes", "").split(";") if note]
        if inventory_rows and not spatial_project_ready and "spatial_inputs_incomplete" not in existing_notes:
            existing_notes.append("spatial_inputs_incomplete")
        project_row["notes"] = ";".join(existing_notes)

    write_tsv(path, rows, ordered_union(existing_fields, spatial_fields))


def replace_path(dst: Path, src: Path, mode: str) -> None:
    if dst.exists() or dst.is_symlink():
        if dst.is_symlink() or dst.is_file():
            dst.unlink()
        else:
            shutil.rmtree(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if mode == "symlink":
        dst.symlink_to(src, target_is_directory=src.is_dir())
    elif mode == "copy":
        if src.is_dir():
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)
    else:
        raise SystemExit("standardize mode must be symlink or copy")


def standardize_spatial(args: argparse.Namespace) -> int:
    rows = read_tsv(args.inventory)
    standardized = 0
    skipped = 0
    for row in rows:
        if row.get("ready") != "true":
            skipped += 1
            continue
        if row.get("ready_status") == "true_pending_intake_python":
            skipped += 1
            continue
        source_root = Path(row["source_path"])
        rel_candidates = ["outs"]
        if source_root.name == "outs":
            rel_candidates = ["."]
        if row.get("bundle_layout") == "spatial_matrix":
            rel_candidates = ["."]
        source = None
        for rel in rel_candidates:
            candidate = source_root if rel == "." else source_root / rel
            if candidate.exists():
                source = candidate
                break
        if source is None:
            skipped += 1
            continue
        replace_path(Path(row["standardized_outs"]), source, args.mode)
        standardized += 1
    print(f"spatial_standardized={standardized}")
    print(f"spatial_skipped={skipped}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Spatial transcriptomics intake helpers.")
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    audit = subparsers.add_parser("audit", help="Audit spatial bundle contracts.")
    audit.add_argument("--samples", required=True, type=Path)
    audit.add_argument("--sections", required=True, type=Path)
    audit.add_argument("--contract", required=True, type=Path)
    audit.add_argument("--output", required=True, type=Path)
    audit.add_argument("--readiness", type=Path)
    audit.add_argument("--data-dir", required=True, type=Path)
    audit.set_defaults(func=audit_spatial)

    standardize = subparsers.add_parser("standardize", help="Link or copy ready spatial bundle inputs.")
    standardize.add_argument("--inventory", required=True, type=Path)
    standardize.add_argument("--mode", choices=["symlink", "copy"], default="symlink")
    standardize.set_defaults(func=standardize_spatial)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
