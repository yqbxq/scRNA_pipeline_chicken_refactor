from __future__ import annotations

import json
import os
import csv
from datetime import datetime, timezone
from pathlib import Path


LIANA_OUTPUT_COLUMNS = [
    "lr_axis_id",
    "condition_value",
    "source",
    "target",
    "ligand",
    "receptor",
    "ligand_human",
    "receptor_human",
    "liana_consensus_score",
    "n_methods_agreed",
    "liana_consensus_hit",
    "status",
    "reason",
]


def resolve_h5ad_path(path: str | Path) -> Path:
    item = Path(path)
    if item.is_file():
        return item
    if item.is_dir():
        files = sorted(item.glob("*.h5ad"), key=lambda value: value.stat().st_mtime)
        if files:
            return files[-1]
    raise FileNotFoundError(f"No H5AD file found at {item}")


def load_scrna_h5ad(path: str | Path):
    try:
        import anndata as ad
    except ImportError as exc:
        raise RuntimeError("Python package anndata is required to read H5AD.") from exc
    return ad.read_h5ad(resolve_h5ad_path(path))


def validate_h5ad_for_liana(adata, cluster_col: str, condition_col: str) -> None:
    missing = [col for col in (cluster_col, condition_col) if col not in adata.obs.columns]
    if missing:
        raise ValueError(f"H5AD obs missing required columns: {', '.join(missing)}")
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("H5AD must contain at least one observation and one variable.")


def read_cluster_eligibility(path: str | Path, allowed_tiers: set[str]) -> set[str]:
    item = Path(path)
    if not item.exists():
        return set()
    with item.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        if not reader.fieldnames or "cluster_id" not in reader.fieldnames or "evidence_tier" not in reader.fieldnames:
            return set()
        return {
            (row.get("cluster_id") or "").strip()
            for row in reader
            if (row.get("evidence_tier") or "").strip() in allowed_tiers and (row.get("cluster_id") or "").strip()
        }


def empty_liana_df(status: str, reason: str) -> list[dict[str, str]]:
    return []


def build_output_entry_py(path: str | Path, fmt: str, module: str, description: str, base_dir: str | Path) -> dict[str, str]:
    item = Path(path)
    base = Path(base_dir)
    try:
        rel = item.relative_to(base)
    except ValueError:
        rel = item
    return {
        "path": str(rel),
        "format": fmt,
        "module": module,
        "description": description,
    }


def write_liana_manifest(
    manifest_path: str | Path,
    output_path: str | Path,
    *,
    base_dir: str | Path,
    status: str,
    n_rows: int,
    inputs: dict[str, str] | None = None,
    version: str | None = None,
) -> None:
    manifest = {
        "module": "07b_liana_consensus",
        "status": status,
        "base_dir": str(base_dir),
        "outputs": {
            "liana_consensus_tsv": build_output_entry_py(
                output_path,
                "tsv",
                "07b_liana_consensus",
                "LIANA multi-method consensus ligand-receptor table",
                base_dir,
            )
        },
        "inputs": inputs or {},
        "depends_on": {},
        "metrics": {"row_n": int(n_rows)},
        "version": version or os.environ.get("MODULE_07B_LIANA_VERSION", "1.0"),
        "created_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    }
    path = Path(manifest_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
