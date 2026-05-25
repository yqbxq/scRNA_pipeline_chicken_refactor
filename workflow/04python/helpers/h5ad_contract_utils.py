from __future__ import annotations

import csv
import os
import warnings
from pathlib import Path
from typing import Any

REQUIRED_COLUMNS = [
    "field_path",
    "field_type",
    "required",
    "dtype",
    "gate",
    "modality",
    "semantics",
    "notes",
]
ALLOWED_FIELD_TYPES = {"scalar", "matrix", "dict", "dataframe"}
ALLOWED_REQUIRED = {"yes", "no", "conditional"}
ALLOWED_GATES = {"primary", "exploratory", "conditional"}
ALLOWED_MODALITIES = {"scrna", "spatial", "both"}
CORE_FIELDS = {"obs.sample_id", "obs.condition", "var.gene_id", "layers.counts"}


def repo_template_path() -> Path:
    return Path(__file__).resolve().parents[3] / "metadata" / "h5ad_export_contract.tsv.template"


def load_h5ad_contract(contract_path: str | Path | None = None, modality: str = "scrna") -> list[dict[str, str]]:
    if modality not in {"scrna", "spatial"}:
        raise ValueError(f"invalid modality: {modality}")
    path = Path(contract_path or os.environ.get("H5AD_CONTRACT_FILE", "metadata/h5ad_export_contract.tsv"))
    if not path.exists():
        path = repo_template_path()
    if not path.exists():
        raise FileNotFoundError(f"missing H5AD contract: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        rows = [dict(row) for row in reader]
    validate_h5ad_contract_schema(rows)
    return [row for row in rows if row["modality"] in {modality, "both"}]


def validate_h5ad_contract_schema(rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError("empty H5AD contract")
    missing_columns = [column for column in REQUIRED_COLUMNS if column not in rows[0]]
    if missing_columns:
        raise ValueError(f"H5AD contract missing columns: {', '.join(missing_columns)}")

    fields = {row["field_path"] for row in rows}
    missing_core = sorted(CORE_FIELDS - fields)
    if missing_core:
        raise ValueError(f"H5AD contract removed core fields: {', '.join(missing_core)}")

    for row in rows:
        if row["field_type"] not in ALLOWED_FIELD_TYPES:
            raise ValueError(f"invalid field_type: {row['field_type']}")
        if row["required"] not in ALLOWED_REQUIRED:
            raise ValueError(f"invalid required: {row['required']}")
        if row["gate"] not in ALLOWED_GATES:
            raise ValueError(f"invalid gate: {row['gate']}")
        if row["modality"] not in ALLOWED_MODALITIES:
            raise ValueError(f"invalid modality: {row['modality']}")


def _split_field(field_path: str) -> tuple[str, str]:
    location, _, field = field_path.partition(".")
    return location, field


def _nested_get(mapping: Any, dotted: str) -> tuple[bool, Any]:
    cursor = mapping
    for part in dotted.split("."):
        if isinstance(cursor, dict) and part in cursor:
            cursor = cursor[part]
        else:
            return False, None
    return True, cursor


def _value_from_mapping(mapping: dict[str, Any], field_path: str) -> tuple[bool, Any]:
    location, field = _split_field(field_path)
    if location not in mapping:
        return False, None
    container = mapping[location]
    if hasattr(container, "columns") and field in container.columns:
        return True, container[field]
    if isinstance(container, dict):
        return _nested_get(container, field)
    return False, None


def _value_from_anndata(adata: Any, field_path: str) -> tuple[bool, Any]:
    location, field = _split_field(field_path)
    if location == "obs":
        if field in getattr(adata, "obs").columns:
            return True, adata.obs[field]
    elif location == "var":
        if field in getattr(adata, "var").columns:
            return True, adata.var[field]
    elif location == "obsm":
        if field in getattr(adata, "obsm").keys():
            return True, adata.obsm[field]
    elif location == "layers":
        if field in getattr(adata, "layers").keys():
            return True, adata.layers[field]
    elif location == "uns":
        return _nested_get(getattr(adata, "uns"), field)
    elif location == "spatial":
        spatial = getattr(adata, "uns", {}).get("spatial", {})
        return _nested_get(spatial, field)
    return False, None


def detect_dtype(value: Any) -> str:
    if value is None:
        return "missing"
    if hasattr(value, "dtypes") and hasattr(value, "columns"):
        return "dataframe"
    if isinstance(value, dict):
        return "dict"
    dtype = getattr(value, "dtype", None)
    if dtype is not None:
        dtype_text = str(dtype)
        if "int" in dtype_text:
            return "int64"
        if "float" in dtype_text:
            return "float32"
        if "bool" in dtype_text:
            return "bool"
        if "object" in dtype_text or "str" in dtype_text:
            return "str"
    if isinstance(value, (list, tuple)) and value:
        if all(isinstance(item, (list, tuple)) for item in value):
            flat = [cell for row in value for cell in row]
            if flat and all(isinstance(item, int) and not isinstance(item, bool) for item in flat):
                return "int64"
            if flat and all(isinstance(item, (int, float)) and not isinstance(item, bool) for item in flat):
                return "float32"
        if all(isinstance(item, int) and not isinstance(item, bool) for item in value):
            return "int64"
        if all(isinstance(item, (int, float)) and not isinstance(item, bool) for item in value):
            return "float32"
        if all(isinstance(item, str) for item in value):
            return "str"
    if isinstance(value, str):
        return "str"
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int):
        return "int64"
    if isinstance(value, float):
        return "float32"
    return type(value).__name__


def dtype_compatible(actual: str, expected: str) -> bool:
    if not expected or expected == "-":
        return True
    if actual == expected:
        return True
    if expected in {"float32", "float64"} and actual in {"float32", "float64", "int64", "int32"}:
        return True
    if expected in {"int64", "int32"} and actual in {"int64", "int32"}:
        return True
    return expected == "str" and actual == "str"


def check_h5ad_contract_anndata(adata: Any, contract: list[dict[str, str]], fail_on: str = "warn") -> dict[str, Any]:
    if fail_on not in {"error", "warn", "none"}:
        raise ValueError(f"invalid fail_on: {fail_on}")
    validate_h5ad_contract_schema(contract)
    is_mapping = isinstance(adata, dict)
    violations: list[dict[str, str]] = []

    for row in contract:
        field_path = row["field_path"]
        found, value = _value_from_mapping(adata, field_path) if is_mapping else _value_from_anndata(adata, field_path)
        severity = "error" if row["gate"] == "primary" else "warn"
        if not found:
            if row["required"] == "yes":
                location, _ = _split_field(field_path)
                violations.append({
                    "field": field_path,
                    "severity": severity,
                    "reason": f"required=yes but field not found in {location}",
                    "expected_dtype": row["dtype"],
                    "actual_dtype": "missing",
                })
            continue

        actual = detect_dtype(value)
        if not dtype_compatible(actual, row["dtype"]):
            violations.append({
                "field": field_path,
                "severity": severity,
                "reason": f"expected dtype={row['dtype']} got {actual}",
                "expected_dtype": row["dtype"],
                "actual_dtype": actual,
            })

    primary_failed = any(item["severity"] == "error" for item in violations)
    if primary_failed and fail_on == "error":
        fields = ", ".join(item["field"] for item in violations if item["severity"] == "error")
        raise ValueError(f"H5AD contract violations: {fields}")
    if violations and fail_on == "warn":
        warnings.warn(f"H5AD contract violations: {len(violations)} field(s)", stacklevel=2)
    return {
        "passed": not primary_failed,
        "violations": violations,
        "summary": {
            "violation_n": len(violations),
            "error_n": sum(item["severity"] == "error" for item in violations),
            "warn_n": sum(item["severity"] == "warn" for item in violations),
        },
    }


def contract_violations_to_manifest(check_result: dict[str, Any]) -> dict[str, Any]:
    violations = check_result.get("violations", [])
    return {
        "contract_passed": bool(check_result.get("passed")),
        "contract_violations_count": len(violations),
        "contract_violations_summary": violations,
    }
