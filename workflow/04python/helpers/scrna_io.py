from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import Callable


def _results_dir(base_dir: str | Path | None = None) -> Path:
    if base_dir is not None:
        return Path(base_dir)
    return Path(os.environ.get("RESULTS_DIR", "results"))


def find_scrna_h5ad(module: str = "03d_panorama", base_dir: str | Path | None = None) -> Path | None:
    h5ad_dir = _results_dir(base_dir) / "90a_export_h5ad" / module
    files = sorted(h5ad_dir.glob("*.h5ad"), key=lambda item: item.stat().st_mtime if item.exists() else 0)
    return files[-1] if files else None


def load_scrna_h5ad(
    module: str = "03d_panorama",
    base_dir: str | Path | None = None,
    fallback: Callable[[], object] | None = None,
):
    path = find_scrna_h5ad(module=module, base_dir=base_dir)
    if path is not None:
        import anndata as ad

        return ad.read_h5ad(path)

    message = f"H5AD not found for scRNA module={module}; expected under {_results_dir(base_dir) / '90a_export_h5ad' / module}"
    if os.environ.get("H5AD_PYTHON_FALLBACK_TO_RDS", "yes").lower() in {"yes", "true", "1", "on"} and fallback is not None:
        warnings.warn(f"{message}; falling back to legacy RDS/explicit input path. DEPRECATED: P-R02-E", stacklevel=2)
        return fallback()
    raise FileNotFoundError(message)


def resolve_scrna_h5ad_path(module: str = "03d_panorama", base_dir: str | Path | None = None, fallback_path: str | Path | None = None) -> Path:
    path = find_scrna_h5ad(module=module, base_dir=base_dir)
    if path is not None:
        return path
    if fallback_path and os.environ.get("H5AD_PYTHON_FALLBACK_TO_RDS", "yes").lower() in {"yes", "true", "1", "on"}:
        warnings.warn("Using explicit fallback H5AD/RDS bridge path. DEPRECATED: P-R02-E", stacklevel=2)
        return Path(fallback_path)
    raise FileNotFoundError(f"H5AD not found for scRNA module={module}")


def resolve_scrna_h5ad_path_strict(module: str = "03d_panorama", base_dir: str | Path | None = None, fallback_path: str | Path | None = None) -> Path:
    path = find_scrna_h5ad(module=module, base_dir=base_dir)
    if path is not None:
        return path
    if fallback_path:
        candidate = Path(fallback_path)
        if candidate.exists() and candidate.suffix.lower() == ".h5ad":
            return candidate
    raise FileNotFoundError(f"H5AD-first scRNA input not found for module={module}; fallback RDS paths are disabled")
