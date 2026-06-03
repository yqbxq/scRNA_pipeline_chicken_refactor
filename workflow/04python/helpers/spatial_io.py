from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import Callable


def _results_dir(base_dir: str | Path | None = None) -> Path:
    if base_dir is not None:
        return Path(base_dir)
    return Path(os.environ.get("RESULTS_DIR", "results"))


def _module_dir(module: str, base_dir: str | Path | None = None) -> Path:
    return _results_dir(base_dir) / "90a_export_h5ad" / module


def find_spatial_h5ad(module: str = "spatial_03_region", section_id: str | None = None, base_dir: str | Path | None = None) -> Path | None:
    h5ad_dir = _module_dir(module, base_dir)
    files = sorted(h5ad_dir.glob("*.h5ad"), key=lambda item: item.stat().st_mtime if item.exists() else 0)
    if section_id:
        section_hits = [item for item in files if section_id in item.stem]
        if section_hits:
            return section_hits[-1]
    return files[-1] if files else None


def load_spatial_h5ad(
    section_id: str | None = None,
    module: str = "spatial_03_region",
    base_dir: str | Path | None = None,
    fallback: Callable[[], object] | None = None,
):
    path = find_spatial_h5ad(module=module, section_id=section_id, base_dir=base_dir)
    if path is not None:
        import anndata as ad

        return ad.read_h5ad(path)

    message = f"H5AD not found for spatial module={module} section_id={section_id or '*'}"
    if os.environ.get("H5AD_PYTHON_FALLBACK_TO_RDS", "yes").lower() in {"yes", "true", "1", "on"} and fallback is not None:
        warnings.warn(f"{message}; falling back to legacy RDS/explicit input path. DEPRECATED: P-R02-E", stacklevel=2)
        return fallback()
    raise FileNotFoundError(message)


def load_all_spatial_h5ad(module: str = "spatial_03_region", base_dir: str | Path | None = None) -> dict[str, object]:
    import anndata as ad

    h5ad_dir = _module_dir(module, base_dir)
    files = sorted(h5ad_dir.glob("*.h5ad"))
    if not files:
        raise FileNotFoundError(f"No spatial H5AD files under {h5ad_dir}")
    return {item.stem: ad.read_h5ad(item) for item in files}


def resolve_spatial_h5ad_path(
    module: str = "spatial_03_region",
    section_id: str | None = None,
    base_dir: str | Path | None = None,
    fallback_path: str | Path | None = None,
) -> Path:
    path = find_spatial_h5ad(module=module, section_id=section_id, base_dir=base_dir)
    if path is not None:
        return path
    if fallback_path and os.environ.get("H5AD_PYTHON_FALLBACK_TO_RDS", "yes").lower() in {"yes", "true", "1", "on"}:
        warnings.warn("Using explicit fallback H5AD/RDS bridge path. DEPRECATED: P-R02-E", stacklevel=2)
        return Path(fallback_path)
    raise FileNotFoundError(f"H5AD not found for spatial module={module} section_id={section_id or '*'}")


def resolve_spatial_h5ad_path_strict(
    module: str = "spatial_03_region",
    section_id: str | None = None,
    base_dir: str | Path | None = None,
    fallback_path: str | Path | None = None,
) -> Path:
    path = find_spatial_h5ad(module=module, section_id=section_id, base_dir=base_dir)
    if path is not None:
        return path
    if fallback_path:
        candidate = Path(fallback_path)
        if candidate.exists() and candidate.suffix.lower() == ".h5ad":
            return candidate
    raise FileNotFoundError(f"H5AD-first spatial input not found for module={module} section_id={section_id or '*'}; fallback RDS paths are disabled")
