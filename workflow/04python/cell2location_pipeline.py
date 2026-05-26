#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys

from helpers.scrna_io import resolve_scrna_h5ad_path
from helpers.spatial_io import resolve_spatial_h5ad_path


def check_environment(use_gpu: bool) -> tuple[str, str]:
    try:
        import cell2location  # noqa: F401
        import anndata  # noqa: F401
        import scanpy  # noqa: F401
        import torch
    except Exception as exc:
        return "skipped_no_python_env", f"cell2location Python stack import failed: {exc}"
    if use_gpu and not torch.cuda.is_available():
        return "skipped_no_gpu", "SPATIAL_C2L_USE_GPU requested but torch.cuda.is_available() is false"
    return "ok", ""


def main() -> int:
    parser = argparse.ArgumentParser(description="cell2location sidecar for ST07.")
    parser.add_argument("--check-env", action="store_true")
    parser.add_argument("--use-gpu", action="store_true")
    parser.add_argument("--ref-h5ad")
    parser.add_argument("--st-h5ad")
    parser.add_argument("--ref-module", default="03d_panorama")
    parser.add_argument("--st-module", default="spatial_03_region")
    parser.add_argument("--section-id", default="")
    parser.add_argument("--annotation-col")
    parser.add_argument("--out-dir")
    parser.add_argument("--ref-epochs", type=int, default=250)
    parser.add_argument("--st-epochs", type=int, default=5000)
    args = parser.parse_args()

    status, reason = check_environment(args.use_gpu)
    if args.check_env:
        if status == "ok":
            print("ok")
            return 0
        sys.stderr.write(f"{status}\t{reason}\n")
        return 20 if status == "skipped_no_python_env" else 21

    if status != "ok":
        sys.stderr.write(f"{status}\t{reason}\n")
        return 20 if status == "skipped_no_python_env" else 21

    try:
        args.ref_h5ad = str(resolve_scrna_h5ad_path(module=args.ref_module, fallback_path=args.ref_h5ad))
        args.st_h5ad = str(resolve_spatial_h5ad_path(module=args.st_module, section_id=args.section_id or None, fallback_path=args.st_h5ad))
    except Exception as exc:
        sys.stderr.write(f"failed_sidecar\tH5AD input resolution failed: {exc}\n")
        return 30

    missing = [name for name in ("ref_h5ad", "st_h5ad", "annotation_col", "out_dir") if not getattr(args, name)]
    if missing:
        sys.stderr.write(f"failed_sidecar\tmissing required runtime args: {','.join(missing)}\n")
        return 30

    sys.stderr.write(
        "failed_sidecar\tformal cell2location training is intentionally gated to the server runtime wrapper\n"
    )
    return 31


if __name__ == "__main__":
    raise SystemExit(main())
