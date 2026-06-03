#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from helpers.spatial_io import resolve_spatial_h5ad_path_strict


def main() -> int:
    parser = argparse.ArgumentParser(description="H5AD-first placeholder for joint Python spatial sidecars.")
    parser.add_argument("--st-h5ad")
    parser.add_argument("--st-module", default="spatial_03_region")
    parser.add_argument("--section-id", default="")
    parser.add_argument("--out-dir", default="results/spatial/tables/joint_h5ad_sidecar")
    parser.add_argument("--results-dir", default=None)
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        h5ad = resolve_spatial_h5ad_path_strict(module=args.st_module, section_id=args.section_id or None, base_dir=args.results_dir, fallback_path=args.st_h5ad)
        status, reason = "ok", ""
    except Exception as exc:
        h5ad, status, reason = "", "skipped_no_h5ad", str(exc)
    try:
        import pandas as pd
        pd.DataFrame([{
            "status": status,
            "reason": reason,
            "input_h5ad_path": str(h5ad),
            "sidecar_role": "joint_python_h5ad_first",
        }]).to_csv(out_dir / "joint_h5ad_sidecar_manifest.tsv", sep="\t", index=False)
    except Exception:
        return 30
    return 0 if status == "ok" else 20


if __name__ == "__main__":
    raise SystemExit(main())
