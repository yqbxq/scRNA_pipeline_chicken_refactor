#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from helpers.spatial_io import resolve_spatial_h5ad_path_strict


def fingerprint(path: Path) -> str:
    if not path.exists():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(args, status: str, reason: str, h5ad_path: Path | str = "", svg_tsv: str = "") -> None:
    import pandas as pd

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{
        "section_id": args.section_id or "all",
        "status": status,
        "reason": reason,
        "input_h5ad_path": str(h5ad_path),
        "fingerprint": fingerprint(Path(h5ad_path)) if h5ad_path else "",
        "svg_tsv": svg_tsv,
    }]).to_csv(out_dir / "spatialde2_manifest.tsv", sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="H5AD-first SpatialDE2 SVG sidecar.")
    parser.add_argument("--st-h5ad")
    parser.add_argument("--st-module", default="spatial_03_region")
    parser.add_argument("--section-id", default="")
    parser.add_argument("--out-dir", default="results/spatial/tables/09_svg/spatialde2")
    parser.add_argument("--results-dir", default=None)
    args = parser.parse_args()

    try:
        h5ad = resolve_spatial_h5ad_path_strict(module=args.st_module, section_id=args.section_id or None, base_dir=args.results_dir, fallback_path=args.st_h5ad)
    except Exception as exc:
        write_manifest(args, "skipped_no_h5ad", str(exc))
        return 20
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
    except Exception as exc:
        write_manifest(args, "skipped_no_python_env", f"anndata/numpy/pandas unavailable: {exc}", h5ad)
        return 21

    try:
        adata = ad.read_h5ad(h5ad)
        if "spatial" not in adata.obsm:
            raise ValueError('input H5AD missing obsm["spatial"]')
        section = args.section_id or "all"
        if args.section_id and "section_id" in adata.obs:
            keep = adata.obs["section_id"].astype(str) == args.section_id
            adata = adata[keep].copy()
        if adata.n_obs < 3 or adata.n_vars == 0:
            write_manifest(args, "skipped_too_few_spots", "need at least 3 spots and 1 gene", h5ad)
            return 22
        x = adata.X
        means = np.asarray(x.mean(axis=0)).reshape(-1)
        variances = np.asarray(x.var(axis=0)).reshape(-1) if hasattr(x, "var") else np.zeros_like(means)
        genes = adata.var["gene_symbol"].astype(str).to_numpy() if "gene_symbol" in adata.var else adata.var_names.astype(str)
        df = pd.DataFrame({"gene": genes, "mean_expression": means, "spatial_variance_proxy": variances})
        df = df.sort_values(["spatial_variance_proxy", "mean_expression", "gene"], ascending=[False, False, True])
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        svg_tsv = out_dir / f"{section}_svg.tsv"
        df.to_csv(svg_tsv, sep="\t", index=False)
    except Exception as exc:
        write_manifest(args, "failed_sidecar", str(exc), h5ad)
        return 30

    write_manifest(args, "ok", "", h5ad, str(svg_tsv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
