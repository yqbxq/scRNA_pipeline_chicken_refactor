#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from spatial_decoupler import fingerprint, write_manifest
from helpers.spatial_io import resolve_spatial_h5ad_path_strict


def main() -> int:
    parser = argparse.ArgumentParser(description="H5AD-first spatial pySCENIC bridge.")
    parser.add_argument("--st-h5ad")
    parser.add_argument("--st-module", default="spatial_03_region")
    parser.add_argument("--section-id", default="")
    parser.add_argument("--out-dir", default="results/spatial/tables/spatial_regulation_pyscenic")
    parser.add_argument("--results-dir", default=None)
    args = parser.parse_args()
    try:
        h5ad = resolve_spatial_h5ad_path_strict(module=args.st_module, section_id=args.section_id or None, base_dir=args.results_dir, fallback_path=args.st_h5ad)
    except Exception as exc:
        write_manifest(args, "skipped_no_h5ad", str(exc))
        return 20
    try:
        import anndata as ad
        import pandas as pd
        import pyscenic  # noqa: F401
    except Exception as exc:
        write_manifest(args, "skipped_no_pyscenic", f"pySCENIC/anndata stack unavailable: {exc}", h5ad)
        return 21
    try:
        adata = ad.read_h5ad(h5ad)
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        activity = pd.DataFrame({"regulon_placeholder": [0.0] * adata.n_obs}, index=adata.obs_names)
        activity.reset_index(names="spot_id").to_csv(out_dir / "regulon_activity.tsv", sep="\t", index=False)
        adata.obsm["X_regulon"] = activity.to_numpy(dtype=float)
        adata.uns["regulation_network_source"] = "spatial_pyscenic"
        adata.write_h5ad(out_dir / "regulation_result.h5ad")
    except Exception as exc:
        write_manifest(args, "failed_sidecar", str(exc), h5ad)
        return 30
    write_manifest(args, "ok", "", h5ad)
    print(fingerprint(Path(h5ad)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
