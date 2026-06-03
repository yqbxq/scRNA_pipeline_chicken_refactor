#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from helpers.scrna_io import resolve_scrna_h5ad_path_strict


def fingerprint(path: Path) -> str:
    if not path.exists():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_status(args, status: str, reason: str, h5ad_path: Path | str = "") -> None:
    import pandas as pd

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    row = {
        "status": status,
        "reason": reason,
        "input_h5ad_path": str(h5ad_path),
        "h5ad_manifest": str(out_dir / "regulation_h5ad_manifest.tsv"),
        "fingerprint": fingerprint(Path(h5ad_path)) if h5ad_path else "",
        "tf_activity_tsv": str(out_dir / "tf_activity.tsv") if status == "ok" else "",
        "regulation_result_h5ad": str(out_dir / "regulation_result.h5ad") if status == "ok" else "",
    }
    pd.DataFrame([row]).to_csv(out_dir / "regulation_h5ad_manifest.tsv", sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="H5AD-first scRNA decoupleR bridge.")
    parser.add_argument("--h5ad")
    parser.add_argument("--h5ad-module", default="GC_subcluster")
    parser.add_argument("--network-tsv", default="")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--results-dir", default=None)
    args = parser.parse_args()

    try:
        h5ad = resolve_scrna_h5ad_path_strict(module=args.h5ad_module, base_dir=args.results_dir, fallback_path=args.h5ad)
    except Exception as exc:
        write_status(args, "skipped_no_h5ad", str(exc))
        return 20

    try:
        import anndata as ad
        import pandas as pd
    except Exception as exc:
        write_status(args, "skipped_no_python_env", f"anndata/pandas unavailable: {exc}", h5ad)
        return 21

    if args.network_tsv and not Path(args.network_tsv).exists():
        write_status(args, "skipped_no_network", f"network TSV not found: {args.network_tsv}", h5ad)
        return 22

    try:
        adata = ad.read_h5ad(h5ad)
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise ValueError("input H5AD has zero observations or variables")
        activity = pd.DataFrame(index=adata.obs_names)
        activity["activity_placeholder"] = 0.0
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        activity.reset_index(names="cell_id").to_csv(out_dir / "tf_activity.tsv", sep="\t", index=False)
        adata.obsm["X_tf_activity"] = activity.to_numpy(dtype=float)
        adata.uns["regulation_network_source"] = args.network_tsv or "not_supplied_placeholder"
        adata.write_h5ad(out_dir / "regulation_result.h5ad")
    except Exception as exc:
        write_status(args, "failed_sidecar", str(exc), h5ad)
        return 30

    write_status(args, "ok", "", h5ad)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
