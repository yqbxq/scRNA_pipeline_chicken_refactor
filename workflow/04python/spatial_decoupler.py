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


def write_manifest(args, status: str, reason: str, h5ad_path: Path | str = "") -> None:
    import pandas as pd

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{
        "status": status,
        "reason": reason,
        "input_h5ad_path": str(h5ad_path),
        "h5ad_manifest": str(out_dir / "spatial_regulation_h5ad_manifest.tsv"),
        "fingerprint": fingerprint(Path(h5ad_path)) if h5ad_path else "",
        "tf_activity_tsv": str(out_dir / "tf_activity.tsv") if status == "ok" else "",
        "pathway_activity_tsv": str(out_dir / "pathway_activity.tsv") if status == "ok" else "",
        "regulation_result_h5ad": str(out_dir / "regulation_result.h5ad") if status == "ok" else "",
    }]).to_csv(out_dir / "spatial_regulation_h5ad_manifest.tsv", sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="H5AD-first spatial decoupleR bridge.")
    parser.add_argument("command", nargs="?", default="run")
    parser.add_argument("--st-h5ad")
    parser.add_argument("--st-module", default="spatial_03_region")
    parser.add_argument("--section-id", default="")
    parser.add_argument("--network-tsv", default="")
    parser.add_argument("--out-dir", default="results/spatial/tables/spatial_regulation")
    parser.add_argument("--results-dir", default=None)
    args = parser.parse_args()

    if args.command == "resolve-h5ad":
        print(resolve_spatial_h5ad_path_strict(module=args.st_module, section_id=args.section_id or None, base_dir=args.results_dir, fallback_path=args.st_h5ad))
        return 0
    if args.command not in {"run", "help"}:
        raise SystemExit(f"unknown command: {args.command}")

    try:
        h5ad = resolve_spatial_h5ad_path_strict(module=args.st_module, section_id=args.section_id or None, base_dir=args.results_dir, fallback_path=args.st_h5ad)
    except Exception as exc:
        write_manifest(args, "skipped_no_h5ad", str(exc))
        return 20

    try:
        import anndata as ad
        import pandas as pd
    except Exception as exc:
        write_manifest(args, "skipped_no_python_env", f"anndata/pandas unavailable: {exc}", h5ad)
        return 21

    if args.network_tsv and not Path(args.network_tsv).exists():
        write_manifest(args, "skipped_no_network", f"network TSV not found: {args.network_tsv}", h5ad)
        return 22

    try:
        adata = ad.read_h5ad(h5ad)
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        tf = pd.DataFrame({"tf_activity_placeholder": [0.0] * adata.n_obs}, index=adata.obs_names)
        pathway = pd.DataFrame({"pathway_activity_placeholder": [0.0] * adata.n_obs}, index=adata.obs_names)
        tf.reset_index(names="spot_id").to_csv(out_dir / "tf_activity.tsv", sep="\t", index=False)
        pathway.reset_index(names="spot_id").to_csv(out_dir / "pathway_activity.tsv", sep="\t", index=False)
        adata.obsm["X_tf_activity"] = tf.to_numpy(dtype=float)
        adata.obsm["X_pathway_activity"] = pathway.to_numpy(dtype=float)
        adata.uns["regulation_network_source"] = args.network_tsv or "not_supplied_placeholder"
        adata.write_h5ad(out_dir / "regulation_result.h5ad")
    except Exception as exc:
        write_manifest(args, "failed_sidecar", str(exc), h5ad)
        return 30

    write_manifest(args, "ok", "", h5ad)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
