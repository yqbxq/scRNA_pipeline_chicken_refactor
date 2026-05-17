#!/usr/bin/env python3
"""STAGATE clustering sidecar for spatial module 03."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def _write_error(path: str | None, message: str) -> None:
    if not path:
        return
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump({"status": "failed", "message": message}, handle, indent=2)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input h5ad with obsm['spatial'].")
    parser.add_argument("--output-cluster-tsv", required=True)
    parser.add_argument("--output-metadata-json", required=True)
    parser.add_argument("--target-clusters", type=int, default=8)
    parser.add_argument("--rad-cutoff", type=float, default=150.0)
    parser.add_argument("--device", default=os.environ.get("STAGATE_DEVICE", "cpu"))
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    try:
        import anndata
        import numpy as np
        import scanpy as sc
        import STAGATE_pyG as st
        import torch
    except Exception as exc:  # pragma: no cover - exercised by smoke via exit code
        message = f"STAGATE bridge dependency import failed: {exc}"
        _write_error(args.output_metadata_json, message)
        print(message, file=sys.stderr)
        return 2

    try:
        adata = anndata.read_h5ad(args.input)
        if "spatial" not in adata.obsm:
            raise ValueError("input h5ad is missing obsm['spatial']")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        st.Cal_Spatial_Net(adata, rad_cutoff=args.rad_cutoff)
        np.random.seed(args.seed)
        torch.manual_seed(args.seed)
        adata = st.train_STAGATE(
            adata,
            hidden_dims=[512, 30],
            n_epochs=300,
            lr=0.001,
            device=args.device,
            random_seed=args.seed,
        )
        sc.pp.neighbors(adata, use_rep="STAGATE")
        chosen_key = "leiden"
        sc.tl.leiden(adata, resolution=0.6, key_added=chosen_key)
        for resolution in [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2]:
            key = f"leiden_{resolution:g}"
            sc.tl.leiden(adata, resolution=resolution, key_added=key)
            chosen_key = key
            if adata.obs[key].nunique() >= args.target_clusters:
                break
        labels = adata.obs[chosen_key].astype(str).to_numpy()
    except Exception as exc:  # pragma: no cover
        message = f"STAGATE bridge run failed: {exc}"
        _write_error(args.output_metadata_json, message)
        print(message, file=sys.stderr)
        return 3

    Path(args.output_cluster_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_cluster_tsv, "w", encoding="utf-8") as handle:
        handle.write("barcode\tcluster\n")
        for barcode, cluster in zip(adata.obs_names, labels):
            handle.write(f"{barcode}\t{cluster}\n")

    with open(args.output_metadata_json, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "status": "ok",
                "n_obs": int(adata.n_obs),
                "n_clusters": int(len(set(labels))),
                "target_clusters": int(args.target_clusters),
                "cluster_key": chosen_key,
                "rad_cutoff": float(args.rad_cutoff),
                "device": args.device,
                "stagate_version": getattr(st, "__version__", "unknown"),
            },
            handle,
            indent=2,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
