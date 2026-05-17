#!/usr/bin/env python3
"""SpaGCN clustering sidecar for spatial module 03.

Input is an h5ad file with ``obsm["spatial"]`` coordinates. The R caller treats
any non-zero exit as a non-blocking ``failed_py_bridge`` backend status.
"""

from __future__ import annotations

import argparse
import json
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
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--beta", type=int, default=49)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    try:
        import anndata
        import numpy as np
        import SpaGCN as spg
        import torch
    except Exception as exc:  # pragma: no cover - exercised by smoke via exit code
        message = f"SpaGCN bridge dependency import failed: {exc}"
        _write_error(args.output_metadata_json, message)
        print(message, file=sys.stderr)
        return 2

    try:
        adata = anndata.read_h5ad(args.input)
        if "spatial" not in adata.obsm:
            raise ValueError("input h5ad is missing obsm['spatial']")
        coords = adata.obsm["spatial"]
        adj = spg.calculate_adj_matrix(x=coords[:, 0], y=coords[:, 1], histology=False)
        l_value = spg.search_l(p=0.5, adj=adj)
        resolution = spg.search_res(
            adata,
            adj,
            l_value,
            args.target_clusters,
            start=0.7,
            step=0.1,
            tol=5e-3,
            lr=0.05,
            max_epochs=20,
            r_seed=args.seed,
            t_seed=args.seed,
            n_seed=args.seed,
        )
        clf = spg.SpaGCN()
        clf.set_l(l_value)
        np.random.seed(args.seed)
        torch.manual_seed(args.seed)
        clf.train(
            adata,
            adj,
            init_spa=True,
            init="louvain",
            res=resolution,
            tol=5e-3,
            lr=0.05,
            max_epochs=200,
        )
        predicted, _ = clf.predict()
    except Exception as exc:  # pragma: no cover
        message = f"SpaGCN bridge run failed: {exc}"
        _write_error(args.output_metadata_json, message)
        print(message, file=sys.stderr)
        return 3

    Path(args.output_cluster_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_cluster_tsv, "w", encoding="utf-8") as handle:
        handle.write("barcode\tcluster\n")
        for barcode, cluster in zip(adata.obs_names, predicted):
            handle.write(f"{barcode}\t{cluster}\n")

    with open(args.output_metadata_json, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "status": "ok",
                "n_obs": int(adata.n_obs),
                "n_clusters": int(len(set(predicted))),
                "target_clusters": int(args.target_clusters),
                "resolution_used": float(resolution),
                "l_used": float(l_value),
                "spagcn_version": getattr(spg, "__version__", "unknown"),
            },
            handle,
            indent=2,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
