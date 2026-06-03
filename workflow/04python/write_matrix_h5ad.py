#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def _fingerprint(path: Path) -> str:
    if not path.exists():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_manifest(path: Path, row: dict[str, object]) -> None:
    import pandas as pd

    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(path, sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Write a simple count-matrix H5AD mirror.")
    parser.add_argument("--mtx", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--obs", required=True)
    parser.add_argument("--obs-id-col", required=True)
    parser.add_argument("--out-h5ad", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--artifact-id", required=True)
    parser.add_argument("--artifact-role", required=True)
    parser.add_argument("--source-object", required=True)
    args = parser.parse_args()

    out_h5ad = Path(args.out_h5ad)
    manifest = Path(args.manifest)
    try:
        import anndata as ad
        import pandas as pd
        import scipy.io
    except Exception as exc:
        _write_manifest(manifest, {
            "artifact_id": args.artifact_id,
            "artifact_role": args.artifact_role,
            "h5ad_path": str(out_h5ad),
            "source_object": args.source_object,
            "n_obs": 0,
            "n_vars": 0,
            "fingerprint": "",
            "status": "skipped_no_python_env",
            "reason": f"anndata/pandas/scipy unavailable: {exc}",
        })
        return 20

    try:
        matrix = scipy.io.mmread(args.mtx).tocsr().T
        obs = pd.read_csv(args.obs, sep="\t", dtype=str).fillna("")
        genes = pd.read_csv(args.genes, sep="\t", dtype=str).fillna("")
        if args.obs_id_col not in obs.columns:
            raise ValueError(f"obs missing id column: {args.obs_id_col}")
        gene_col = "gene_id" if "gene_id" in genes.columns else genes.columns[0]
        obs.index = obs[args.obs_id_col].astype(str)
        genes.index = genes[gene_col].astype(str)
        if matrix.shape != (obs.shape[0], genes.shape[0]):
            raise ValueError(f"matrix shape {matrix.shape} does not match obs={obs.shape[0]} var={genes.shape[0]}")
        adata = ad.AnnData(X=matrix, obs=obs, var=genes)
        adata.layers["counts"] = matrix.copy()
        adata.uns["module"] = args.source_object
        adata.uns["export_timestamp"] = pd.Timestamp.utcnow().isoformat()
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out_h5ad)
    except Exception as exc:
        _write_manifest(manifest, {
            "artifact_id": args.artifact_id,
            "artifact_role": args.artifact_role,
            "h5ad_path": str(out_h5ad),
            "source_object": args.source_object,
            "n_obs": 0,
            "n_vars": 0,
            "fingerprint": "",
            "status": "failed_write_h5ad",
            "reason": str(exc),
        })
        return 30

    _write_manifest(manifest, {
        "artifact_id": args.artifact_id,
        "artifact_role": args.artifact_role,
        "h5ad_path": str(out_h5ad),
        "source_object": args.source_object,
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "fingerprint": _fingerprint(out_h5ad),
        "status": "ok",
        "reason": "",
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
