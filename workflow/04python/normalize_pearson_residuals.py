#!/usr/bin/env python3
"""Compute Pearson-residual normalized spatial counts for the R ST pipeline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scanpy as sc
from scipy import io as scipy_io
from scipy import sparse


def read_lines(path: Path) -> list[str]:
    return [line.rstrip("\n") for line in path.read_text(encoding="utf-8").splitlines() if line.rstrip("\n")]


def build_anndata_from_mtx(counts_mtx: Path, features: Path, barcodes: Path):
    import anndata as ad

    counts = scipy_io.mmread(str(counts_mtx))
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    gene_names = read_lines(features)
    spot_names = read_lines(barcodes)
    adata = ad.AnnData(X=counts.T.tocsr())
    adata.var_names = gene_names
    adata.obs_names = spot_names
    return adata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="Optional input h5ad path.")
    parser.add_argument("--output", help="Optional output h5ad path.")
    parser.add_argument("--counts-mtx", help="Gene x spot Matrix Market counts path.")
    parser.add_argument("--features", help="Feature name file, one feature per line.")
    parser.add_argument("--barcodes", help="Barcode file, one barcode per line.")
    parser.add_argument("--output-mtx", help="Gene x spot Matrix Market residual output path.")
    parser.add_argument("--metadata-json", help="Optional metadata JSON output path.")
    parser.add_argument("--clip", type=float, default=None)
    args = parser.parse_args()

    if args.input:
        adata = sc.read_h5ad(args.input)
    else:
        required = [args.counts_mtx, args.features, args.barcodes]
        if any(value is None for value in required):
            parser.error("provide either --input or --counts-mtx/--features/--barcodes")
        adata = build_anndata_from_mtx(Path(args.counts_mtx), Path(args.features), Path(args.barcodes))

    kwargs = {}
    if args.clip is not None:
        kwargs["clip"] = args.clip
    sc.experimental.pp.normalize_pearson_residuals(adata, inplace=True, **kwargs)
    if sparse.issparse(adata.X):
        residuals = adata.X.T.tocsr()
    else:
        residuals = sparse.csr_matrix(np.asarray(adata.X).T)

    if args.output_mtx:
        out_mtx = Path(args.output_mtx)
        out_mtx.parent.mkdir(parents=True, exist_ok=True)
        scipy_io.mmwrite(str(out_mtx), residuals)
    if args.output:
        out_h5ad = Path(args.output)
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out_h5ad)
    if args.metadata_json:
        meta = {
            "scanpy_version": sc.__version__,
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "clip": args.clip,
            "input_mode": "h5ad" if args.input else "matrix_market",
        }
        Path(args.metadata_json).write_text(json.dumps(meta, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
