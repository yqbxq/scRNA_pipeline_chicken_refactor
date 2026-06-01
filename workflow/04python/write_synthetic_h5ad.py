#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys


def _read_table(path: Path):
    import pandas as pd

    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _fingerprint(path: Path) -> str:
    import hashlib

    if not path.exists():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_manifest(rows: list[dict[str, object]], path: Path) -> None:
    import pandas as pd

    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _load_mtx(mtx: Path, genes: Path, obs: Path, obs_id_col: str):
    import pandas as pd
    import scipy.io

    matrix = scipy.io.mmread(mtx).tocsr().T
    gene_df = pd.read_csv(genes, sep="\t", dtype=str).fillna("")
    obs_df = pd.read_csv(obs, sep="\t", dtype=str).fillna("")
    gene_col = "gene_id" if "gene_id" in gene_df.columns else gene_df.columns[0]
    if obs_id_col not in obs_df.columns:
        raise ValueError(f"{obs} missing required column {obs_id_col}")
    var_names = gene_df[gene_col].astype(str).tolist()
    obs_names = obs_df[obs_id_col].astype(str).tolist()
    if matrix.shape != (len(obs_names), len(var_names)):
        raise ValueError(f"matrix shape {matrix.shape} does not match obs={len(obs_names)} var={len(var_names)}")
    obs_df.index = obs_names
    var_df = gene_df.copy()
    var_df.index = var_names
    return matrix, obs_df, var_df


def main() -> int:
    parser = argparse.ArgumentParser(description="Write ST07f synthetic reference/spatial H5AD artifacts.")
    parser.add_argument("--reference-mtx", required=True)
    parser.add_argument("--reference-genes", required=True)
    parser.add_argument("--reference-obs", required=True)
    parser.add_argument("--spatial-mtx", required=True)
    parser.add_argument("--spatial-genes", required=True)
    parser.add_argument("--spatial-obs", required=True)
    parser.add_argument("--spatial-coords", required=True)
    parser.add_argument("--truth-tsv", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--validation-id", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "synthetic_h5ad_manifest.tsv"
    ref_h5ad = out_dir / "synthetic_reference.h5ad"
    st_h5ad = out_dir / "synthetic_spatial.h5ad"

    try:
        import anndata as ad
        import numpy as np
    except Exception as exc:
        rows = []
        for role, h5ad in (("synthetic_reference", ref_h5ad), ("synthetic_spatial", st_h5ad)):
            rows.append(
                {
                    "artifact_id": f"{args.validation_id}:{role}",
                    "artifact_role": role,
                    "h5ad_path": str(h5ad),
                    "source_object": "scDesign3_synthetic",
                    "n_obs": 0,
                    "n_vars": 0,
                    "obs_required_cols": "synthetic_cell_id,cell_type" if role == "synthetic_reference" else "spot_id,section_id,validation_id",
                    "obsm_required_keys": "" if role == "synthetic_reference" else "spatial",
                    "var_required_cols": "gene_id,gene_name",
                    "fingerprint": "",
                    "status": "skipped_no_python_env",
                    "reason": f"anndata/scipy stack unavailable: {exc}",
                }
            )
        _write_manifest(rows, manifest_path)
        return 20

    try:
        ref_x, ref_obs, ref_var = _load_mtx(Path(args.reference_mtx), Path(args.reference_genes), Path(args.reference_obs), "synthetic_cell_id")
        st_x, st_obs, st_var = _load_mtx(Path(args.spatial_mtx), Path(args.spatial_genes), Path(args.spatial_obs), "spot_id")
        coords = _read_table(Path(args.spatial_coords))
        if not {"spot_id", "x", "y"}.issubset(coords.columns):
            raise ValueError("spatial coordinates table requires spot_id, x, y")
        coords = coords.set_index("spot_id").loc[st_obs.index, ["x", "y"]].astype(float)
        ref_obs["validation_id"] = args.validation_id
        st_obs["validation_id"] = args.validation_id
        truth_tsv = str(Path(args.truth_tsv))

        ref = ad.AnnData(X=ref_x, obs=ref_obs, var=ref_var)
        ref.uns["source"] = "scDesign3"
        ref.uns["validation_id"] = args.validation_id
        ref.write_h5ad(ref_h5ad)

        st = ad.AnnData(X=st_x, obs=st_obs, var=st_var)
        st.obsm["spatial"] = coords.to_numpy(dtype=float)
        st.uns["source"] = "scDesign3_mixed_spots"
        st.uns["validation_id"] = args.validation_id
        st.uns["synthetic_truth_tsv"] = truth_tsv
        st.write_h5ad(st_h5ad)
    except Exception as exc:
        rows = []
        for role, h5ad in (("synthetic_reference", ref_h5ad), ("synthetic_spatial", st_h5ad)):
            rows.append(
                {
                    "artifact_id": f"{args.validation_id}:{role}",
                    "artifact_role": role,
                    "h5ad_path": str(h5ad),
                    "source_object": "scDesign3_synthetic",
                    "n_obs": 0,
                    "n_vars": 0,
                    "obs_required_cols": "synthetic_cell_id,cell_type" if role == "synthetic_reference" else "spot_id,section_id,validation_id",
                    "obsm_required_keys": "" if role == "synthetic_reference" else "spatial",
                    "var_required_cols": "gene_id,gene_name",
                    "fingerprint": "",
                    "status": "failed_write_h5ad",
                    "reason": str(exc),
                }
            )
        _write_manifest(rows, manifest_path)
        return 30

    rows = [
        {
            "artifact_id": f"{args.validation_id}:synthetic_reference",
            "artifact_role": "synthetic_reference",
            "h5ad_path": str(ref_h5ad),
            "source_object": "scDesign3_synthetic_cells",
            "n_obs": ref.n_obs,
            "n_vars": ref.n_vars,
            "obs_required_cols": "synthetic_cell_id,cell_type,validation_id",
            "obsm_required_keys": "",
            "var_required_cols": "gene_id,gene_name",
            "fingerprint": _fingerprint(ref_h5ad),
            "status": "ok",
            "reason": "",
        },
        {
            "artifact_id": f"{args.validation_id}:synthetic_spatial",
            "artifact_role": "synthetic_spatial",
            "h5ad_path": str(st_h5ad),
            "source_object": "scDesign3_mixed_spots",
            "n_obs": st.n_obs,
            "n_vars": st.n_vars,
            "obs_required_cols": "spot_id,section_id,validation_id",
            "obsm_required_keys": "spatial",
            "var_required_cols": "gene_id,gene_name",
            "fingerprint": _fingerprint(st_h5ad),
            "status": "ok",
            "reason": "",
        },
    ]
    _write_manifest(rows, manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
