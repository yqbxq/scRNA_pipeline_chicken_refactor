#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

from helpers.scrna_io import resolve_scrna_h5ad_path
from helpers.spatial_io import resolve_spatial_h5ad_path


def check_environment(use_gpu: bool) -> tuple[str, str]:
    try:
        import cell2location  # noqa: F401
        import anndata  # noqa: F401
        import scanpy  # noqa: F401
        import torch
    except Exception as exc:
        return "skipped_no_python_env", f"cell2location Python stack import failed: {exc}"
    if use_gpu and not torch.cuda.is_available():
        return "skipped_no_gpu", "SPATIAL_C2L_USE_GPU requested but torch.cuda.is_available() is false"
    return "ok", ""


def _train_model(model, max_epochs: int, use_gpu: bool) -> None:
    try:
        model.train(max_epochs=max_epochs, use_gpu=use_gpu)
        return
    except TypeError:
        pass
    accelerator = "gpu" if use_gpu else "cpu"
    model.train(max_epochs=max_epochs, accelerator=accelerator)


def _export_posterior(model, adata, use_gpu: bool):
    kwargs = {"num_samples": 1000, "batch_size": 2500}
    try:
        exported = model.export_posterior(adata, sample_kwargs={**kwargs, "use_gpu": use_gpu})
    except TypeError:
        exported = model.export_posterior(adata, sample_kwargs=kwargs)
    return adata if exported is None else exported


def _cell_state_df(adata_ref) -> pd.DataFrame:
    import pandas as pd

    key = "means_per_cluster_mu_fg"
    if key not in adata_ref.varm:
        raise KeyError(f"reference posterior missing adata_ref.varm[{key!r}]")
    df = pd.DataFrame(adata_ref.varm[key], index=adata_ref.var_names)
    factor_names = []
    try:
        factor_names = list(adata_ref.uns.get("mod", {}).get("factor_names", []))
    except Exception:
        factor_names = []
    if factor_names:
        cols = [f"{key}_{name}" for name in factor_names]
        present = [col for col in cols if col in df.columns]
        if present:
            df = df[present]
            df.columns = factor_names[: len(present)]
        elif len(factor_names) == df.shape[1]:
            df.columns = factor_names
    df = df.loc[:, df.sum(axis=0) > 0]
    if df.shape[1] == 0:
        raise ValueError("reference posterior has zero nonzero cell-state factors")
    return df


def _abundance_frame(adata_st) -> tuple[str, pd.DataFrame]:
    import pandas as pd

    preferred = ["q05_cell_abundance_w_sf", "means_cell_abundance_w_sf"]
    keys = [key for key in preferred if key in adata_st.obsm]
    keys.extend([key for key in adata_st.obsm.keys() if "cell_abundance" in key and key not in keys])
    if not keys:
        raise KeyError("spatial posterior missing cell abundance matrix in adata_st.obsm")
    key = keys[0]
    value = adata_st.obsm[key]
    if isinstance(value, pd.DataFrame):
        frame = value.copy()
    else:
        frame = pd.DataFrame(value, index=adata_st.obs_names)
    if frame.shape[1] == 0:
        raise ValueError(f"{key} has zero columns")
    frame.index = adata_st.obs_names
    frame.columns = [str(col).replace("q05cell_abundance_w_sf_", "").replace("meanscell_abundance_w_sf_", "") for col in frame.columns]
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    frame[frame < 0] = 0.0
    row_sum = frame.sum(axis=1)
    keep = row_sum > 0
    frame.loc[keep, :] = frame.loc[keep, :].div(row_sum[keep], axis=0)
    frame = frame.loc[keep, frame.sum(axis=0) > 0]
    if frame.empty:
        raise ValueError(f"{key} normalized to an empty spot-celltype matrix")
    return key, frame


def _write_outputs(prop: pd.DataFrame, adata_st, out_dir: str | Path, abundance_key: str, runtime_sec: float | None = None) -> None:
    import numpy as np
    import pandas as pd

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    wide = prop.reset_index().rename(columns={"index": "spot_id"})
    wide.to_csv(out / "spot_celltype_proportions_wide.tsv", sep="\t", index=False)

    long = wide.melt(id_vars=["spot_id"], var_name="cell_type", value_name="proportion")
    long.to_csv(out / "spot_celltype_proportions.tsv", sep="\t", index=False)

    arr = prop.to_numpy(dtype=float)
    entropy = -np.sum(np.where(arr > 0, arr * np.log(arr), 0.0), axis=1)
    dominant = prop.columns[np.argmax(arr, axis=1)]
    spot_meta = pd.DataFrame(
        {
            "spot_id": prop.index,
            "dominant_celltype": dominant,
            "mixing_entropy": entropy,
        }
    )
    for col in ("section_id", "sample_id"):
        if col in adata_st.obs.columns:
            spot_meta[col] = adata_st.obs.loc[prop.index, col].astype(str).to_numpy()
    spot_meta.to_csv(out / "spot_metadata.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "tool": "cell2location",
                "status": "ok",
                "reason": "",
                "n_spots": prop.shape[0],
                "n_celltypes": prop.shape[1],
                "abundance_key": abundance_key,
                "runtime_sec": "" if runtime_sec is None else runtime_sec,
            }
        ]
    )
    summary.to_csv(out / "method_summary.tsv", sep="\t", index=False)


def run_cell2location(args) -> tuple[int, str]:
    import anndata as ad
    from cell2location.models import Cell2location, RegressionModel

    ref_h5ad = resolve_scrna_h5ad_path(module=args.ref_module, fallback_path=args.ref_h5ad)
    st_h5ad = resolve_spatial_h5ad_path(module=args.st_module, section_id=args.section_id or None, fallback_path=args.st_h5ad)
    adata_ref = ad.read_h5ad(ref_h5ad)
    adata_st = ad.read_h5ad(st_h5ad)

    if args.annotation_col not in adata_ref.obs.columns:
        raise ValueError(f"reference H5AD obs missing annotation column: {args.annotation_col}")
    if args.section_id and "section_id" in adata_st.obs.columns:
        keep = adata_st.obs["section_id"].astype(str) == str(args.section_id)
        if int(keep.sum()) == 0:
            raise ValueError(f"spatial H5AD has no spots for section_id={args.section_id}")
        adata_st = adata_st[keep].copy()
    if adata_ref.n_obs < 2 or adata_st.n_obs < 2:
        raise ValueError("cell2location needs at least two reference cells and two spatial spots")

    RegressionModel.setup_anndata(adata_ref, labels_key=args.annotation_col)
    ref_model = RegressionModel(adata_ref)
    _train_model(ref_model, args.ref_epochs, args.use_gpu)
    adata_ref = _export_posterior(ref_model, adata_ref, args.use_gpu)
    cell_state = _cell_state_df(adata_ref)

    shared = sorted(set(adata_st.var_names).intersection(cell_state.index))
    if len(shared) < args.min_shared_genes:
        raise ValueError(f"only {len(shared)} shared genes between reference and spatial H5AD; minimum is {args.min_shared_genes}")
    adata_st = adata_st[:, shared].copy()
    cell_state = cell_state.loc[shared, :]

    batch_key = args.batch_key if args.batch_key and args.batch_key in adata_st.obs.columns else None
    Cell2location.setup_anndata(adata_st, batch_key=batch_key)
    st_model = Cell2location(
        adata_st,
        cell_state_df=cell_state,
        N_cells_per_location=args.n_cells_per_location,
        detection_alpha=args.detection_alpha,
    )
    _train_model(st_model, args.st_epochs, args.use_gpu)
    adata_st = _export_posterior(st_model, adata_st, args.use_gpu)
    abundance_key, prop = _abundance_frame(adata_st)
    _write_outputs(prop, adata_st, args.out_dir, abundance_key)
    try:
        adata_st.write_h5ad(Path(args.out_dir) / "cell2location_result.h5ad")
    except Exception as exc:
        sys.stderr.write(f"warning\tcell2location_result.h5ad write failed: {exc}\n")
    return 0, f"ok\tspots={prop.shape[0]}\tcelltypes={prop.shape[1]}\tabundance_key={abundance_key}"


def main() -> int:
    parser = argparse.ArgumentParser(description="cell2location sidecar for ST07.")
    parser.add_argument("--check-env", action="store_true")
    parser.add_argument("--use-gpu", action="store_true")
    parser.add_argument("--ref-h5ad")
    parser.add_argument("--st-h5ad")
    parser.add_argument("--ref-module", default="03d_panorama")
    parser.add_argument("--st-module", default="spatial_03_region")
    parser.add_argument("--section-id", default="")
    parser.add_argument("--annotation-col")
    parser.add_argument("--out-dir")
    parser.add_argument("--ref-epochs", type=int, default=250)
    parser.add_argument("--st-epochs", type=int, default=5000)
    parser.add_argument("--batch-key", default="section_id")
    parser.add_argument("--n-cells-per-location", type=float, default=30.0)
    parser.add_argument("--detection-alpha", type=float, default=20.0)
    parser.add_argument("--min-shared-genes", type=int, default=100)
    args = parser.parse_args()

    status, reason = check_environment(args.use_gpu)
    if args.check_env:
        if status == "ok":
            print("ok")
            return 0
        sys.stderr.write(f"{status}\t{reason}\n")
        return 20 if status == "skipped_no_python_env" else 21

    if status != "ok":
        sys.stderr.write(f"{status}\t{reason}\n")
        return 20 if status == "skipped_no_python_env" else 21

    missing = [name for name in ("ref_h5ad", "st_h5ad", "annotation_col", "out_dir") if not getattr(args, name)]
    if missing:
        optional_resolvers = {"ref_h5ad", "st_h5ad"}
        required_missing = [name for name in missing if name not in optional_resolvers]
        if required_missing:
            sys.stderr.write(f"failed_sidecar\tmissing required runtime args: {','.join(required_missing)}\n")
            return 30

    try:
        code, message = run_cell2location(args)
    except FileNotFoundError as exc:
        sys.stderr.write(f"failed_sidecar\tH5AD input resolution failed: {exc}\n")
        return 30
    except Exception as exc:
        sys.stderr.write(f"failed_sidecar\tcell2location training failed: {exc}\n")
        return 31
    print(message)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
