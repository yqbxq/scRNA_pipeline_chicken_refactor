#!/usr/bin/env python3
import argparse
import sys


def fail(code, status, reason):
    sys.stderr.write(f"{status}\t{reason}\n")
    return code


def import_runtime():
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import squidpy as sq
    except Exception as exc:
        raise RuntimeError(f"anndata/numpy/pandas/squidpy import failed: {exc}") from exc
    return ad, np, pd, sq


def read_spots_h5ad(path, group_by, adata_mod, np_mod, pd_mod):
    adata = adata_mod.read_h5ad(path)
    if group_by not in adata.obs:
        raise ValueError(f"h5ad obs missing group_by column: {group_by}")
    if "spatial" not in adata.obsm:
        raise ValueError('h5ad missing adata.obsm["spatial"]')
    coords = np_mod.asarray(adata.obsm["spatial"], dtype=float)
    labels = adata.obs[group_by].astype(str)
    keep = pd_mod.Series(labels.to_numpy(), index=adata.obs_names).str.len().gt(0).to_numpy()
    keep = keep & np_mod.isfinite(coords).all(axis=1)
    adata = adata[keep].copy()
    adata.obs[group_by] = adata.obs[group_by].astype(str).astype("category")
    return adata


def read_spots_tsv(path, adata_mod, np_mod, pd_mod):
    spots = pd_mod.read_csv(path, sep="\t")
    required = {"spot_id", "x", "y", "label"}
    missing = required.difference(spots.columns)
    if missing:
        raise ValueError(f"spots TSV missing columns: {','.join(sorted(missing))}")
    if "section_id" not in spots.columns:
        spots["section_id"] = ""
    spots = spots.dropna(subset=["spot_id", "x", "y", "label"]).copy()
    spots["label"] = spots["label"].astype(str)
    spots = spots[spots["label"].str.len() > 0].copy()
    coords = spots[["x", "y"]].to_numpy(dtype=float)
    keep = np_mod.isfinite(coords).all(axis=1)
    spots = spots.loc[keep].copy()
    coords = coords[keep]
    obs = spots[["label", "section_id"]].copy()
    obs.index = spots["spot_id"].astype(str)
    obs["label"] = obs["label"].astype("category")
    return adata_mod.AnnData(
        X=np_mod.zeros((len(spots), 1), dtype=float),
        obs=obs,
        var=pd_mod.DataFrame(index=["placeholder_gene"]),
        obsm={"spatial": coords},
    )


def label_categories(adata, group_by):
    values = adata.obs[group_by]
    if hasattr(values, "cat"):
        return [str(x) for x in values.cat.categories]
    return sorted({str(x) for x in values})


def write_matrix(path, matrix, labels, pd_mod, np_mod):
    arr = np_mod.asarray(matrix)
    pd_mod.DataFrame(arr, index=labels, columns=labels).to_csv(path, sep="\t")


def write_co_occurrence(path, co_occurrence_uns, labels, pd_mod, np_mod):
    occ = np_mod.asarray(co_occurrence_uns.get("occ"))
    intervals = co_occurrence_uns.get("interval")
    if intervals is None:
        intervals = list(range(occ.shape[-1] + 1 if occ.ndim > 1 else occ.shape[0] + 1))
    intervals = list(np_mod.asarray(intervals).reshape(-1))

    rows = []
    if occ.ndim == 3:
        for i, region_a in enumerate(labels):
            for j, region_b in enumerate(labels):
                for k in range(occ.shape[2]):
                    rows.append(
                        {
                            "region_a": region_a,
                            "region_b": region_b,
                            "interval_idx": k,
                            "interval_start": intervals[k] if k < len(intervals) else "",
                            "interval_end": intervals[k + 1] if k + 1 < len(intervals) else "",
                            "co_occurrence": occ[i, j, k],
                        }
                    )
    elif occ.ndim == 2:
        for i, region in enumerate(labels[: occ.shape[0]]):
            for k in range(occ.shape[1]):
                rows.append(
                    {
                        "region": region,
                        "interval_idx": k,
                        "interval_start": intervals[k] if k < len(intervals) else "",
                        "interval_end": intervals[k + 1] if k + 1 < len(intervals) else "",
                        "co_occurrence": occ[i, k],
                    }
                )
    else:
        rows.append({"co_occurrence": occ.reshape(-1)[0] if occ.size else ""})
    pd_mod.DataFrame(rows).to_csv(path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="Squidpy spatial neighborhood sidecar for ST06.")
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--h5ad", help='AnnData input with obsm["spatial"] and obs group labels')
    input_group.add_argument("--spots", help="Legacy TSV with spot_id, x, y, label columns")
    parser.add_argument("--group-by", default="label", help="AnnData obs column used as Squidpy cluster_key")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--radius", type=float, default=200.0)
    parser.add_argument("--perms", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--coord-type", default="generic")
    parser.add_argument("--interval", type=int, default=10)
    args = parser.parse_args()

    try:
        adata_mod, np_mod, pd_mod, sq = import_runtime()
    except RuntimeError as exc:
        return fail(20, "skipped_no_packages", str(exc))

    try:
        if args.h5ad:
            adata = read_spots_h5ad(args.h5ad, args.group_by, adata_mod, np_mod, pd_mod)
            group_by = args.group_by
        else:
            adata = read_spots_tsv(args.spots, adata_mod, np_mod, pd_mod)
            group_by = "label"
    except Exception as exc:
        return fail(30, "failed_sidecar", str(exc))

    labels = label_categories(adata, group_by)
    labels = [label for label in labels if label in set(adata.obs[group_by].astype(str))]
    if adata.n_obs < 3 or len(labels) < 2:
        return fail(22, "skipped_too_few_spots", "need at least 3 spots and 2 labels")

    try:
        sq.gr.spatial_neighbors(adata, coord_type=args.coord_type, radius=args.radius)
        sq.gr.nhood_enrichment(adata, cluster_key=group_by, n_perms=max(1, args.perms), seed=args.seed)
        sq.gr.co_occurrence(adata, cluster_key=group_by, interval=args.interval)
        interaction_matrix = sq.gr.interaction_matrix(adata, cluster_key=group_by, copy=True)
    except Exception as exc:
        return fail(31, "failed_sidecar", f"squidpy graph calculation failed: {exc}")

    try:
        nhood_key = f"{group_by}_nhood_enrichment"
        co_occurrence_key = f"{group_by}_co_occurrence"
        zscore = adata.uns[nhood_key]["zscore"]
        co_occurrence_uns = adata.uns[co_occurrence_key]
        if interaction_matrix is None:
            interaction_matrix = adata.uns.get(f"{group_by}_interactions")
        if interaction_matrix is None:
            raise ValueError("Squidpy interaction_matrix output was not found")

        write_matrix(f"{args.out_prefix}_interaction_matrix.tsv", interaction_matrix, labels, pd_mod, np_mod)
        write_matrix(f"{args.out_prefix}_nhood_enrichment_zscore.tsv", zscore, labels, pd_mod, np_mod)
        write_co_occurrence(f"{args.out_prefix}_co_occurrence.tsv", co_occurrence_uns, labels, pd_mod, np_mod)
        pd_mod.DataFrame(
            [
                {
                    "status": "ok",
                    "spot_n": adata.n_obs,
                    "label_n": len(labels),
                    "radius": args.radius,
                    "perms": max(1, args.perms),
                    "coord_type": args.coord_type,
                    "interval": args.interval,
                }
            ]
        ).to_csv(f"{args.out_prefix}_summary.tsv", sep="\t", index=False)
    except Exception as exc:
        return fail(32, "failed_sidecar", f"output serialization failed: {exc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
