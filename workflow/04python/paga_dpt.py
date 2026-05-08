#!/usr/bin/env python3

import argparse
from pathlib import Path


def read_lines(path):
    p = Path(path)
    if not p.exists():
        return []
    return [x.strip() for x in p.read_text(encoding="utf-8").splitlines() if x.strip()]


def empty_outputs(job, status, reason):
    import pandas as pd

    pd.DataFrame(columns=["cell_id", "pseudotime", "label"]).to_csv(job["pseudotime_csv"], index=False)
    pd.DataFrame(columns=["source", "target", "connectivity"]).to_csv(job["connectivity_tsv"], sep="\t", index=False)
    Path(job["status"]).write_text(f"{status}\t{reason}\n", encoding="utf-8")


def run_job(job):
    try:
        import numpy as np
        import pandas as pd
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import scanpy as sc
        from scipy import io
        from scipy import sparse
    except Exception as exc:
        empty_outputs(job, "skipped_package_missing", str(exc))
        return {"status": "skipped_package_missing", "reason": str(exc)}

    input_dir = Path(job["input_dir"])
    try:
        counts = io.mmread(input_dir / "counts.mtx")
        if not sparse.issparse(counts):
            counts = sparse.csr_matrix(counts)
        counts = counts.T.tocsr()
        features = read_lines(input_dir / "features.tsv")
        barcodes = read_lines(input_dir / "barcodes.tsv")
        meta = pd.read_csv(input_dir / "metadata.tsv", sep="\t")
        if "cell_id" not in meta.columns:
            meta["cell_id"] = barcodes
        meta = meta.set_index("cell_id").reindex(barcodes)
        adata = sc.AnnData(X=counts, obs=meta)
        adata.var_names = features
        adata.obs_names = barcodes

        umap_path = input_dir / "umap.tsv"
        if umap_path.exists():
            umap = pd.read_csv(umap_path, sep="\t").set_index("cell_id").reindex(barcodes)
            adata.obsm["X_umap"] = umap[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)

        label_var = job["label_var"]
        if label_var not in adata.obs.columns:
            raise ValueError(f"label_var missing from metadata: {label_var}")
        adata.obs[label_var] = adata.obs[label_var].astype("category")
        if adata.obs[label_var].nunique() < 2:
            raise ValueError(f"label_var has fewer than 2 categories: {label_var}")

        if "X_umap" in adata.obsm:
            sc.pp.neighbors(adata, use_rep="X_umap", n_neighbors=min(30, max(2, adata.n_obs - 1)))
        else:
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            sc.pp.pca(adata, n_comps=min(30, max(2, adata.n_obs - 1)))
            sc.pp.neighbors(adata, n_pcs=min(30, max(2, adata.n_obs - 1)))
            sc.tl.umap(adata)

        sc.tl.paga(adata, groups=label_var)
        sc.tl.diffmap(adata)
        root_cells = read_lines(input_dir / "root_cells.txt")
        root_cells = [x for x in root_cells if x in adata.obs_names]
        if root_cells:
            adata.uns["iroot"] = int(np.where(adata.obs_names == root_cells[0])[0][0])
        else:
            adata.uns["iroot"] = 0
        sc.tl.dpt(adata)

        pst = np.asarray(adata.obs["dpt_pseudotime"], dtype=float)
        pst[~np.isfinite(pst)] = np.nan
        out = pd.DataFrame({
            "cell_id": adata.obs_names,
            "pseudotime": pst,
            "label": adata.obs[label_var].astype(str).to_numpy(),
        })
        out.to_csv(job["pseudotime_csv"], index=False)

        conn = adata.uns["paga"]["connectivities"]
        cats = list(adata.obs[label_var].cat.categories)
        coo = conn.tocoo()
        conn_df = pd.DataFrame({
            "source": [cats[i] for i in coo.row],
            "target": [cats[i] for i in coo.col],
            "connectivity": coo.data,
        })
        conn_df.to_csv(job["connectivity_tsv"], sep="\t", index=False)
        adata.write_h5ad(job["h5ad_path"])

        fig_path = Path(job["figure_png"])
        fig_path.parent.mkdir(parents=True, exist_ok=True)
        if "X_umap" in adata.obsm:
            xy = adata.obsm["X_umap"]
            plt.figure(figsize=(6, 4.8))
            plt.scatter(xy[:, 0], xy[:, 1], c=pst, s=2, cmap="viridis")
            plt.colorbar(label="DPT pseudotime")
            plt.xlabel("UMAP 1")
            plt.ylabel("UMAP 2")
            plt.title(f"PAGA-DPT {job['pair_id']} {job['split_value']}")
            plt.tight_layout()
            plt.savefig(fig_path, dpi=300)
            plt.close()

        Path(job["status"]).write_text("ok\t\n", encoding="utf-8")
        return {"status": "ok", "reason": ""}
    except Exception as exc:
        empty_outputs(job, "failed", str(exc))
        return {"status": "failed", "reason": str(exc)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", required=True)
    parser.add_argument("--index", required=True)
    args = parser.parse_args()

    import pandas as pd

    jobs = pd.read_csv(args.jobs, sep="\t")
    rows = []
    for _, job_row in jobs.iterrows():
        job = {k: "" if pd.isna(v) else str(v) for k, v in job_row.to_dict().items()}
        result = run_job(job)
        rows.append({
            "pair_id": job["pair_id"],
            "split_value": job["split_value"],
            "method": "paga_dpt",
            "methods_enabled": "yes",
            "input_rds": job["input_rds"],
            "output_path": job["pseudotime_csv"],
            "extra_path": job["connectivity_tsv"],
            "figure_path": job["figure_png"],
            "n_cells": job.get("n_cells", ""),
            "status": result["status"],
            "reason": result["reason"],
            "runtime_s": "",
            "h5ad_path": job["h5ad_path"],
        })
    pd.DataFrame(rows).to_csv(args.index, sep="\t", index=False)


if __name__ == "__main__":
    main()
