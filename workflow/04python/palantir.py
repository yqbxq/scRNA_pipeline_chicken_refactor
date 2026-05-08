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

    pd.DataFrame(columns=["cell_id", "pseudotime"]).to_csv(job["pseudotime_csv"], index=False)
    pd.DataFrame(columns=["cell_id", "entropy"]).to_csv(job["entropy_csv"], index=False)
    pd.DataFrame(columns=["cell_id", "terminal_state", "probability"]).to_csv(job["fate_csv"], index=False)
    Path(job["status"]).write_text(f"{status}\t{reason}\n", encoding="utf-8")


def entropy_from_probs(prob):
    import numpy as np

    prob = np.asarray(prob, dtype=float)
    prob = np.clip(prob, 1e-12, 1.0)
    return -np.sum(prob * np.log(prob), axis=1)


def run_scanpy_fallback(job, package_status, package_reason):
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
        empty_outputs(job, "skipped_package_missing", f"{package_reason}; scanpy fallback unavailable: {exc}")
        return {"status": "skipped_package_missing", "reason": str(exc)}

    input_dir = Path(job["input_dir"])
    counts = io.mmread(input_dir / "counts.mtx")
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    counts = counts.T.tocsr()
    features = read_lines(input_dir / "features.tsv")
    barcodes = read_lines(input_dir / "barcodes.tsv")
    meta = pd.read_csv(input_dir / "metadata.tsv", sep="\t").set_index("cell_id").reindex(barcodes)
    adata = sc.AnnData(X=counts, obs=meta)
    adata.var_names = features
    adata.obs_names = barcodes
    if adata.n_obs < 200:
        empty_outputs(job, "skipped_low_data", f"Palantir requires >=200 cells by policy; observed {adata.n_obs}")
        return {"status": "skipped_low_data", "reason": f"n_cells={adata.n_obs}"}

    terminal_cells = read_lines(input_dir / "terminal_cells.txt")
    terminal_cells = [x for x in terminal_cells if x in adata.obs_names]
    if not terminal_cells:
        empty_outputs(job, "skipped_no_terminal_states", "no terminal cells resolved from terminal_group")
        return {"status": "skipped_no_terminal_states", "reason": "no terminal cells"}

    umap_path = input_dir / "umap.tsv"
    if umap_path.exists():
        umap = pd.read_csv(umap_path, sep="\t").set_index("cell_id").reindex(barcodes)
        adata.obsm["X_umap"] = umap[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)
        sc.pp.neighbors(adata, use_rep="X_umap", n_neighbors=min(30, max(2, adata.n_obs - 1)))
    else:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata, n_comps=min(30, max(2, adata.n_obs - 1)))
        sc.pp.neighbors(adata, n_pcs=min(30, max(2, adata.n_obs - 1)))
        sc.tl.umap(adata)

    sc.tl.diffmap(adata)
    root_cells = read_lines(input_dir / "root_cells.txt")
    root_cells = [x for x in root_cells if x in adata.obs_names]
    adata.uns["iroot"] = int(np.where(adata.obs_names == (root_cells[0] if root_cells else adata.obs_names[0]))[0][0])
    sc.tl.dpt(adata)
    pst = np.asarray(adata.obs["dpt_pseudotime"], dtype=float)
    pd.DataFrame({"cell_id": adata.obs_names, "pseudotime": pst}).to_csv(job["pseudotime_csv"], index=False)

    label_var = job["label_var"]
    if label_var in adata.obs.columns:
        terminal_labels = sorted(set(adata.obs.loc[terminal_cells, label_var].astype(str)))
    else:
        terminal_labels = ["terminal"]
    probs = np.zeros((adata.n_obs, len(terminal_labels)), dtype=float)
    for j, label in enumerate(terminal_labels):
        if label_var in adata.obs.columns:
            target_idx = np.where(adata.obs[label_var].astype(str).to_numpy() == label)[0]
        else:
            target_idx = np.array([adata.obs_names.get_loc(x) for x in terminal_cells])
        target_pt = np.nanmedian(pst[target_idx]) if len(target_idx) else np.nanmax(pst)
        probs[:, j] = 1.0 / (np.abs(pst - target_pt) + 1e-3)
    probs = probs / probs.sum(axis=1, keepdims=True)
    fate = []
    for i, cell in enumerate(adata.obs_names):
        for j, label in enumerate(terminal_labels):
            fate.append({"cell_id": cell, "terminal_state": label, "probability": probs[i, j]})
    pd.DataFrame(fate).to_csv(job["fate_csv"], index=False)
    pd.DataFrame({"cell_id": adata.obs_names, "entropy": entropy_from_probs(probs)}).to_csv(job["entropy_csv"], index=False)

    fig_path = Path(job["figure_png"])
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    xy = adata.obsm["X_umap"]
    plt.figure(figsize=(6, 4.8))
    plt.scatter(xy[:, 0], xy[:, 1], c=pst, s=2, cmap="viridis")
    plt.colorbar(label="Palantir/DPT pseudotime")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title(f"Palantir fallback {job['pair_id']} {job['split_value']}")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close()
    Path(job["status"]).write_text(f"{package_status}\t{package_reason}\n", encoding="utf-8")
    return {"status": package_status, "reason": package_reason}


def run_job(job):
    try:
        import palantir  # noqa: F401
    except Exception as exc:
        return run_scanpy_fallback(job, "skipped_package_missing", f"palantir package unavailable: {exc}")

    # Palantir APIs have changed across releases. Until a project env pins a
    # known-compatible version, use the scanpy diffusion fallback but mark the
    # status so downstream reports do not over-interpret it as a full Palantir run.
    return run_scanpy_fallback(job, "fallback_scanpy_dpt", "palantir import succeeded; stable wrapper not pinned")


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
        try:
            result = run_job(job)
        except Exception as exc:
            empty_outputs(job, "failed", str(exc))
            result = {"status": "failed", "reason": str(exc)}
        rows.append({
            "pair_id": job["pair_id"],
            "split_value": job["split_value"],
            "method": "palantir",
            "methods_enabled": "yes",
            "input_rds": job["input_rds"],
            "output_path": job["pseudotime_csv"],
            "extra_path": job["entropy_csv"],
            "figure_path": job["figure_png"],
            "n_cells": job.get("n_cells", ""),
            "status": result["status"],
            "reason": result["reason"],
            "runtime_s": "",
            "fate_path": job["fate_csv"],
        })
    pd.DataFrame(rows).to_csv(args.index, sep="\t", index=False)


if __name__ == "__main__":
    main()
