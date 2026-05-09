#!/usr/bin/env python3
import argparse
import csv
import json
import math
import os
import re
import time
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D


CELLTYPE_COLOR_POOL = [
    "#79BB7D",
    "#C4A9D6",
    "#FDBE6F",
    "#2F6CB3",
    "#E30073",
    "#A54F37",
    "#69706A",
    "#0E9B84",
    "#D95F02",
    "#6B63B7",
]

CLUSTER_COLORS = {
    "1": "#74BA73",
    "2": "#C8C6E4",
    "3": "#FDBE6F",
    "4": "#2E63AF",
    "5": "#E30073",
    "6": "#A54F37",
    "7": "#69706A",
    "8": "#0E9B84",
    "9": "#D95F02",
    "10": "#6B63B7",
    "11": "#D90445",
    "12": "#1893D1",
    "13": "#FDBA12",
    "14": "#4C2A7A",
    "15": "#D6457A",
}

RAINBOW_CMAP = LinearSegmentedColormap.from_list(
    "blue_to_red_rainbow",
    [
        "#2C7BB6",
        "#00A6CA",
        "#00CCBC",
        "#90EB9D",
        "#FFFF8C",
        "#F9D057",
        "#F29E2E",
        "#E76818",
        "#D7191C",
    ],
)

INDEX_COLUMNS = [
    "pair_id",
    "split_value",
    "method",
    "methods_enabled",
    "input_path",
    "output_path",
    "extra_path",
    "figure_path",
    "n_cells",
    "status",
    "reason",
    "runtime_s",
    "qc_tsv",
    "vector_tsv",
    "stochastic_status",
]

QC_COLUMNS = [
    "pair_id",
    "split_value",
    "n_cells",
    "n_genes",
    "spliced_total",
    "unspliced_total",
    "unspliced_spliced_ratio",
    "fit_likelihood_mean",
    "fit_likelihood_median",
    "fit_likelihood_max",
    "velocity_gene_n",
    "velocity_confidence_mean",
    "velocity_confidence_median",
    "velocity_confidence_q05",
    "velocity_confidence_q95",
    "velocity_length_mean",
    "velocity_confidence_stochastic_mean",
    "velocity_length_stochastic_mean",
    "status",
    "reason",
]

VECTOR_COLUMNS = [
    "pair_id",
    "split_value",
    "cell_id",
    "UMAP_1",
    "UMAP_2",
    "velocity_umap_1",
    "velocity_umap_2",
    "velocity_umap_1_stochastic",
    "velocity_umap_2_stochastic",
    "latent_time",
    "velocity_pseudotime",
    "velocity_length",
    "velocity_confidence",
    "velocity_length_stochastic",
    "velocity_confidence_stochastic",
    "cell_type",
    "cluster",
]


def safe_id(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    return value.strip("_") or "NA"


def display_split(value):
    value = str(value or "").strip()
    return value if value and value != "NA" else "pooled"


def unit_id(pair_id, split_value):
    split_value = display_split(split_value)
    if split_value == "pooled":
        return safe_id(pair_id)
    return f"{safe_id(pair_id)}__{safe_id(split_value)}"


def method_key(value):
    return str(value or "").strip().lower().replace("_", "-")


def token_list(value):
    value = str(value or "").strip()
    if not value:
        return []
    return [part for part in re.split(r"[,;\s]+", value) if part]


def token_matches(tokens, method):
    keys = {method_key(x) for x in tokens}
    key = method_key(method)
    aliases = {key, key.replace("-", "_")}
    return bool(keys.intersection({method_key(x) for x in aliases}))


def method_enabled(row, method, default=True):
    extra = token_list(row.get("methods_extra", ""))
    extra_keys = {method_key(x) for x in extra}
    key = method_key(method)
    if f"-{key}" in extra_keys or f"no-{key}" in extra_keys:
        return False
    if f"+{key}" in extra_keys or key in extra_keys:
        return True
    tools = token_list(row.get("tools_to_run", ""))
    if tools:
        return token_matches(tools, method)
    return bool(default)


def disabled_reason(row, method):
    if method_enabled(row, method):
        return ""
    tools = str(row.get("tools_to_run", "") or "")
    extra = str(row.get("methods_extra", "") or "")
    if extra:
        return f"methods_extra disables {method}"
    if tools:
        return f"tools_to_run does not include {method}"
    return f"{method} is disabled"


def read_tsv(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment="#")


def write_tsv(rows, path, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def normalize_obs_name(obs_name, fallback_sample=None):
    obs_name = str(obs_name)
    if ":" in obs_name:
        sample_id, barcode = obs_name.split(":", 1)
    else:
        sample_id, barcode = fallback_sample, obs_name

    barcode = barcode.strip()
    if barcode.endswith("x"):
        barcode = f"{barcode[:-1]}-1"
    elif not re.search(r"-\d+$", barcode):
        barcode = f"{barcode}-1"

    if sample_id is None:
        raise ValueError(f"Cannot determine sample id for barcode: {obs_name}")
    return f"{sample_id}:{barcode}"


def load_pairs(pairs_path):
    pairs = read_tsv(pairs_path)
    if pairs.empty:
        return pairs
    for col in ["trajectory_id", "method", "enabled", "tools_to_run", "methods_extra"]:
        if col not in pairs:
            pairs[col] = ""
    pairs["method"] = pairs["method"].str.lower().replace("", "trajectory")
    pairs["enabled"] = pairs["enabled"].str.lower().replace("", "yes")
    pairs = pairs[(pairs["method"] == "velocity") & (pairs["enabled"] != "no")].copy()
    pairs = pairs.rename(columns={"trajectory_id": "pair_id"})
    return pairs


def load_units(reference_index, pairs_path):
    ref = read_tsv(reference_index)
    if ref.empty:
        return []
    for col in ["pair_id", "split_value", "metadata_csv", "umap_csv", "status"]:
        if col not in ref:
            ref[col] = ""
    ref = ref[ref["status"] == "ok"].copy()
    pairs = load_pairs(pairs_path)
    if not pairs.empty:
        ref = ref.merge(pairs, on="pair_id", how="left", suffixes=("", ".pair"))
        for col in ["tools_to_run", "methods_extra"]:
            pair_col = f"{col}.pair"
            if pair_col in ref:
                ref[col] = ref[col].where(ref[col].astype(str).str.len() > 0, ref[pair_col])
    return ref.fillna("").to_dict(orient="records")


def load_loom_table(loom_index):
    idx = read_tsv(loom_index)
    if idx.empty:
        return pd.DataFrame(columns=["sample_id", "loom_path", "status"])
    for col in ["sample_id", "loom_path", "status"]:
        if col not in idx:
            idx[col] = ""
    return idx[idx["status"].isin(["ok", "ok_existing"])].copy()


def load_looms(loom_table, loom_dir, allowed_samples):
    loom_paths = []
    if not loom_table.empty:
        for _, row in loom_table.iterrows():
            sample_id = str(row.get("sample_id", "") or "")
            if allowed_samples and sample_id not in allowed_samples:
                continue
            path = Path(str(row.get("loom_path", "") or ""))
            if path.exists():
                loom_paths.append(path)
    if not loom_paths:
        for path in sorted(Path(loom_dir).glob("*.loom")):
            sample_id = path.stem
            if allowed_samples and sample_id not in allowed_samples:
                continue
            loom_paths.append(path)
    if not loom_paths:
        raise FileNotFoundError("No .loom files found for this velocity unit.")

    adata_list = []
    for loom_path in sorted(set(loom_paths)):
        sample_name = loom_path.stem
        adata_i = scv.read(str(loom_path), cache=True)
        adata_i.obs_names = [normalize_obs_name(barcode, fallback_sample=sample_name) for barcode in adata_i.obs_names]
        adata_i.var_names_make_unique()
        adata_i.obs["sample_id"] = sample_name
        adata_list.append(adata_i)

    return ad.concat(adata_list, join="outer", merge="same")


def read_unit_metadata(unit):
    meta_df = pd.read_csv(unit["metadata_csv"], index_col=0)
    umap_df = pd.read_csv(unit["umap_csv"], index_col=0)
    meta_df.index = meta_df.index.astype(str)
    umap_df.index = umap_df.index.astype(str)
    return meta_df, umap_df


def integrate_metadata(adata, meta_df, umap_df):
    common_cells = adata.obs_names.intersection(meta_df.index).intersection(umap_df.index)
    if len(common_cells) == 0:
        raise ValueError("Loom cells do not overlap exported velocity metadata/UMAP.")

    adata = adata[common_cells].copy()
    adata.obsm["X_umap"] = umap_df.loc[common_cells, ["UMAP_1", "UMAP_2"]].to_numpy()

    for column in meta_df.columns:
        adata.obs[column] = meta_df.loc[common_cells, column].values

    if "cell_type" not in adata.obs and "cell_subtype" in adata.obs:
        adata.obs["cell_type"] = adata.obs["cell_subtype"].astype(str)

    if "cell_type" in adata.obs:
        observed = [value for value in pd.unique(adata.obs["cell_type"].astype(str)) if value and value.lower() != "nan"]
        adata.obs["cell_type"] = pd.Categorical(adata.obs["cell_type"].astype(str), categories=observed, ordered=True)
        colors = assign_discrete_colors(observed)
        adata.uns["cell_type_colors"] = [colors[x] for x in adata.obs["cell_type"].cat.categories]

    cluster_col = None
    for candidate in ("seurat_clusters", "cluster", "predicted_cluster"):
        if candidate in adata.obs:
            cluster_col = candidate
            break
    if cluster_col is not None:
        categories = sorted(pd.unique(adata.obs[cluster_col].astype(str)), key=lambda x: (not x.isdigit(), x))
        adata.obs["cluster_for_plot"] = pd.Categorical(adata.obs[cluster_col].astype(str), categories=categories, ordered=True)
        colors = assign_discrete_colors(categories, list(CLUSTER_COLORS.values()))
        adata.uns["cluster_for_plot_colors"] = [colors[x] for x in adata.obs["cluster_for_plot"].cat.categories]

    return adata


def assign_discrete_colors(categories, color_pool=None):
    clean = [str(x) for x in categories if pd.notna(x) and str(x)]
    if color_pool is None:
        color_pool = CELLTYPE_COLOR_POOL
    if len(clean) > len(color_pool):
        color_pool = color_pool + [plt.cm.tab20(i) for i in range(len(clean) - len(color_pool))]
    return {category: color_pool[idx] for idx, category in enumerate(clean)}


def annotate_cluster_numbers(ax, adata):
    if "cluster_for_plot" not in adata.obs or "X_umap" not in adata.obsm:
        return
    coords = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"], index=adata.obs_names)
    coords["cluster"] = adata.obs["cluster_for_plot"].astype(str).values
    centers = coords.dropna().groupby("cluster", as_index=False)[["UMAP_1", "UMAP_2"]].median()
    for _, row in centers.iterrows():
        ax.text(
            row["UMAP_1"],
            row["UMAP_2"],
            row["cluster"],
            fontsize=10,
            fontweight="bold",
            color="white",
            ha="center",
            va="center",
            bbox=dict(boxstyle="circle,pad=0.25", fc="black", ec="none", alpha=0.65),
        )


def placeholder_figure(path, title, message):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=13)
    ax.set_axis_off()
    ax.set_title(title)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_stream(adata, color_key, path, title, palette=None, annotate_clusters=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if color_key not in adata.obs:
        placeholder_figure(path, title, f"Missing {color_key}")
        return
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color=color_key,
        palette=palette,
        legend_loc="right margin",
        title=title,
        frameon=False,
        ax=ax,
        show=False,
    )
    if color_key == "cell_type" and "cell_type_colors" in adata.uns:
        categories = list(adata.obs["cell_type"].cat.categories)
        colors = list(adata.uns.get("cell_type_colors", []))
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
        handles = [
            Line2D([0], [0], marker="o", linestyle="", markerfacecolor=colors[idx], markeredgecolor="none", markersize=8, label=label)
            for idx, label in enumerate(categories)
        ]
        ax.legend(handles=handles, labels=[x.get_label() for x in handles], loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    if annotate_clusters:
        annotate_cluster_numbers(ax, adata)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_feature(adata, color_key, path, title, cmap):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if color_key not in adata.obs:
        placeholder_figure(path, title, f"Missing {color_key}")
        return
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    scv.pl.scatter(adata, basis="umap", color=color_key, color_map=cmap, size=30, title=title, frameon=False, ax=ax, show=False)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def compute_moments_compat(adata, n_pcs=30, n_neighbors=30):
    n_pcs = max(2, min(int(n_pcs), adata.n_obs - 1, adata.n_vars - 1))
    n_neighbors = max(2, min(int(n_neighbors), adata.n_obs - 1))
    try:
        scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
        return
    except TypeError as exc:
        if "write_knn_indices" not in str(exc):
            raise
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, method="gauss")
    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)


def numeric_summary(values):
    values = pd.to_numeric(pd.Series(values), errors="coerce")
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {"mean": math.nan, "median": math.nan, "q05": math.nan, "q95": math.nan, "max": math.nan}
    return {
        "mean": float(values.mean()),
        "median": float(values.median()),
        "q05": float(values.quantile(0.05)),
        "q95": float(values.quantile(0.95)),
        "max": float(values.max()),
    }


def matrix_total(matrix):
    if matrix is None:
        return math.nan
    total = matrix.sum()
    if hasattr(total, "item"):
        return float(total.item())
    return float(total)


def copy_stochastic_outputs(adata):
    copied = False
    if "velocity_umap" in adata.obsm:
        adata.obsm["velocity_umap_stochastic"] = adata.obsm["velocity_umap"].copy()
        copied = True
    if "velocity" in adata.layers:
        adata.layers["velocity_stochastic"] = adata.layers["velocity"].copy()
        copied = True
    for key in ("velocity_length", "velocity_confidence"):
        if key in adata.obs:
            adata.obs[f"{key}_stochastic"] = adata.obs[key].values
            copied = True
    for key in ("velocity_graph", "velocity_graph_neg"):
        if key in adata.uns:
            adata.uns[f"{key}_stochastic"] = adata.uns[key].copy()
            copied = True
    return copied


def compute_velocity_embedding_safe(adata):
    try:
        scv.tl.velocity_embedding(adata, basis="umap")
        return "ok"
    except Exception as exc:
        return str(exc)


def obs_values(adata, key, default=np.nan):
    if key in adata.obs:
        return pd.to_numeric(adata.obs[key], errors="coerce")
    return pd.Series([default] * adata.n_obs, index=adata.obs_names)


def export_velocity_vectors(pair_id, split_value, adata, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    umap = adata.obsm.get("X_umap")
    dyn = adata.obsm.get("velocity_umap")
    stoch = adata.obsm.get("velocity_umap_stochastic")
    if umap is None:
        umap = np.full((adata.n_obs, 2), np.nan)
    if dyn is None:
        dyn = np.full((adata.n_obs, 2), np.nan)
    if stoch is None:
        stoch = np.full((adata.n_obs, 2), np.nan)
    cluster = adata.obs["cluster_for_plot"].astype(str).values if "cluster_for_plot" in adata.obs else np.repeat("", adata.n_obs)
    cell_type = adata.obs["cell_type"].astype(str).values if "cell_type" in adata.obs else np.repeat("", adata.n_obs)
    df = pd.DataFrame(
        {
            "pair_id": pair_id,
            "split_value": display_split(split_value),
            "cell_id": adata.obs_names.astype(str),
            "UMAP_1": umap[:, 0],
            "UMAP_2": umap[:, 1],
            "velocity_umap_1": dyn[:, 0],
            "velocity_umap_2": dyn[:, 1],
            "velocity_umap_1_stochastic": stoch[:, 0],
            "velocity_umap_2_stochastic": stoch[:, 1],
            "latent_time": obs_values(adata, "latent_time").values,
            "velocity_pseudotime": obs_values(adata, "velocity_pseudotime").values,
            "velocity_length": obs_values(adata, "velocity_length").values,
            "velocity_confidence": obs_values(adata, "velocity_confidence").values,
            "velocity_length_stochastic": obs_values(adata, "velocity_length_stochastic").values,
            "velocity_confidence_stochastic": obs_values(adata, "velocity_confidence_stochastic").values,
            "cell_type": cell_type,
            "cluster": cluster,
        }
    )
    df.to_csv(path, sep="\t", index=False)
    return df


def build_qc_row(pair_id, split_value, adata, status="ok", reason=""):
    spliced_total = matrix_total(adata.layers.get("spliced"))
    unspliced_total = matrix_total(adata.layers.get("unspliced"))
    fit = numeric_summary(adata.var.get("fit_likelihood", pd.Series(dtype=float)))
    confidence = numeric_summary(adata.obs.get("velocity_confidence", pd.Series(dtype=float)))
    length = numeric_summary(adata.obs.get("velocity_length", pd.Series(dtype=float)))
    confidence_stoch = numeric_summary(adata.obs.get("velocity_confidence_stochastic", pd.Series(dtype=float)))
    length_stoch = numeric_summary(adata.obs.get("velocity_length_stochastic", pd.Series(dtype=float)))
    velocity_gene_n = int(pd.Series(adata.var.get("velocity_genes", False)).astype(bool).sum()) if adata.n_vars else 0
    return {
        "pair_id": pair_id,
        "split_value": display_split(split_value),
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "spliced_total": spliced_total,
        "unspliced_total": unspliced_total,
        "unspliced_spliced_ratio": unspliced_total / spliced_total if spliced_total and np.isfinite(spliced_total) else math.nan,
        "fit_likelihood_mean": fit["mean"],
        "fit_likelihood_median": fit["median"],
        "fit_likelihood_max": fit["max"],
        "velocity_gene_n": velocity_gene_n,
        "velocity_confidence_mean": confidence["mean"],
        "velocity_confidence_median": confidence["median"],
        "velocity_confidence_q05": confidence["q05"],
        "velocity_confidence_q95": confidence["q95"],
        "velocity_length_mean": length["mean"],
        "velocity_confidence_stochastic_mean": confidence_stoch["mean"],
        "velocity_length_stochastic_mean": length_stoch["mean"],
        "status": status,
        "reason": reason,
    }


def run_unit(unit, args, loom_table):
    pair_id = unit["pair_id"]
    split_value = display_split(unit.get("split_value", ""))
    unit = dict(unit)
    start = time.time()
    uid = unit_id(pair_id, split_value)
    out_dir = Path(args.velocity_output_dir) / "scvelo" / uid
    figure_dir = Path(args.figure_dir) / uid
    out_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    h5ad_path = out_dir / f"scvelo_result_{uid}.h5ad"
    qc_path = out_dir / f"velocity_qc_{uid}.tsv"
    vector_path = out_dir / f"velocity_vectors_{uid}.tsv"
    figures = {
        "stream_clusters": figure_dir / f"Figure_Velocity_StreamClusters_{uid}.png",
        "stream_celltypes": figure_dir / f"Figure_Velocity_StreamCellTypes_{uid}.png",
        "latent_time": figure_dir / f"Figure_Velocity_LatentTime_{uid}.png",
        "velocity_length": figure_dir / f"Figure_Velocity_Length_{uid}.png",
        "velocity_confidence": figure_dir / f"Figure_Velocity_Confidence_{uid}.png",
    }

    if not method_enabled(unit, "scvelo_dynamical", default=True):
        reason = disabled_reason(unit, "scvelo_dynamical")
        row = index_row(pair_id, split_value, unit, h5ad_path, qc_path, vector_path, figures, "no", "skipped_disabled", reason, time.time() - start, "not_run")
        return row, [], {}

    try:
        meta_df, umap_df = read_unit_metadata(unit)
        samples = {str(x).split(":", 1)[0] for x in meta_df.index if ":" in str(x)}
        adata = load_looms(loom_table, args.loom_dir, samples)
        adata = integrate_metadata(adata, meta_df, umap_df)
        if adata.n_obs < 3 or adata.n_vars < 3:
            raise ValueError("Fewer than 3 cells or genes after metadata/loom intersection.")

        scv.pp.filter_and_normalize(adata, min_shared_counts=args.min_shared_counts, n_top_genes=args.top_genes)
        compute_moments_compat(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)
        scv.tl.recover_dynamics(adata, n_jobs=args.threads)

        stochastic_status = "skipped_disabled"
        if method_enabled(unit, "scvelo_stochastic", default=True):
            scv.tl.velocity(adata, mode="stochastic")
            scv.tl.velocity_graph(adata)
            compute_velocity_embedding_safe(adata)
            scv.tl.velocity_confidence(adata)
            stochastic_status = "ok" if copy_stochastic_outputs(adata) else "not_available"

        scv.tl.velocity(adata, mode="dynamical")
        scv.tl.velocity_graph(adata)
        compute_velocity_embedding_safe(adata)
        scv.tl.latent_time(adata)
        scv.tl.velocity_confidence(adata)

        plot_stream(adata, "cluster_for_plot", figures["stream_clusters"], "RNA velocity by clusters", annotate_clusters=True)
        plot_stream(adata, "cell_type", figures["stream_celltypes"], "RNA velocity by cell types", palette=list(adata.uns.get("cell_type_colors", [])))
        plot_feature(adata, "latent_time", figures["latent_time"], "latent time", RAINBOW_CMAP)
        plot_feature(adata, "velocity_length", figures["velocity_length"], "dynamical velocity length", RAINBOW_CMAP)
        plot_feature(adata, "velocity_confidence", figures["velocity_confidence"], "dynamical velocity confidence", RAINBOW_CMAP)

        adata.write(str(h5ad_path))
        qc_row = build_qc_row(pair_id, split_value, adata)
        write_tsv([qc_row], qc_path, QC_COLUMNS)
        export_velocity_vectors(pair_id, split_value, adata, vector_path)
        row = index_row(pair_id, split_value, unit, h5ad_path, qc_path, vector_path, figures, "yes", "ok", "", time.time() - start, stochastic_status, adata.n_obs)
        dynamic = dynamic_outputs(args.project_root, pair_id, split_value, h5ad_path, qc_path, vector_path, figures)
        return row, [qc_row], dynamic
    except Exception as exc:
        reason = str(exc)
        qc_row = {key: "" for key in QC_COLUMNS}
        qc_row.update({"pair_id": pair_id, "split_value": split_value, "status": "failed", "reason": reason})
        write_tsv([qc_row], qc_path, QC_COLUMNS)
        row = index_row(pair_id, split_value, unit, h5ad_path, qc_path, vector_path, figures, "yes", "failed", reason, time.time() - start, "not_evaluated")
        return row, [qc_row], {}


def index_row(pair_id, split_value, unit, h5ad_path, qc_path, vector_path, figures, enabled, status, reason, runtime_s, stochastic_status, n_cells=0):
    return {
        "pair_id": pair_id,
        "split_value": display_split(split_value),
        "method": "scvelo_dynamical",
        "methods_enabled": enabled,
        "input_path": unit.get("metadata_csv", ""),
        "output_path": str(h5ad_path) if status == "ok" else "",
        "extra_path": unit.get("umap_csv", ""),
        "figure_path": ";".join(str(path) for path in figures.values() if Path(path).exists()),
        "n_cells": n_cells,
        "status": status,
        "reason": reason,
        "runtime_s": f"{runtime_s:.3f}",
        "qc_tsv": str(qc_path),
        "vector_tsv": str(vector_path),
        "stochastic_status": stochastic_status,
    }


def relative_path(path, base_dir):
    path = Path(path)
    base_dir = Path(base_dir)
    try:
        return str(path.resolve().relative_to(base_dir.resolve()))
    except Exception:
        return str(path)


def output_entry(path, typ, module, semantics, base_dir, schema=None):
    entry = {
        "path": relative_path(path, base_dir),
        "type": typ,
        "produced_by": module,
        "row_semantics": semantics,
    }
    if schema:
        entry["schema"] = schema
    return entry


def dynamic_outputs(project_root, pair_id, split_value, h5ad_path, qc_path, vector_path, figures):
    uid = unit_id(pair_id, split_value)
    outputs = {
        f"scvelo_result__{uid}": output_entry(h5ad_path, "h5ad", "10c_scvelo_dynamical", "scVelo dynamical result for one velocity unit", project_root),
        f"velocity_qc__{uid}": output_entry(qc_path, "tsv", "10c_scvelo_dynamical", "scVelo QC metrics for one velocity unit", project_root),
        f"velocity_vectors__{uid}": output_entry(vector_path, "tsv", "10c_scvelo_dynamical", "per-cell scVelo UMAP velocity vectors", project_root),
    }
    for name, path in figures.items():
        if Path(path).exists():
            outputs[f"velocity_figure__{uid}__{name}"] = output_entry(path, "png", "10c_scvelo_dynamical", "scVelo velocity figure", project_root)
    return outputs


def write_manifest(args, outputs):
    manifest = {
        "module": "10c_scvelo_dynamical",
        "version": args.module_version,
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "base_dir": str(args.project_root),
        "inputs": {
            "velocity_reference_index": str(args.reference_index),
            "velocity_loom_index": str(args.loom_index),
            "trajectory_pairs": str(args.pairs),
            "loom_dir": str(args.loom_dir),
        },
        "outputs": outputs,
        "depends_on": {
            "module_10a": str(args.module_10a_manifest),
            "module_10b": str(args.module_10b_manifest),
        },
    }
    manifest_path = Path(args.manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(description="Run split-aware scVelo dynamical velocity.")
    env = os.environ
    project_root = Path(env.get("PROJECT_ROOT", Path.cwd()))
    results_dir = Path(env.get("RESULTS_DIR", project_root / "results"))
    table_dir = Path(env.get("TABLE_DIR", results_dir / "tables"))
    manifest_dir = Path(env.get("MANIFEST_DIR", results_dir / "manifests"))
    velocity_dir = Path(env.get("VELOCITY_DIR", results_dir / "velocity"))
    parser.add_argument("--project-root", default=str(project_root))
    parser.add_argument("--reference-index", default=str(table_dir / "velocity/inputs/velocity_reference_index.tsv"))
    parser.add_argument("--loom-index", default=str(table_dir / "velocity/inputs/velocity_loom_index.tsv"))
    parser.add_argument("--pairs", default=str(project_root / "metadata/trajectory_pairs.tsv"))
    parser.add_argument("--loom-dir", default=str(env.get("VELOCITY_LOOM_DIR", velocity_dir / "loom")))
    parser.add_argument("--velocity-output-dir", default=str(env.get("VELOCITY_OUTPUT_DIR", velocity_dir / "output")))
    parser.add_argument("--figure-dir", default=str(Path(env.get("FIGURE_DIR", results_dir / "figures")) / "velocity/scvelo"))
    parser.add_argument("--out-index", default=str(table_dir / "velocity/methods/scvelo/scvelo_index.tsv"))
    parser.add_argument("--qc-out", default=str(table_dir / "velocity/methods/scvelo/velocity_qc.tsv"))
    parser.add_argument("--manifest", default=str(manifest_dir / "10c_scvelo_dynamical/_manifest.json"))
    parser.add_argument("--module-10a-manifest", default=str(manifest_dir / "10a_run_velocyto/_manifest.json"))
    parser.add_argument("--module-10b-manifest", default=str(manifest_dir / "10b_prepare_velocity_reference/_manifest.json"))
    parser.add_argument("--module-version", default=env.get("MODULE_10C_VERSION", env.get("MODULE_10_VERSION", "1.0")))
    parser.add_argument("--threads", type=int, default=int(env.get("SCVELO_THREADS", env.get("MAIN_THREADS", "8"))))
    parser.add_argument("--min-shared-counts", type=int, default=int(env.get("SCVELO_MIN_SHARED_COUNTS", "20")))
    parser.add_argument("--top-genes", type=int, default=int(env.get("SCVELO_TOP_GENES", "2000")))
    parser.add_argument("--n-pcs", type=int, default=int(env.get("SCVELO_N_PCS", "30")))
    parser.add_argument("--n-neighbors", type=int, default=int(env.get("SCVELO_N_NEIGHBORS", "30")))
    return parser.parse_args()


def main():
    args = parse_args()
    scv.settings.verbosity = 3
    scv.settings.presenter_view = True
    scv.set_figure_params("scvelo")

    units = load_units(args.reference_index, args.pairs)
    loom_table = load_loom_table(args.loom_index)
    index_rows = []
    qc_rows = []
    outputs = {
        "scvelo_index_tsv": output_entry(args.out_index, "tsv", "10c_scvelo_dynamical", "scVelo status by velocity unit", args.project_root),
        "velocity_qc_tsv": output_entry(args.qc_out, "tsv", "10c_scvelo_dynamical", "scVelo QC metrics by velocity unit", args.project_root),
    }

    for unit in units:
        row, qc, dynamic = run_unit(unit, args, loom_table)
        index_rows.append(row)
        qc_rows.extend(qc)
        outputs.update(dynamic)

    write_tsv(index_rows, args.out_index, INDEX_COLUMNS)
    write_tsv(qc_rows, args.qc_out, QC_COLUMNS)
    write_manifest(args, outputs)
    print(f"10c completed. scVelo index: {args.out_index}", flush=True)


if __name__ == "__main__":
    main()
