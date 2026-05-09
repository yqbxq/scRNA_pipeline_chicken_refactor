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
    "fate_csv",
    "macrostates_tsv",
    "terminal_figure_path",
]

MACRO_COLUMNS = [
    "pair_id",
    "split_value",
    "cell_id",
    "macrostate",
    "terminal_state",
    "initial_state",
    "max_fate_state",
    "max_fate_probability",
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
    extra = str(row.get("methods_extra", "") or "")
    tools = str(row.get("tools_to_run", "") or "")
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
    return pairs.rename(columns={"trajectory_id": "pair_id"})


def load_units(scvelo_index, pairs_path, qc_path):
    idx = read_tsv(scvelo_index)
    if idx.empty:
        return []
    for col in ["pair_id", "split_value", "output_path", "status", "n_cells"]:
        if col not in idx:
            idx[col] = ""
    idx = idx[idx["status"] == "ok"].copy()
    pairs = load_pairs(pairs_path)
    if not pairs.empty:
        idx = idx.merge(pairs, on="pair_id", how="left", suffixes=("", ".pair"))
        for col in ["tools_to_run", "methods_extra"]:
            pair_col = f"{col}.pair"
            if pair_col in idx:
                idx[col] = idx[col].where(idx[col].astype(str).str.len() > 0, idx[pair_col])
    qc = read_tsv(qc_path)
    if not qc.empty:
        for col in ["pair_id", "split_value", "velocity_confidence_mean"]:
            if col not in qc:
                qc[col] = ""
        idx = idx.merge(qc[["pair_id", "split_value", "velocity_confidence_mean"]], on=["pair_id", "split_value"], how="left")
    return idx.fillna("").to_dict(orient="records")


def relative_path(path, base_dir):
    path = Path(path)
    base_dir = Path(base_dir)
    try:
        return str(path.resolve().relative_to(base_dir.resolve()))
    except Exception:
        return str(path)


def output_entry(path, typ, module, semantics, base_dir):
    return {
        "path": relative_path(path, base_dir),
        "type": typ,
        "produced_by": module,
        "row_semantics": semantics,
    }


def lineage_to_frame(obj, index):
    if obj is None:
        return pd.DataFrame(index=index)
    if hasattr(obj, "to_df"):
        df = obj.to_df()
        df.index = df.index.astype(str)
        return df
    arr = np.asarray(obj)
    if arr.ndim == 1:
        arr = arr.reshape((-1, 1))
    names = getattr(obj, "names", None)
    if names is None or len(names) != arr.shape[1]:
        names = [f"fate_{i + 1}" for i in range(arr.shape[1])]
    return pd.DataFrame(arr, index=index, columns=[str(x) for x in names])


def obs_col(adata, col):
    if col not in adata.obs:
        return pd.Series([""] * adata.n_obs, index=adata.obs_names)
    return adata.obs[col].astype(str)


def run_cellrank(adata, n_states):
    import cellrank as cr

    if "velocity_graph" not in adata.uns:
        raise ValueError("scVelo result lacks velocity_graph")
    kernel = cr.kernels.VelocityKernel(adata)
    kernel.compute_transition_matrix()
    estimator = cr.estimators.GPCCA(kernel)
    n_states = max(2, min(int(n_states), adata.n_obs - 1))
    cluster_key = "cell_type" if "cell_type" in adata.obs else ("cluster_for_plot" if "cluster_for_plot" in adata.obs else None)
    try:
        if cluster_key:
            estimator.compute_macrostates(n_states=n_states, cluster_key=cluster_key)
        else:
            estimator.compute_macrostates(n_states=n_states)
    except TypeError:
        estimator.compute_macrostates(n_states=n_states)
    try:
        estimator.compute_terminal_states()
    except Exception:
        pass
    try:
        estimator.compute_initial_states()
    except Exception:
        pass
    estimator.compute_fate_probabilities()
    return estimator


def build_macro_df(pair_id, split_value, adata, fate_df):
    max_state = pd.Series([""] * adata.n_obs, index=adata.obs_names)
    max_prob = pd.Series([math.nan] * adata.n_obs, index=adata.obs_names)
    if not fate_df.empty:
        fate_numeric = fate_df.apply(pd.to_numeric, errors="coerce")
        max_state = fate_numeric.idxmax(axis=1).astype(str)
        max_prob = fate_numeric.max(axis=1)
    out = pd.DataFrame(
        {
            "pair_id": pair_id,
            "split_value": display_split(split_value),
            "cell_id": adata.obs_names.astype(str),
            "macrostate": obs_col(adata, "macrostates").values,
            "terminal_state": obs_col(adata, "terminal_states").values,
            "initial_state": obs_col(adata, "initial_states").values,
            "max_fate_state": max_state.reindex(adata.obs_names).values,
            "max_fate_probability": max_prob.reindex(adata.obs_names).values,
            "cell_type": obs_col(adata, "cell_type").values,
            "cluster": obs_col(adata, "cluster_for_plot").values,
        }
    )
    return out


def plot_fate(pair_id, split_value, adata, macro_df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    coords = adata.obsm.get("X_umap")
    fig, ax = plt.subplots(figsize=(7, 5.5))
    if coords is None or macro_df.empty:
        ax.text(0.5, 0.5, "No CellRank fate data", ha="center", va="center")
        ax.set_axis_off()
    else:
        values = pd.to_numeric(macro_df["max_fate_probability"], errors="coerce")
        sc = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=8, cmap="viridis", linewidths=0)
        fig.colorbar(sc, ax=ax, label="max fate probability")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
    ax.set_title(f"CellRank fate {pair_id} {display_split(split_value)}")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_terminal(pair_id, split_value, adata, macro_df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    coords = adata.obsm.get("X_umap")
    fig, ax = plt.subplots(figsize=(7, 5.5))
    if coords is None or macro_df.empty:
        ax.text(0.5, 0.5, "No CellRank terminal data", ha="center", va="center")
        ax.set_axis_off()
    else:
        labels = macro_df["max_fate_state"].astype(str)
        cats = sorted(x for x in labels.unique() if x and x.lower() != "nan")
        cmap = {cat: plt.cm.tab20(i % 20) for i, cat in enumerate(cats)}
        colors = [cmap.get(x, (0.7, 0.7, 0.7, 1)) for x in labels]
        ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=8, linewidths=0)
        handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=cmap[cat], label=cat) for cat in cats]
        if handles:
            ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
    ax.set_title(f"CellRank terminals {pair_id} {display_split(split_value)}")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def empty_outputs(fate_path, macro_path):
    fate_path = Path(fate_path)
    fate_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(columns=["cell_id"]).to_csv(fate_path, index=False)
    write_tsv([], macro_path, MACRO_COLUMNS)


def index_row(unit, fate_path, macro_path, fate_fig, terminal_fig, status, reason, runtime_s, n_cells=0, enabled="yes"):
    return {
        "pair_id": unit.get("pair_id", ""),
        "split_value": display_split(unit.get("split_value", "")),
        "method": "cellrank",
        "methods_enabled": enabled,
        "input_path": unit.get("output_path", ""),
        "output_path": str(fate_path) if status == "ok" else "",
        "extra_path": str(macro_path),
        "figure_path": str(fate_fig) if Path(fate_fig).exists() else "",
        "n_cells": n_cells,
        "status": status,
        "reason": reason,
        "runtime_s": f"{runtime_s:.3f}",
        "fate_csv": str(fate_path),
        "macrostates_tsv": str(macro_path),
        "terminal_figure_path": str(terminal_fig) if Path(terminal_fig).exists() else "",
    }


def run_unit(unit, args):
    start = time.time()
    pair_id = unit["pair_id"]
    split_value = display_split(unit.get("split_value", ""))
    uid = unit_id(pair_id, split_value)
    fate_path = Path(args.cellrank_dir) / f"cellrank_fate_{uid}.csv"
    macro_path = Path(args.cellrank_dir) / f"cellrank_macrostates_{uid}.tsv"
    fate_fig = Path(args.figure_dir) / uid / f"Figure_CR_fate_{uid}.png"
    terminal_fig = Path(args.figure_dir) / uid / f"Figure_CR_terminal_{uid}.png"

    if not method_enabled(unit, "cellrank", default=True):
        empty_outputs(fate_path, macro_path)
        reason = disabled_reason(unit, "cellrank")
        return index_row(unit, fate_path, macro_path, fate_fig, terminal_fig, "skipped_disabled", reason, time.time() - start, enabled="no"), {}

    n_cells = int(float(unit.get("n_cells") or 0))
    confidence = pd.to_numeric(pd.Series([unit.get("velocity_confidence_mean", "")]), errors="coerce").iloc[0]
    if n_cells < args.min_cells or (np.isfinite(confidence) and confidence < args.min_velocity_confidence):
        empty_outputs(fate_path, macro_path)
        reason = f"n_cells={n_cells}, velocity_confidence_mean={confidence}"
        return index_row(unit, fate_path, macro_path, fate_fig, terminal_fig, "skipped_low_quality", reason, time.time() - start, n_cells=n_cells), {}

    try:
        adata = ad.read_h5ad(unit["output_path"])
        estimator = run_cellrank(adata, args.n_states)
        fate_df = lineage_to_frame(getattr(estimator, "fate_probabilities", None), adata.obs_names.astype(str))
        fate_df.insert(0, "cell_id", fate_df.index.astype(str))
        fate_df.insert(0, "split_value", split_value)
        fate_df.insert(0, "pair_id", pair_id)
        fate_path.parent.mkdir(parents=True, exist_ok=True)
        fate_df.to_csv(fate_path, index=False)

        macro_df = build_macro_df(pair_id, split_value, adata, fate_df.set_index("cell_id").drop(columns=["pair_id", "split_value"], errors="ignore"))
        write_tsv(macro_df.to_dict(orient="records"), macro_path, MACRO_COLUMNS)
        plot_fate(pair_id, split_value, adata, macro_df, fate_fig)
        plot_terminal(pair_id, split_value, adata, macro_df, terminal_fig)
        dynamic = {
            f"cellrank_fate__{uid}": output_entry(fate_path, "csv", "10f_cellrank_fate", "CellRank fate probabilities by cell", args.project_root),
            f"cellrank_macrostates__{uid}": output_entry(macro_path, "tsv", "10f_cellrank_fate", "CellRank macrostates by cell", args.project_root),
            f"cellrank_fate_figure__{uid}": output_entry(fate_fig, "png", "10f_cellrank_fate", "CellRank fate probability UMAP", args.project_root),
            f"cellrank_terminal_figure__{uid}": output_entry(terminal_fig, "png", "10f_cellrank_fate", "CellRank terminal state UMAP", args.project_root),
        }
        return index_row(unit, fate_path, macro_path, fate_fig, terminal_fig, "ok", "", time.time() - start, n_cells=adata.n_obs), dynamic
    except ImportError as exc:
        empty_outputs(fate_path, macro_path)
        return index_row(unit, fate_path, macro_path, fate_fig, terminal_fig, "failed_no_package", str(exc), time.time() - start, n_cells=n_cells), {}
    except Exception as exc:
        empty_outputs(fate_path, macro_path)
        return index_row(unit, fate_path, macro_path, fate_fig, terminal_fig, "failed", str(exc), time.time() - start, n_cells=n_cells), {}


def write_manifest(args, outputs):
    manifest = {
        "module": "10f_cellrank_fate",
        "version": args.module_version,
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "base_dir": str(args.project_root),
        "inputs": {
            "module_10c_index": str(args.scvelo_index),
            "velocity_qc": str(args.qc),
            "trajectory_pairs": str(args.pairs),
        },
        "outputs": outputs,
        "depends_on": {"module_10c": str(args.module_10c_manifest)},
    }
    manifest_path = Path(args.manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(description="Run CellRank fate probabilities for velocity units.")
    env = os.environ
    project_root = Path(env.get("PROJECT_ROOT", Path.cwd()))
    results_dir = Path(env.get("RESULTS_DIR", project_root / "results"))
    table_dir = Path(env.get("TABLE_DIR", results_dir / "tables"))
    manifest_dir = Path(env.get("MANIFEST_DIR", results_dir / "manifests"))
    figure_dir = Path(env.get("FIGURE_DIR", results_dir / "figures")) / "velocity/cellrank"
    parser.add_argument("--project-root", default=str(project_root))
    parser.add_argument("--scvelo-index", default=str(table_dir / "velocity/methods/scvelo/scvelo_index.tsv"))
    parser.add_argument("--qc", default=str(table_dir / "velocity/methods/scvelo/velocity_qc.tsv"))
    parser.add_argument("--pairs", default=str(project_root / "metadata/trajectory_pairs.tsv"))
    parser.add_argument("--out-index", default=str(table_dir / "velocity/methods/cellrank/cellrank_index.tsv"))
    parser.add_argument("--cellrank-dir", default=str(table_dir / "velocity/methods/cellrank"))
    parser.add_argument("--figure-dir", default=str(figure_dir))
    parser.add_argument("--manifest", default=str(manifest_dir / "10f_cellrank_fate/_manifest.json"))
    parser.add_argument("--module-10c-manifest", default=str(manifest_dir / "10c_scvelo_dynamical/_manifest.json"))
    parser.add_argument("--module-version", default=env.get("MODULE_10F_VERSION", env.get("MODULE_10_VERSION", "1.0")))
    parser.add_argument("--min-cells", type=int, default=int(env.get("CELLRANK_MIN_CELLS", "200")))
    parser.add_argument("--min-velocity-confidence", type=float, default=float(env.get("CELLRANK_MIN_VELOCITY_CONFIDENCE", "0.05")))
    parser.add_argument("--n-states", type=int, default=int(env.get("CELLRANK_N_STATES", "6")))
    return parser.parse_args()


def main():
    args = parse_args()
    units = load_units(args.scvelo_index, args.pairs, args.qc)
    rows = []
    outputs = {
        "cellrank_index_tsv": output_entry(args.out_index, "tsv", "10f_cellrank_fate", "CellRank status by velocity unit", args.project_root)
    }
    for unit in units:
        row, dynamic = run_unit(unit, args)
        rows.append(row)
        outputs.update(dynamic)
    write_tsv(rows, args.out_index, INDEX_COLUMNS)
    write_manifest(args, outputs)
    print(f"10f completed. CellRank index: {args.out_index}", flush=True)


if __name__ == "__main__":
    main()
