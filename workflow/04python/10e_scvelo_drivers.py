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
    "driver_tsv",
    "driver_n",
    "velocity_gene_n",
]

DRIVER_COLUMNS = [
    "pair_id",
    "split_value",
    "rank",
    "gene",
    "velocity_gene",
    "fit_likelihood",
    "fit_alpha",
    "fit_beta",
    "fit_gamma",
    "fit_t_",
]

OVERLAP_COLUMNS = [
    "pair_id",
    "split_a",
    "split_b",
    "driver_n_a",
    "driver_n_b",
    "overlap_n",
    "jaccard",
    "overlap_genes",
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


def read_tsv(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_tsv(rows, path, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


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


def read_scvelo_units(index_path):
    idx = read_tsv(index_path)
    if idx.empty:
        return []
    for col in ["pair_id", "split_value", "output_path", "n_cells", "status"]:
        if col not in idx:
            idx[col] = ""
    idx = idx[idx["status"] == "ok"].copy()
    return idx.fillna("").to_dict(orient="records")


def numeric_col(df, col):
    if col not in df:
        return pd.Series([math.nan] * len(df), index=df.index)
    return pd.to_numeric(df[col], errors="coerce")


def bool_col(df, col):
    if col not in df:
        return pd.Series([False] * len(df), index=df.index)
    values = df[col]
    if values.dtype == bool:
        return values.fillna(False)
    return values.astype(str).str.lower().isin(["true", "1", "yes"])


def driver_rows_for_unit(unit, top_n, driver_path):
    pair_id = unit["pair_id"]
    split_value = display_split(unit.get("split_value", ""))
    h5ad_path = Path(unit.get("output_path", ""))
    start = time.time()
    if not h5ad_path.exists():
        return index_row(unit, driver_path, "failed", f"missing scVelo h5ad: {h5ad_path}", time.time() - start), []

    adata = ad.read_h5ad(h5ad_path)
    var = adata.var.copy()
    var["gene"] = var.index.astype(str)
    var["fit_likelihood_numeric"] = numeric_col(var, "fit_likelihood")
    var["velocity_gene_bool"] = bool_col(var, "velocity_genes")
    var = var.sort_values(["velocity_gene_bool", "fit_likelihood_numeric", "gene"], ascending=[False, False, True])
    var = var.head(top_n)

    rows = []
    for rank, (_, row) in enumerate(var.iterrows(), start=1):
        out = {
            "pair_id": pair_id,
            "split_value": split_value,
            "rank": rank,
            "gene": row.get("gene", ""),
            "velocity_gene": str(bool(row.get("velocity_gene_bool", False))).lower(),
            "fit_likelihood": row.get("fit_likelihood_numeric", math.nan),
        }
        for col in ["fit_alpha", "fit_beta", "fit_gamma", "fit_t_"]:
            out[col] = row.get(col, "")
        rows.append(out)

    write_tsv(rows, driver_path, DRIVER_COLUMNS)
    velocity_gene_n = int(bool_col(adata.var, "velocity_genes").sum())
    return index_row(unit, driver_path, "ok", "", time.time() - start, len(rows), velocity_gene_n), rows


def index_row(unit, driver_path, status, reason, runtime_s, driver_n=0, velocity_gene_n=0):
    return {
        "pair_id": unit.get("pair_id", ""),
        "split_value": display_split(unit.get("split_value", "")),
        "method": "scvelo_drivers",
        "methods_enabled": "yes",
        "input_path": unit.get("output_path", ""),
        "output_path": str(driver_path) if status == "ok" else "",
        "extra_path": "",
        "figure_path": "",
        "n_cells": unit.get("n_cells", ""),
        "status": status,
        "reason": reason,
        "runtime_s": f"{runtime_s:.3f}",
        "driver_tsv": str(driver_path),
        "driver_n": driver_n,
        "velocity_gene_n": velocity_gene_n,
    }


def build_overlap(driver_rows):
    by_pair = {}
    for row in driver_rows:
        by_pair.setdefault(row["pair_id"], {}).setdefault(row["split_value"], set()).add(row["gene"])

    overlaps = []
    for pair_id, split_map in by_pair.items():
        splits = sorted(split_map)
        for i, split_a in enumerate(splits):
            for split_b in splits[i + 1 :]:
                genes_a = split_map[split_a]
                genes_b = split_map[split_b]
                overlap = sorted(genes_a.intersection(genes_b))
                union_n = len(genes_a.union(genes_b))
                overlaps.append(
                    {
                        "pair_id": pair_id,
                        "split_a": split_a,
                        "split_b": split_b,
                        "driver_n_a": len(genes_a),
                        "driver_n_b": len(genes_b),
                        "overlap_n": len(overlap),
                        "jaccard": (len(overlap) / union_n) if union_n else math.nan,
                        "overlap_genes": ",".join(overlap),
                    }
                )
    return overlaps


def write_manifest(args, outputs):
    manifest = {
        "module": "10e_scvelo_drivers",
        "version": args.module_version,
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "base_dir": str(args.project_root),
        "inputs": {"module_10c_index": str(args.scvelo_index)},
        "outputs": outputs,
        "depends_on": {"module_10c": str(args.module_10c_manifest)},
    }
    manifest_path = Path(args.manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(description="Extract scVelo driver genes and split overlap.")
    env = os.environ
    project_root = Path(env.get("PROJECT_ROOT", Path.cwd()))
    results_dir = Path(env.get("RESULTS_DIR", project_root / "results"))
    table_dir = Path(env.get("TABLE_DIR", results_dir / "tables"))
    manifest_dir = Path(env.get("MANIFEST_DIR", results_dir / "manifests"))
    parser.add_argument("--project-root", default=str(project_root))
    parser.add_argument("--scvelo-index", default=str(table_dir / "velocity/methods/scvelo/scvelo_index.tsv"))
    parser.add_argument("--out-index", default=str(table_dir / "velocity/methods/drivers/velocity_driver_index.tsv"))
    parser.add_argument("--overlap-out", default=str(table_dir / "velocity/methods/drivers/velocity_driver_overlap.tsv"))
    parser.add_argument("--driver-dir", default=str(table_dir / "velocity/methods/drivers"))
    parser.add_argument("--manifest", default=str(manifest_dir / "10e_scvelo_drivers/_manifest.json"))
    parser.add_argument("--module-10c-manifest", default=str(manifest_dir / "10c_scvelo_dynamical/_manifest.json"))
    parser.add_argument("--module-version", default=env.get("MODULE_10E_VERSION", env.get("MODULE_10_VERSION", "1.0")))
    parser.add_argument("--top-n", type=int, default=int(env.get("VELOCITY_DRIVER_TOP_N", "200")))
    return parser.parse_args()


def main():
    args = parse_args()
    units = read_scvelo_units(args.scvelo_index)
    index_rows = []
    driver_rows = []
    outputs = {
        "velocity_driver_index_tsv": output_entry(args.out_index, "tsv", "10e_scvelo_drivers", "driver extraction status by velocity unit", args.project_root),
        "velocity_driver_overlap_tsv": output_entry(args.overlap_out, "tsv", "10e_scvelo_drivers", "top driver overlap by split", args.project_root),
    }

    for unit in units:
        uid = unit_id(unit["pair_id"], unit.get("split_value", ""))
        driver_path = Path(args.driver_dir) / f"velocity_drivers_{uid}_top{args.top_n}.tsv"
        row, rows = driver_rows_for_unit(unit, args.top_n, driver_path)
        index_rows.append(row)
        driver_rows.extend(rows)
        if row["status"] == "ok":
            outputs[f"velocity_drivers__{uid}"] = output_entry(driver_path, "tsv", "10e_scvelo_drivers", "top scVelo driver genes for one velocity unit", args.project_root)

    overlap_rows = build_overlap(driver_rows)
    write_tsv(index_rows, args.out_index, INDEX_COLUMNS)
    write_tsv(overlap_rows, args.overlap_out, OVERLAP_COLUMNS)
    write_manifest(args, outputs)
    print(f"10e completed. driver index: {args.out_index}", flush=True)


if __name__ == "__main__":
    main()
