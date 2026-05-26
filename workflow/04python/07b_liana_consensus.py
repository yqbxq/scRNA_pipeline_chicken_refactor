#!/usr/bin/env python3
"""LIANA+ multi-method consensus communication analysis."""

from __future__ import annotations

import argparse
import os
import sys
import csv
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from helpers.chicken_human_ortholog import OrthologLookup
from helpers.liana_io_utils import (
    LIANA_OUTPUT_COLUMNS,
    empty_liana_df,
    load_scrna_h5ad,
    read_cluster_eligibility,
    validate_h5ad_for_liana,
    write_liana_manifest,
)
from helpers.lr_axis_id_utils import standardize_lr_axis_id_df
from helpers.lr_axis_id_utils import standardize_lr_axis_id


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad-path", default=os.environ.get("LIANA_H5AD_PATH", "results/90a_export_h5ad/03d_panorama"))
    parser.add_argument("--cluster-eligibility", default=os.environ.get("LIANA_CLUSTER_ELIGIBILITY", ""))
    parser.add_argument("--cluster-col", default=os.environ.get("LIANA_CLUSTER_COL", os.environ.get("COMMUNICATION_CELL_TYPE_COL", "cell_subtype")))
    parser.add_argument("--condition-col", default=os.environ.get("LIANA_CONDITION_COL", "condition"))
    parser.add_argument("--methods", default=os.environ.get("LIANA_METHODS_LIST", "cellphonedb,connectome,sca,natmi,logfc,rank_aggregate"))
    parser.add_argument("--aggregate-by", default=os.environ.get("LIANA_CONSENSUS_AGGREGATE", "rank_aggregate"))
    parser.add_argument("--min-methods-agreed", type=int, default=int(os.environ.get("LIANA_MIN_METHODS_AGREED_INSIDE", "3")))
    parser.add_argument("--min-cells-per-group", type=int, default=int(os.environ.get("LIANA_MIN_CELLS_PER_GROUP", "100")))
    parser.add_argument("--resource", default=os.environ.get("LIANA_RESOURCE_DB", "consensus"))
    parser.add_argument("--output", default=os.environ.get("LIANA_CONSENSUS_OUTPUT", "results/tables/communication/liana_consensus/liana_consensus_lr.tsv"))
    parser.add_argument("--manifest", default=os.environ.get("MODULE_07B_MANIFEST", "results/manifests/07b_liana_consensus/_manifest.json"))
    parser.add_argument("--ortholog-lut", default=os.environ.get("ORTHOLOG_CHICKEN_HUMAN_TSV", "metadata/ortholog_chicken_human.tsv"))
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--mock-run", action="store_true", help="Write deterministic fixture output without importing LIANA/anndata.")
    return parser.parse_args()


def output_schema_cols() -> list[str]:
    return list(LIANA_OUTPUT_COLUMNS)


def row_count(rows) -> int:
    if hasattr(rows, "__len__"):
        return len(rows)
    return 0


def write_output(rows, path: str | Path) -> None:
    item = Path(path)
    item.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(rows, "columns") and hasattr(rows, "to_csv"):
        out = rows.copy()
        for col in output_schema_cols():
            if col not in out.columns:
                out[col] = ""
        out = out[output_schema_cols() + [col for col in out.columns if col not in output_schema_cols()]]
        out.to_csv(item, sep="\t", index=False)
        return

    rows = list(rows)
    fieldnames = list(output_schema_cols())
    extras = []
    for row in rows:
        for key in row:
            if key not in fieldnames and key not in extras:
                extras.append(key)
    with item.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames + extras)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames + extras})


def mock_liana_output(args: argparse.Namespace) -> list[dict[str, object]]:
    rows = [
        {
            "condition_value": "mock_condition",
            "source": "pGC",
            "target": "rgGC",
            "ligand": "WNT5A_WNT5B",
            "receptor": "FZD3+LRP6",
            "liana_consensus_score": 0.01,
            "n_methods_agreed": max(args.min_methods_agreed, 3),
            "liana_consensus_hit": 1,
            "status": "ok",
            "reason": "mock_run",
        },
        {
            "condition_value": "mock_condition",
            "source": "TC",
            "target": "GC",
            "ligand": "TGFB1",
            "receptor": "TGFBR1/TGFBR2",
            "liana_consensus_score": 0.2,
            "n_methods_agreed": 1,
            "liana_consensus_hit": 0,
            "status": "below_consensus",
            "reason": "mock_run",
        },
    ]
    ortholog = OrthologLookup(args.ortholog_lut)
    for row in rows:
        row["ligand_human"] = ortholog.chicken_to_human_complex(row["ligand"])
        row["receptor_human"] = ortholog.chicken_to_human_complex(row["receptor"])
        row["lr_axis_id"] = standardize_lr_axis_id(
            row["ligand_human"],
            row["receptor_human"],
            row["source"],
            row["target"],
        )
    return rows


def liana_result_from_adata(args: argparse.Namespace):
    import pandas as pd

    adata = load_scrna_h5ad(args.h5ad_path)
    validate_h5ad_for_liana(adata, args.cluster_col, args.condition_col)

    if args.check_only:
        print(f"[07b] check-only: H5AD ok (n_obs={adata.n_obs}, n_vars={adata.n_vars})")
        return [], "check_only"

    allowed_tiers = set(os.environ.get("INVENTORY_GATE_ALLOWED_TIERS_COMMUNICATION", "primary,exploratory,primary_merged").split(","))
    allowed_tiers = {item.strip() for item in allowed_tiers if item.strip()}
    eligible = read_cluster_eligibility(args.cluster_eligibility, allowed_tiers) if args.cluster_eligibility else set()
    if eligible:
        mask = adata.obs[args.cluster_col].astype(str).isin(eligible)
        adata = adata[mask].copy()
    if adata.n_obs == 0:
        return [], "empty_no_eligible_clusters"

    try:
        import liana as li
    except ImportError as exc:
        raise RuntimeError("Python package liana is required unless --mock-run is used.") from exc

    result_rows: list[pd.DataFrame] = []
    methods = [item.strip() for item in args.methods.split(",") if item.strip() and item.strip() != "rank_aggregate"]
    ortholog = OrthologLookup(args.ortholog_lut)

    for condition_value in sorted(adata.obs[args.condition_col].astype(str).unique()):
        subset = adata[adata.obs[args.condition_col].astype(str) == condition_value].copy()
        if subset.n_obs < args.min_cells_per_group:
            continue
        li.mt.rank_aggregate(
            subset,
            groupby=args.cluster_col,
            resource_name=args.resource,
            methods=methods,
            verbose=False,
            inplace=True,
        )
        liana_df = subset.uns.get("liana_res")
        if liana_df is None or len(liana_df) == 0:
            continue
        df = pd.DataFrame(liana_df).copy()
        rename = {"ligand_complex": "ligand", "receptor_complex": "receptor"}
        df = df.rename(columns={key: value for key, value in rename.items() if key in df.columns})
        for col in ("source", "target", "ligand", "receptor"):
            if col not in df.columns:
                raise ValueError(f"LIANA output missing required column: {col}")
        df["condition_value"] = condition_value
        df["ligand_human"] = df["ligand"].map(ortholog.chicken_to_human_complex)
        df["receptor_human"] = df["receptor"].map(ortholog.chicken_to_human_complex)
        df["n_methods_agreed"] = count_methods_agreed(df)
        score_col = "magnitude_rank" if "magnitude_rank" in df.columns else args.aggregate_by
        df["liana_consensus_score"] = pd.to_numeric(df.get(score_col, pd.Series([float("nan")] * len(df))), errors="coerce")
        df["liana_consensus_hit"] = (df["n_methods_agreed"] >= args.min_methods_agreed).astype(int)
        df["status"] = "ok"
        df["reason"] = ""
        result_rows.append(standardize_lr_axis_id_df(df))

    if not result_rows:
        return [], "empty_insufficient_cells"
    out = pd.concat(result_rows, ignore_index=True)
    return out.drop_duplicates(subset=["lr_axis_id", "condition_value"]), "ok"


def count_methods_agreed(df: pd.DataFrame) -> pd.Series:
    import pandas as pd

    agreed = pd.Series(0, index=df.index)
    p_threshold = float(os.environ.get("LIANA_CELLPHONEDB_PVAL_THRESHOLD", "0.05"))
    for p_col in [col for col in df.columns if col.endswith("_pvals")]:
        agreed += (pd.to_numeric(df[p_col], errors="coerce") < p_threshold).fillna(False).astype(int)
    for score_col in [col for col in df.columns if col.endswith("_score")]:
        agreed += (pd.to_numeric(df[score_col], errors="coerce") > 0).fillna(False).astype(int)
    return agreed


def main() -> int:
    args = parse_args()
    base_dir = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()

    if args.mock_run:
        out = mock_liana_output(args)
        status = "ok_mock"
    else:
        try:
            out, status = liana_result_from_adata(args)
        except Exception as exc:
            if os.environ.get("COMMUNICATION_REQUIRE_FULL_PIPELINE", "no").lower() in {"yes", "true", "1", "on"}:
                raise
            reason = str(exc)
            print(f"[07b] {reason}; writing empty LIANA output", file=sys.stderr)
            out = empty_liana_df("empty_runtime_unavailable", reason)
            status = "empty_runtime_unavailable"

    write_output(out, args.output)
    write_liana_manifest(
        args.manifest,
        args.output,
        base_dir=base_dir,
        status=status,
            n_rows=row_count(out),
        inputs={
            "h5ad_path": args.h5ad_path,
            "cluster_eligibility": args.cluster_eligibility,
            "ortholog_lut": args.ortholog_lut,
            "methods": args.methods,
        },
    )
    print(f"[07b] wrote {args.output} rows={row_count(out)} status={status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
