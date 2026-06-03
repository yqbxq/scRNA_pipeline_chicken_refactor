#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import importlib.util
from pathlib import Path

from helpers.spatial_io import resolve_spatial_h5ad_path_strict


SVG_COLUMNS = [
    "section_id", "condition", "method", "gene_id", "gene_symbol", "score", "pvalue", "qvalue", "rank",
    "mean_expression", "n_spots", "n_genes_tested", "status", "reason",
    "is_ligand_candidate", "is_receptor_candidate", "is_receiver_target_candidate",
    "is_region_marker_candidate", "is_regulator_target_candidate", "linked_question_ids",
]


def fingerprint(path: Path) -> str:
    if not path.exists():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(args, status: str, reason: str, h5ad_path: Path | str = "", svg_tsv: str = "", method: str = "SpatialDE2") -> None:
    import pandas as pd

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{
        "section_id": args.section_id or "all",
        "status": status,
        "reason": reason,
        "method": method,
        "proxy_allowed": "yes" if args.allow_proxy else "no",
        "input_h5ad_path": str(h5ad_path),
        "fingerprint": fingerprint(Path(h5ad_path)) if h5ad_path else "",
        "svg_tsv": svg_tsv,
    }]).to_csv(out_dir / "spatialde2_manifest.tsv", sep="\t", index=False)


def empty_svg_df(section: str, condition: str, status: str, reason: str, method: str):
    import pandas as pd

    return pd.DataFrame([{
        "section_id": section,
        "condition": condition,
        "method": method,
        "gene_id": "",
        "gene_symbol": "",
        "score": "",
        "pvalue": "",
        "qvalue": "",
        "rank": "",
        "mean_expression": "",
        "n_spots": 0,
        "n_genes_tested": 0,
        "status": status,
        "reason": reason,
        "is_ligand_candidate": "no",
        "is_receptor_candidate": "no",
        "is_receiver_target_candidate": "no",
        "is_region_marker_candidate": "no",
        "is_regulator_target_candidate": "no",
        "linked_question_ids": "",
    }], columns=SVG_COLUMNS)


def spatialde2_available() -> bool:
    return importlib.util.find_spec("spatialde2") is not None or importlib.util.find_spec("SpatialDE") is not None


def formal_spatialde2(adata, section: str, condition: str):
    import numpy as np
    import pandas as pd

    coords = np.asarray(adata.obsm["spatial"])
    if coords.shape[1] < 2:
        raise ValueError('obsm["spatial"] must have at least two columns')
    x = adata.X
    expr = pd.DataFrame(np.asarray(x.toarray() if hasattr(x, "toarray") else x), index=adata.obs_names, columns=adata.var_names)
    sample_info = pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1]}, index=adata.obs_names)

    result = None
    if importlib.util.find_spec("SpatialDE") is not None:
        import SpatialDE
        if hasattr(SpatialDE, "run"):
            result = SpatialDE.run(sample_info, expr)
    elif importlib.util.find_spec("spatialde2") is not None:
        import spatialde2
        for fn_name in ("run", "run_spatialde2", "spatialde2", "fit"):
            fn = getattr(spatialde2, fn_name, None)
            if callable(fn):
                try:
                    result = fn(sample_info=sample_info, counts=expr)
                except TypeError:
                    result = fn(sample_info, expr)
                break

    if result is None:
        raise RuntimeError("SpatialDE2/SpatialDE package detected but no supported run API was found")
    df = pd.DataFrame(result)
    gene_col = next((c for c in ("g", "gene", "gene_id", "gene_symbol") if c in df.columns), None)
    score_col = next((c for c in ("FSV", "score", "LLR", "stat", "spatial_score") if c in df.columns), None)
    p_col = next((c for c in ("pval", "pvalue", "p_value", "p") if c in df.columns), None)
    q_col = next((c for c in ("qval", "qvalue", "q_value", "padj") if c in df.columns), None)
    if gene_col is None:
        raise RuntimeError("SpatialDE2 result is missing a gene column")
    if score_col is None:
        score_col = q_col or p_col
    gene_values = df[gene_col].astype(str)
    gene_symbol_by_id = adata.var["gene_symbol"].astype(str).to_dict()
    gene_id_by_symbol = adata.var["gene_id"].astype(str).to_dict()
    means = expr.mean(axis=0)
    out = pd.DataFrame({
        "section_id": section,
        "condition": condition,
        "method": "SpatialDE2",
        "gene_id": [gene_id_by_symbol.get(g, g) for g in gene_values],
        "gene_symbol": [gene_symbol_by_id.get(g, g) for g in gene_values],
        "score": pd.to_numeric(df[score_col], errors="coerce") if score_col else np.nan,
        "pvalue": pd.to_numeric(df[p_col], errors="coerce") if p_col else np.nan,
        "qvalue": pd.to_numeric(df[q_col], errors="coerce") if q_col else np.nan,
        "mean_expression": [float(means.get(g, np.nan)) for g in gene_values],
        "n_spots": adata.n_obs,
        "n_genes_tested": adata.n_vars,
        "status": "ok",
        "reason": "",
        "is_ligand_candidate": "no",
        "is_receptor_candidate": "no",
        "is_receiver_target_candidate": "no",
        "is_region_marker_candidate": "no",
        "is_regulator_target_candidate": "no",
        "linked_question_ids": "",
    })
    sort_col = "qvalue" if out["qvalue"].notna().any() else "pvalue" if out["pvalue"].notna().any() else "score"
    ascending = sort_col != "score"
    out = out.sort_values([sort_col, "gene_symbol"], ascending=[ascending, True])
    out["rank"] = range(1, len(out) + 1)
    return out[SVG_COLUMNS]


def expression_proxy(adata, section: str, condition: str):
    import numpy as np
    import pandas as pd

    x = adata.X
    means = np.asarray(x.mean(axis=0)).reshape(-1)
    variances = np.asarray(x.var(axis=0)).reshape(-1) if hasattr(x, "var") else np.zeros_like(means)
    gene_symbols = adata.var["gene_symbol"].astype(str).to_numpy() if "gene_symbol" in adata.var else adata.var_names.astype(str)
    gene_ids = adata.var["gene_id"].astype(str).to_numpy() if "gene_id" in adata.var else adata.var_names.astype(str)
    df = pd.DataFrame({
        "section_id": section,
        "condition": condition,
        "method": "SpatialDE2_proxy",
        "gene_id": gene_ids,
        "gene_symbol": gene_symbols,
        "score": variances,
        "pvalue": 1.0,
        "qvalue": 1.0,
        "mean_expression": means,
        "n_spots": adata.n_obs,
        "n_genes_tested": adata.n_vars,
        "status": "ok_proxy",
        "reason": "SPATIALDE2_ALLOW_PROXY enabled; ranking by expression variance proxy",
        "is_ligand_candidate": "no",
        "is_receptor_candidate": "no",
        "is_receiver_target_candidate": "no",
        "is_region_marker_candidate": "no",
        "is_regulator_target_candidate": "no",
        "linked_question_ids": "",
    })
    df = df.sort_values(["score", "mean_expression", "gene_symbol"], ascending=[False, False, True])
    df["rank"] = range(1, len(df) + 1)
    return df[SVG_COLUMNS]


def main() -> int:
    parser = argparse.ArgumentParser(description="H5AD-first SpatialDE2 SVG sidecar.")
    parser.add_argument("--st-h5ad")
    parser.add_argument("--st-module", default="spatial_03_region")
    parser.add_argument("--section-id", default="")
    parser.add_argument("--out-dir", default="results/spatial/tables/09_svg/spatialde2")
    parser.add_argument("--results-dir", default=None)
    parser.add_argument("--allow-proxy", action="store_true", default=False)
    args = parser.parse_args()
    import os
    args.allow_proxy = args.allow_proxy or os.environ.get("SPATIALDE2_ALLOW_PROXY", "").lower() in {"yes", "true", "1", "on"}

    try:
        h5ad = resolve_spatial_h5ad_path_strict(module=args.st_module, section_id=args.section_id or None, base_dir=args.results_dir, fallback_path=args.st_h5ad)
    except Exception as exc:
        write_manifest(args, "skipped_no_h5ad", str(exc))
        return 20
    try:
        import anndata as ad
    except Exception as exc:
        write_manifest(args, "skipped_no_python_env", f"anndata unavailable: {exc}", h5ad)
        return 21

    try:
        adata = ad.read_h5ad(h5ad)
        if "spatial" not in adata.obsm:
            raise ValueError('input H5AD missing obsm["spatial"]')
        section = args.section_id or "all"
        for col in ("section_id", "condition"):
            if col not in adata.obs:
                raise ValueError(f"input H5AD missing obs.{col}")
        for col in ("gene_symbol", "gene_id"):
            if col not in adata.var:
                raise ValueError(f"input H5AD missing var.{col}")
        if args.section_id and "section_id" in adata.obs:
            keep = adata.obs["section_id"].astype(str) == args.section_id
            adata = adata[keep].copy()
        condition_values = sorted({str(x) for x in adata.obs["condition"].dropna().astype(str)})
        condition = condition_values[0] if len(condition_values) == 1 else ("mixed" if condition_values else "unknown")
        if adata.n_obs < 3 or adata.n_vars == 0:
            out_dir = Path(args.out_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            svg_tsv = out_dir / f"{section}_svg.tsv"
            empty_svg_df(section, condition, "skipped_too_few_spots", "need at least 3 spots and 1 gene", "SpatialDE2").to_csv(svg_tsv, sep="\t", index=False)
            write_manifest(args, "skipped_too_few_spots", "need at least 3 spots and 1 gene", h5ad, str(svg_tsv))
            return 22

        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        svg_tsv = out_dir / f"{section}_svg.tsv"
        if not spatialde2_available():
            if not args.allow_proxy:
                empty_svg_df(section, condition, "skipped_no_spatialde2", "SpatialDE2 Python package unavailable; set SPATIALDE2_ALLOW_PROXY=yes for exploratory proxy", "SpatialDE2").to_csv(svg_tsv, sep="\t", index=False)
                write_manifest(args, "skipped_no_spatialde2", "SpatialDE2 Python package unavailable", h5ad, str(svg_tsv), "SpatialDE2")
                return 23
            df = expression_proxy(adata, section, condition)
            df.to_csv(svg_tsv, sep="\t", index=False)
            write_manifest(args, "ok_proxy", "SpatialDE2 unavailable; proxy enabled", h5ad, str(svg_tsv), "SpatialDE2_proxy")
            return 0

        df = formal_spatialde2(adata, section, condition)
        df.to_csv(svg_tsv, sep="\t", index=False)
    except Exception as exc:
        write_manifest(args, "failed_sidecar", str(exc), h5ad)
        return 30

    write_manifest(args, "ok", "", h5ad, str(svg_tsv), "SpatialDE2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
