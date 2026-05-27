#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import sys
import csv
from datetime import datetime, timezone
from pathlib import Path

try:
    import pandas as pd
except Exception:  # pragma: no cover - local smoke may not have pandas
    pd = None

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow" / "04python"))
from helpers.lr_axis_id_utils import standardize_lr_axis_id  # noqa: E402


def env_path(name: str, default: str) -> Path:
    return Path(os.environ.get(name) or default)


def read_tsv(path: Path) -> pd.DataFrame:
    if pd is None:
        raise RuntimeError("pandas is unavailable")
    if not path.exists() or path.stat().st_size == 0:
      return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def read_tsv_records(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv_records(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(empty_output().columns) if pd is not None else [
        "section_id", "lr_axis_id", "pair_id", "condition_value", "ligand",
        "receptor", "sender", "receiver", "signal_score", "spatial_support",
        "status", "reason",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def relative(path: Path, base: Path) -> str:
    try:
        return str(path.resolve().relative_to(base.resolve()))
    except ValueError:
        return str(path)


def write_manifest(manifest_path: Path, output_path: Path, base_dir: Path, input_manifest: Path, lr_candidates: Path) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "module": "spatial_08f_commot_spatial",
        "version": os.environ.get("MODULE_SPATIAL_07_VERSION", "1.0"),
        "timestamp": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "base_dir": str(base_dir),
        "inputs": {
            "commot_input_manifest": str(input_manifest),
            "commot_lr_candidates_tsv": str(lr_candidates),
            "commot_distance_threshold": os.environ.get("COMMOT_DISTANCE_THRESHOLD", "500"),
            "commot_permutations": os.environ.get("COMMOT_PERMUTATIONS", "100"),
        },
        "outputs": {
            "commot_lr_tsv": {
                "path": relative(output_path, base_dir),
                "type": "tsv",
                "produced_by": "spatial_08f_commot_spatial",
                "row_semantics": "one row per section and LR axis evaluated by COMMOT",
            }
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def empty_output() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "section_id",
            "lr_axis_id",
            "pair_id",
            "condition_value",
            "ligand",
            "receptor",
            "sender",
            "receiver",
            "signal_score",
            "spatial_support",
            "status",
            "reason",
        ]
    )


def candidate_rows_for_status(input_df: pd.DataFrame, candidates: pd.DataFrame, status: str, reason: str) -> pd.DataFrame:
    rows = []
    if candidates.empty:
        candidates = pd.DataFrame([{"lr_axis_id": "", "ligand": "", "receptor": "", "sender": "", "receiver": "", "pair_id": "", "condition_value": ""}])
    if input_df.empty:
        input_df = pd.DataFrame([{"section_id": "", "status": status, "reason": reason}])
    for _, section in input_df.iterrows():
        for _, cand in candidates.iterrows():
            rows.append(
                {
                    "section_id": section.get("section_id", ""),
                    "lr_axis_id": cand.get("lr_axis_id", ""),
                    "pair_id": cand.get("pair_id", ""),
                    "condition_value": cand.get("condition_value", ""),
                    "ligand": cand.get("ligand", ""),
                    "receptor": cand.get("receptor", ""),
                    "sender": cand.get("sender", ""),
                    "receiver": cand.get("receiver", ""),
                    "signal_score": "",
                    "spatial_support": "no",
                    "status": status,
                    "reason": reason,
                }
            )
    return pd.DataFrame(rows, columns=empty_output().columns)


def import_commot_runtime(require_runtime: bool):
    try:
        import anndata as ad  # noqa: F401
        import commot as ct  # noqa: F401
        return ad, ct, ""
    except Exception as exc:  # pragma: no cover - depends on optional server env
        if require_runtime:
            raise
        return None, None, f"COMMOT/anndata runtime unavailable: {exc}"


def score_from_commot_uns(adata, ligand: str, receptor: str, sender: str, receiver: str) -> float | None:
    needle = f"{ligand}-{receptor}".lower()
    scores = []
    for key, value in getattr(adata, "uns", {}).items():
        if "commot" not in str(key).lower() or needle not in str(key).lower():
            continue
        try:
            if hasattr(value, "loc") and sender in value.index and receiver in value.columns:
                scores.append(float(value.loc[sender, receiver]))
            elif hasattr(value, "to_numpy"):
                arr = value.to_numpy()
                if arr.size:
                    scores.append(float(pd.Series(arr.ravel()).dropna().max()))
        except Exception:
            continue
    return max(scores) if scores else None


def run_commot(input_df: pd.DataFrame, candidates: pd.DataFrame, distance_threshold: float, permutations: int, require_runtime: bool) -> pd.DataFrame:
    if candidates.empty:
        return candidate_rows_for_status(input_df, candidates, "skipped_no_lr_candidates", "No LR candidates were provided.")

    ad, ct, runtime_reason = import_commot_runtime(require_runtime)
    if ad is None or ct is None:
        return candidate_rows_for_status(input_df, candidates, "skipped_missing_commot", runtime_reason)

    rows = []
    ready = input_df[input_df.get("status", "") == "ready"].copy()
    if ready.empty:
        return candidate_rows_for_status(input_df, candidates, "skipped_no_ready_h5ad", "No ready section H5AD input was available.")

    df_ligrec = candidates[["ligand", "receptor"]].copy()
    df_ligrec["pathway"] = candidates.get("pathway", "user")
    df_ligrec = df_ligrec.drop_duplicates()

    for _, section in ready.iterrows():  # pragma: no cover - requires optional COMMOT runtime
        h5ad_path = Path(section.get("h5ad_path", ""))
        label_col = section.get("label_col", "cell_type_main") or "cell_type_main"
        section_id = section.get("section_id", "")
        if not h5ad_path.exists():
            rows.append(candidate_rows_for_status(pd.DataFrame([section]), candidates, "skipped_missing_h5ad", f"Missing H5AD: {h5ad_path}"))
            continue
        adata = ad.read_h5ad(h5ad_path)
        if label_col not in adata.obs:
            rows.append(candidate_rows_for_status(pd.DataFrame([section]), candidates, "skipped_missing_cluster_col", f"Missing adata.obs column: {label_col}"))
            continue
        ct.tl.cluster_communication_spatial_permutation(
            adata,
            df_ligrec=df_ligrec,
            database_name="user",
            dis_thr=distance_threshold,
            clustering=label_col,
            n_permutations=permutations,
            random_seed=int(os.environ.get("RANDOM_SEED", "42")),
            verbose=False,
        )
        for _, cand in candidates.iterrows():
            lr_axis_id = cand.get("lr_axis_id") or standardize_lr_axis_id(cand.get("ligand", ""), cand.get("receptor", ""), cand.get("sender", ""), cand.get("receiver", ""))
            score = score_from_commot_uns(adata, cand.get("ligand", ""), cand.get("receptor", ""), cand.get("sender", ""), cand.get("receiver", ""))
            rows.append(
                {
                    "section_id": section_id,
                    "lr_axis_id": lr_axis_id,
                    "pair_id": cand.get("pair_id", ""),
                    "condition_value": cand.get("condition_value", ""),
                    "ligand": cand.get("ligand", ""),
                    "receptor": cand.get("receptor", ""),
                    "sender": cand.get("sender", ""),
                    "receiver": cand.get("receiver", ""),
                    "signal_score": "" if score is None else score,
                    "spatial_support": "no",
                    "status": "ok_commot_run",
                    "reason": "" if score is not None else "COMMOT completed but no score matrix was found for this LR axis.",
                }
            )
    if not rows:
        return empty_output()
    return pd.concat([item if isinstance(item, pd.DataFrame) else pd.DataFrame([item]) for item in rows], ignore_index=True)[empty_output().columns]


def main() -> int:
    project_root = Path(os.environ.get("PROJECT_ROOT") or os.getcwd())
    results_dir = Path(os.environ.get("RESULTS_DIR") or project_root / "results")
    input_manifest = env_path("COMMOT_INPUT_MANIFEST", str(results_dir / "spatial" / "tables" / "08_commot" / "commot_input_manifest.tsv"))
    lr_candidates = env_path("COMMOT_LR_CANDIDATES_PREPARED_TSV", str(results_dir / "spatial" / "tables" / "08_commot" / "commot_lr_candidates.tsv"))
    output_tsv = env_path("COMMOT_OUTPUT_TSV", str(results_dir / "spatial" / "tables" / "08_commot" / "commot_lr.tsv"))
    manifest_path = env_path("COMMOT_MANIFEST", str(results_dir / "manifests" / "spatial_08f_commot_spatial" / "_manifest.json"))
    distance_threshold = float(os.environ.get("COMMOT_DISTANCE_THRESHOLD", "500") or "500")
    permutations = int(os.environ.get("COMMOT_PERMUTATIONS", "100") or "100")
    require_runtime = (os.environ.get("COMMOT_REQUIRE_RUNTIME", "no").lower() in {"yes", "true", "1", "on"})

    if pd is None:
        input_records = read_tsv_records(input_manifest) or [{"section_id": "", "status": "skipped_missing_pandas", "reason": "pandas is unavailable"}]
        candidate_records = read_tsv_records(lr_candidates) or [{"lr_axis_id": "", "ligand": "", "receptor": "", "sender": "", "receiver": "", "pair_id": "", "condition_value": ""}]
        rows = []
        for section in input_records:
            for cand in candidate_records:
                rows.append({
                    "section_id": section.get("section_id", ""),
                    "lr_axis_id": cand.get("lr_axis_id", ""),
                    "pair_id": cand.get("pair_id", ""),
                    "condition_value": cand.get("condition_value", ""),
                    "ligand": cand.get("ligand", ""),
                    "receptor": cand.get("receptor", ""),
                    "sender": cand.get("sender", ""),
                    "receiver": cand.get("receiver", ""),
                    "signal_score": "",
                    "spatial_support": "no",
                    "status": "skipped_missing_pandas",
                    "reason": "Python pandas is unavailable; COMMOT runtime was not attempted.",
                })
        write_tsv_records(rows, output_tsv)
        row_n = len(rows)
    else:
        input_df = read_tsv(input_manifest)
        candidates = read_tsv(lr_candidates)
        out = run_commot(input_df, candidates, distance_threshold, permutations, require_runtime)
        write_tsv(out, output_tsv)
        row_n = len(out)
    write_manifest(manifest_path, output_tsv, project_root, input_manifest, lr_candidates)
    print(f"08f COMMOT completed rows={row_n} output={output_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
