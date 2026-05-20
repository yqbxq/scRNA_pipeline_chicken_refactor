#!/usr/bin/env python3
import argparse
import itertools
import math
import sys


def fail(code, status, reason):
    sys.stderr.write(f"{status}\t{reason}\n")
    return code


def main():
    parser = argparse.ArgumentParser(description="Spatial neighborhood sidecar for ST06.")
    parser.add_argument("--spots", required=True, help="TSV with spot_id, x, y, label columns")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--radius", type=float, default=200.0)
    parser.add_argument("--perms", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    try:
        import squidpy  # noqa: F401
    except Exception as exc:
        return fail(20, "skipped_no_packages", f"squidpy import failed: {exc}")

    try:
        import numpy as np
        import pandas as pd
    except Exception as exc:
        return fail(21, "skipped_no_python_env", f"numpy/pandas import failed: {exc}")

    spots = pd.read_csv(args.spots, sep="\t")
    required = {"spot_id", "x", "y", "label"}
    missing = required.difference(spots.columns)
    if missing:
        return fail(30, "failed_sidecar", f"spots TSV missing columns: {','.join(sorted(missing))}")
    spots = spots.dropna(subset=["x", "y", "label"]).copy()
    spots["label"] = spots["label"].astype(str)
    labels = sorted([x for x in spots["label"].unique() if x])
    if len(spots) < 3 or len(labels) < 2:
        return fail(22, "skipped_too_few_spots", "need at least 3 spots and 2 labels")

    coords = spots[["x", "y"]].to_numpy(dtype=float)
    label_vec = spots["label"].to_numpy()
    n = len(spots)
    observed = pd.DataFrame(0, index=labels, columns=labels, dtype=float)
    distance_sum = pd.DataFrame(0.0, index=labels, columns=labels, dtype=float)

    for i, j in itertools.combinations(range(n), 2):
        dist = math.dist(coords[i], coords[j])
        if dist <= args.radius:
            a = label_vec[i]
            b = label_vec[j]
            observed.loc[a, b] += 1
            observed.loc[b, a] += 1
            distance_sum.loc[a, b] += dist
            distance_sum.loc[b, a] += dist

    rng = np.random.default_rng(args.seed)
    expected_sum = pd.DataFrame(0.0, index=labels, columns=labels)
    expected_sq = pd.DataFrame(0.0, index=labels, columns=labels)
    perms = max(1, int(args.perms))
    for _ in range(perms):
        perm = rng.permutation(label_vec)
        mat = pd.DataFrame(0, index=labels, columns=labels, dtype=float)
        for i, j in itertools.combinations(range(n), 2):
            if math.dist(coords[i], coords[j]) <= args.radius:
                a = perm[i]
                b = perm[j]
                mat.loc[a, b] += 1
                mat.loc[b, a] += 1
        expected_sum += mat
        expected_sq += mat * mat

    expected = expected_sum / perms
    variance = (expected_sq / perms) - (expected * expected)
    zscore = (observed - expected) / np.sqrt(variance.replace(0, np.nan))
    zscore = zscore.replace([np.inf, -np.inf], np.nan).fillna(0)
    mean_distance = distance_sum / observed.replace(0, np.nan)
    mean_distance = mean_distance.replace([np.inf, -np.inf], np.nan).fillna(0)

    observed.to_csv(f"{args.out_prefix}_interaction_matrix.tsv", sep="\t")
    zscore.to_csv(f"{args.out_prefix}_nhood_enrichment_zscore.tsv", sep="\t")
    mean_distance.to_csv(f"{args.out_prefix}_co_occurrence.tsv", sep="\t")
    pd.DataFrame(
        [{"status": "ok", "spot_n": n, "label_n": len(labels), "radius": args.radius, "perms": perms}]
    ).to_csv(f"{args.out_prefix}_summary.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
