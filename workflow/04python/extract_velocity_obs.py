#!/usr/bin/env python3

import argparse

import anndata as ad
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    adata = ad.read_h5ad(args.h5ad)
    obs = adata.obs.copy()
    cols = [
        col
        for col in ["latent_time", "velocity_pseudotime", "velocity_length", "velocity_confidence"]
        if col in obs.columns
    ]
    df = obs[cols].copy() if cols else pd.DataFrame(index=obs.index)
    df.insert(0, "cell_id", obs.index.astype(str))
    df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
