"""Permutation test for O/E enrichment of dark-region mutations.

Null: mutations fall on random coding positions with the tumor's own context
spectrum. Within each tumor, for each SBS96 channel c, the number of dark hits
is redrawn as Binomial(n_c, f_dark(c)) — i.e. positions are reshuffled among
coding sites while mutation types are preserved. The empirical p-value compares
the observed statistic to the permutation distribution.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", required=True, help="annotated mutations with context_96 and is_dark")
    parser.add_argument("--background", required=True, help="CSV: context_96,dark_fraction")
    parser.add_argument("--perms", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mut = pd.read_csv(args.mutations, low_memory=False)
    bg = pd.read_csv(args.background)
    bg_map = dict(zip(bg["context_96"], bg["dark_fraction"]))

    n_all = len(mut)
    mut = mut[mut["is_dark"].notna()].copy()
    mut["p_dark_bg"] = mut["context_96"].map(bg_map)
    n_ctx = int(mut["p_dark_bg"].isna().sum())
    mut = mut[mut["p_dark_bg"].notna()]
    print(f"Mutations: {n_all} total, {len(mut)} used, {n_all - len(mut)} dropped (NaN darkness or unmapped context)")
    if n_ctx:
        print(f"WARNING: {n_ctx} with context_96 absent from background")

    counts = (mut.groupby(["Tumor_Sample_Barcode", "context_96"])
              .agg(n=("is_dark", "size"), o=("is_dark", "sum"))
              .reset_index())
    counts["p"] = counts["context_96"].map(bg_map)

    obs_oe = counts["o"].sum() / (counts["n"] * counts["p"]).sum()
    rng = np.random.default_rng(args.seed)
    null = np.empty(args.perms)
    for i in range(args.perms):
        o_perm = rng.binomial(counts["n"].values, counts["p"].values).sum()
        null[i] = o_perm / (counts["n"] * counts["p"]).sum()

    p_hi = (1 + np.sum(null >= obs_oe)) / (args.perms + 1)
    p_lo = (1 + np.sum(null <= obs_oe)) / (args.perms + 1)
    print(f"Observed global O/E: {obs_oe:.4f}")
    print(f"Null mean {null.mean():.4f} sd {null.std():.4f}")
    print(f"Empirical p (enrichment): {p_hi:.4g}, p (depletion): {p_lo:.4g}")

    out = pd.DataFrame({"perm_oe": null})
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
