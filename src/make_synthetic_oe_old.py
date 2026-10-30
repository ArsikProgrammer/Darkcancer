"""OLD synthetic generator (for bias reproduction)."""
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--odds_ratio", type=float, default=1.5)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    mut = pd.read_csv(ROOT / "data/processed/george_coding_snps.csv", low_memory=False)
    ctx = pd.read_csv(ROOT / "data/processed/george_coding_snps_context.csv", low_memory=False)
    sig = pd.read_csv(ROOT / "data/processed/tumor_signature_activities.csv", low_memory=False)

    sig_map = dict(zip(sig.Tumor_Sample_Barcode, sig.SBS4))
    vals = pd.Series(list(sig_map.values()))
    std = (vals - vals.mean()) / vals.std()
    std_map = dict(zip(sig.Tumor_Sample_Barcode, std))

    base = 0.24
    beta = np.log(args.odds_ratio)

    mut["Start_Position"] = pd.to_numeric(mut["Start_Position"], errors="coerce")
    ctx["Start_Position"] = pd.to_numeric(ctx["Start_Position"], errors="coerce")
    mut = mut.merge(ctx[["Tumor_Sample_Barcode", "Start_Position", "context_96"]], on=["Tumor_Sample_Barcode", "Start_Position"], how="left")
    mut = mut.dropna(subset=["context_96"])

    rng = np.random.default_rng(args.seed)
    probs = []
    for _, r in mut.iterrows():
        sbs4_std = std_map.get(r["Tumor_Sample_Barcode"], 0.0)
        p = 1 / (1 + np.exp(-(np.log(base/(1-base)) + beta * sbs4_std)))
        probs.append(p)
    mut["is_dark"] = rng.binomial(1, probs).astype(int)

    contexts = mut.context_96.unique()
    bg = pd.DataFrame({"context_96": contexts, "dark_fraction": np.random.default_rng(7).uniform(0.15, 0.35, size=len(contexts))})
    bg.to_csv(ROOT / "data/processed/_synth_bg.csv", index=False)

    mut.to_csv(ROOT / "data/processed/_synth_mutations_dark.csv", index=False)
    print(f"OLD OR={args.odds_ratio}, seed={args.seed}: mutations={len(mut)}, dark={mut.is_dark.sum()}")

if __name__ == "__main__":
    main()
