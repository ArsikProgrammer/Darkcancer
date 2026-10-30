"""Raw (unadjusted) association between smoke signature and dark-region share.

Descriptive analysis without the 96-context correction: per-sample fraction of
mutations in dark regions vs SBS4 fraction. Pearson and Spearman correlations
are reported for each cohort. Not the pre-registered test (see ANALYSIS_PLAN).
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

COHORTS = {
    "SCLC": (["annotated_george.csv", "annotated_jhu.csv", "annotated_clcgp.csv"],
             "tumor_signature_activities_refit_dedup.csv"),
    "LUAD": (["annotated_luad.csv"], "luad_signature_activities.csv"),
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed", default="data/processed")
    parser.add_argument("--output", default=None, help="save the report to this file")
    args = parser.parse_args()
    proc = Path(args.processed)

    lines = []
    def log(s=""):
        print(s)
        lines.append(s)

    log("Raw (unadjusted) association: dark mutation share vs SBS4 fraction")
    for name, (mut_files, sig_file) in COHORTS.items():
        frames = []
        for f in mut_files:
            d = pd.read_csv(proc / f, low_memory=False)
            d["Tumor_Sample_Barcode"] = d["Tumor_Sample_Barcode"].astype(str)
            frames.append(d)
        mut = pd.concat(frames, ignore_index=True)
        mut = mut[mut["is_dark"].notna()]
        per = mut.groupby("Tumor_Sample_Barcode")["is_dark"].mean().rename("dark_fraction").reset_index()
        sig = pd.read_csv(proc / sig_file, low_memory=False)
        sig["Tumor_Sample_Barcode"] = sig["Tumor_Sample_Barcode"].astype(str)
        df = per.merge(sig, on="Tumor_Sample_Barcode", how="inner")
        r, p = pearsonr(df["SBS4_fraction"], df["dark_fraction"])
        rs, ps = spearmanr(df["SBS4_fraction"], df["dark_fraction"])
        log(f"{name}: n = {len(df)} samples")
        log(f"  Pearson  r = {r:+.3f} (p = {p:.3f})")
        log(f"  Spearman rho = {rs:+.3f} (p = {ps:.3f})")

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        print("Saved:", args.output)


if __name__ == "__main__":
    main()
