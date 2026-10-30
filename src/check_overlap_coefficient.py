"""Check sample overlap by overlap coefficient."""
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def load(name):
    df = pd.read_csv(ROOT / "data" / "processed" / name, low_memory=False)
    df["key"] = df["Chromosome"].astype(str) + ":" + df["Start_Position"].astype(str) + ":" + df["Reference_Allele"].astype(str) + ":" + df["Tumor_Seq_Allele2"].astype(str)
    return df

george = load("george_coding_snps.csv")
jhu = load("sclc_jhu_coding_snps.csv")
clcgp = load("sclc_clcgp_coding_snps.csv")

def overlap_coeff(a, b):
    ag = a.groupby("Tumor_Sample_Barcode")["key"].apply(set)
    bg = b.groupby("Tumor_Sample_Barcode")["key"].apply(set)
    best = 0
    pair = None
    for sa, va in ag.items():
        for sb, vb in bg.items():
            inter = len(va & vb)
            denom = min(len(va), len(vb))
            if denom == 0:
                continue
            coef = inter / denom
            if coef > best:
                best = coef
                pair = (sa, sb, inter, len(va), len(vb))
    return best, pair

for label, da, db in [("George vs JHU", george, jhu),
                      ("George vs CLCGP", george, clcgp),
                      ("JHU vs CLCGP", jhu, clcgp)]:
    coef, pair = overlap_coeff(da, db)
    print(f"{label}: overlap coeff = {coef:.4f}, pair = {pair}")
    if coef > 0.5:
        print("  WARNING: possible duplicate sample!")
