"""Check sample overlap between cohorts by mutation profiles."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

def load(name):
    df = pd.read_csv(ROOT / "data" / "processed" / name, low_memory=False)
    df["key"] = df["Chromosome"].astype(str) + ":" + df["Start_Position"].astype(str) + ":" + df["Reference_Allele"].astype(str) + ":" + df["Tumor_Seq_Allele2"].astype(str)
    return df

george = load("george_coding_snps.csv")
jhu = load("sclc_jhu_coding_snps.csv")
clcgp = load("sclc_clcgp_coding_snps.csv")

def profile_overlap(a, b):
    # Jaccard of variant keys per sample pair
    a_groups = a.groupby("Tumor_Sample_Barcode")["key"].apply(set)
    b_groups = b.groupby("Tumor_Sample_Barcode")["key"].apply(set)
    max_jaccard = 0
    best_pair = None
    for sa, va in a_groups.items():
        for sb, vb in b_groups.items():
            inter = len(va & vb)
            union = len(va | vb)
            if union == 0:
                continue
            jacc = inter / union
            if jacc > max_jaccard:
                max_jaccard = jacc
                best_pair = (sa, sb, inter, union)
    return max_jaccard, best_pair

for name_a, da, name_b, db in [("George","george_coding_snps.csv","JHU","sclc_jhu_coding_snps.csv"),
                                ("George","george_coding_snps.csv","CLCGP","sclc_clcgp_coding_snps.csv"),
                                ("JHU","sclc_jhu_coding_snps.csv","CLCGP","sclc_clcgp_coding_snps.csv")]:
    # reload properly
    pass

# simpler: use loaded
for label, da, db in [("George vs JHU", george, jhu), ("George vs CLCGP", george, clcgp), ("JHU vs CLCGP", jhu, clcgp)]:
    mj, pair = profile_overlap(da, db)
    print(f"{label}: max Jaccard = {mj:.4f}, pair = {pair}")
    if mj > 0.5:
        print("  WARNING: possible duplicate sample!")
