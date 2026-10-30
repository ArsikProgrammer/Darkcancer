"""Build per-tumor 96-context counts from the three context tables."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed" / "tumor_context_counts.csv"
OUT.parent.mkdir(parents=True, exist_ok=True)

frames = []
for name, dataset in [
    ("george", "george"),
    ("sclc_clcgp", "clcgp"),
    ("sclc_jhu", "jhu"),
]:
    p = ROOT / "data" / "processed" / f"{name}_coding_snps_context.csv"
    d = pd.read_csv(p, low_memory=False)
    counts = (
        d.groupby(["Tumor_Sample_Barcode", "context_96"])
        .size()
        .rename("count")
        .reset_index()
    )
    counts["dataset"] = dataset
    frames.append(counts)

all_counts = pd.concat(frames, ignore_index=True)
all_counts.to_csv(OUT, index=False)
print("Rows:", len(all_counts))
print("Samples:", all_counts.Tumor_Sample_Barcode.nunique())
print("Datasets:", all_counts.dataset.value_counts().to_dict())
print("Saved:", OUT)
