"""Estimate COSMIC SBS activities per tumor using non-negative least squares.

Restricted to the lung-relevant signature set (no 79-signature overfit) and
opportunity-adjusted to the coding observation space (see refit_signatures).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import nnls

from refit_signatures import opportunity_adjust, SIGS

ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / "data" / "processed" / "tumor_context_counts.csv"
MATRIX = ROOT / "data" / "raw" / "COSMIC_v3.3.1_SBS_GRCh37.txt"
OPP = ROOT / "data" / "processed" / "opportunity_96.csv"
OUT = ROOT / "data" / "processed" / "tumor_signature_activities.csv"
OUT.parent.mkdir(parents=True, exist_ok=True)

counts = pd.read_csv(COUNTS)
mat = pd.read_csv(MATRIX, sep="\t")
mat = mat.set_index("Type")
missing = [s for s in SIGS if s not in mat.columns]
if missing:
    raise ValueError(f"COSMIC matrix lacks signatures: {missing}")
mat = opportunity_adjust(mat[SIGS], OPP)

# Align contexts exactly.
contexts = list(mat.index)
pivot = counts.pivot_table(
    index="Tumor_Sample_Barcode", columns="context_96", values="count",
    aggfunc="sum", fill_value=0,
)
pivot = pivot.reindex(columns=contexts, fill_value=0)

# Keep tumor metadata from counts table.
meta = counts.drop_duplicates("Tumor_Sample_Barcode").set_index("Tumor_Sample_Barcode")

A = mat.to_numpy(dtype=float)
rows = []
for sample, vec in pivot.iterrows():
    b = vec.to_numpy(dtype=float)
    coef, _ = nnls(A, b)
    total = float(b.sum())
    act = dict(zip(mat.columns, coef))
    rows.append({
        "Tumor_Sample_Barcode": sample,
        "dataset": meta.loc[sample, "dataset"],
        "total_mutations": total,
        "SBS4": act.get("SBS4", 0.0),
        "SBS40": act.get("SBS40", 0.0),
        "SBS5": act.get("SBS5", 0.0),
        "SBS1": act.get("SBS1", 0.0),
        "SBS2": act.get("SBS2", 0.0),
        "SBS13": act.get("SBS13", 0.0),
        "SBS31": act.get("SBS31", 0.0),
        "SBS35": act.get("SBS35", 0.0),
        "SBS4_fraction": (act.get("SBS4", 0.0) / total) if total else 0.0,
    })

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
print("Samples:", len(out))
print("SBS4 > 0:", (out.SBS4 > 0).sum())
print("SBS4 median:", out.SBS4.median())
print("Saved:", OUT)
