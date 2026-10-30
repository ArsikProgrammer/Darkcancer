"""Estimate SBS4 fraction from a context CSV."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "data" / "raw" / "COSMIC_v3.3.1_SBS_GRCh37.txt"
SIGS = ["SBS1", "SBS2", "SBS4", "SBS5", "SBS13", "SBS17a", "SBS40"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, low_memory=False)
    mat = pd.read_csv(MATRIX, sep="\t").set_index("Type")[SIGS]
    contexts = list(mat.index)
    A = mat.to_numpy(dtype=float)

    counts = df.groupby(["Tumor_Sample_Barcode", "context_96"]).size().rename("count").reset_index()
    pivot = counts.pivot_table(index="Tumor_Sample_Barcode", columns="context_96", values="count", aggfunc="sum", fill_value=0)
    pivot = pivot.reindex(columns=contexts, fill_value=0)

    rows = []
    for sample, vec in pivot.iterrows():
        b = vec.to_numpy(dtype=float)
        coef, _ = nnls(A, b)
        total = b.sum()
        act = dict(zip(SIGS, coef))
        rows.append({
            "Tumor_Sample_Barcode": sample,
            "total_mutations": total,
            "SBS4": act["SBS4"],
            "SBS4_fraction": (act["SBS4"] / total) if total else 0.0,
        })
    out = pd.DataFrame(rows)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Samples:", len(out))
    print("Saved:", args.output)

if __name__ == "__main__":
    main()
