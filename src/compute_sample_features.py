"""Compute per-sample dark-mutation fraction and merge with signature activities."""
from pathlib import Path
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", required=True, help="annotated mutations CSV with is_dark")
    parser.add_argument("--signatures", required=True, help="tumor_signature_activities.csv")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mut = pd.read_csv(args.mutations, low_memory=False)
    sig = pd.read_csv(args.signatures, low_memory=False)

    if "is_dark" not in mut.columns:
        raise ValueError("Mutations file must contain is_dark column")

    n_all = len(mut)
    mut = mut[mut["is_dark"].notna()].copy()
    print(f"Mutations: {n_all} total, {len(mut)} with darkness label, {n_all - len(mut)} dropped (is_dark NaN)")

    agg = (
        mut.groupby("Tumor_Sample_Barcode")
        .agg(
            n_mutations=("is_dark", "size"),
            n_dark=("is_dark", "sum"),
        )
        .reset_index()
    )
    agg["dark_fraction"] = agg["n_dark"] / agg["n_mutations"].replace(0, pd.NA)

    out = agg.merge(sig, on="Tumor_Sample_Barcode", how="left", validate="one_to_one")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Samples:", len(out))
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
