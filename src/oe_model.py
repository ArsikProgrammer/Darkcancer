"""Poisson GEE for dark-mutation counts with offset=log(E), clustered by patient.

Robustness variant of the primary model (see oe_analysis.py): same terms
(SBS4_std + TMB_std + wes, cancer_type/interaction only when >=2 types),
Exchangeable working correlation. exp(beta_SBS4_std) = O/E ratio per 1 SD.
Cluster = patient_id when available in --signatures, else sample.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.genmod.families import Poisson
from statsmodels.genmod.cov_struct import Exchangeable


def cancer_type_of(dataset):
    d = str(dataset).lower()
    if "luad" in d:
        return "LUAD"
    if "lusc" in d:
        return "LUSC"
    return "SCLC"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", required=True)
    parser.add_argument("--background", required=True)
    parser.add_argument("--signatures", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mut = pd.read_csv(args.mutations, low_memory=False)
    bg = pd.read_csv(args.background)
    sig = pd.read_csv(args.signatures, low_memory=False)

    bg_map = dict(zip(bg["context_96"], bg["dark_fraction"]))
    if bg["dark_fraction"].isna().any() or ((bg["dark_fraction"] <= 0) | (bg["dark_fraction"] >= 1)).any():
        raise ValueError("background dark_fraction must be in (0,1) for all contexts")

    n_all = len(mut)
    mut = mut[mut["is_dark"].notna()].copy()
    print(f"Mutations: {n_all} total, {len(mut)} with darkness label, {n_all - len(mut)} dropped (is_dark NaN)")

    mut["p_dark_bg"] = mut["context_96"].map(bg_map)
    n_ctx = int(mut["p_dark_bg"].isna().sum())
    if n_ctx:
        print(f"WARNING: dropping {n_ctx} mutations with context_96 absent from background")
        mut = mut[mut["p_dark_bg"].notna()]

    agg = (
        mut.groupby("Tumor_Sample_Barcode")
        .agg(n_mutations=("is_dark", "size"), O=("is_dark", "sum"), E=("p_dark_bg", "sum"))
        .reset_index()
    )

    sig_cols = [c for c in ["Tumor_Sample_Barcode", "SBS4_fraction", "SBS4", "total_mutations",
                            "dataset", "patient_id"] if c in sig.columns]
    out = agg.merge(sig[sig_cols], on="Tumor_Sample_Barcode", how="left", validate="one_to_one")
    if "SBS4_fraction" not in out.columns:
        out["SBS4_fraction"] = out["SBS4"] / out["total_mutations"].replace(0, np.nan)
    out["SBS4_std"] = (out["SBS4_fraction"] - out["SBS4_fraction"].mean()) / out["SBS4_fraction"].std()
    out["cancer_type"] = out["dataset"].map(cancer_type_of) if "dataset" in out.columns else "SCLC"
    out["TMB_std"] = (out["total_mutations"] - out["total_mutations"].mean()) / out["total_mutations"].std()
    out["wes"] = out["dataset"].astype(str).str.lower().isin(["jhu", "clcgp"]).astype(float)

    df = out.dropna(subset=["O", "E", "SBS4_std", "total_mutations"]).copy()
    groups = df["patient_id"] if "patient_id" in df.columns else df["Tumor_Sample_Barcode"]
    n_clusters = groups.nunique()
    print(f"Clusters: {n_clusters} ({'patient' if 'patient_id' in df.columns else 'sample'})")

    terms = ["SBS4_std"]
    if df["cancer_type"].nunique() > 1:
        terms.append("C(cancer_type, Treatment('SCLC'))")
        terms.append("SBS4_std:C(cancer_type, Treatment('SCLC'))")
    terms.append("TMB_std")
    if 0 < df["wes"].sum() < len(df):
        terms.append("wes")
    formula = "O ~ " + " + ".join(terms)
    model = smf.gee(formula, groups=groups, data=df, family=Poisson(),
                    cov_struct=Exchangeable(), offset=np.log(df["E"])).fit()
    print(model.summary().as_text()[:1500])

    coefs = pd.DataFrame({
        "term": model.params.index,
        "coef": model.params.values,
        "oe_ratio": np.exp(model.params.values),
        "se": model.bse.values,
        "pvalue": model.pvalues.values,
    })
    if "SBS4_std" in model.params.index:
        b, se = model.params["SBS4_std"], model.bse["SBS4_std"]
        lo, hi = np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)
        print(f"\nSBS4 (per 1 SD): O/E ratio = {np.exp(b):.4f} (95% CI: {lo:.4f}, {hi:.4f})")
        eq_lo, eq_hi = 0.833, 1.2
        verdict = ("within equivalence bounds" if lo >= eq_lo and hi <= eq_hi
                   else "excludes 1" if (lo > 1 or hi < 1)
                   else "inconclusive vs equivalence bounds")
        print(f"Equivalence bounds ({eq_lo}, {eq_hi}): {verdict}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    coefs.to_csv(args.output, index=False)
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
