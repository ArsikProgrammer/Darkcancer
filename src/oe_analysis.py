"""Compute per-sample O/E for mutations in IDR and fit the primary model.

Primary (SCLC-only data): O ~ offset(log E) + SBS4_std + TMB_std + wes,
Poisson GLM, HC1 robust SE; SBS4_std = z-scored SBS4_fraction; wes = JHU+CLCGP
(WES sources) vs George. Primary measure: exp(beta_SBS4_std) per 1 SD,
equivalence bounds (0.833, 1.2). cancer_type and SBS4*cancer_type enter only
when >=2 cancer types are present (after TCGA inclusion).
Samples with unknown darkness (is_dark NaN) are excluded from both O and E.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import t as tdist


def cancer_type_of(dataset):
    d = str(dataset).lower()
    if "luad" in d:
        return "LUAD"
    if "lusc" in d:
        return "LUSC"
    return "SCLC"


def build_frame(agg, sig):
    out = agg.merge(sig, on="Tumor_Sample_Barcode", how="left", validate="one_to_one")
    out["cancer_type"] = out["dataset"].map(cancer_type_of) if "dataset" in out.columns else "SCLC"
    out["SBS4_fraction"] = out.get("SBS4_fraction", out["SBS4"] / out["total_mutations"].replace(0, np.nan))
    out["SBS4_std"] = (out["SBS4_fraction"] - out["SBS4_fraction"].mean()) / out["SBS4_fraction"].std()
    out["TMB_std"] = (out["total_mutations"] - out["total_mutations"].mean()) / out["total_mutations"].std()
    out["wes"] = out["dataset"].astype(str).str.lower().isin(["jhu", "clcgp"]).astype(float)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", required=True, help="annotated mutations with context_96 and is_dark")
    parser.add_argument("--background", required=True, help="CSV: context_96,dark_fraction")
    parser.add_argument("--signatures", required=True, help="tumor signature activities")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mut = pd.read_csv(args.mutations, low_memory=False)
    bg = pd.read_csv(args.background)
    sig = pd.read_csv(args.signatures, low_memory=False)

    if "context_96" not in mut.columns or "is_dark" not in mut.columns:
        raise ValueError("Mutations must contain context_96 and is_dark")

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
        .agg(
            n_mutations=("is_dark", "size"),
            O=("is_dark", "sum"),
            E=("p_dark_bg", "sum"),
        )
        .reset_index()
    )
    agg["OE"] = agg["O"] / agg["E"]

    out = build_frame(agg, sig)

    df = out.dropna(subset=["O", "E", "SBS4_std", "total_mutations"])
    terms = ["SBS4_std"]
    if df["cancer_type"].nunique() > 1:
        terms.append("C(cancer_type, Treatment('SCLC'))")
        terms.append("SBS4_std:C(cancer_type, Treatment('SCLC'))")
    terms.append("TMB_std")
    if 0 < df["wes"].sum() < len(df):
        terms.append("wes")
    formula = "O ~ " + " + ".join(terms)
    # primary inference: NB2 with t CI (dispersion-adjusted, Var = mu + theta*mu^2;
    # calibrated in sim_primary_model.py at theta = 0.15 and 0.30)
    model = smf.negativebinomial(formula, data=df, offset=np.log(df["E"])).fit(disp=0)
    crit = tdist.ppf(0.975, model.nobs - len(model.params))
    print(model.summary().as_text()[:1500])
    if "SBS4_std" in model.params.index:
        b, se = model.params["SBS4_std"], model.bse["SBS4_std"]
        lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
        print(f"\nSBS4 (per 1 SD, NB2-t primary): O/E ratio = {np.exp(b):.4f} (95% CI: {lo:.4f}, {hi:.4f})")
        for tag, kwargs in [("quasi-Poisson", {"scale": "x2"}), ("HC1", {"cov_type": "HC1"})]:
            sens = smf.glm(formula, data=df, family=sm.families.Poisson(),
                           offset=np.log(df["E"])).fit(**kwargs)
            b2, s2 = sens.params["SBS4_std"], sens.bse["SBS4_std"]
            print(f"  sensitivity {tag}: O/E = {np.exp(b2):.4f} "
                  f"(CI {np.exp(b2 - 1.96 * s2):.4f}, {np.exp(b2 + 1.96 * s2):.4f})")
        eq_lo, eq_hi = 0.833, 1.2
        verdict = ("within equivalence bounds" if lo >= eq_lo and hi <= eq_hi
                   else "excludes 1" if (lo > 1 or hi < 1)
                   else "inconclusive vs equivalence bounds")
        print(f"Equivalence bounds ({eq_lo}, {eq_hi}): {verdict}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Samples:", len(out))
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
