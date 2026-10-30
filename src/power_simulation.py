"""Power simulation for Dark Proteome x SBS4 effect.

Simulates per-sample dark counts with a logit model and analyzes them with a
Binomial GLM (same logit link as the DGP). Random effects are on the logit
scale; the cohort effect is drawn per cohort, not per sample. Power is the
fraction of significant fits among converged fits.
"""
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.genmod.families import Binomial

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/processed/sample_features_50_annotated.csv")
    parser.add_argument("--odds_ratio", type=float, default=1.2)
    parser.add_argument("--sims", type=int, default=300)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--coverage", type=float, default=1.0)
    parser.add_argument("--random_effects", action="store_true", default=False)
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    df = df.dropna(subset=["dark_fraction", "SBS4", "n_mutations"])
    sbs4_std = (df["SBS4"] - df["SBS4"].mean()) / df["SBS4"].std()
    base = df["dark_fraction"].mean()
    n = len(df)
    beta = np.log(args.odds_ratio)
    rng = np.random.default_rng(args.seed)

    cohorts = df["dataset"].astype(str).values if "dataset" in df.columns else np.zeros(n, dtype=int)
    uniq = np.unique(cohorts)
    # fraction-scale SD -> logit-scale SD (delta method at base)
    sample_sd = df["dark_fraction"].std() / (base * (1 - base))
    cohort_sd = 0.02 / (base * (1 - base))

    sig = 0
    attempted = 0
    var_fracs = []
    for _ in range(args.sims):
        if args.random_effects:
            sample_re = rng.normal(0, sample_sd, size=n)
            cohort_re_by_level = rng.normal(0, cohort_sd, size=len(uniq))
            cohort_re = pd.Series(cohorts).map(dict(zip(uniq, cohort_re_by_level))).values
            linpred = np.log(base/(1-base)) + beta * sbs4_std.values + sample_re + cohort_re
        else:
            linpred = np.log(base/(1-base)) + beta * sbs4_std.values
        p = 1 / (1 + np.exp(-linpred))
        n_mut = (df["n_mutations"].values * args.coverage).astype(int)
        n_mut = np.clip(n_mut, 1, None)
        n_dark = rng.binomial(n_mut, p)
        frac = n_dark / n_mut
        var_fracs.append(frac.var())
        X = sm.add_constant(sbs4_std.values)
        try:
            model = sm.GLM(frac, X, family=Binomial(), var_weights=n_mut.astype(float)).fit()
            attempted += 1
            if model.pvalues[1] < 0.05:
                sig += 1
        except Exception:
            continue
    power = sig / attempted if attempted else float("nan")
    print(f"OR={args.odds_ratio}, sims={args.sims}, converged={attempted}, n={n}, "
          f"coverage={args.coverage}, RE={args.random_effects}")
    print(f"Power (p<0.05): {power:.3f}")
    print(f"Base dark fraction: {base:.3f}")
    print(f"Mean simulated frac variance: {np.mean(var_fracs):.4f}")
    if power >= 0.999:
        print("WARNING: power is suspiciously 1.000; check simulation design and coverage")


if __name__ == "__main__":
    main()
