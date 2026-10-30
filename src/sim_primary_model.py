"""Simulate the primary model (Poisson, offset=log E) to assess bias and coverage.

DGP deliberately differs from the fitted model: per-sample Gamma frailty adds
overdispersion the GLM does not model, and O is generated as a Poisson count
around E * exp(beta1*SBS4_std + beta2*TMB_std + beta3*wes) * u, u ~ Gamma.
Two scenarios per run set: null (beta1 = 0) and effect (beta1 = log 1.2).
Metrics: bias of beta1_hat, empirical SD vs mean model SE, 95% CI coverage,
rejection rate at p < 0.05, and equivalence-bound classification.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats import t as tdist

EQ_LO, EQ_HI = 0.833, 1.2


def simulate(rng, n, beta1, beta2, beta3, disp, base_e, frailty="gamma"):
    sbs4 = rng.normal(0, 1, n)
    tmb = rng.normal(0, 1, n)
    wes = rng.binomial(1, 0.25, n)
    E = base_e * rng.lognormal(0, 0.5, n)
    lin = np.exp(beta1 * sbs4 + beta2 * tmb + beta3 * wes)
    if disp <= 0:
        lam = E * lin
    elif frailty == "gamma":
        # NB2-type: Var(O) = mu + disp*mu^2
        lam = E * lin * rng.gamma(1 / disp, disp, n)
    elif frailty == "lognormal":
        # Poisson-lognormal: mean-1 lognormal factor, heavier tails than gamma
        s2 = np.log1p(disp)
        lam = E * lin * np.exp(-s2 / 2 + np.sqrt(s2) * rng.normal(0, 1, n))
    elif frailty == "nb1":
        # NB1-type: Var(O) = mu + disp*mu (non-quadratic structure)
        lam = rng.gamma(np.maximum(E * lin, 1e-8) / disp, disp)
    else:
        raise ValueError(frailty)
    O = rng.poisson(lam)
    return pd.DataFrame({"O": O, "E": E, "SBS4_std": sbs4, "TMB_std": tmb, "wes": wes})


def fit(df, cov):
    if cov == "nb2_t":
        # primary: NB2 (Var = mu + theta*mu^2) with t-based CI
        m = sm.NegativeBinomial(df["O"], sm.add_constant(df[["SBS4_std", "TMB_std", "wes"]]),
                                offset=np.log(df["E"])).fit(disp=0)
        return m.params["SBS4_std"], m.bse["SBS4_std"], tdist.ppf(0.975, len(df) - len(m.params))
    if cov == "quasi":
        # quasi-Poisson: model-based SE scaled by sqrt(Pearson chi2 / df)
        model = smf.glm("O ~ SBS4_std + TMB_std + wes", data=df,
                        family=sm.families.Poisson(), offset=np.log(df["E"])).fit(scale="x2")
    else:
        model = smf.glm("O ~ SBS4_std + TMB_std + wes", data=df,
                        family=sm.families.Poisson(), offset=np.log(df["E"])).fit(cov_type=cov)
    return model.params["SBS4_std"], model.bse["SBS4_std"], 1.96


def run_scenario(rng, args, beta1):
    rows = []
    for _ in range(args.sims):
        df = simulate(rng, args.n, beta1, args.beta2, args.beta3, args.disp, args.base_e, args.frailty)
        try:
            b, se, crit = fit(df, args.cov)
        except Exception:
            continue
        if not (np.isfinite(b) and np.isfinite(se)):
            continue
        lo, hi = b - crit * se, b + crit * se
        rows.append({"beta1_hat": b, "se": se, "lo": lo, "hi": hi,
                     "covers": lo <= beta1 <= hi, "rejects": (lo > 0) or (hi < 0),
                     "equiv": lo >= np.log(EQ_LO) and hi <= np.log(EQ_HI)})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sims", type=int, default=500)
    parser.add_argument("--n", type=int, default=143)
    parser.add_argument("--beta2", type=float, default=0.05)
    parser.add_argument("--beta3", type=float, default=0.1)
    parser.add_argument("--disp", type=float, default=0.15, help="overdispersion parameter theta")
    parser.add_argument("--frailty", default="gamma", choices=["gamma", "lognormal", "nb1"])
    parser.add_argument("--cov", default="nb2_t", help="SE method: nb2_t (primary), quasi, HC1/HC2/HC3")
    parser.add_argument("--base_e", type=float, default=60.0)
    parser.add_argument("--or_effect", type=float, default=1.2)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", default=None, help="per-run CSV (synthetic; keep outside data/)")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    all_rows = []
    for name, beta1 in [("H0 (beta1=0)", 0.0), (f"effect (O/E={args.or_effect})", np.log(args.or_effect))]:
        res = run_scenario(rng, args, beta1)
        res["scenario"] = name
        all_rows.append(res)
        if len(res) < args.sims:
            print(f"WARNING: only {len(res)}/{args.sims} converged fits in {name}")
        b_true = beta1
        print(f"\n=== {name}: {len(res)} runs, n={args.n}, disp={args.disp} ===")
        print(f"bias:   mean beta1_hat = {res.beta1_hat.mean():.4f} (true {b_true:.4f}), "
              f"bias = {res.beta1_hat.mean() - b_true:+.4f}")
        print(f"spread: empirical SD = {res.beta1_hat.std():.4f}, mean model SE = {res.se.mean():.4f}")
        print(f"coverage of 95% CI: {res.covers.mean():.3f}")
        print(f"rejection rate (p<0.05): {res.rejects.mean():.3f}")
        print(f"CI fully within equivalence bounds (log scale): {res.equiv.mean():.3f}")

    out = pd.concat(all_rows, ignore_index=True)
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.output, index=False)
        print("\nSaved:", args.output)


if __name__ == "__main__":
    main()
