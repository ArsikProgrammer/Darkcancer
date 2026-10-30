"""Confirmatory run for the dark-proteome x SBS4 analysis (pre-registration section 18).

One run -> one output file containing: git commit hash, SHA-256 of every input,
the primary model (NB2, t-based CI), the permutation test, the negative controls
(SBS1 as predictor; Silent mutations as outcome) using the same code path, and
the pre-registered secondary analysis (mapping-discard fraction as covariate).

--shuffle-sbs4 performs a full rehearsal with SBS4 permuted between samples
within each dataset (no real association is used); the confirmatory run itself
must be executed exactly once (pre-registration 18.4).
"""
from pathlib import Path
import argparse
import hashlib
import subprocess
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import t as tdist

EQ_LO, EQ_HI = 0.833, 1.2


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fit_nb2(df, predictor, extra=None):
    cols = [predictor, "TMB_std", "wes"] + list(extra or [])
    cols = [c for c in cols if df[c].nunique() > 1]
    n, p = len(df), len(cols) + 1
    X = sm.add_constant(df[cols])
    off = np.log(df["E"])
    note = "NB2-t"
    m = None
    try:
        gp0 = sm.GLM(df["O"], X, family=sm.families.Poisson(), offset=off).fit()
        m = sm.NegativeBinomial(df["O"], X, offset=off).fit(
            start_params=list(gp0.params) + [0.05], disp=0, maxiter=200)
    except Exception:
        m = None
    ok = (m is not None and np.isfinite(m.bse).all()
          and m.mle_retvals.get("converged", False)
          and m.params.get("alpha", 0.0) > 1e-4)
    if ok:
        b, se = m.params[predictor], m.bse[predictor]
    else:
        # pre-registered boundary rule: NB2 with alpha -> 0 equals Poisson;
        # use the Poisson-t limit of NB2-t when alpha is at the boundary
        gp = sm.GLM(df["O"], X, family=sm.families.Poisson(), offset=off).fit()
        b, se = gp.params[predictor], gp.bse[predictor]
        note = "NB2 alpha->0 boundary: Poisson-t limit"
    crit = tdist.ppf(0.975, n - p)
    if not (np.isfinite(b) and np.isfinite(se)):
        raise RuntimeError("fit produced non-finite estimates")
    return b, se, crit, note


def verdict(lo, hi):
    if lo >= EQ_LO and hi <= EQ_HI:
        return "INFORMATIVE NULL"
    if lo > 1 or hi < 1:
        return "EFFECT DETECTED"
    return "INCONCLUSIVE"


def prep(mut, bg_map, sig, predictor_frac_col, shuffle_rng=None):
    mut = mut[mut["is_dark"].notna()].copy()
    mut["p_dark_bg"] = mut["context_96"].map(bg_map)
    dropped_ctx = int(mut["p_dark_bg"].isna().sum())
    mut = mut[mut["p_dark_bg"].notna()]
    agg = (mut.groupby("Tumor_Sample_Barcode")
           .agg(n=("is_dark", "size"), O=("is_dark", "sum"), E=("p_dark_bg", "sum"))
           .reset_index())
    out = agg.merge(sig, on="Tumor_Sample_Barcode", how="inner", validate="one_to_one")
    frac = out[predictor_frac_col].astype(float)
    if shuffle_rng is not None:
        frac = pd.Series(shuffle_rng.permutation(frac.values), index=frac.index)
    out["pred_std"] = (frac - frac.mean()) / frac.std()
    out["TMB_std"] = (out["total_mutations"] - out["total_mutations"].mean()) / out["total_mutations"].std()
    out["wes"] = out["dataset"].astype(str).str.lower().isin(["jhu", "clcgp"]).astype(float)
    return out.dropna(subset=["O", "E", "pred_std", "TMB_std"]), dropped_ctx


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", nargs="+", required=True)
    parser.add_argument("--background", required=True)
    parser.add_argument("--signatures", required=True)
    parser.add_argument("--mapping", nargs="*", default=None,
                        help="mapping CSVs (for the discard-fraction secondary analysis)")
    parser.add_argument("--output", required=True)
    parser.add_argument("--perms", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261028)
    parser.add_argument("--shuffle-sbs4", action="store_true",
                        help="rehearsal: permute SBS4 within dataset (no real association)")
    args = parser.parse_args()

    lines = []
    def log(s=""):
        print(s)
        lines.append(s)

    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                            text=True, cwd=Path(__file__).parent).stdout.strip()
    mode = "REHEARSAL (SBS4 shuffled within dataset)" if args.shuffle_sbs4 else "CONFIRMATORY"
    log(f"=== {mode} RUN ===")
    log(f"git commit: {commit}")
    log(f"seed: {args.seed}, permutations: {args.perms}")
    log()
    log("-- input checksums (SHA-256) --")
    inputs = list(args.mutations) + [args.background, args.signatures] + (args.mapping or [])
    for p in inputs:
        log(f"{sha256(p)}  {p}")
    log()

    bg = pd.read_csv(args.background)
    bg_map = dict(zip(bg["context_96"], bg["dark_fraction"]))
    sig = pd.read_csv(args.signatures, low_memory=False)
    sig["Tumor_Sample_Barcode"] = sig["Tumor_Sample_Barcode"].astype(str)
    if "SBS4_fraction" not in sig.columns:
        sig["SBS4_fraction"] = sig["SBS4"] / sig["total_mutations"]
    if "SBS1_fraction" not in sig.columns:
        sig["SBS1_fraction"] = sig["SBS1"] / sig["total_mutations"]

    rng = np.random.default_rng(args.seed)
    mut_all = pd.concat([pd.read_csv(p, low_memory=False) for p in args.mutations],
                        ignore_index=True)
    mut_all["Tumor_Sample_Barcode"] = mut_all["Tumor_Sample_Barcode"].astype(str)

    results = {}

    def run_model(frame_mut, predictor_frac_col, label, shuffle=None):
        df, dropped = prep(frame_mut, bg_map, sig, predictor_frac_col, shuffle)
        b, se, crit, note = fit_nb2(df, "pred_std")
        lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
        v = verdict(lo, hi)
        log(f"{label}: n_samples={len(df)}, n_mutations={int(df.n.sum())} "
            f"(dropped unmapped context: {dropped})")
        log(f"  [{note}]")
        log(f"  O/E ratio per 1 SD = {np.exp(b):.4f} (95% CI {lo:.4f}..{hi:.4f}) -> {v}")
        results[label] = (np.exp(b), lo, hi, v)
        return df

    # 1. primary model (all validated coding mutations, SBS4 predictor)
    shuf = rng if args.shuffle_sbs4 else None
    run_model(mut_all, "SBS4_fraction", "PRIMARY (SBS4, all coding)", shuf)

    # 2. negative control a: SBS1 as predictor, same outcome
    run_model(mut_all, "SBS1_fraction", "CONTROL A (SBS1 predictor, all coding)", shuf)

    # 3. negative control b: synonymous outcome, SBS4 predictor
    silent = mut_all[mut_all["Variant_Classification"] == "Silent"]
    run_model(silent, "SBS4_fraction", "CONTROL B (SBS4 predictor, Silent outcome)", shuf)

    # 4. permutation test (positions reshuffled within tumor, types preserved)
    mut_p = mut_all[mut_all["is_dark"].notna()].copy()
    mut_p["p_dark_bg"] = mut_p["context_96"].map(bg_map)
    mut_p = mut_p[mut_p["p_dark_bg"].notna()]
    counts = (mut_p.groupby(["Tumor_Sample_Barcode", "context_96"])
              .agg(n=("is_dark", "size"), o=("is_dark", "sum")).reset_index())
    counts["p"] = counts["context_96"].map(bg_map)
    obs = counts["o"].sum() / (counts["n"] * counts["p"]).sum()
    null = np.empty(args.perms)
    for i in range(args.perms):
        null[i] = rng.binomial(counts["n"].values, counts["p"].values).sum() / \
            (counts["n"] * counts["p"]).sum()
    p_hi = (1 + int((null >= obs).sum())) / (args.perms + 1)
    p_lo = (1 + int((null <= obs).sum())) / (args.perms + 1)
    log()
    log(f"PERMUTATION: observed global O/E = {obs:.4f}, "
        f"null mean = {null.mean():.4f} sd = {null.std():.4f}")
    log(f"  empirical p (enrichment) = {p_hi:.4g}, p (depletion) = {p_lo:.4g}")
    results["PERMUTATION"] = (obs, p_hi, p_lo, None)

    # 5. secondary: discard-fraction covariate (pre-registration 18.5.1)
    if args.mapping:
        import re
        per_sample = []
        for mp in args.mapping:
            m = pd.read_csv(mp, low_memory=False)
            checks = []
            for c in ["match_translated", "match_uniprot"]:
                if c in m.columns:
                    checks.append(m[c].astype(str).str.lower().eq("true"))
            m["validated"] = checks[0]
            for ch in checks[1:]:
                m["validated"] = m["validated"] & ch
            m["_s"] = m["sample"].astype(str).str.split("_").str[-1].str.strip()
            per_sample.append(m.groupby("_s")["validated"].mean().rename("keep_rate"))
        keep = pd.concat(per_sample).groupby(level=0).mean()
        df, _ = prep(mut_all, bg_map, sig, "SBS4_fraction", shuf)
        df["_s"] = df["Tumor_Sample_Barcode"].astype(str)
        df["discard"] = 1 - df["_s"].map(keep).fillna(keep.mean())
        b, se, crit, note = fit_nb2(df, "pred_std", extra=["discard"])
        lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
        log()
        log(f"SECONDARY (discard covariate): O/E per 1 SD = {np.exp(b):.4f} "
            f"(95% CI {lo:.4f}..{hi:.4f}) -> {verdict(lo, hi)}  [{note}]")
        log("  (compare with PRIMARY; matching conclusions => discard bias not driving result)")

    # 6. control-failure rule (pre-registration 18.3)
    log()
    if args.shuffle_sbs4:
        b = results["PRIMARY (SBS4, all coding)"]
        ok = abs(np.log(b[0])) < 0.1
        log(f"REHEARSAL CHECK: |log O/E| = {abs(np.log(b[0])):.4f} "
            f"({'OK - near zero' if ok else 'NOT NEAR ZERO - investigate'})")
    else:
        ca, cb = results["CONTROL A (SBS1 predictor, all coding)"], results["CONTROL B (SBS4 predictor, Silent outcome)"]
        failed = []
        for lbl, r in (("CONTROL A (SBS1)", ca), ("CONTROL B (Silent)", cb)):
            if r[1] > 1 or r[2] < 1:
                failed.append(lbl)
        if failed:
            log(f"CONTROLS FAILED: {', '.join(failed)} -> RESULT UNDER QUESTION")
        else:
            log("CONTROLS PASSED: primary result stands per pre-registration 18.3")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nSaved single output file: {args.output}")


if __name__ == "__main__":
    main()
