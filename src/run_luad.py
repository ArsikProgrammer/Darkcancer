"""LUAD dose-response run (pre-registration section 19). One run -> one file.

Same functions as the SCLC confirmatory run (run_confirmatory: sha256, verdict,
fit_nb2 with the pre-registered alpha->0 Poisson-t boundary rule). Adds the
pre-registered LUAD secondaries: pack-years, smoking-status contrast, and the
hypermutated-sample exclusion. --shuffle-predictor performs a rehearsal only.
"""
from pathlib import Path
import argparse
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_confirmatory import sha256, verdict, fit_nb2

EQ_LO, EQ_HI = 0.833, 1.2


def prep(mut, bg_map, sig, frac_col, rng=None):
    mut = mut[mut["is_dark"].notna()].copy()
    mut["p_dark_bg"] = mut["context_96"].map(bg_map)
    dropped = int(mut["p_dark_bg"].isna().sum())
    mut = mut[mut["p_dark_bg"].notna()]
    agg = (mut.groupby("Tumor_Sample_Barcode")
           .agg(n=("is_dark", "size"), O=("is_dark", "sum"), E=("p_dark_bg", "sum"))
           .reset_index())
    out = agg.merge(sig, on="Tumor_Sample_Barcode", how="inner", validate="one_to_one")
    frac = out[frac_col].astype(float)
    if rng is not None:
        frac = pd.Series(rng.permutation(frac.values), index=frac.index)
    out["pred_std"] = (frac - frac.mean()) / frac.std()
    out["TMB_std"] = (out["total_mutations"] - out["total_mutations"].mean()) / out["total_mutations"].std()
    out["wes"] = 0.0
    return out.dropna(subset=["O", "E", "pred_std", "TMB_std"]), dropped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", required=True)
    parser.add_argument("--background", required=True)
    parser.add_argument("--signatures", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--smoking", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--perms", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261028)
    parser.add_argument("--shuffle-predictor", action="store_true")
    args = parser.parse_args()

    lines = []
    def log(s=""):
        print(s)
        lines.append(s)

    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                            text=True, cwd=Path(__file__).parent).stdout.strip()
    mode = "REHEARSAL (predictor shuffled)" if args.shuffle_predictor else "CONFIRMATORY (LUAD dose-response)"
    log(f"=== {mode} RUN ===")
    log(f"git commit: {commit}")
    log(f"seed: {args.seed}, permutations: {args.perms}")
    log()
    log("-- input checksums (SHA-256) --")
    for p in [args.mutations, args.background, args.signatures, args.mapping, args.smoking]:
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

    mut_all = pd.read_csv(args.mutations, low_memory=False)
    mut_all["Tumor_Sample_Barcode"] = mut_all["Tumor_Sample_Barcode"].astype(str)

    # pre-registered exclusion: hypermutated samples (> 10x median coding SNPs)
    cnt = mut_all.groupby("Tumor_Sample_Barcode").size()
    keep = cnt[cnt <= 10 * cnt.median()].index
    dropped_hm = sorted(set(cnt.index) - set(keep))
    log(f"Exclusion (coding SNPs > 10x median {int(cnt.median())}): "
        f"{len(dropped_hm)} samples {dropped_hm}")
    mut_all = mut_all[mut_all["Tumor_Sample_Barcode"].isin(keep)]
    log()

    rng = np.random.default_rng(args.seed)
    shuf = rng if args.shuffle_predictor else None
    results = {}

    def run_model(frame, frac_col, label):
        df, dropped = prep(frame, bg_map, sig, frac_col, shuf)
        b, se, crit, note = fit_nb2(df, "pred_std")
        lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
        v = verdict(lo, hi)
        log(f"{label}: n_samples={len(df)}, n_mutations={int(df.n.sum())} "
            f"(dropped unmapped context: {dropped})")
        log(f"  [{note}]")
        log(f"  O/E ratio per 1 SD = {np.exp(b):.4f} (95% CI {lo:.4f}..{hi:.4f}) -> {v}")
        results[label] = (np.exp(b), lo, hi, v)

    run_model(mut_all, "SBS4_fraction", "PRIMARY (SBS4 dose, all coding)")
    run_model(mut_all, "SBS1_fraction", "CONTROL A (SBS1 predictor)")
    silent = mut_all[mut_all["Variant_Classification"] == "Silent"]
    run_model(silent, "SBS4_fraction", "CONTROL B (SBS4, Silent outcome)")

    # permutation
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
    log(f"PERMUTATION: observed global O/E = {obs:.4f}, null mean = {null.mean():.4f} "
        f"sd = {null.std():.4f}; p(enrich) = {p_hi:.4g}, p(deplete) = {p_lo:.4g}")

    # secondary 1: pack-years as independent smoking measure
    sm = pd.read_csv(args.smoking, sep="\t", comment="#", dtype=str)
    sm["patient"] = sm["PATIENT_ID"].astype(str)
    sm["pack_years"] = pd.to_numeric(sm["SMOKING_PACK_YEARS"].replace("[Not Available]", None),
                                     errors="coerce")
    sm["status"] = pd.to_numeric(sm["TOBACCO_SMOKING_HISTORY_INDICATOR"].replace("[Not Available]", None),
                                 errors="coerce")
    sig["patient"] = sig["Tumor_Sample_Barcode"].str[:12]
    sig_py = sig.merge(sm[["patient", "pack_years", "status"]], on="patient", how="inner")
    sig_py = sig_py.drop_duplicates("Tumor_Sample_Barcode")
    df_py, _ = prep(mut_all, bg_map, sig_py, "pack_years", shuf)
    df_py = df_py.dropna(subset=["pred_std"])
    if len(df_py) >= 30:
        b, se, crit, note = fit_nb2(df_py, "pred_std")
        lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
        log()
        log(f"SECONDARY 1 (pack-years, n={len(df_py)}): O/E per 1 SD = {np.exp(b):.4f} "
            f"(95% CI {lo:.4f}..{hi:.4f}) -> {verdict(lo, hi)}  [{note}]")

    # secondary 2: never (1) vs current (2) smokers, SBS4 predictor
    sig_nc = sig_py[sig_py["status"].isin([1, 2])]
    df_nc, _ = prep(mut_all[mut_all["Tumor_Sample_Barcode"].isin(sig_nc["Tumor_Sample_Barcode"])],
                    bg_map, sig_nc, "SBS4_fraction", shuf)
    b, se, crit, note = fit_nb2(df_nc, "pred_std")
    lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
    log()
    log(f"SECONDARY 2 (never vs current, n={len(df_nc)}): O/E per 1 SD = {np.exp(b):.4f} "
        f"(95% CI {lo:.4f}..{hi:.4f}) -> {verdict(lo, hi)}  [{note}]")

    # secondary 3: discard covariate
    m = pd.read_csv(args.mapping, low_memory=False)
    checks = [m[c].astype(str).str.lower().eq("true")
              for c in ["match_translated", "match_uniprot"] if c in m.columns]
    m["validated"] = checks[0]
    for ch in checks[1:]:
        m["validated"] = m["validated"] & ch
    m["_s"] = m["sample"].astype(str)
    keep_rate = m.groupby("_s")["validated"].mean()
    df_d, _ = prep(mut_all, bg_map, sig, "SBS4_fraction", shuf)
    df_d["discard"] = 1 - df_d["Tumor_Sample_Barcode"].map(keep_rate).fillna(keep_rate.mean())
    b, se, crit, note = fit_nb2(df_d, "pred_std", extra=["discard"])
    lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
    log()
    log(f"SECONDARY 3 (discard covariate): O/E per 1 SD = {np.exp(b):.4f} "
        f"(95% CI {lo:.4f}..{hi:.4f}) -> {verdict(lo, hi)}  [{note}]")

    # control rule
    log()
    if args.shuffle_predictor:
        b0 = results["PRIMARY (SBS4 dose, all coding)"]
        ok = abs(np.log(b0[0])) < 0.1
        log(f"REHEARSAL CHECK: |log O/E| = {abs(np.log(b0[0])):.4f} "
            f"({'OK - near zero' if ok else 'NOT NEAR ZERO - investigate'})")
    else:
        failed = [lbl for lbl in ("CONTROL A (SBS1 predictor)", "CONTROL B (SBS4, Silent outcome)")
                  if results[lbl][1] > 1 or results[lbl][2] < 1]
        if failed:
            log(f"CONTROLS FAILED: {', '.join(failed)} -> RESULT UNDER QUESTION")
        else:
            log("CONTROLS PASSED: primary result stands per pre-registration 18.3/19.4")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nSaved single output file: {args.output}")


if __name__ == "__main__":
    main()
