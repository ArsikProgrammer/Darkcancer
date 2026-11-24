"""LUAD variant of the disorder-transition run (plan sections 21 + 19.5 secondaries).

Same code path as run_disorder_transitions with the pre-registered LUAD
additions: hypermutated-sample exclusion, pack-years and never-vs-current
secondaries on the disorder-changing model and the delta regressions.

Run:
  python src/run_disorder_transitions_luad.py --mutations data/processed/disorder_transitions_luad.csv ...
"""
from pathlib import Path
import argparse
import subprocess
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import t as tdist

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_confirmatory import sha256, fit_nb2, verdict  # noqa: E402
from run_disorder_transitions import (CHANGING, class_frames, delta_models, fit_and_log,
                                      make_sig, prep_class, resampling_test)  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", nargs="+", required=True)
    parser.add_argument("--transition-background", required=True)
    parser.add_argument("--darkness-background", required=True)
    parser.add_argument("--signatures", required=True)
    parser.add_argument("--smoking", required=True)
    parser.add_argument("--mapping", nargs="*", default=None)
    parser.add_argument("--output", required=True)
    parser.add_argument("--perms", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261124)
    parser.add_argument("--shuffle-sbs4", action="store_true")
    args = parser.parse_args()

    lines = []
    def log(s=""):
        print(s)
        lines.append(s)

    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                            text=True, cwd=Path(__file__).parent).stdout.strip()
    mode = "REHEARSAL (SBS4 shuffled)" if args.shuffle_sbs4 else "CONFIRMATORY (LUAD transitions)"
    log(f"=== {mode} RUN ===")
    log(f"git commit: {commit}")
    log(f"seed: {args.seed}, permutations: {args.perms}")
    log()
    log("-- input checksums (SHA-256) --")
    inputs = list(args.mutations) + [args.transition_background, args.darkness_background,
                                     args.signatures, args.smoking] + (args.mapping or [])
    for p in inputs:
        log(f"{sha256(p)}  {p}")
    log()

    trans_bg = pd.read_csv(args.transition_background)
    static_bg = pd.read_csv(args.darkness_background)
    static_map = dict(zip(static_bg["context_96"], static_bg["dark_fraction"]))
    sig = make_sig(args.signatures)
    rng = np.random.default_rng(args.seed)

    mut_all = pd.concat([pd.read_csv(p, low_memory=False) for p in args.mutations],
                        ignore_index=True)
    mut_all["Tumor_Sample_Barcode"] = mut_all["Tumor_Sample_Barcode"].astype(str)

    # pre-registered hypermutator exclusion (19.3): coding SNPs > 10x cohort median
    cnt = mut_all.groupby("Tumor_Sample_Barcode").size()
    keep = cnt[cnt <= 10 * cnt.median()].index
    dropped_hm = sorted(set(cnt.index) - set(keep))
    log(f"Exclusion (coding SNPs > 10x median {int(cnt.median())}): "
        f"{len(dropped_hm)} samples {dropped_hm}")
    mut_all = mut_all[mut_all["Tumor_Sample_Barcode"].isin(keep)]

    silent = mut_all[mut_all["Variant_Classification"] == "Silent"]
    nonsyn = mut_all[mut_all["Variant_Classification"] != "Silent"]
    ok_delta = bool((silent["delta_pos"] == 0).all() and (silent["delta_window10"] == 0).all())
    ok_cls = bool(silent["transition"].isin(["O-O", "D-D"]).all())
    log(f"SILENT INVARIANTS (21.2): deltas zero = {ok_delta}, classes in O-O/D-D = {ok_cls} "
        f"(silent rows: {len(silent)})")
    if not (ok_delta and ok_cls):
        raise SystemExit("silent invariants violated - pipeline bug, run invalid")
    log()
    log(f"Non-synonymous rows: {len(nonsyn)}; transition counts: "
        f"{nonsyn['transition'].value_counts().to_dict()}")
    log()

    shuf = rng if args.shuffle_sbs4 else None
    results = {}

    frames = class_frames(nonsyn, trans_bg, sig, "SBS4_fraction", shuf)
    for c in frames:
        df, dropped, zero_e = frames[c]
        log(f"[class {c}] dropped unmapped context rows: {dropped}, zero-E samples dropped: {zero_e}")
        fit_and_log(df, f"{c} (SBS4)", log, results)
    log()
    for c in frames:
        df, _, _ = class_frames(nonsyn, trans_bg, sig, "SBS1_fraction", shuf)[c]
        fit_and_log(df, f"{c} (CONTROL A: SBS1)", log, results)

    log()
    s = silent.copy()
    s["O"] = (s["state_wt"] == "D").astype(int)
    df_b, dropped_b, zero_b = prep_class(s, static_map, sig, "SBS4_fraction", shuf)
    log(f"[CONTROL B] dropped unmapped context rows: {dropped_b}, zero-E samples: {zero_b}")
    fit_and_log(df_b, "CONTROL B (Silent outcome, static background)", log, results)

    resampling_test(nonsyn, trans_bg, rng, args.perms, log)
    delta_models(nonsyn, sig, log)

    # 19.5 secondaries on the disorder-changing model: pack-years, never-vs-current
    sm_df = pd.read_csv(args.smoking, sep="\t", comment="#", dtype=str)
    sm_df["patient"] = sm_df["PATIENT_ID"].astype(str)
    sm_df["pack_years"] = pd.to_numeric(sm_df["SMOKING_PACK_YEARS"].replace("[Not Available]", None),
                                        errors="coerce")
    sm_df["status"] = pd.to_numeric(sm_df["TOBACCO_SMOKING_HISTORY_INDICATOR"].replace("[Not Available]", None),
                                    errors="coerce")
    sig["patient"] = sig["Tumor_Sample_Barcode"].str[:12]
    sig_py = sig.merge(sm_df[["patient", "pack_years", "status"]], on="patient", how="inner")
    sig_py = sig_py.drop_duplicates("Tumor_Sample_Barcode")

    df_py, _, _ = class_frames(nonsyn, trans_bg, sig_py, "pack_years", shuf)[CHANGING]
    df_py = df_py.dropna(subset=["pred_std"])
    if len(df_py) >= 30:
        b, se, crit, note = fit_nb2(df_py, "pred_std")
        lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
        log()
        log(f"SECONDARY 1 (pack-years, {CHANGING}, n={len(df_py)}): O/E per 1 SD = {np.exp(b):.4f} "
            f"(95% CI {lo:.4f}..{hi:.4f}) -> {verdict(lo, hi)}  [{note}]")

    sig_nc = sig_py[sig_py["status"].isin([1, 2])]
    df_nc, _, _ = class_frames(
        nonsyn[nonsyn["Tumor_Sample_Barcode"].isin(sig_nc["Tumor_Sample_Barcode"])],
        trans_bg, sig_nc, "SBS4_fraction", shuf)[CHANGING]
    b, se, crit, note = fit_nb2(df_nc, "pred_std")
    lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
    log()
    log(f"SECONDARY 2 (never vs current, {CHANGING}, n={len(df_nc)}): O/E per 1 SD = {np.exp(b):.4f} "
        f"(95% CI {lo:.4f}..{hi:.4f}) -> {verdict(lo, hi)}  [{note}]")

    if args.mapping:
        per_sample = []
        for mp in args.mapping:
            m = pd.read_csv(mp, low_memory=False)
            checks = [m[c].astype(str).str.lower().eq("true")
                      for c in ["match_translated", "match_uniprot"] if c in m.columns]
            m["validated"] = checks[0]
            for ch in checks[1:]:
                m["validated"] = m["validated"] & ch
            m["_s"] = m["sample"].astype(str)
            per_sample.append(m.groupby("_s")["validated"].mean().rename("keep_rate"))
        keep_rate = pd.concat(per_sample).groupby(level=0).mean()
        df_d, _, _ = class_frames(nonsyn, trans_bg, sig, "SBS4_fraction", shuf)[CHANGING]
        df_d = df_d.copy()
        df_d["discard"] = 1 - df_d["Tumor_Sample_Barcode"].astype(str).map(keep_rate).fillna(keep_rate.mean())
        b, se, crit, note = fit_nb2(df_d, "pred_std", extra=["discard"])
        lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
        log()
        log(f"SECONDARY 3 (discard covariate, {CHANGING}): O/E per 1 SD = {np.exp(b):.4f} "
            f"(95% CI {lo:.4f}..{hi:.4f}) -> {verdict(lo, hi)}  [{note}]")

    log()
    if args.shuffle_sbs4:
        r0 = results[f"{CHANGING} (SBS4)"]
        ok = abs(np.log(r0[0])) < 0.1
        log(f"REHEARSAL CHECK: |log O/E| = {abs(np.log(r0[0])):.4f} "
            f"({'OK - near zero' if ok else 'NOT NEAR ZERO - investigate'})")
    else:
        failed = []
        for lbl in results:
            if "CONTROL A" in lbl or "CONTROL B" in lbl:
                r = results[lbl]
                if r[1] > 1 or r[2] < 1:
                    failed.append(lbl)
        if failed:
            log(f"CONTROLS FAILED: {'; '.join(failed)} -> RESULT UNDER QUESTION")
        else:
            log("CONTROLS PASSED: class results stand per pre-registration 18.3/21.4")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nSaved single output file: {args.output}")


if __name__ == "__main__":
    main()
