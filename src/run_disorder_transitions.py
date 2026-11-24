"""Transition-class run for the disorder analysis (pre-registration section 21).

One run -> one output file: git commit, SHA-256 of every input, per-class NB2
models with the section-21.4 rules, the combined disorder-changing model
(O-D + D-O), negative controls (SBS1 predictor; Silent outcome on the static
background), the multinomial background-resampling test, the descriptive
delta regressions, and the mapping-discard secondary (18.5.1).

--shuffle-sbs4 performs a rehearsal only; the confirmatory run is executed
exactly once (18.4).

Run:
  python src/run_disorder_transitions.py --mutations data/processed/disorder_transitions_george.csv ...
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

EQ_LO, EQ_HI = 0.833, 1.2
CLASSES = ["O-O", "O-D", "D-O", "D-D"]
CHANGING = "O-D+D-O"


def make_sig(path):
    sig = pd.read_csv(path, low_memory=False)
    sig["Tumor_Sample_Barcode"] = sig["Tumor_Sample_Barcode"].astype(str)
    if "SBS4_fraction" not in sig.columns:
        sig["SBS4_fraction"] = sig["SBS4"] / sig["total_mutations"]
    if "SBS1_fraction" not in sig.columns:
        sig["SBS1_fraction"] = sig["SBS1"] / sig["total_mutations"]
    return sig


def prep_class(mut, pmap, sig, frac_col, shuffle_rng=None):
    """Per-sample O/E frame for one outcome column, mirroring run_confirmatory.prep."""
    mut = mut.copy()
    mut["p"] = mut["context_96"].map(pmap)
    dropped_ctx = int(mut["p"].isna().sum())
    mut = mut[mut["p"].notna()]
    agg = (mut.groupby("Tumor_Sample_Barcode")
           .agg(n=("p", "size"), O=("O", "sum"), E=("p", "sum"))
           .reset_index())
    zero_e = int((agg["E"] <= 0).sum())
    agg = agg[agg["E"] > 0]
    out = agg.merge(sig, on="Tumor_Sample_Barcode", how="inner", validate="one_to_one")
    if frac_col not in out.columns and frac_col.endswith("_fraction"):
        base = frac_col[:-len("_fraction")]
        if base in out.columns and "total_mutations" in out.columns:
            out[frac_col] = out[base] / out["total_mutations"]
    frac = out[frac_col].astype(float)
    if shuffle_rng is not None:
        frac = pd.Series(shuffle_rng.permutation(frac.values), index=frac.index)
    out["pred_std"] = (frac - frac.mean()) / frac.std()
    out["TMB_std"] = (out["total_mutations"] - out["total_mutations"].mean()) / \
        out["total_mutations"].std()
    out["wes"] = out["dataset"].astype(str).str.lower().isin(["jhu", "clcgp"]).astype(float) \
        if "dataset" in out.columns else 0.0
    return out.dropna(subset=["O", "E", "pred_std", "TMB_std"]), dropped_ctx, zero_e


def fit_and_log(df, label, log, results):
    b, se, crit, note = fit_nb2(df, "pred_std")
    lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
    v = verdict(lo, hi)
    log(f"{label}: n_samples={len(df)}, n_mutations={int(df.n.sum())}")
    log(f"  [{note}]")
    log(f"  O/E ratio per 1 SD = {np.exp(b):.4f} (95% CI {lo:.4f}..{hi:.4f}) -> {v}")
    results[label] = (np.exp(b), lo, hi, v)
    return b, se, crit


def class_frames(nonsyn, trans_bg, sig, frac_col, shuffle=None):
    """Frames per class and the combined disorder-changing class."""
    pmap = {c: dict(zip(trans_bg["context_96"], trans_bg[f"p_{c.replace('-', '').lower()}"]))
            for c in CLASSES}
    frames = {}
    for c in CLASSES + [CHANGING]:
        m = nonsyn.copy()
        if c == CHANGING:
            m["O"] = m["transition"].isin(["O-D", "D-O"]).astype(int)
            pm = {ctx: pmap["O-D"].get(ctx, np.nan) + pmap["D-O"].get(ctx, np.nan)
                  if pd.notna(pmap["O-D"].get(ctx)) and pd.notna(pmap["D-O"].get(ctx))
                  else np.nan for ctx in m["context_96"].unique()}
        else:
            m["O"] = (m["transition"] == c).astype(int)
            pm = pmap[c]
        frames[c] = prep_class(m, pm, sig, frac_col, shuffle)
    return frames


def resampling_test(nonsyn, trans_bg, rng, perms, log):
    p_cols = [f"p_{c.replace('-', '').lower()}" for c in CLASSES]
    pmap = trans_bg.set_index("context_96")[p_cols + ["n_sampled"]]
    nonsyn = nonsyn[nonsyn["context_96"].isin(pmap.index)]
    counts = (nonsyn.groupby(["Tumor_Sample_Barcode", "context_96"])
              .size().rename("n").reset_index())
    n_c = counts.groupby("context_96")["n"].sum().reindex(pmap.index).to_numpy()
    P = pmap[p_cols].to_numpy(dtype=float)
    E = (n_c[:, None] * P).sum(axis=0)
    obs_counts = nonsyn["transition"].value_counts()
    O = np.array([obs_counts.get(c, 0) for c in CLASSES], dtype=float)
    obs = O / E
    null = np.empty((perms, len(CLASSES)))
    for i in range(perms):
        draws = np.zeros(len(CLASSES))
        for j in range(len(P)):
            draws += rng.multinomial(int(n_c[j]), P[j])
        null[i] = draws / E
    log()
    log("BACKGROUND RESAMPLING (multinomial, parametric null):")
    labels = CLASSES + [CHANGING]
    obs_all = np.append(obs, (O[1] + O[2]) / (E[1] + E[2]))
    null_all = np.column_stack([null, (null[:, 1] + null[:, 2]) / (E[1] + E[2])])
    for k, lbl in enumerate(labels):
        p_hi = (1 + int((null_all[:, k] >= obs_all[k]).sum())) / (perms + 1)
        p_lo = (1 + int((null_all[:, k] <= obs_all[k]).sum())) / (perms + 1)
        log(f"  {lbl}: observed global O/E = {obs_all[k]:.4f}, "
            f"null mean = {null_all[:, k].mean():.4f} sd = {null_all[:, k].std():.4f}; "
            f"p(enrich) = {p_hi:.4g}, p(deplete) = {p_lo:.4g}")
    log(f"  (minimum at {perms} runs: {1 / (perms + 1):.1e}; read as bound, not estimate)")


def delta_models(nonsyn, sig, log, extra_note=""):
    agg = (nonsyn.groupby("Tumor_Sample_Barcode")
           .agg(mean_delta_pos=("delta_pos", "mean"),
                mean_delta_window=("delta_window10", "mean")).reset_index())
    df = agg.merge(sig, on="Tumor_Sample_Barcode", how="inner", validate="one_to_one")
    frac = df["SBS4_fraction"].astype(float)
    df["pred_std"] = (frac - frac.mean()) / frac.std()
    df["TMB_std"] = (df["total_mutations"] - df["total_mutations"].mean()) / \
        df["total_mutations"].std()
    df["wes"] = df["dataset"].astype(str).str.lower().isin(["jhu", "clcgp"]).astype(float) \
        if "dataset" in df.columns else 0.0
    log()
    log(f"DELTA REGRESSIONS (descriptive, HC1 t-CI){extra_note}:")
    for col in ["mean_delta_pos", "mean_delta_window"]:
        d = df.dropna(subset=[col, "pred_std", "TMB_std"])
        cols = [c for c in ["pred_std", "TMB_std", "wes"] if d[c].nunique() > 1]
        X = sm.add_constant(d[cols])
        m = sm.OLS(d[col], X).fit(cov_type="HC1")
        b, se = m.params["pred_std"], m.bse["pred_std"]
        crit = tdist.ppf(0.975, len(d) - len(cols) - 1)
        log(f"  {col} ~ SBS4_std: b = {b:.5f} (95% CI {b - crit * se:.5f}..{b + crit * se:.5f}), "
            f"n = {len(d)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", nargs="+", required=True,
                        help="disorder_transitions_*.csv (annotate_disorder_transitions output)")
    parser.add_argument("--transition-background", required=True)
    parser.add_argument("--darkness-background", required=True,
                        help="static build_background output on the disorder track (Silent control)")
    parser.add_argument("--signatures", required=True)
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
    mode = "REHEARSAL (SBS4 shuffled across samples)" if args.shuffle_sbs4 else "CONFIRMATORY"
    log(f"=== {mode} RUN (disorder transitions) ===")
    log(f"git commit: {commit}")
    log(f"seed: {args.seed}, permutations: {args.perms}")
    log()
    log("-- input checksums (SHA-256) --")
    inputs = list(args.mutations) + [args.transition_background, args.darkness_background,
                                     args.signatures] + (args.mapping or [])
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

    # 21.2 invariants
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

    # 1. class models + combined, SBS4 and SBS1 predictors
    frames = class_frames(nonsyn, trans_bg, sig, "SBS4_fraction", shuf)
    for c in CLASSES + [CHANGING]:
        df, dropped, zero_e = frames[c]
        log(f"[class {c}] dropped unmapped context rows: {dropped}, zero-E samples dropped: {zero_e}")
        fit_and_log(df, f"{c} (SBS4)", log, results)
    log()
    for c in CLASSES + [CHANGING]:
        df, _, _ = class_frames(nonsyn, trans_bg, sig, "SBS1_fraction", shuf)[c]
        fit_and_log(df, f"{c} (CONTROL A: SBS1)", log, results)

    # 2. CONTROL B: silent outcome on the static disorder background
    log()
    s = silent.copy()
    s["O"] = (s["state_wt"] == "D").astype(int)
    df_b, dropped_b, zero_b = prep_class(s, static_map, sig, "SBS4_fraction", shuf)
    log(f"[CONTROL B] dropped unmapped context rows: {dropped_b}, zero-E samples: {zero_b}")
    fit_and_log(df_b, "CONTROL B (Silent outcome, static background)", log, results)

    # 3. background resampling
    resampling_test(nonsyn, trans_bg, rng, args.perms, log)

    # 4. descriptive deltas
    delta_models(nonsyn, sig, log)

    # 5. discard-fraction secondary (18.5.1)
    if args.mapping:
        per_sample = []
        for mp in args.mapping:
            m = pd.read_csv(mp, low_memory=False)
            checks = [m[c].astype(str).str.lower().eq("true")
                      for c in ["match_translated", "match_uniprot"] if c in m.columns]
            m["validated"] = checks[0]
            for ch in checks[1:]:
                m["validated"] = m["validated"] & ch
            m["_s"] = m["sample"].astype(str).str.split("_").str[-1].str.strip()
            per_sample.append(m.groupby("_s")["validated"].mean().rename("keep_rate"))
        keep = pd.concat(per_sample).groupby(level=0).mean()
        df_d, _, _ = class_frames(nonsyn, trans_bg, sig, "SBS4_fraction", shuf)[CHANGING]
        df_d = df_d.copy()
        df_d["discard"] = 1 - df_d["Tumor_Sample_Barcode"].astype(str).map(keep).fillna(keep.mean())
        b, se, crit, note = fit_nb2(df_d, "pred_std", extra=["discard"])
        lo, hi = np.exp(b - crit * se), np.exp(b + crit * se)
        log()
        log(f"SECONDARY (discard covariate, {CHANGING}): O/E per 1 SD = {np.exp(b):.4f} "
            f"(95% CI {lo:.4f}..{hi:.4f}) -> {verdict(lo, hi)}  [{note}]")

    # 6. control rule (18.3)
    log()
    if args.shuffle_sbs4:
        r0 = results[f"{CHANGING} (SBS4)"]
        ok = abs(np.log(r0[0])) < 0.1
        log(f"REHEARSAL CHECK: |log O/E| = {abs(np.log(r0[0])):.4f} "
            f"({'OK - near zero' if ok else 'NOT NEAR ZERO - investigate'})")
    else:
        failed = []
        for lbl in [f"{c} (CONTROL A: SBS1)" for c in CLASSES + [CHANGING]] + \
                   ["CONTROL B (Silent outcome, static background)"]:
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
