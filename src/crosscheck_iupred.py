"""Cross-predictor robustness check with IUPred2A (plan section 21.5).

Pre-fixed rule (21.5): this check is descriptive and cannot alter the primary
results. Stratified subset of mutations (strata: transition class x cohort x
|wt_score - 0.5| band), reported as transition-state concordance and score
correlations.

Two adapters:
  --iupred-dir DIR   IUPred2A standalone (full check: WT and mutant sequences)
  --mode rest        IUPred2A web API by UniProt accession (partial check:
                     WT scores only; mutant sequences are not accepted by the API)

Run:
  python src/crosscheck_iupred.py --mutations data/processed/disorder_transitions_*.csv \
      --output results/iupred_crosscheck_20261124.txt --mode rest --n 300
"""
from pathlib import Path
import argparse
import re
import subprocess
import time
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from score_disorder_effects import load_uniprot_seqs, load_wt_scores, window_delta

REST_URL = "https://iupred2a.elte.hu/iupred2a/long/{acc}"
BANDS = [0.1, 0.3]  # |wt - 0.5| band edges


def band(v):
    d = abs(v - 0.5)
    return 0 if d < BANDS[0] else (1 if d < BANDS[1] else 2)


def stratified_sample(mut, n, rng):
    mut = mut.copy()
    mut["band"] = mut["wt_score"].map(band)
    mut["stratum"] = (mut["transition"].astype(str) + "|" + mut["file"] + "|" +
                      mut["band"].astype(str))
    out = []
    for _, grp in mut.groupby("stratum"):
        k = max(1, int(round(n * len(grp) / len(mut))))
        out.append(grp.sample(min(k, len(grp)), random_state=rng))
    return pd.concat(out, ignore_index=True)


class IUPredREST:
    def __init__(self):
        self.cache = {}

    def wt_scores(self, acc):
        if acc in self.cache:
            return self.cache[acc]
        try:
            with urllib.request.urlopen(REST_URL.format(acc=acc), timeout=60) as r:
                text = r.read().decode("utf-8", "ignore")
        except Exception as e:
            print(f"  REST failed for {acc}: {e}")
            self.cache[acc] = None
            return None
        scores = {}
        for line in text.splitlines():
            m = re.match(r"^(\d+)\s+\S+\s+([0-9.]+)$", line.strip())
            if m:
                scores[int(m.group(1))] = float(m.group(2))
        self.cache[acc] = scores if scores else None
        time.sleep(0.5)
        return self.cache[acc]


class IUPredStandalone:
    def __init__(self, iupred_dir):
        self.dir = Path(iupred_dir)

    def score_sequence(self, seq):
        fasta = Path(f"_tmp_{abs(hash(seq))}.fasta")
        fasta.write_text(f">tmp\n{seq}\n")
        try:
            r = subprocess.run(["python", str(self.dir / "iupred2a.py"), str(fasta), "long"],
                               capture_output=True, text=True, cwd=self.dir, timeout=300)
            scores = {}
            for line in r.stdout.splitlines():
                m = re.match(r"^(\d+)\s+\S+\s+([0-9.]+)$", line.strip())
                if m:
                    scores[int(m.group(1))] = float(m.group(2))
            return scores
        finally:
            fasta.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--mode", choices=["auto", "rest", "standalone"], default="auto")
    parser.add_argument("--iupred-dir", default=None)
    parser.add_argument("--n", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=20261124)
    args = parser.parse_args()

    frames = []
    for p in args.mutations:
        d = pd.read_csv(p, low_memory=False)
        d["file"] = Path(p).name
        frames.append(d)
    mut = pd.concat(frames, ignore_index=True)
    mut = mut[(mut["Variant_Classification"] == "Missense_Mutation") & mut["wt_score"].notna()]
    print("Missense rows available:", len(mut))

    rng = np.random.default_rng(args.seed)
    n = args.n
    mode = args.mode
    if mode == "auto":
        mode = "standalone" if args.iupred_dir else "rest"
    if mode == "rest":
        n = min(n, 300)

    sub = stratified_sample(mut, n, rng)
    print(f"Mode: {mode}; subset: {len(sub)} mutations")

    lines = []
    def log(s=""):
        print(s)
        lines.append(s)

    log("=== CROSS-PREDICTOR CHECK (IUPred2A, plan 21.5; descriptive, cannot alter primary) ===")
    log(f"mode: {mode}, seed: {args.seed}, subset n: {len(sub)}")
    log()

    if mode == "rest":
        client = IUPredREST()
        accs = sorted(sub["uniprot_accession"].unique())[:120]
        log(f"REST: fetching WT scores for {len(accs)} accessions (capped at 120)")
        wt_cache = load_wt_scores()
        rows = []
        for acc in accs:
            scores = client.wt_scores(acc)
            if not scores:
                continue
            meta = sub[sub["uniprot_accession"] == acc]
            mp_wt = wt_cache.get(acc)
            for _, r in meta.iterrows():
                pos = int(r["uniprot_pos"])
                if pos not in scores or mp_wt is None:
                    continue
                rows.append({
                    "uniprot_accession": acc, "uniprot_pos": pos,
                    "iupred_wt": scores[pos], "mp_wt": float(mp_wt[pos - 1]),
                    "state_iupred": "D" if scores[pos] >= 0.5 else "O",
                    "state_mp": r["state_wt"], "transition_mp": r["transition"],
                })
        if not rows:
            log("REST: no comparable rows - SKIPPED (see mode=standalone for the full check)")
        else:
            res = pd.DataFrame(rows)
            agree = (res["state_iupred"] == res["state_mp"]).mean()
            rho, p = spearmanr(res["iupred_wt"], res["mp_wt"])
            log(f"WT state concordance at mutated positions: {agree:.4f} (n = {len(res)})")
            log(f"WT score Spearman rho = {rho:.4f} (p = {p:.2g})")
            log(f"mean shift (IUPred - metapredict) = {(res['iupred_wt'] - res['mp_wt']).mean():+.4f}")
            by_cls = res.groupby("transition_mp").apply(
                lambda g: (g["state_iupred"] == g["state_mp"]).mean(), include_groups=False)
            log(f"concordance by metapredict class: {by_cls.round(4).to_dict()}")
            log()
            log("NOTE: mutant-side cross-check (delta concordance) requires IUPred2A standalone; "
                "the web API scores UniProt accessions only.")
    else:
        client = IUPredStandalone(args.iupred_dir)
        seqs = load_uniprot_seqs()
        wt_cache = load_wt_scores()
        rows = []
        for _, r in sub.iterrows():
            acc, pos, alt = r["uniprot_accession"], int(r["uniprot_pos"]), r["alt_aa"]
            seq = seqs.get(acc)
            if seq is None:
                continue
            mut_seq = seq[:pos - 1] + alt + seq[pos:]
            sw = client.score_sequence(seq)
            sm = client.score_sequence(mut_seq)
            if not sw or not sm or pos not in sw:
                continue
            dw, wlen = window_delta([sw[i] for i in sorted(sw)], [sm[i] for i in sorted(sm)],
                                    pos - 1)
            rows.append({
                "uniprot_accession": acc, "uniprot_pos": pos,
                "iupred_wt": sw[pos], "iupred_mut": sm.get(pos, np.nan),
                "iupred_delta": sm.get(pos, np.nan) - sw[pos], "iupred_delta_window": dw,
                "mp_delta": r["delta_pos"], "mp_delta_window": r["delta_window10"],
                "transition_mp": r["transition"],
                "transition_iupred": ("D" if sw[pos] >= 0.5 else "O") + "-" +
                                     ("D" if sm.get(pos, 0) >= 0.5 else "O"),
            })
        res = pd.DataFrame(rows)
        if len(res):
            agree = (res["transition_mp"] == res["transition_iupred"]).mean()
            rho1, _ = spearmanr(res["iupred_delta"], res["mp_delta"])
            rho2, _ = spearmanr(res["iupred_delta_window"], res["mp_delta_window"])
            log(f"Transition-class concordance: {agree:.4f} (n = {len(res)})")
            log(f"Spearman rho (delta_pos) = {rho1:.4f}; rho (delta_window10) = {rho2:.4f}")
        else:
            log("Standalone produced no rows - SKIPPED")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()
