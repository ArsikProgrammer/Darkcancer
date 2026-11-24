"""Score disorder effects of single amino-acid substitutions (plan section 21).

For each substitution (uniprot_accession, uniprot_pos, ref_aa, alt_aa) the WT
UniProt sequence and the single-substitution mutant sequence are scored with
metapredict; per-residue deltas are computed at the mutated position and over
the clipped window [pos-10, pos+10]. Scores are cached per substitution so the
observed-mutation set and the transition background share one scoring pass.

Also exports the pure helpers used by the transition pipeline:
parse_hgvsp, transition_class, window_delta.

Run:
  python src/score_disorder_effects.py --substitutions data/processed/_subs_observed.csv \
      --output data/processed/disorder_substitution_scores.csv
"""
from pathlib import Path
import argparse
import pickle
import re

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "data" / "processed" / "_uniprot_human_reviewed.fasta"
WT_CACHE = ROOT / "data" / "processed" / "_disorder_wt_scores.pkl"
MUT_CACHE = ROOT / "data" / "processed" / "_disorder_mutant_scores.pkl"

WINDOW_HALF = 10
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")


def is_scoreable(seq):
    """metapredict one-hot encodes the 20 standard residues only."""
    return set(seq) <= AA_ALPHABET

_HGVSP_RE = re.compile(r"^p?\.?([A-Z*])(\d+)([A-Z*=])?$")


def parse_hgvsp(s):
    """Parse a short protein change: p.R175H / R175H / R175= / R175.

    Returns (ref_aa, pos, alt_aa) with alt_aa == ref_aa for synonymous
    notation ("=" or omitted), or None if the string is not a simple
    amino-acid substitution notation (nonsense "*", indels, "?", ...).
    """
    m = _HGVSP_RE.match(str(s).strip())
    if not m:
        return None
    ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
    if alt is None or alt == "=":
        alt = ref
    return ref, pos, alt


def transition_class(wt_score, mut_score, threshold):
    """O-O / O-D / D-O / D-D from two disorder scores (>= threshold = D)."""
    s_wt = "D" if wt_score >= threshold else "O"
    s_mut = "D" if mut_score >= threshold else "O"
    return f"{s_wt}-{s_mut}"


def window_delta(wt_scores, mut_scores, pos0, half=WINDOW_HALF):
    """Mean (mut - wt) over the clipped window around 0-based pos0.

    Returns (delta_window, window_len).
    """
    lo = max(0, pos0 - half)
    hi = min(len(wt_scores), pos0 + half + 1)
    d = np.asarray(mut_scores[lo:hi]) - np.asarray(wt_scores[lo:hi])
    return float(d.mean()), hi - lo


def load_wt_scores():
    with open(WT_CACHE, "rb") as f:
        return pickle.load(f)


def batch_scores(pred):
    """Normalize one predict_disorder_batch element to a score array.

    metapredict returns [sequence, scores] pairs per input.
    """
    if isinstance(pred, (list, tuple)) and len(pred) == 2 and isinstance(pred[0], str):
        return pred[1]
    return pred


def load_mut_scores():
    if MUT_CACHE.exists():
        with open(MUT_CACHE, "rb") as f:
            return pickle.load(f)
    return {}


def save_mut_scores(cache):
    MUT_CACHE.parent.mkdir(parents=True, exist_ok=True)
    with open(MUT_CACHE, "wb") as f:
        pickle.dump(cache, f)


def load_uniprot_seqs():
    from map_mutations_to_residues import load_uniprot_fasta
    return load_uniprot_fasta()


def score_substitutions(subs, wt_cache=None, batch_size=2000, window_half=WINDOW_HALF):
    """Score substitutions; subs is a DataFrame with columns
    uniprot_accession, uniprot_pos (1-based), ref_aa, alt_aa.

    Returns a DataFrame with wt_score, mut_score, delta_pos, delta_window10,
    window_len. Results are cached in MUT_CACHE across calls.
    """
    import metapredict as mp

    if wt_cache is None:
        wt_cache = load_wt_scores()
    seqs = load_uniprot_seqs()
    mut_cache = load_mut_scores()

    keys = list(subs[["uniprot_accession", "uniprot_pos", "alt_aa"]]
                .itertuples(index=False, name=None))
    need = []
    seen = set()
    skipped = 0
    for acc, pos, alt in keys:
        k = (str(acc), int(pos), str(alt))
        if k not in mut_cache and k not in seen:
            if k[0] in seqs and 1 <= k[1] <= len(seqs[k[0]]):
                if not is_scoreable(seqs[k[0]]):
                    skipped += 1
                    continue
                seen.add(k)
                need.append(k)
    if skipped:
        print(f"Skipped {skipped} substitutions in proteins with non-standard residues")

    by_acc = {}
    for acc, pos, alt in need:
        by_acc.setdefault(acc, []).append((pos, alt))

    n_done = 0
    for acc, items in sorted(by_acc.items()):
        seq = seqs[acc]
        wt = wt_cache.get(acc)
        mutants = []
        for pos, alt in items:
            mut_seq = seq[:pos - 1] + alt + seq[pos:]
            mutants.append(mut_seq)
        for i in range(0, len(mutants), batch_size):
            chunk = mutants[i:i + batch_size]
            preds = mp.predict_disorder_batch(chunk, show_progress_bar=False)
            for (pos, alt), pred in zip(items[i:i + batch_size], preds):
                pred = np.asarray(batch_scores(pred), dtype=np.float64)
                p0 = pos - 1
                mut_score = float(pred[p0])
                if wt is not None:
                    dw, wlen = window_delta(wt, pred, p0, window_half)
                else:
                    dw, wlen = float("nan"), 0
                mut_cache[(acc, pos, alt)] = (mut_score, dw, wlen)
                n_done += 1
        print(f"  {acc}: scored {len(items)} mutants", flush=True)

    if need:
        save_mut_scores(mut_cache)
    print(f"Substitutions scored this run: {n_done}; cache size: {len(mut_cache)}")

    rows = []
    for acc, pos, alt in keys:
        k = (str(acc), int(pos), str(alt))
        if k not in mut_cache:
            rows.append((np.nan, np.nan, np.nan, np.nan, 0))
            continue
        mut_score, dw, wlen = mut_cache[k]
        wt = wt_cache.get(str(acc))
        wt_score = float(wt[int(pos) - 1]) if wt is not None else np.nan
        rows.append((wt_score, mut_score, mut_score - wt_score if wt is not None else np.nan,
                     dw, wlen))
    out = subs.copy()
    out[["wt_score", "mut_score", "delta_pos", "delta_window10", "window_len"]] = pd.DataFrame(
        rows, index=out.index)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--substitutions", required=True,
                        help="CSV: uniprot_accession,uniprot_pos,ref_aa,alt_aa")
    parser.add_argument("--output", required=True)
    parser.add_argument("--batch-size", type=int, default=2000)
    args = parser.parse_args()

    subs = pd.read_csv(args.substitutions, low_memory=False)
    need = {"uniprot_accession", "uniprot_pos", "ref_aa", "alt_aa"}
    if not need.issubset(subs.columns):
        raise ValueError(f"substitution file must contain {need}")
    subs = subs.drop_duplicates(["uniprot_accession", "uniprot_pos", "ref_aa", "alt_aa"])
    print("Unique substitutions:", len(subs))

    out = score_substitutions(subs, batch_size=args.batch_size)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
