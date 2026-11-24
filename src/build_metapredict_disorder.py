"""Build a disorder track from metapredict scores on reviewed human UniProt sequences.

Plan section 20: is_disordered = score >= threshold (primary 0.5; sensitivity
0.4 / 0.6). Scores are computed once and cached in _disorder_wt_scores.pkl;
thresholds are applied to the cached scores. Track schema matches
build_alphafold_darkness: Hugo_Symbol,Protein_position,is_dark,<score>,source.

Run:
  python src/build_metapredict_disorder.py --threshold 0.5 \
      --output data/processed/disorder_proteome_metapredict_t0.5.csv
"""
from pathlib import Path
import argparse
import pickle

import numpy as np
import pandas as pd

from build_alphafold_darkness import fetch_uniprot_mapping
from score_disorder_effects import load_uniprot_seqs, batch_scores, is_scoreable, WT_CACHE

ROOT = Path(__file__).resolve().parents[1]
MAPPING = ROOT / "data" / "processed" / "_uniprot_gene_accession.csv"


def build_wt_cache(batch_size=2000):
    """Score all reviewed human proteins once; return {accession: np.float32 array}."""
    import metapredict as mp

    if WT_CACHE.exists():
        with open(WT_CACHE, "rb") as f:
            cache = pickle.load(f)
        print("WT cache loaded:", len(cache), "proteins")
        return cache

    seqs = load_uniprot_seqs()
    print("UniProt sequences:", len(seqs))
    items = sorted((a, s) for a, s in seqs.items() if is_scoreable(s))
    skipped = len(seqs) - len(items)
    if skipped:
        print(f"Skipped {skipped} proteins with non-standard residues (not encodable by metapredict)")
    cache = {}
    for i in range(0, len(items), batch_size):
        chunk = items[i:i + batch_size]
        preds = mp.predict_disorder_batch([s for _, s in chunk], show_progress_bar=False)
        for (acc, _), pred in zip(chunk, preds):
            cache[acc] = np.asarray(batch_scores(pred), dtype=np.float32)
        print(f"  scored {min(i + batch_size, len(items))}/{len(items)}", flush=True)
    WT_CACHE.parent.mkdir(parents=True, exist_ok=True)
    with open(WT_CACHE, "wb") as f:
        pickle.dump(cache, f)
    print("WT cache saved:", len(cache), "proteins")
    return cache


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--output", required=True)
    parser.add_argument("--batch-size", type=int, default=2000)
    args = parser.parse_args()

    cache = build_wt_cache(args.batch_size)
    mapping = fetch_uniprot_mapping(MAPPING)
    print("Gene-accession entries:", len(mapping))

    src = f"metapredict>={args.threshold:g}"
    written = 0
    done = set()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        import csv
        w = csv.writer(f)
        w.writerow(["Hugo_Symbol", "Protein_position", "is_dark", "disorder_score", "source"])
        for gene, grp in mapping.groupby("gene", sort=False):
            if gene in done:
                continue
            for acc in grp["accession"]:
                scores = cache.get(acc)
                if scores is None:
                    continue
                for pos, val in enumerate(scores, start=1):
                    w.writerow([gene, pos, int(val >= args.threshold), float(val), src])
                done.add(gene)
                written += len(scores)
                break
    print("Disorder track rows:", written, "genes:", len(done))
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
