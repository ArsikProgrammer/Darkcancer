"""Build the transition background: per SBS96 channel, the distribution of
disorder-transition classes over non-synonymous coding opportunities (plan 21.3).

Universe: representative GENCODE v19 transcripts (longest CDS per gene, as in
build_background) whose translated protein is identical to the gene's UniProt
sequence (exact coordinates). At every coding site each of the 3 SNV alleles is
translated; synonymous and stop outcomes are excluded. Each remaining
opportunity yields (context_96, alt residue) and, after mutant scoring, a
transition class at the mutated position (threshold rule).

Per channel, class fractions are estimated from K sampled opportunities
(reservoir sampling, seeded). Output columns: context_96, n_sites, n_sampled,
n_oo, n_od, n_do, n_dd, p_oo..p_dd, se_oo..se_dd.

Run:
  python src/build_transition_background.py --threshold 0.5 \
      --output data/processed/transition_background_metapredict_t0.5.csv
"""
from pathlib import Path
import argparse
import gzip
import re
from collections import defaultdict

import numpy as np
import pandas as pd

from build_background import load_reference, transcript_symbols
from build_alphafold_darkness import fetch_uniprot_mapping
from context_utils import context_for_position, rc, COMP_MAP
from map_mutations_to_residues import translate
from score_disorder_effects import load_wt_scores, score_substitutions, transition_class

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "data" / "raw" / "reference" / "hg19.fa"
GTF = ROOT / "data" / "raw" / "gencode.v19.annotation.gtf.gz"
CDS = ROOT / "data" / "processed" / "gencode_v19_cds.csv"
MAPPING = ROOT / "data" / "processed" / "_uniprot_gene_accession.csv"

CLASSES = ["O-O", "O-D", "D-O", "D-D"]


def rep_transcripts(cds):
    cds = cds.copy()
    cds["tlen"] = cds["end"] - cds["start"] + 1
    cds_total = cds.groupby("transcript_id")["tlen"].sum()
    t2g = cds.drop_duplicates("transcript_id").set_index("transcript_id")["gene_id"]
    keep = cds_total.groupby(t2g.reindex(cds_total.index)).idxmax()
    return cds[cds["transcript_id"].isin(set(keep.values))]


def transcript_coding_walk(grp, strand, chrom_seq):
    """Yield (coding_base, genomic_pos) in CDS order."""
    exons = grp.sort_values("start", ascending=(strand == "+"))
    for _, ex in exons.iterrows():
        positions = range(int(ex["start"]), int(ex["end"]) + 1)
        if strand == "-":
            positions = range(int(ex["end"]), int(ex["start"]) - 1, -1)
        for p in positions:
            if p - 1 < 0 or p > len(chrom_seq):
                return
            base = chrom_seq[p - 1]
            yield (COMP_MAP.get(base, "N") if strand == "-" else base), p


def collect_opportunities(cds_rep, sym, gene_accs, seqs, ref, rng, k_per_channel):
    """Walk all coding sites; reservoir-sample K non-synonymous opportunities per channel."""
    reservoir = defaultdict(list)   # context -> list of (acc, pos, ref_aa, alt_aa)
    counts = defaultdict(int)       # context -> opportunities seen
    n_sites = defaultdict(int)
    stats = {"genes_identical": 0, "genes_skipped_nonidentical": 0, "genes_no_symbol": 0,
             "genes_no_uniprot": 0, "ref_mismatch": 0, "stop_or_syn": 0}

    for tid, grp in cds_rep.groupby("transcript_id"):
        symbol = sym.get(tid) or sym.get(tid.split(".")[0])
        if symbol is None:
            stats["genes_no_symbol"] += 1
            continue
        accs = gene_accs.get(symbol)
        if not accs:
            stats["genes_no_uniprot"] += 1
            continue
        strand = grp.iloc[0]["strand"]
        chrom = grp.iloc[0]["chrom"]
        chrom_seq = ref.get(chrom)
        if chrom_seq is None:
            continue
        walk = list(transcript_coding_walk(grp, strand, chrom_seq))
        coding = "".join(b for b, _ in walk)
        protein = translate(coding)
        if protein.endswith("*"):
            protein = protein[:-1]  # CDS features include the stop codon; UniProt does not
        acc = next((a for a in accs if seqs.get(a) == protein), None)
        if acc is None:
            stats["genes_skipped_nonidentical"] += 1
            continue
        stats["genes_identical"] += 1

        n_codons = len(protein)
        for cidx in range(n_codons):
            pos = cidx + 1
            ref_aa = protein[cidx]
            if ref_aa not in "ACDEFGHIKLMNPQRSTVWY":
                continue
            codon = coding[3 * cidx:3 * cidx + 3]
            for j in range(3):
                _, gpos = walk[3 * cidx + j]
                gref = chrom_seq[gpos - 1]
                for galt in "ACGT":
                    if galt == gref:
                        continue
                    label = context_for_position(chrom_seq, gpos, gref, galt)
                    if label is None:
                        continue
                    ctx = label[2]
                    calt = COMP_MAP.get(galt, "N") if strand == "-" else galt
                    new_codon = codon[:j] + calt + codon[j + 1:]
                    alt_aa = translate(new_codon)
                    if alt_aa == ref_aa or alt_aa not in "ACDEFGHIKLMNPQRSTVWY":
                        stats["stop_or_syn"] += 1
                        continue
                    n_sites[ctx] += 1
                    counts[ctx] += 1
                    res = reservoir[ctx]
                    if len(res) < k_per_channel:
                        res.append((acc, pos, ref_aa, alt_aa))
                    else:
                        idx = int(rng.integers(0, counts[ctx]))
                        if idx < k_per_channel:
                            res[idx] = (acc, pos, ref_aa, alt_aa)
    return reservoir, n_sites, stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--sites-per-channel", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260115)
    parser.add_argument("--reservoir-cache", default=str(ROOT / "data" / "processed" / "_transition_reservoir.pkl"))
    parser.add_argument("--output", required=True)
    parser.add_argument("--reference", default=str(REF))
    parser.add_argument("--gtf", default=str(GTF))
    parser.add_argument("--cds", default=str(CDS))
    args = parser.parse_args()

    cds = pd.read_csv(args.cds, low_memory=False)
    cds_rep = rep_transcripts(cds)
    sym = transcript_symbols(args.gtf)
    print("Representative transcripts:", cds_rep["transcript_id"].nunique())

    mapping = fetch_uniprot_mapping(MAPPING)
    gene_accs = {g: list(grp["accession"]) for g, grp in mapping.groupby("gene", sort=False)}
    from score_disorder_effects import load_uniprot_seqs
    seqs = load_uniprot_seqs()
    ref = load_reference(args.reference)
    print("Reference chromosomes:", len(ref), "gene-accession entries:", len(gene_accs))

    rng = np.random.default_rng(args.seed)
    cache_path = Path(args.reservoir_cache)
    if cache_path.exists():
        import pickle
        with open(cache_path, "rb") as f:
            cached = pickle.load(f)
        reservoir, n_sites, stats = cached["reservoir"], cached["n_sites"], cached["stats"]
        print("Reservoir cache loaded:", cache_path)
    else:
        reservoir, n_sites, stats = collect_opportunities(
            cds_rep, sym, gene_accs, seqs, ref, rng, args.sites_per_channel)
        import pickle
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "wb") as f:
            pickle.dump({"reservoir": reservoir, "n_sites": n_sites, "stats": stats}, f)
        print("Reservoir cache saved:", cache_path)
    print("Walk stats:", stats)
    print("Channels with opportunities:", len(n_sites))

    rows = []
    for ctx in sorted(reservoir):
        rows.extend(reservoir[ctx])
    subs = pd.DataFrame(rows, columns=["uniprot_accession", "uniprot_pos", "ref_aa", "alt_aa"])
    subs = subs.drop_duplicates()
    print("Sampled substitutions to score:", len(subs))
    wt_cache = load_wt_scores()
    scored = score_substitutions(subs, wt_cache=wt_cache)
    scored["transition"] = [transition_class(w, m, args.threshold)
                            for w, m in zip(scored["wt_score"], scored["mut_score"])]
    key_cls = {(a, p, alt): tr for a, p, alt, tr in
               scored[["uniprot_accession", "uniprot_pos", "alt_aa", "transition"]]
               .itertuples(index=False, name=None)}

    out_rows = []
    for ctx in sorted(n_sites):
        res = reservoir[ctx]
        n = len(res)
        cnt = {c: 0 for c in CLASSES}
        for acc, pos, _ref, alt in res:
            cnt[key_cls[(acc, pos, alt)]] += 1
        row = {"context_96": ctx, "n_sites": n_sites[ctx], "n_sampled": n}
        for c in CLASSES:
            k = c.replace("-", "").lower()
            p = cnt[c] / n if n else float("nan")
            row[f"n_{k}"] = cnt[c]
            row[f"p_{k}"] = p
            row[f"se_{k}"] = float(np.sqrt(p * (1 - p) / n)) if n else float("nan")
        out_rows.append(row)

    out = pd.DataFrame(out_rows)
    if len(out) != 96:
        raise RuntimeError(f"expected 96 channels, got {len(out)}")
    if (out["n_sampled"] == 0).any():
        raise RuntimeError("some channels have no sampled opportunities")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Overall class fractions:",
          {c: round(out[f"n_{c.replace('-', '').lower()}"].sum() / out["n_sampled"].sum(), 4)
           for c in CLASSES})
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
