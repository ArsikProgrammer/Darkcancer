"""Build SBS96 opportunity counts for the coding space and the whole genome.

Used to reweight COSMIC (genome-normalized) signatures to the coding-only
observation space: expected coding spectrum of a signature is proportional to
P_genome(c) * n_coding(c) / n_genome(c).

Output CSV: context_96, n_coding, n_genome (sites whose pyrimidine-normalized
trinucleotide supports that channel).
"""
from pathlib import Path
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

from context_utils import rc
from build_background import load_reference, REF, CDS

BASES = "ACGT"


def genome_trimer_counts(ref):
    """Counts of ACGT-only trinucleotides over the whole genome, index = 16i+4j+k."""
    lut = np.full(256, 255, dtype=np.uint8)
    for b, v in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
        lut[ord(b)] = v
    counts = np.zeros(64, dtype=np.int64)
    for chrom, seq in ref.items():
        m = lut[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
        step = 50_000_000
        for start in range(0, max(len(m) - 2, 0), step):
            seg = m[start:start + step + 2]
            a, b, c = seg[:-2], seg[1:-1], seg[2:]
            code = (a.astype(np.uint16) << 4) | (b.astype(np.uint16) << 2) | c.astype(np.uint16)
            valid = (a < 4) & (b < 4) & (c < 4)
            u, n = np.unique(code[valid], return_counts=True)
            counts[u] += n
    return counts


def coding_trimer_counts(cds, ref):
    """Counts of normalized trinucleotides over CDS positions (one transcript per gene)."""
    counts = defaultdict(int)
    for tid, grp in cds.groupby("transcript_id"):
        chrom = grp.iloc[0]["chrom"]
        seq = ref.get(chrom)
        if seq is None:
            continue
        strand = grp.iloc[0]["strand"]
        exons = grp.sort_values("start", ascending=(strand == "+"))
        for _, ex in exons.iterrows():
            positions = range(int(ex["start"]), int(ex["end"]) + 1)
            if strand == "-":
                positions = range(int(ex["end"]), int(ex["start"]) - 1, -1)
            for p in positions:
                if p - 2 < 0 or p + 1 > len(seq):
                    continue
                tri = seq[p - 2:p + 1]
                if len(tri) != 3 or tri[1] not in BASES:
                    continue
                if tri[1] in "AG":
                    tri = rc(tri)
                counts[tri] += 1
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--reference", default=str(REF))
    parser.add_argument("--cds", default=str(CDS))
    args = parser.parse_args()

    ref = load_reference(args.reference)
    print("Reference chromosomes:", len(ref))

    cds = pd.read_csv(args.cds, low_memory=False)
    cds["tlen"] = cds["end"] - cds["start"] + 1
    cds_total = cds.groupby("transcript_id")["tlen"].sum()
    t2g = cds.drop_duplicates("transcript_id").set_index("transcript_id")["gene_id"]
    keep = cds_total.groupby(t2g.reindex(cds_total.index)).idxmax()
    cds = cds[cds["transcript_id"].isin(set(keep.values))]
    print("Representative transcripts:", cds["transcript_id"].nunique())

    print("Counting genome trinucleotides...")
    g_counts = genome_trimer_counts(ref)
    genome = defaultdict(int)
    for code in range(64):
        tri = BASES[code >> 4] + BASES[(code >> 2) & 3] + BASES[code & 3]
        ntri = rc(tri) if tri[1] in "AG" else tri
        genome[ntri] += int(g_counts[code])

    print("Counting coding trinucleotides...")
    coding = coding_trimer_counts(cds, ref)

    rows = []
    for tri in sorted(genome):
        for alt in BASES:
            if alt == tri[1]:
                continue
            rows.append({
                "context_96": f"{tri[0]}[{tri[1]}>{alt}]{tri[2]}",
                "n_coding": coding.get(tri, 0),
                "n_genome": genome[tri],
            })
    out = pd.DataFrame(rows)
    if len(out) != 96:
        raise RuntimeError(f"expected 96 channels, got {len(out)}")
    if (out["n_genome"] == 0).any() or (out["n_coding"] == 0).any():
        raise RuntimeError("some channels have zero opportunities")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
