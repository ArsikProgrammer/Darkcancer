"""Build the O/E background: dark_fraction per SBS96 channel.

f_dark(c) = (# coding sites with channel c in dark residues)
          / (# coding sites with channel c with known darkness).

Computed from GENCODE v19 CDS x hg19 x a darkness track. Conventions match the
mutation pipeline: the trinucleotide is genomic (pyrimidine-normalized), the
darkness key is (gene symbol, protein position). One representative transcript
per gene (longest CDS). Sites without pLDDT (absent from the darkness track)
are excluded from both numerator and denominator.
"""
import gzip
import re
from pathlib import Path
import argparse
from collections import defaultdict

import pandas as pd

from context_utils import rc

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "data" / "raw" / "reference" / "hg19.fa"
GTF = ROOT / "data" / "raw" / "gencode.v19.annotation.gtf.gz"
CDS = ROOT / "data" / "processed" / "gencode_v19_cds.csv"


def load_reference(path):
    ref = {}
    chrom = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if chrom is not None:
                    ref[chrom] = "".join(parts).upper()
                chrom = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if chrom is not None:
        ref[chrom] = "".join(parts).upper()
    return ref


def transcript_symbols(gtf_path):
    """transcript_id(+version) -> gene_name, from CDS features."""
    out = {}
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parts[8]
            tid = re.search(r'transcript_id "([^"]+)"', attrs)
            ver = re.search(r'transcript_version "([^"]+)"', attrs)
            name = re.search(r'gene_name "([^"]+)"', attrs)
            if tid and name:
                key = tid.group(1) + ("." + ver.group(1) if ver else "")
                out[key] = name.group(1)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--darkness", required=True, help="darkness CSV: Hugo_Symbol,Protein_position,is_dark")
    parser.add_argument("--genes", default=None,
                        help="optional file with gene symbols to restrict the universe to")
    parser.add_argument("--output", required=True)
    parser.add_argument("--reference", default=str(REF))
    parser.add_argument("--gtf", default=str(GTF))
    parser.add_argument("--cds", default=str(CDS))
    args = parser.parse_args()

    dark_df = pd.read_csv(args.darkness, low_memory=False)
    dark_by_gene = defaultdict(dict)
    for g, p, d in dark_df[["Hugo_Symbol", "Protein_position", "is_dark"]].itertuples(index=False):
        if pd.notna(p) and pd.notna(d):
            dark_by_gene[str(g)][int(p)] = int(d)
    print("Darkness track genes:", len(dark_by_gene), "sites:", len(dark_df))

    sym = transcript_symbols(args.gtf)
    print("Transcript->symbol entries:", len(sym))

    cds = pd.read_csv(args.cds, low_memory=False)
    cds["tlen"] = cds["end"] - cds["start"] + 1
    cds_total = cds.groupby("transcript_id")["tlen"].sum()
    t2g = cds.drop_duplicates("transcript_id").set_index("transcript_id")["gene_id"]
    keep = cds_total.groupby(t2g.reindex(cds_total.index)).idxmax()
    rep = set(keep.values)
    cds = cds[cds["transcript_id"].isin(rep)]
    print("Representative transcripts:", len(rep))

    ref = load_reference(args.reference)
    print("Reference chromosomes:", len(ref))
    genes_ok = None
    if args.genes:
        genes_ok = {ln.strip() for ln in open(args.genes) if ln.strip()}
        print("Gene universe restricted to", len(genes_ok), "genes")

    total = defaultdict(int)
    dark = defaultdict(int)
    skipped_no_symbol = 0
    n_sites = 0

    for tid, grp in cds.groupby("transcript_id"):
        symbol = sym.get(tid) or sym.get(tid.split(".")[0])
        if symbol is None:
            skipped_no_symbol += len(grp)
            continue
        if genes_ok is not None and symbol not in genes_ok:
            continue
        gene_dark = dark_by_gene.get(symbol)
        if gene_dark is None:
            continue
        chrom = grp.iloc[0]["chrom"]
        seq = ref.get(chrom)
        if seq is None:
            continue
        strand = grp.iloc[0]["strand"]
        exons = grp.sort_values("start", ascending=(strand == "+"))
        cds_index = 0
        for _, ex in exons.iterrows():
            positions = range(int(ex["start"]), int(ex["end"]) + 1)
            if strand == "-":
                positions = range(int(ex["end"]), int(ex["start"]) - 1, -1)
            for p in positions:
                prot_pos = cds_index // 3 + 1
                cds_index += 1
                if p - 2 < 0 or p + 1 > len(seq):
                    continue
                tri = seq[p - 2:p + 1]
                if len(tri) != 3 or tri[1] not in "ACGT":
                    continue
                if tri[1] in "AG":
                    tri = rc(tri)
                d = gene_dark.get(prot_pos)
                if d is None:
                    continue
                total[tri] += 1
                dark[tri] += d
                n_sites += 1

    if skipped_no_symbol:
        print(f"WARNING: {skipped_no_symbol} CDS exons skipped (no gene_name for transcript)")
    print("Counted sites with known darkness:", n_sites)

    rows = []
    for tri in sorted(total):
        if total[tri] == 0:
            continue
        frac = dark[tri] / total[tri]
        for alt in "ACGT":
            if alt == tri[1]:
                continue
            rows.append({
                "context_96": f"{tri[0]}[{tri[1]}>{alt}]{tri[2]}",
                "dark_fraction": frac,
                "n_dark": dark[tri],
                "n_sites": total[tri],
            })
    out = pd.DataFrame(rows)
    if len(out) != 96:
        raise RuntimeError(f"expected 96 channels, got {len(out)}")
    if ((out["dark_fraction"] <= 0) | (out["dark_fraction"] >= 1)).any():
        raise RuntimeError("some channels have dark_fraction outside (0,1)")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Overall dark fraction:", round(sum(dark.values()) / sum(total.values()), 4))
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
