"""Extract trinucleotide context from a reference FASTA for a coding-SNP table."""
from pathlib import Path
import argparse
import pandas as pd

from context_utils import context_for_position, normalize_chrom

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "data" / "raw" / "reference" / "hg19.fa"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--reference", default=str(REF))
    args = parser.parse_args()

    inp = Path(args.input)
    outp = Path(args.output)
    ref = Path(args.reference)
    df = pd.read_csv(inp, low_memory=False)
    print("Input rows:", len(df), "samples:", df.Tumor_Sample_Barcode.nunique())

    groups = {}
    for idx, row in df.iterrows():
        chrom = "chr" + normalize_chrom(row["Chromosome"])
        groups.setdefault(chrom, []).append(idx)

    raw_results = [None] * len(df)
    norm_results = [None] * len(df)
    label_results = [None] * len(df)
    bad = 0

    current_chrom = None
    seq_parts = []

    def process_chrom(chrom, parts):
        nonlocal bad
        if chrom not in groups:
            return
        s = "".join(parts).upper()
        for idx in groups[chrom]:
            row = df.loc[idx]
            try:
                pos = int(row["Start_Position"])
            except (ValueError, TypeError):
                bad += 1
                continue
            ref_allele = str(row["Reference_Allele"]).upper()
            alt_allele = str(row["Tumor_Seq_Allele2"]).upper()
            result = context_for_position(s, pos, ref_allele, alt_allele)
            if result is None:
                bad += 1
                continue
            raw, norm, label = result
            raw_results[idx] = raw
            norm_results[idx] = norm
            label_results[idx] = label

    with ref.open("r") as f:
        for line in f:
            if line.startswith(">"):
                if current_chrom is not None:
                    process_chrom(current_chrom, seq_parts)
                current_chrom = line[1:].split()[0]
                seq_parts = []
            else:
                if current_chrom in groups:
                    seq_parts.append(line.strip())
        if current_chrom is not None:
            process_chrom(current_chrom, seq_parts)

    mask = [x is not None for x in label_results]
    ctx = pd.DataFrame({
        "raw_trinucleotide": [raw_results[i] for i in range(len(df)) if mask[i]],
        "norm_trinucleotide": [norm_results[i] for i in range(len(df)) if mask[i]],
        "context_96": [label_results[i] for i in range(len(df)) if mask[i]],
    })
    out = pd.concat([df.loc[mask].reset_index(drop=True), ctx.reset_index(drop=True)], axis=1)
    out.to_csv(outp, index=False)
    print("Rows with context:", len(out), "bad:", bad)
    print("Saved:", outp)


if __name__ == "__main__":
    main()
