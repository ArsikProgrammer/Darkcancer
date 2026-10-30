"""Extract primary chromosomes 1-22, X, Y from Ensembl GRCh37 toplevel FASTA."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "raw" / "reference" / "GRCh37.toplevel.fa"
DST = ROOT / "data" / "raw" / "reference" / "GRCh37.primary.fa"

primary = {str(i) for i in range(1, 23)} | {"X", "Y"}
keep = False
written = 0
with SRC.open("r") as fin, DST.open("w") as fout:
    for line in fin:
        if line.startswith(">"):
            chrom = line[1:].split()[0]
            keep = chrom in primary
            if keep:
                written += 1
        if keep:
            fout.write(line)
print("Primary chromosome records written:", written)
print("Saved:", DST)
