"""Parse GENCODE v19 GTF and keep protein-coding CDS."""
from pathlib import Path
import gzip
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GTF = ROOT / "data" / "raw" / "gencode.v19.annotation.gtf.gz"
OUT = ROOT / "data" / "processed" / "gencode_v19_cds.csv"

rows = []
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) != 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = cols
        if feature != "CDS":
            continue
        if 'gene_type "protein_coding"' not in attrs:
            continue
        if 'transcript_type "protein_coding"' not in attrs:
            continue
        gene_id = attrs.split('gene_id "')[1].split('"')[0]
        transcript_id = attrs.split('transcript_id "')[1].split('"')[0]
        rows.append({
            "chrom": chrom,
            "start": int(start),
            "end": int(end),
            "strand": strand,
            "frame": frame,
            "gene_id": gene_id,
            "transcript_id": transcript_id,
        })

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print("CDS rows:", len(df), "transcripts:", df.transcript_id.nunique())
print("Saved:", OUT)
