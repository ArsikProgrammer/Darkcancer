"""Step 02: convert cBioPortal JHU API mutations to a flat table."""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

rows = json.load(open(RAW / "sclc_jhu_api.json"))
df = pd.DataFrame([{
    "Hugo_Symbol": r.get("gene", {}).get("hugoGeneSymbol"),
    "Tumor_Sample_Barcode": r.get("sampleId"),
    "Chromosome": r.get("chr"),
    "Start_Position": r.get("startPosition"),
    "End_Position": r.get("endPosition"),
    "Variant_Classification": r.get("mutationType"),
    "Variant_Type": r.get("variantType"),
    "Reference_Allele": r.get("referenceAllele"),
    "Tumor_Seq_Allele2": r.get("variantAllele"),
    "Transcript_ID": r.get("refseqMrnaId"),
    "HGVSp_Short": r.get("proteinChange"),
    "Protein_position": r.get("proteinPosStart"),
    "NCBI_Build": r.get("ncbiBuild"),
} for r in rows])

print("JHU rows:", len(df))
print(df.Variant_Type.value_counts().to_dict())
print("Samples:", df.Tumor_Sample_Barcode.nunique())
df.to_csv(OUT / "sclc_jhu_mutations.csv", index=False)
print("Saved:", OUT / "sclc_jhu_mutations.csv")
