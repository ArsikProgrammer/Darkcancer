"""Filter TCGA LUAD mutations to coding SNPs in primary tumors."""
from pathlib import Path
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", default="data/raw/luad_data_mutations.txt")
parser.add_argument("--output", default="data/processed/luad_coding_snps.csv")
args = parser.parse_args()

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / args.input
OUT = ROOT / args.output
OUT.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(RAW, sep="\t", comment="#", low_memory=False)
print("Input rows:", len(df), "samples:", df.Tumor_Sample_Barcode.nunique())

if not (df["NCBI_Build"].astype(str) == "GRCh37").all():
    bad = df["NCBI_Build"].astype(str) != "GRCh37"
    print(f"Dropping {int(bad.sum())} non-GRCh37 rows")
    df = df[~bad]

# primary tumors only (TCGA sample type code -01)
n_before = df["Tumor_Sample_Barcode"].nunique()
df = df[df["Tumor_Sample_Barcode"].astype(str).str.split("-").str[3].str.startswith("01")]
print(f"Primary tumor filter: {n_before} -> {df.Tumor_Sample_Barcode.nunique()} samples")

coding = {
    "Missense_Mutation", "Silent", "Nonsense_Mutation",
    "Splice_Site", "Splice_Region", "Translation_Start_Site",
    "Nonstop_Mutation", "NON_SYNONYMOUS_CODING",
}
out = df[(df.Variant_Type == "SNP") & (df.Variant_Classification.isin(coding))
         & (df.Reference_Allele != df.Tumor_Seq_Allele2)].copy()
keep = [
    "Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Start_Position",
    "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
    "Variant_Classification", "Variant_Type", "Protein_position",
    "HGVSp_Short", "Transcript_ID",
]
out[keep].to_csv(OUT, index=False)
print("LUAD coding SNPs:", len(out), "samples:", out.Tumor_Sample_Barcode.nunique())
print("Saved:", OUT)
