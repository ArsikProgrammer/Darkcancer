"""Filter CLCGP MAF to coding SNPs."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "sclc_clcgp_maf.txt"
OUT = ROOT / "data" / "processed" / "sclc_clcgp_coding_snps.csv"
OUT.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(RAW, sep="\t", low_memory=False)
coding = {
    "Missense_Mutation", "Silent", "Nonsense_Mutation",
    "Splice_Site", "Splice_Region", "Translation_Start_Site",
    "Nonstop_Mutation",
}
out = df[(df.Variant_Type == "SNP") & (df.Variant_Classification.isin(coding))
         & (df.Reference_Allele != df.Tumor_Seq_Allele2)].copy()
keep = [
    "Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Start_Position",
    "End_Position", "Strand", "Variant_Classification", "Variant_Type",
    "Reference_Allele", "Tumor_Seq_Allele2", "Transcript_ID", "HGVSc",
    "HGVSp", "HGVSp_Short", "Protein_position", "Codons",
]
out[keep].to_csv(OUT, index=False)
print("CLCGP coding SNPs:", len(out), "samples:", out.Tumor_Sample_Barcode.nunique())
print("Saved:", OUT)
