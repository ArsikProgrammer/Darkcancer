"""Step 01: inspect George et al. 2015 SCLC MAF and clinical metadata."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

maf = pd.read_csv(RAW / "sclc_ucologne_2015_maf.txt", sep="\t", low_memory=False)
clinical = pd.read_csv(RAW / "sclc_ucologne_clinical_sample.txt", sep="\t", comment="#", low_memory=False)

print("MAF rows:", len(maf))
print("Unique samples:", maf.Tumor_Sample_Barcode.nunique())
print("\nVariant types:")
print(maf.Variant_Type.value_counts())
print("\nTop variant classifications:")
print(maf.Variant_Classification.value_counts().head(15))

coding_classes = {
    "Missense_Mutation",
    "Silent",
    "Nonsense_Mutation",
    "Splice_Site",
    "Splice_Region",
    "Translation_Start_Site",
}

snps = maf[
    (maf.Variant_Type == "SNP")
    & (maf.Variant_Classification.isin(coding_classes))
    & (maf.Reference_Allele != maf.Tumor_Seq_Allele2)
].copy()

snps = snps[snps.Tumor_Sample_Barcode.isin(
    clinical.loc[clinical.SAMPLE_TYPE == "Primary", "SAMPLE_ID"]
)].copy()
print("Primary-tumor coding SNPs:", len(snps))
print("Primary samples:", snps.Tumor_Sample_Barcode.nunique())

keep = [
    "Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Start_Position",
    "End_Position", "Strand", "Variant_Classification", "Variant_Type",
    "Reference_Allele", "Tumor_Seq_Allele2", "Transcript_ID", "HGVSc",
    "HGVSp", "HGVSp_Short", "Protein_position", "Codons",
]
snps[keep].to_csv(OUT / "george_coding_snps.csv", index=False)

merged = snps.merge(
    clinical[["SAMPLE_ID", "SAMPLE_TYPE", "TMB_NONSYNONYMOUS", "ONCOTREE_CODE"]],
    left_on="Tumor_Sample_Barcode",
    right_on="SAMPLE_ID",
    how="left",
)
print("\nSample types in merged SNP table:")
print(merged.SAMPLE_TYPE.value_counts(dropna=False))
print("\nSaved:", OUT / "george_coding_snps.csv")
