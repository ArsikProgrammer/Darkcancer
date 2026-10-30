"""Filter JHU mutations to coding SNPs."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "processed" / "sclc_jhu_mutations.csv"
OUT = ROOT / "data" / "processed" / "sclc_jhu_coding_snps.csv"
OUT.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(RAW, low_memory=False)
coding = {
    "Missense_Mutation", "Silent", "Nonsense_Mutation",
    "Splice_Site", "Splice_Region", "Translation_Start_Site",
    "Nonstop_Mutation", "NON_SYNONYMOUS_CODING",
}
if "NCBI_Build" in df.columns:
    bad_build = df["NCBI_Build"].astype(str) != "GRCh37"
    if bad_build.any():
        print(f"Dropping {int(bad_build.sum())} non-GRCh37 rows: "
              f"{df.loc[bad_build, 'NCBI_Build'].value_counts().to_dict()}")
        df = df[~bad_build]

# exclude cell lines (Rudin 2012: 30 tumors + 12 cell lines; pre-registration: tumors only)
clin = pd.read_csv(ROOT / "data" / "raw" / "sclc_jhu_clinical_sample.txt",
                   sep="\t", comment="#", dtype=str)
tumors = set(clin.loc[clin["SAMPLE_CLASS"].str.strip().str.lower() == "tumor", "SAMPLE_ID"])
n_before = df["Tumor_Sample_Barcode"].nunique()
df = df[df["Tumor_Sample_Barcode"].astype(str).isin(tumors)]
print(f"Sample filter (SAMPLE_CLASS == Tumor): {n_before} -> "
      f"{df['Tumor_Sample_Barcode'].nunique()} samples, "
      f"excluded {n_before - df['Tumor_Sample_Barcode'].nunique()} (cell lines etc.)")
if df["Tumor_Sample_Barcode"].nunique() != 30:
    print(f"WARNING: expected 30 JHU tumors, got {df['Tumor_Sample_Barcode'].nunique()}")

out = df[(df.Variant_Type == "SNP") & (df.Variant_Classification.isin(coding))
         & (df.Reference_Allele != df.Tumor_Seq_Allele2)].copy()
keep = [
    "Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Start_Position",
    "End_Position", "Variant_Classification", "Variant_Type",
    "Reference_Allele", "Tumor_Seq_Allele2", "Transcript_ID", "HGVSp_Short",
    "Protein_position",
]
out[keep].to_csv(OUT, index=False)
print("JHU coding SNPs:", len(out), "samples:", out.Tumor_Sample_Barcode.nunique())
print("Saved:", OUT)
