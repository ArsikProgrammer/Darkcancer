"""Annotate coding SNPs with a darkness label from an external darkness track.

With --mapping (outputs of map_mutations_to_residues / map_jhu_coordinates):
only mutations whose translated residue matches the UniProt residue are kept
(PRE_REGISTRATION), and the darkness lookup uses the UniProt position
(Protein_position + offset). Without --mapping the transcript position is used
directly as before. Unknown darkness stays NaN (never treated as ordered).
"""
import re
from pathlib import Path
import argparse
import pandas as pd


def parse_position(s):
    pm = re.match(r"\d+", str(s))
    return int(pm.group()) if pm else None


def load_mapping(paths):
    mp = pd.concat([pd.read_csv(p, low_memory=False) for p in paths], ignore_index=True)
    mp = mp.rename(columns={"gene": "Hugo_Symbol", "protein_position": "Protein_position",
                            "sample": "Tumor_Sample_Barcode"})
    need = {"Hugo_Symbol", "Protein_position"}
    if not need.issubset(mp.columns):
        raise ValueError(f"Mapping file must contain columns: {need}")
    checks = []
    for c in ["match_translated", "match_uniprot"]:
        if c in mp.columns:
            checks.append(mp[c].astype(str).str.lower().map({"true": True, "false": False}))
    if not checks:
        raise ValueError("Mapping file must contain match_translated and/or match_uniprot")
    mp["validated"] = checks[0]
    for ch in checks[1:]:
        mp["validated"] = mp["validated"] & ch
    mp["validated"] = mp["validated"].fillna(False)
    mp["Protein_position"] = mp["Protein_position"].map(parse_position)
    offset = mp["offset"].fillna(0) if "offset" in mp.columns else 0
    mp["uniprot_pos"] = mp["Protein_position"] + offset
    return mp.dropna(subset=["Protein_position"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", required=True, help="coding SNPs CSV")
    parser.add_argument("--darkness", required=True, help="darkness CSV: Hugo_Symbol,Protein_position,is_dark")
    parser.add_argument("--mapping", nargs="*", default=None,
                        help="validated mapping CSVs (map_mutations_to_residues / map_jhu_coordinates)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mut = pd.read_csv(args.mutations, low_memory=False)
    dark = pd.read_csv(args.darkness, low_memory=False)

    required = {"Hugo_Symbol", "Protein_position", "is_dark"}
    if not required.issubset(dark.columns):
        raise ValueError(f"Darkness file must contain columns: {required}")

    mut["Protein_position"] = mut["Protein_position"].map(parse_position)

    if args.mapping:
        mp = load_mapping(args.mapping)
        key = ["Hugo_Symbol", "Protein_position"]
        if "Tumor_Sample_Barcode" in mp.columns and \
                set(mut["Tumor_Sample_Barcode"].astype(str)) & set(mp["Tumor_Sample_Barcode"].astype(str)):
            key = ["Tumor_Sample_Barcode"] + key
        else:
            print("WARNING: mapping sample ids do not match mutations; validating at (gene, position) level")
        mp = mp.sort_values("validated", ascending=False).drop_duplicates(key)
        n0 = len(mut)
        mut = mut.merge(mp[key + ["validated", "uniprot_pos"]], on=key, how="left", validate="many_to_one")
        n_invalid = int((mut["validated"] == False).sum())
        n_unmapped = int(mut["validated"].isna().sum())
        mut = mut[mut["validated"] == True].copy()
        print(f"Validation: {n0} mutations, kept {len(mut)}, "
              f"residue mismatch {n_invalid}, no mapping {n_unmapped}")
    else:
        mut["uniprot_pos"] = mut["Protein_position"]

    dark = dark[["Hugo_Symbol", "Protein_position", "is_dark"]].drop_duplicates()
    dark["Protein_position"] = dark["Protein_position"].map(parse_position)
    dark = dark.dropna(subset=["Protein_position"]).rename(
        columns={"Protein_position": "uniprot_pos"})

    out = mut.merge(dark, on=["Hugo_Symbol", "uniprot_pos"], how="left", validate="many_to_one")
    out["is_dark"] = out["is_dark"]  # keep unknown as NaN; do not treat missing pLDDT as ordered
    out["has_plddt"] = out["is_dark"].notna()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Rows:", len(out), "with pLDDT:", int(out["has_plddt"].sum()),
          "missing:", int((~out["has_plddt"]).sum()))

    print("Saved:", args.output)


if __name__ == "__main__":
    main()
