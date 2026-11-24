"""Annotate validated mutations with disorder-transition quantities (plan section 21).

Join order matches annotate_mutations_darkness: mapping rows are validated
(match flags) and keyed at (sample, gene, protein position); the disorder
scores are joined at (uniprot_accession, uniprot_pos, alt_aa). Transition
class = state_wt -> state_mut at the mutated position (threshold rule).

Silent mutations are the negative control (21.2): the protein sequence is
unchanged, so delta_pos = delta_window10 = 0 exactly and the class is O-O or
D-D; these invariants are asserted before writing.

Modes:
  --export-subs FILE   write the unique substitution list for scoring and exit
  --scores FILE        full annotation (after scoring) to --output

Run:
  python src/annotate_disorder_transitions.py --mutations data/processed/george_coding_snps_context.csv \
      --mapping data/processed/george_mapping.csv --export-subs data/processed/_subs_george.csv
  python src/annotate_disorder_transitions.py --mutations data/processed/george_coding_snps_context.csv \
      --mapping data/processed/george_mapping.csv --scores data/processed/disorder_substitution_scores.csv \
      --threshold 0.5 --output data/processed/disorder_transitions_george.csv
"""
from pathlib import Path
import argparse

import numpy as np
import pandas as pd

from annotate_mutations_darkness import load_mapping, parse_position
from score_disorder_effects import (AA_ALPHABET, MUT_CACHE, WT_CACHE, load_wt_scores,
                                    parse_hgvsp, transition_class)

MISSNSE = "Missense_Mutation"
SILENT = "Silent"


def join_validated(mut, mapping_paths):
    """Validate mutations via the mapping and add uniprot_accession/uniprot_pos/uniprot_aa."""
    mp = load_mapping(mapping_paths)
    if "Tumor_Sample_Barcode" in mp.columns:
        mp["Tumor_Sample_Barcode"] = mp["Tumor_Sample_Barcode"].astype(str)
    extra = [c for c in ["uniprot_accession", "uniprot_aa"] if c in mp.columns]
    key = ["Hugo_Symbol", "Protein_position"]
    if "Tumor_Sample_Barcode" in mp.columns and \
            set(mut["Tumor_Sample_Barcode"].astype(str)) & set(mp["Tumor_Sample_Barcode"].astype(str)):
        key = ["Tumor_Sample_Barcode"] + key
    else:
        print("WARNING: mapping sample ids do not match mutations; validating at (gene, position) level")
    mp = mp.sort_values("validated", ascending=False).drop_duplicates(key)
    n0 = len(mut)
    out = mut.merge(mp[key + ["validated", "uniprot_pos"] + extra], on=key, how="left",
                    validate="many_to_one")
    n_invalid = int((out["validated"] == False).sum())
    n_unmapped = int(out["validated"].isna().sum())
    out = out[out["validated"] == True].copy()
    print(f"Validation: {n0} mutations, kept {len(out)}, "
          f"residue mismatch {n_invalid}, no mapping {n_unmapped}")
    return out


def classify_mutations(mut):
    """Split validated mutations into (missense scored rows, silent rows, exclusions)."""
    counts = {"missense": 0, "silent": 0, "not_missense_or_silent": 0,
              "unparsable_hgvsp": 0, "ref_mismatch": 0, "pos_mismatch": 0}
    rows_mis, rows_sil = [], []
    for _, r in mut.iterrows():
        vc = str(r["Variant_Classification"])
        if vc == SILENT:
            rows_sil.append(r)
            counts["silent"] += 1
            continue
        if vc != MISSNSE:
            counts["not_missense_or_silent"] += 1
            continue
        parsed = parse_hgvsp(r.get("HGVSp_Short"))
        if parsed is None:
            counts["unparsable_hgvsp"] += 1
            continue
        ref, pos, alt = parsed
        if alt == ref or alt not in AA_ALPHABET:
            counts["unparsable_hgvsp"] += 1
            continue
        if "uniprot_aa" in r.index and pd.notna(r["uniprot_aa"]) and ref != str(r["uniprot_aa"]):
            counts["ref_mismatch"] += 1
            continue
        if pos != parse_position(r["Protein_position"]):
            counts["pos_mismatch"] += 1
            continue
        r = r.copy()
        r["ref_aa"], r["alt_aa"] = ref, alt
        rows_mis.append(r)
        counts["missense"] += 1
    mis = pd.DataFrame(rows_mis)
    sil = pd.DataFrame(rows_sil)
    print("Classification:", counts)
    return mis, sil, counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", required=True, help="context CSV with HGVSp_Short")
    parser.add_argument("--mapping", nargs="*", required=True)
    parser.add_argument("--export-subs", default=None,
                        help="write unique substitutions (uniprot_accession,uniprot_pos,ref_aa,alt_aa) and exit")
    parser.add_argument("--scores", default=None,
                        help="disorder_substitution_scores.csv (score_disorder_effects output)")
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    mut = pd.read_csv(args.mutations, low_memory=False)
    mut["Tumor_Sample_Barcode"] = mut["Tumor_Sample_Barcode"].astype(str)
    mut["Protein_position"] = mut["Protein_position"].map(parse_position)
    joined = join_validated(mut, args.mapping)
    mis, sil, counts = classify_mutations(joined)

    if args.export_subs:
        subs = mis[["uniprot_accession", "uniprot_pos", "ref_aa", "alt_aa"]].drop_duplicates()
        Path(args.export_subs).parent.mkdir(parents=True, exist_ok=True)
        subs.to_csv(args.export_subs, index=False)
        print("Unique substitutions:", len(subs), "-> saved:", args.export_subs)
        return

    if not args.scores:
        raise SystemExit("either --export-subs or --scores is required")
    if not args.output:
        raise SystemExit("--output is required with --scores")

    scores = pd.read_csv(args.scores, low_memory=False)
    scores = scores.drop_duplicates(["uniprot_accession", "uniprot_pos", "alt_aa"])
    mis = mis.merge(
        scores[["uniprot_accession", "uniprot_pos", "alt_aa", "wt_score", "mut_score",
                "delta_pos", "delta_window10", "window_len"]],
        on=["uniprot_accession", "uniprot_pos", "alt_aa"], how="left", validate="many_to_one")
    n_scored = int(mis["mut_score"].notna().sum())
    mis = mis[mis["mut_score"].notna()].copy()
    print(f"Missense with scores: {n_scored}/{counts['missense']} "
          f"(unscored dropped: {counts['missense'] - n_scored})")

    # silent rows: sequence unchanged -> zero deltas, class from WT state only
    wt_cache = load_wt_scores()
    sil["wt_score"] = [float(wt_cache[a][int(p) - 1])
                       for a, p in zip(sil["uniprot_accession"], sil["uniprot_pos"])]
    sil["mut_score"] = sil["wt_score"]
    sil["delta_pos"] = 0.0
    sil["delta_window10"] = 0.0
    sil["window_len"] = 0
    sil["ref_aa"] = sil.get("uniprot_aa")
    sil["alt_aa"] = sil["ref_aa"]

    thr = args.threshold
    for df in (mis, sil):
        df["state_wt"] = np.where(df["wt_score"] >= thr, "D", "O")
        df["state_mut"] = np.where(df["mut_score"] >= thr, "D", "O")
        df["transition"] = df["state_wt"] + "-" + df["state_mut"]

    # 21.2 invariants: silent mutations cannot change the protein
    assert (sil["delta_pos"] == 0).all() and (sil["delta_window10"] == 0).all()
    assert sil["transition"].isin(["O-O", "D-D"]).all()
    assert (mis["state_wt"] + "-" + mis["state_mut"] == mis["transition"]).all()

    keep_cols = ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification",
                 "context_96", "uniprot_accession", "uniprot_pos", "ref_aa", "alt_aa",
                 "wt_score", "mut_score", "delta_pos", "delta_window10", "window_len",
                 "state_wt", "state_mut", "transition"]
    out = pd.concat([mis, sil], ignore_index=True)
    out = out[[c for c in keep_cols if c in out.columns]]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Rows:", len(out), "missense:", len(mis), "silent:", len(sil))
    print("Transition counts:", out["transition"].value_counts().to_dict())
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
