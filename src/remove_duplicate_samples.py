"""Remove duplicate samples across cohorts using overlap coefficient >=0.5.

Duplicates form connected components (A~B, B~C => A,B,C one component); exactly
one sample is kept per component: George (WGS) first, then more mutations, then
name. Pairwise log goes to duplicate_removal_log_full.csv, the authoritative
per-component decisions to duplicate_removal_log.csv.
"""
from pathlib import Path
from collections import defaultdict
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"

def load(name):
    df = pd.read_csv(PROC / name, low_memory=False)
    df["key"] = (df["Chromosome"].astype(str) + ":" + df["Start_Position"].astype(str)
                 + ":" + df["Reference_Allele"].astype(str) + ":" + df["Tumor_Seq_Allele2"].astype(str))
    return df

def overlap_pairs(a, b):
    ag = a.groupby("Tumor_Sample_Barcode")["key"].apply(set)
    bg = b.groupby("Tumor_Sample_Barcode")["key"].apply(set)
    out = []
    for sa, va in ag.items():
        for sb, vb in bg.items():
            denom = min(len(va), len(vb))
            if denom == 0:
                continue
            coef = len(va & vb) / denom
            if coef >= 0.5:
                out.append((sa, sb, len(va & vb), len(va), len(vb), coef))
    return out

def main():
    george = load("george_coding_snps.csv")
    jhu = load("sclc_jhu_coding_snps.csv")
    clcgp = load("sclc_clcgp_coding_snps.csv")

    cohorts = {"George": george, "JHU": jhu, "CLCGP": clcgp}
    pairs = []
    pairs += [("George", "JHU", *p) for p in overlap_pairs(george, jhu)]
    pairs += [("George", "CLCGP", *p) for p in overlap_pairs(george, clcgp)]
    pairs += [("JHU", "CLCGP", *p) for p in overlap_pairs(jhu, clcgp)]

    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    nmut = {}
    for cohort, df in cohorts.items():
        for s, keys in df.groupby("Tumor_Sample_Barcode")["key"].apply(set).items():
            nmut[(cohort, s)] = len(keys)
            find((cohort, s))

    log_rows = []
    for cohort_a, cohort_b, sa, sb, inter, na, nb, coef in pairs:
        union((cohort_a, sa), (cohort_b, sb))
        log_rows.append({
            "cohort_a": cohort_a, "sample_a": sa,
            "cohort_b": cohort_b, "sample_b": sb,
            "overlap_coeff": coef,
        })

    # within-cohort duplicates
    for cohort, df in cohorts.items():
        g = df.groupby("Tumor_Sample_Barcode")["key"].apply(set)
        names = list(g.index)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                va, vb = g[names[i]], g[names[j]]
                denom = min(len(va), len(vb))
                if denom == 0:
                    continue
                coef = len(va & vb) / denom
                if coef >= 0.5:
                    union((cohort, names[i]), (cohort, names[j]))
                    log_rows.append({
                        "cohort_a": cohort, "sample_a": names[i],
                        "cohort_b": cohort, "sample_b": names[j],
                        "overlap_coeff": coef,
                    })

    # one keeper per connected component: George first, then more mutations, then name
    comps = defaultdict(list)
    for node in parent:
        comps[find(node)].append(node)
    removed = set()
    decision_rows = []
    for nodes in comps.values():
        if len(nodes) == 1:
            continue
        keeper = max(nodes, key=lambda x: (x[0] == "George", nmut.get(x, 0), str(x[1])))
        for node in sorted(nodes):
            if node != keeper:
                removed.add(node)
                decision_rows.append({
                    "cohort": node[0], "sample": node[1],
                    "kept": f"{keeper[0]} {keeper[1]}",
                    "component_size": len(nodes),
                })

    pd.DataFrame(log_rows).to_csv(PROC / "duplicate_removal_log_full.csv", index=False)
    pd.DataFrame(decision_rows).to_csv(PROC / "duplicate_removal_log.csv", index=False)
    removed_samples = {s for _, s in removed}
    print("Removed samples:", len(removed_samples))

    # filter context counts
    for inp, outp in [("tumor_context_counts.csv", "tumor_context_counts_dedup.csv"),
                      ("sample_features.csv", "sample_features_dedup.csv")]:
        p = PROC / inp
        if not p.exists():
            continue
        df = pd.read_csv(p, low_memory=False)
        before = df["Tumor_Sample_Barcode"].nunique() if "Tumor_Sample_Barcode" in df.columns else len(df)
        df2 = df[~df["Tumor_Sample_Barcode"].isin(removed_samples)].copy()
        after = df2["Tumor_Sample_Barcode"].nunique() if "Tumor_Sample_Barcode" in df2.columns else len(df2)
        df2.to_csv(PROC / outp, index=False)
        print(f"{inp}: samples {before} -> {after}, rows {len(df)} -> {len(df2)}")

    print("Saved duplicate_removal_log_full.csv with", len(log_rows), "pairs")


if __name__ == "__main__":
    main()
