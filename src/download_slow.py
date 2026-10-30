"""Rate-limited AlphaFold pLDDT download for missing accessions."""
import time
import urllib.request
import pandas as pd
from pathlib import Path

CACHE = Path("data/processed/_plddt_needed.csv")
MAPPING = Path("data/processed/_uniprot_gene_accession.csv")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", nargs="+", required=True)
    args = parser.parse_args()

    mapping = pd.read_csv(MAPPING)
    cache = pd.read_csv(CACHE) if CACHE.exists() else pd.DataFrame(columns=["accession", "position", "plddt"])

    # needed positions
    frames = [pd.read_csv(m, low_memory=False)[["Hugo_Symbol", "Protein_position"]] for m in args.mutations]
    mut = pd.concat(frames, ignore_index=True)
    mut["Protein_position"] = pd.to_numeric(mut["Protein_position"], errors="coerce").dropna().astype(int)
    uniq = mut.drop_duplicates()

    # accessions needed
    needed_accs = set()
    acc_positions = {}
    for gene, grp in uniq.groupby("Hugo_Symbol"):
        accs = mapping[mapping["gene"] == gene]["accession"].tolist()
        positions = set(grp["Protein_position"].tolist())
        for acc in accs:
            needed_accs.add(acc)
            acc_positions.setdefault(acc, set()).update(positions)

    cached_accs = set(cache["accession"].unique()) if len(cache) else set()
    to_fetch = [a for a in needed_accs if a not in cached_accs]
    print(f"Needed: {len(needed_accs)}, cached: {len(cached_accs)}, to fetch: {len(to_fetch)}", flush=True)

    new_rows = []
    saved_count = 0
    for i, acc in enumerate(to_fetch, 1):
        positions = acc_positions[acc]
        url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.pdb"
        parsed = {}
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=30) as r:
                    text = r.read().decode("utf-8", "ignore")
                for line in text.splitlines():
                    if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                        continue
                    try:
                        rn = int(line[22:26])
                        bf = float(line[60:66])
                    except ValueError:
                        continue
                    if rn in positions and rn not in parsed:
                        parsed[rn] = bf
                break
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    time.sleep(2 ** attempt + 1)
                    continue
                break
            except Exception:
                time.sleep(1)
                break
        if parsed:
            for pos, val in parsed.items():
                new_rows.append({"accession": acc, "position": pos, "plddt": val})
        # save incrementally every 50
        if i % 50 == 0:
            if new_rows:
                new_df = pd.DataFrame(new_rows)
                cache = pd.concat([cache, new_df], ignore_index=True).drop_duplicates(["accession", "position"])
                cache.to_csv(CACHE, index=False)
                new_rows = []
            print(f"Progress: {i}/{len(to_fetch)}, cache rows: {len(cache)}", flush=True)
        time.sleep(0.5)
    # final save
    if new_rows:
        new_df = pd.DataFrame(new_rows)
        cache = pd.concat([cache, new_df], ignore_index=True).drop_duplicates(["accession", "position"])
    cache.to_csv(CACHE, index=False)
    print(f"Done. Cache rows: {len(cache)}", flush=True)


if __name__ == "__main__":
    main()
