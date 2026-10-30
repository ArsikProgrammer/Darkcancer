"""Download AlphaFold pLDDT only for residues that are actually mutated."""
import concurrent.futures
import gzip
import io
import json
import urllib.request
from pathlib import Path
import pandas as pd

CACHE = Path("data/processed/_plddt_needed.csv")
MAPPING = Path("data/processed/_uniprot_gene_accession.csv")


def fetch_mapping():
    if MAPPING.exists():
        return pd.read_csv(MAPPING)
    # fallback: fetch all reviewed human
    import urllib.request
    url = "https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606%20AND%20reviewed:true&format=json&fields=accession,gene_names&size=500"
    rows = []
    while url:
        with urllib.request.urlopen(url, timeout=60) as r:
            data = json.load(r)
            link = r.headers.get("Link", "")
        for e in data.get("results", []):
            acc = e.get("primaryAccession")
            for g in e.get("genes", []):
                name = g.get("geneName", {}).get("value")
                if acc and name:
                    rows.append({"accession": acc, "gene": name})
        import re
        m = re.search(r"<([^>]+)>", link)
        url = m.group(1) if (m and 'rel="next"' in link) else None
    df = pd.DataFrame(rows).drop_duplicates()
    MAPPING.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(MAPPING, index=False)
    return df


def parse_positions(pdb_text, positions):
    res = {}
    wanted = set(positions)
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            rn = int(line[22:26])
            bf = float(line[60:66])
        except ValueError:
            continue
        if rn in wanted and rn not in res:
            res[rn] = bf
    return res


def download_one(acc, positions):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.pdb"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            text = r.read().decode("utf-8", "ignore")
        parsed = parse_positions(text, positions)
        return acc, parsed, None
    except Exception as e:
        return acc, {}, str(e)[:120]


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", nargs="+", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    mapping = fetch_mapping()

    # needed positions per gene
    needed = {}
    for mfile in args.mutations:
        df = pd.read_csv(mfile, low_memory=False)
        for gene, grp in df.groupby("Hugo_Symbol"):
            needed.setdefault(str(gene), set()).update(
                int(x) for x in grp["Protein_position"].dropna().tolist()
                if pd.notna(x)
            )

    # map genes to accessions
    acc_positions = {}
    acc_gene = {}
    for gene, positions in needed.items():
        accs = mapping[mapping["gene"] == gene]["accession"].tolist()
        for acc in accs:
            acc_positions.setdefault(acc, set()).update(positions)
            acc_gene[acc] = gene

    print("Accessions to fetch:", len(acc_positions))

    # load cache
    if CACHE.exists():
        cache = pd.read_csv(CACHE)
    else:
        cache = pd.DataFrame(columns=["accession", "position", "plddt"])
    cached_accs = set(cache["accession"].unique()) if len(cache) else set()
    to_fetch = {a: p for a, p in acc_positions.items() if a not in cached_accs}
    print("Already cached:", len(cached_accs), "to fetch:", len(to_fetch))

    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(download_one, acc, positions): acc for acc, positions in to_fetch.items()}
        done = 0
        for fut in concurrent.futures.as_completed(futures):
            acc, parsed, err = fut.result()
            done += 1
            if parsed:
                for pos, val in parsed.items():
                    rows.append({"accession": acc, "position": pos, "plddt": val})
            if done % 100 == 0:
                print(f"Downloaded {done}/{len(to_fetch)}")
    if rows:
        new_df = pd.DataFrame(rows)
        cache = pd.concat([cache, new_df], ignore_index=True).drop_duplicates(["accession", "position"])
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        cache.to_csv(CACHE, index=False)
    print("Cache rows:", len(cache))
    print("Saved:", CACHE)


if __name__ == "__main__":
    main()
