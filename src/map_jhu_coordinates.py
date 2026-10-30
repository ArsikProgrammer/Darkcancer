"""Map JHU mutations via genomic coordinates to GENCODE CDS and UniProt."""
from pathlib import Path
import argparse
import json
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CDS = ROOT / "data" / "processed" / "gencode_v19_cds.csv"
REF_PKL = ROOT / "data" / "processed" / "_hg19_ref.pkl"
UNIPROT = ROOT / "data" / "processed" / "_uniprot_human_reviewed.fasta"
MAPPING = ROOT / "data" / "processed" / "_uniprot_gene_accession.csv"

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def translate_codon(c):
    return CODON_TABLE.get(c.upper(), "X")

def translate(seq):
    aa = []
    for i in range(0, len(seq)-2, 3):
        aa.append(translate_codon(seq[i:i+3]))
    return "".join(aa)

import pickle
def get_ref():
    with open(REF_PKL, "rb") as f:
        return pickle.load(f)

def load_uniprot():
    seqs = {}
    cur = None
    buf = []
    with open(UNIPROT) as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0].split("|")[1] if "|" in line else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs

_T_CACHE = {}
_SUB_INDEX = {}

def _sub_for(base, cds_df):
    if not _SUB_INDEX:
        tmp = cds_df.copy()
        tmp["_base"] = tmp["transcript_id"].str.split(".").str[0]
        for b, grp in tmp.groupby("_base"):
            _SUB_INDEX[b] = grp
    return _SUB_INDEX.get(base)

def coord_to_protein(transcript_id, chrom, pos, ref, cds_df):
    base = transcript_id.split(".")[0]
    ent = _T_CACHE.get(base)
    if ent is None:
        sub = _sub_for(base, cds_df)
        if sub is None or sub.empty:
            _T_CACHE[base] = False
            return None, None
        strand = sub["strand"].iloc[0]
        sub = sub.sort_values("start", ascending=(strand == "+"))
        parts = []
        for _, r2 in sub.iterrows():
            chrom2 = r2["chrom"] if r2["chrom"].startswith("chr") else "chr" + r2["chrom"]
            seq = ref[chrom2][int(r2["start"])-1:int(r2["end"])]
            if r2["strand"] == "-":
                comp = str.maketrans("ACGT", "TGCA")
                seq = seq.translate(comp)[::-1]
            parts.append(seq)
        cds = "".join(parts)
        prot = translate(cds)
        if "*" in prot:
            prot = prot.split("*")[0]
        exons, cum, c = [], [], 0
        for _, r2 in sub.iterrows():
            s, e = int(r2["start"]), int(r2["end"])
            exons.append((s, e))
            cum.append(c)
            c += e - s + 1
        ent = (strand, list(zip(exons, cum)), prot)
        _T_CACHE[base] = ent
    if ent is False:
        return None, None
    strand, exons_cum, prot = ent
    for (s, e), cum in exons_cum:
        if s <= pos <= e:
            offset = pos - s if strand == "+" else e - pos
            cds_pos = cum + offset
            prot_pos = cds_pos // 3 + 1
            if 1 <= prot_pos <= len(prot):
                return prot_pos, prot
            return None, None
    return None, None

mapping_df = None
def get_mapping():
    global mapping_df
    if mapping_df is None:
        mapping_df = pd.read_csv(MAPPING)
    return mapping_df

CDS_DF = None
def _load_cds():
    global CDS_DF
    if CDS_DF is None:
        CDS_DF = pd.read_csv(CDS, low_memory=False)
    return CDS_DF

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--output", default="data/processed/jhu_mapping_coord.csv")
    parser.add_argument("--mutations-csv", default=None,
                        help="MAF-like CSV (Hugo_Symbol, Tumor_Sample_Barcode, Chromosome, "
                             "Start_Position, Reference_Allele, Tumor_Seq_Allele2, "
                             "Protein_position) instead of the JHU API JSON")
    args = parser.parse_args()

    if args.mutations_csv:
        src = pd.read_csv(args.mutations_csv, low_memory=False)
        df = pd.DataFrame({
            "Hugo_Symbol": src["Hugo_Symbol"],
            "Tumor_Sample_Barcode": src["Tumor_Sample_Barcode"],
            "chrom": src["Chromosome"].astype(str).map(
                lambda c: c if c.startswith("chr") else "chr" + c),
            "pos": src["Start_Position"],
            "ref": src["Reference_Allele"],
            "alt": src["Tumor_Seq_Allele2"],
            "protein_pos_api": pd.to_numeric(src["Protein_position"], errors="coerce"),
            "mutationType": src.get("Variant_Classification"),
        })
    else:
        rows = json.load(open(ROOT / "data" / "raw" / "sclc_jhu_api.json"))
        df = pd.DataFrame([{
            "Hugo_Symbol": r.get("gene", {}).get("hugoGeneSymbol"),
            "Tumor_Sample_Barcode": r.get("sampleId"),
            "chrom": "chr" + str(r.get("chr")),
            "pos": r.get("startPosition"),
            "ref": r.get("referenceAllele"),
            "alt": r.get("variantAllele"),
            "protein_pos_api": r.get("proteinPosStart"),
            "mutationType": r.get("mutationType"),
        } for r in rows])

    # keep SNPs only
    df = df[df["ref"].str.len() == 1]
    df = df[df["alt"].str.len() == 1]
    if args.limit:
        df = df.head(args.limit)

    ref = get_ref()
    cds_df = pd.read_csv(CDS, low_memory=False)
    uniprot_seqs = None
    # lazy load
    def get_uniprot(acc):
        nonlocal uniprot_seqs
        if uniprot_seqs is None:
            seqs = {}
            cur = None
            buf = []
            with open(UNIPROT) as f:
                for line in f:
                    if line.startswith(">"):
                        if cur is not None:
                            seqs[cur] = "".join(buf)
                        cur = line[1:].split()[0].split("|")[1] if "|" in line else line[1:].split()[0]
                        buf = []
                    else:
                        buf.append(line.strip())
            if cur is not None:
                seqs[cur] = "".join(buf)
            uniprot_seqs = seqs
        return uniprot_seqs.get(acc)

    mapping = get_mapping()

    from build_background import transcript_symbols, GTF
    sym = transcript_symbols(GTF)
    print("Transcript->symbol entries:", len(sym))

    from bisect import bisect_right
    gene_acc = mapping.groupby("gene")["accession"].apply(list).to_dict()
    ex_index = {}
    for chrom2, grp in cds_df.groupby("chrom"):
        grp = grp.sort_values("start")
        ex_index[chrom2] = (grp["start"].tolist(), grp["end"].tolist(),
                            grp["transcript_id"].tolist())

    rows_out = []
    for _, r in df.iterrows():
        gene = r["Hugo_Symbol"]
        chrom = r["chrom"]
        if chrom == "chr23": chrom = "chrX"
        if chrom == "chr24": chrom = "chrY"
        pos = int(r["pos"])
        # find transcripts overlapping pos, preferring transcripts of the same gene
        starts, ends, tids = ex_index.get(chrom, ([], [], []))
        j = bisect_right(starts, pos)
        trans_all = []
        seen = set()
        k = j - 1
        while k >= 0 and pos - starts[k] <= 500000:
            if ends[k] >= pos and tids[k] not in seen:
                seen.add(tids[k])
                trans_all.append(tids[k])
            k -= 1
        trans_all.reverse()
        trans = [t for t in trans_all
                 if sym.get(t) == gene or sym.get(t.split(".")[0]) == gene] or trans_all

        api_pos = r.get("protein_pos_api")
        api_pos = float(api_pos) if pd.notna(api_pos) else None
        prot_pos = None
        prot = None
        transcript_used = None
        # prefer a transcript whose computed position agrees with the source
        # annotation (within +-5); fall back to the first that yields a position
        for pass_api in ([True, False] if api_pos is not None else [False]):
            for t in trans:
                pp, p = coord_to_protein(t, chrom, pos, ref, cds_df)
                if pp is None:
                    continue
                if pass_api and abs(pp - api_pos) > 5:
                    continue
                prot_pos = pp
                prot = p
                transcript_used = t
                break
            if transcript_used is not None:
                break

        accs = gene_acc.get(gene, [])
        uniprot_seq = None
        acc_used = None
        for acc in accs:
            s = get_uniprot(acc)
            if s:
                uniprot_seq = s
                acc_used = acc
                break

        translated_aa = None
        if prot is not None and prot_pos is not None and 1 <= prot_pos <= len(prot):
            translated_aa = prot[prot_pos-1]
        uniprot_aa = None
        if uniprot_seq is not None and prot_pos is not None and 1 <= prot_pos <= len(uniprot_seq):
            uniprot_aa = uniprot_seq[prot_pos-1]

        ref_aa = str(r["ref"]).upper()
        # check ref allele matches hg19
        ref_base = ref[chrom][pos-1:pos].upper()
        ref_match = (ref_base == ref_aa)

        match_uniprot = (uniprot_aa is not None and translated_aa is not None and uniprot_aa == translated_aa)

        rows_out.append({
            "sample": r["Tumor_Sample_Barcode"],
            "gene": gene,
            "chrom": chrom,
            "pos": pos,
            "ref_allele": ref_aa,
            "ref_match_hg19": ref_match,
            "transcript_used": transcript_used,
            "protein_position": prot_pos,
            "translated_aa": translated_aa,
            "uniprot_aa": uniprot_aa,
            "uniprot_accession": acc_used,
            "cds_found": prot is not None,
            "uniprot_found": uniprot_seq is not None,
            "match_uniprot": match_uniprot,
        })

    out = pd.DataFrame(rows_out)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Rows:", len(out))
    print("CDS found:", int(out.cds_found.sum()), "/", len(out))
    print("UniProt found:", int(out.uniprot_found.sum()), "/", len(out))
    print("Ref match hg19:", int(out.ref_match_hg19.sum()), "/", len(out))
    print("Match UniProt:", int(out.match_uniprot.sum()), "/", int(out.uniprot_found.sum()))
    print("Saved:", args.output)

if __name__ == "__main__":
    main()
