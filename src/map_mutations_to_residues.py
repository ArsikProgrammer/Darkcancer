"""Map coding mutations to protein residues and validate against UniProt."""
from pathlib import Path
import argparse
import urllib.request
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CDS = ROOT / "data" / "processed" / "gencode_v19_cds.csv"
REF = ROOT / "data" / "raw" / "reference" / "hg19.fa"

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

def load_ref(path):
    seqs = {}
    cur = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf).upper()
    return seqs

REF_SEQS = None
import pickle
def get_ref():
    global REF_SEQS
    if REF_SEQS is None:
        pkl = ROOT / "data" / "processed" / "_hg19_ref.pkl"
        if pkl.exists():
            with open(pkl, "rb") as f:
                REF_SEQS = pickle.load(f)
        else:
            REF_SEQS = load_ref(REF)
    return REF_SEQS

_uniprot_fasta_cache = None
def load_uniprot_fasta():
    global _uniprot_fasta_cache
    if _uniprot_fasta_cache is not None:
        return _uniprot_fasta_cache
    path = ROOT / "data" / "processed" / "_uniprot_human_reviewed.fasta"
    seqs = {}
    cur = None
    buf = []
    with open(path) as f:
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
    _uniprot_fasta_cache = seqs
    return seqs

def fetch_uniprot_seq(acc):
    seqs = load_uniprot_fasta()
    return seqs.get(acc)

cds_df = None
def get_cds_exons():
    global cds_df
    if cds_df is None:
        cds_df = pd.read_csv(CDS, low_memory=False)
        cds_df["transcript_base"] = cds_df["transcript_id"].str.split(".").str[0]
    return cds_df

def get_cds(transcript_id):
    ex = get_cds_exons()
    base = transcript_id.split(".")[0]
    sub = ex[ex["transcript_base"] == base].copy()
    if sub.empty:
        return None
    # CDS exons sorted by genomic position; handle strand
    if sub["strand"].iloc[0] == "+":
        sub = sub.sort_values("start")
    else:
        sub = sub.sort_values("start", ascending=False)
    ref = get_ref()
    cds_parts = []
    for _, row in sub.iterrows():
        chrom = row["chrom"]
        key = chrom if chrom.startswith("chr") else "chr" + chrom
        if key not in ref:
            if chrom in ref:
                key = chrom
            else:
                return None
        seq = ref[key][int(row["start"])-1:int(row["end"])]
        if row["strand"] == "-":
            comp = str.maketrans("ACGT", "TGCA")
            seq = seq.translate(comp)[::-1]
        cds_parts.append(seq)
    cds = "".join(cds_parts)
    return cds

transcript_cache = {}
def get_protein(transcript_id):
    if transcript_id in transcript_cache:
        return transcript_cache[transcript_id]
    cds = get_cds(transcript_id)
    if cds is None:
        transcript_cache[transcript_id] = None
        return None
    prot = translate(cds)
    if "*" in prot:
        prot = prot.split("*")[0]
    transcript_cache[transcript_id] = prot
    return prot


mapping_cache = None
def get_mapping():
    global mapping_cache
    if mapping_cache is None:
        mapping_cache = pd.read_csv(ROOT / "data" / "processed" / "_uniprot_gene_accession.csv")
    return mapping_cache

def align_offset(prot, uniprot_seq):
    if prot == uniprot_seq:
        return 0
    if prot in uniprot_seq:
        return uniprot_seq.find(prot)
    if uniprot_seq in prot:
        return -prot.find(uniprot_seq)
    window = prot[:30]
    idx = uniprot_seq.find(window)
    if idx != -1:
        return idx
    for L in [20, 15, 10]:
        w = prot[:L]
        idx = uniprot_seq.find(w)
        if idx != -1:
            return idx
    return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    mut = pd.read_csv(args.mutations, low_memory=False)
    if args.limit:
        mut = mut.head(args.limit)

    mapping = get_mapping()

    rows = []
    for tid, grp in mut.groupby("Transcript_ID"):
        is_refseq = str(tid).startswith("NM_")
        if is_refseq:
            prot = None
        else:
            prot = get_protein(str(tid))
        gene = grp["Hugo_Symbol"].iloc[0]
        accs = mapping[mapping["gene"] == gene]["accession"].tolist()
        uniprot_seq = None
        acc_used = None
        for acc in accs:
            s = fetch_uniprot_seq(acc)
            if s:
                uniprot_seq = s
                acc_used = acc
                break
        offset = None
        if prot is not None and uniprot_seq is not None:
            offset = align_offset(prot, uniprot_seq)
        if is_refseq and uniprot_seq is not None:
            offset = 0
        for _, r in grp.iterrows():
            pos = r.get("Protein_position")
            try:
                pos = int(pos)
            except (ValueError, TypeError):
                pos = None
            ref_aa_codon = None
            codons = r.get("Codons")
            if isinstance(codons, str) and "/" in codons:
                ref_codon = codons.split("/")[0]
                ref_aa_codon = translate_codon(ref_codon)
            # JHU RefSeq: HGVSp_Short like G722C, M1?, X70_splice
            if ref_aa_codon is None:
                hgvsp = r.get("HGVSp_Short")
                if isinstance(hgvsp, str):
                    import re
                    mm = re.match(r"^([A-Z*])(\d+)", hgvsp)
                    if mm:
                        ref_aa_codon = mm.group(1)
            translated_aa = None
            if prot is not None and pos is not None and 1 <= pos <= len(prot):
                translated_aa = prot[pos-1]
            uniprot_aa = None
            if uniprot_seq is not None and offset is not None and pos is not None:
                upos = pos - 1 + offset
                if 0 <= upos < len(uniprot_seq):
                    uniprot_aa = uniprot_seq[upos]
            match_translated = (translated_aa is not None and ref_aa_codon is not None and translated_aa == ref_aa_codon)
            match_uniprot = (uniprot_aa is not None and ref_aa_codon is not None and uniprot_aa == ref_aa_codon)
            rows.append({
                "transcript_id": tid,
                "gene": gene,
                "sample": r.get("Tumor_Sample_Barcode"),
                "uniprot_accession": acc_used,
                "protein_position": pos,
                "ref_aa_codon": ref_aa_codon,
                "translated_aa": translated_aa,
                "uniprot_aa": uniprot_aa,
                "offset": offset,
                "match_translated": match_translated,
                "match_uniprot": match_uniprot,
                "cds_found": (not is_refseq) and (prot is not None),
                "uniprot_found": uniprot_seq is not None,
            })

    out = pd.DataFrame(rows)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Rows:", len(out))
    print("CDS found:", int(out.cds_found.sum()), "/", len(out))
    print("UniProt found:", int(out.uniprot_found.sum()), "/", len(out))
    print("Match translated:", int(out.match_translated.sum()), "/", int(out.cds_found.sum()))
    print("Match UniProt:", int(out.match_uniprot.sum()), "/", int(out.uniprot_found.sum()))
    print("Saved:", args.output)

if __name__ == "__main__":
    main()
