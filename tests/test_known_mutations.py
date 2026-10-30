import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd
from map_jhu_coordinates import coord_to_protein, get_ref, load_uniprot, CDS


def _check_hotspot(chrom, pos, expected_pos, expected_aa, uniprot_acc, expected_ref):
    ref = get_ref()
    assert ref[chrom][pos - 1] == expected_ref, \
        f"unexpected ref base at {chrom}:{pos}: {ref[chrom][pos - 1]} (check coordinate)"
    cds_df = pd.read_csv(CDS, low_memory=False)
    trans = cds_df[(cds_df["chrom"] == chrom) & (cds_df["start"] <= pos)
                   & (cds_df["end"] >= pos)]["transcript_id"].unique()
    assert len(trans) > 0, f"no CDS exon covers {chrom}:{pos}"
    hits = set()
    for t in trans:
        pp, prot = coord_to_protein(t, chrom, pos, ref, cds_df)
        if pp is not None:
            hits.add((pp, prot[pp - 1]))
    assert (expected_pos, expected_aa) in hits, \
        f"{chrom}:{pos} mapped to {sorted(hits)[:10]}, expected codon {expected_pos} ({expected_aa})"
    uni = load_uniprot()[uniprot_acc]
    assert uni[expected_pos - 1] == expected_aa, \
        f"UniProt {uniprot_acc} position {expected_pos} is {uni[expected_pos - 1]}, not {expected_aa}"


def test_tp53_r175h():
    # c.524G>A (p.R175H), minus strand -> genomic C>T at hg19 chr17:7578406
    _check_hotspot("chr17", 7578406, 175, "R", "P04637", "C")


def test_kras_g12c():
    # c.34G>T (p.G12C), minus strand -> genomic C>A at hg19 chr12:25398285 (GRCh37)
    _check_hotspot("chr12", 25398285, 12, "G", "P01116", "C")
