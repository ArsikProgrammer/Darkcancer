"""Tests for the disorder-transition pipeline (plan section 21) pure logic."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from annotate_disorder_transitions import classify_mutations  # noqa: E402
from build_transition_background import collect_opportunities  # noqa: E402
from run_disorder_transitions import CHANGING, CLASSES, class_frames, prep_class  # noqa: E402
from score_disorder_effects import (parse_hgvsp, transition_class, window_delta)  # noqa: E402


def test_parse_hgvsp():
    assert parse_hgvsp("p.R175H") == ("R", 175, "H")
    assert parse_hgvsp("R175H") == ("R", 175, "H")
    assert parse_hgvsp("R175=") == ("R", 175, "R")
    assert parse_hgvsp("R175") == ("R", 175, "R")
    assert parse_hgvsp("p.R175*") == ("R", 175, "*")
    assert parse_hgvsp("M1?") is None
    assert parse_hgvsp("p.R175_H176del") is None
    assert parse_hgvsp(float("nan")) is None


def test_transition_class_threshold():
    assert transition_class(0.4, 0.6, 0.5) == "O-D"
    assert transition_class(0.6, 0.4, 0.5) == "D-O"
    assert transition_class(0.9, 0.8, 0.5) == "D-D"
    assert transition_class(0.1, 0.2, 0.5) == "O-O"
    # silent invariant: equal scores can never give O-D / D-O at any threshold
    for s in [0.0, 0.3, 0.5, 0.7, 1.0]:
        for thr in [0.4, 0.5, 0.6]:
            assert transition_class(s, s, thr) in ("O-O", "D-D")


def test_window_delta_clipping():
    wt = np.zeros(5)
    mut = np.array([0.0, 0.0, 1.0, 0.0, 0.0])
    d, wlen = window_delta(wt, mut, 0, half=10)
    assert wlen == 5
    assert d == pytest.approx(0.2)
    d, wlen = window_delta(wt, mut, 2, half=1)
    assert wlen == 3
    assert d == pytest.approx(1 / 3)


def _mut_df():
    return pd.DataFrame({
        "Tumor_Sample_Barcode": ["s1"] * 6 + ["s2"] * 4,
        "Hugo_Symbol": ["G"] * 10,
        "Protein_position": [1] * 10,
        "Variant_Classification": ["Missense_Mutation"] * 8 + ["Silent"] * 2,
        "HGVSp_Short": ["M1V"] * 4 + ["M1L"] * 4 + ["M1="] * 2,
        "context_96": ["A[C>T]G"] * 5 + ["A[C>G]G"] * 5,
        "uniprot_accession": ["A1"] * 10,
        "uniprot_pos": [1] * 10,
        "uniprot_aa": ["M"] * 10,
        "validated": [True] * 10,
    })


def test_classify_mutations():
    mis, sil, counts = classify_mutations(_mut_df())
    assert counts["missense"] == 8
    assert counts["silent"] == 2
    assert (mis["alt_aa"] == "V").sum() == 4 and (mis["alt_aa"] == "L").sum() == 4
    assert len(sil) == 2


def test_prep_class_oe_sums():
    mut = pd.DataFrame({
        "Tumor_Sample_Barcode": ["s1", "s1", "s1", "s2"],
        "context_96": ["c1", "c1", "c2", "c1"],
        "O": [1, 0, 1, 0],
    })
    sig = pd.DataFrame({
        "Tumor_Sample_Barcode": ["s1", "s2"],
        "total_mutations": [10, 20],
        "SBS4": [1.0, 5.0],
        "dataset": ["george", "jhu"],
    })
    pmap = {"c1": 0.25, "c2": 0.5}
    df, dropped, zero_e = prep_class(mut, pmap, sig, "SBS4_fraction")
    assert dropped == 0 and zero_e == 0
    s1 = df[df["Tumor_Sample_Barcode"] == "s1"].iloc[0]
    assert s1["O"] == 2 and s1["E"] == pytest.approx(1.0) and s1["n"] == 3
    # E = 0 cells are dropped
    df2, _, zero2 = prep_class(mut.assign(context_96="cX"), {"cX": 0.0}, sig, "SBS4_fraction")
    assert zero2 == 2 and len(df2) == 0


def test_class_frames_changing():
    mut = pd.DataFrame({
        "Tumor_Sample_Barcode": ["s1"] * 4 + ["s2"] * 2,
        "context_96": ["c1"] * 6,
        "transition": ["O-O", "O-D", "D-O", "D-D", "O-D", "D-D"],
    })
    bg = pd.DataFrame({
        "context_96": ["c1"],
        "p_oo": [0.4], "p_od": [0.1], "p_do": [0.1], "p_dd": [0.4],
        "n_sampled": [100],
    })
    sig = pd.DataFrame({
        "Tumor_Sample_Barcode": ["s1", "s2"],
        "total_mutations": [10, 12],
        "SBS4": [1.0, 4.0],
        "dataset": ["george", "jhu"],
    })
    frames = class_frames(mut, bg, sig, "SBS4_fraction")
    df_ch, _, _ = frames[CHANGING]
    s1 = df_ch[df_ch["Tumor_Sample_Barcode"] == "s1"].iloc[0]
    s2 = df_ch[df_ch["Tumor_Sample_Barcode"] == "s2"].iloc[0]
    assert s1["O"] == 2 and s1["E"] == pytest.approx(0.8)
    assert s2["O"] == 1 and s2["E"] == pytest.approx(0.4)
    df_oo, _, _ = frames["O-O"]
    s1oo = df_oo[df_oo["Tumor_Sample_Barcode"] == "s1"].iloc[0]
    assert s1oo["O"] == 1 and s1oo["E"] == pytest.approx(1.6)
    assert set(CLASSES) == {"O-O", "O-D", "D-O", "D-D"}


def test_collect_opportunities_synthetic():
    # chr1: 6 flanking C + ATGGCTTAA + 6 flanking C; exon covers the CDS
    chrom = "CCCCCC" + "ATGGCTTAA" + "CCCCCC"
    cds = pd.DataFrame([{"transcript_id": "T1", "gene_id": "G1", "chrom": "chr1",
                         "start": 7, "end": 15, "strand": "+"}])
    rng = np.random.default_rng(1)
    reservoir, n_sites, stats = collect_opportunities(
        cds, {"T1": "GENE1"}, {"GENE1": ["A1"]}, {"A1": "MA"}, {"chr1": chrom},
        rng, k_per_channel=5)
    opps = [o for res in reservoir.values() for o in res]
    keys = {(a, p, r, alt) for a, p, r, alt in opps}
    # M1 via A>C at chr1:7 -> CTG = L
    assert ("A1", 1, "M", "L") in keys
    # synonymous GCT->GCT-ish and stop TAA outcomes must never appear
    assert not any(alt == ref for _, _, ref, alt in keys)
    assert not any(p == 3 for _, p, _, _ in keys)  # stop codon never processed
    assert all(n <= 5 for n in (len(r) for r in reservoir.values()))
    assert stats["genes_identical"] == 1
    # determinism: same seed -> same reservoir
    rng2 = np.random.default_rng(1)
    reservoir2, _, _ = collect_opportunities(
        cds, {"T1": "GENE1"}, {"GENE1": ["A1"]}, {"A1": "MA"}, {"chr1": chrom},
        rng2, k_per_channel=5)
    assert {k: sorted(v) for k, v in reservoir.items()} == \
           {k: sorted(v) for k, v in reservoir2.items()}
