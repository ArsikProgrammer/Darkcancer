"""Statistics tests on synthetic data only (no local data required)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from run_confirmatory import fit_nb2, verdict, prep
from refit_signatures import opportunity_adjust
from annotate_mutations_darkness import parse_position


def _synthetic_df(n=400, beta=0.3, seed=7):
    rng = np.random.default_rng(seed)
    pred = rng.normal(0, 1, n)
    E = 60 * rng.lognormal(0, 0.5, n)
    O = rng.poisson(E * np.exp(beta * pred))
    return pd.DataFrame({
        "O": O, "E": E, "pred_std": pred,
        "TMB_std": rng.normal(0, 1, n), "wes": rng.binomial(1, 0.3, n),
    })


def test_fit_nb2_recovers_coefficient():
    beta_true = 0.3
    df = _synthetic_df(beta=beta_true)
    b, se, crit, note = fit_nb2(df, "pred_std")
    assert np.isfinite(b) and np.isfinite(se)
    assert abs(b - beta_true) < 0.15
    assert crit > 1.9


def test_fit_nb2_drops_constant_columns():
    df = _synthetic_df()
    df["wes"] = 0.0  # constant column must be dropped, not crash the fit
    import pytest
    with pytest.warns(UserWarning, match="[Hh]essian"):
        b, se, crit, note = fit_nb2(df, "pred_std")
    assert np.isfinite(b)  # NB2 hits the alpha boundary -> Poisson-t fallback is expected


def test_verdict_rules():
    assert verdict(0.9, 1.1) == "INFORMATIVE NULL"
    assert verdict(0.8, 0.95) == "EFFECT DETECTED"
    assert verdict(1.05, 1.5) == "EFFECT DETECTED"
    assert verdict(0.9, 1.5) == "INCONCLUSIVE"


def test_opportunity_adjust_shapes():
    idx = [f"{a}[{b}>{c}]{d}" for a in "ACGT" for b in "CT" for c in "ACGT"
           if c != b for d in "ACGT"]
    mat = pd.DataFrame({"SBS1": np.random.default_rng(1).uniform(0.001, 0.02, len(idx))},
                       index=idx)
    opp = pd.DataFrame({"context_96": idx,
                        "n_coding": 1000, "n_genome": 10000})
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "opp.csv"
        opp.to_csv(p, index=False)
        adj = opportunity_adjust(mat, p)
    assert adj.shape == mat.shape
    assert np.allclose(adj.sum(axis=0), 1.0)


def test_prep_oe_counts_same_mutation_set():
    """O and E must come from the same rows; NaN darkness excluded from both."""
    mut = pd.DataFrame({
        "Tumor_Sample_Barcode": ["s1"] * 4 + ["s2"] * 2 + ["s3"] * 2,
        "context_96": ["A[C>A]A"] * 3 + ["A[C>T]A", "A[C>A]A", "A[C>T]A",
                                        "A[C>A]A", "A[C>A]A"],
        "is_dark": [1, 1, 0, 1, 0, 0, np.nan, 1],
    })
    bg_map = {"A[C>A]A": 0.25, "A[C>T]A": 0.5}
    sig = pd.DataFrame({"Tumor_Sample_Barcode": ["s1", "s2", "s3"],
                        "SBS4_fraction": [0.4, 0.1, 0.2],
                        "total_mutations": [100, 50, 80],
                        "dataset": ["george", "jhu", "george"]})
    df, dropped = prep(mut, bg_map, sig, "SBS4_fraction")
    assert dropped == 0
    r1 = df[df.Tumor_Sample_Barcode == "s1"].iloc[0]
    assert r1["O"] == 3 and abs(r1["E"] - (0.25 * 3 + 0.5)) < 1e-9
    r2 = df[df.Tumor_Sample_Barcode == "s2"].iloc[0]
    assert r2["O"] == 0 and abs(r2["E"] - 0.75) < 1e-9
    # s3: the NaN row is excluded from BOTH O and E; the remaining row counts
    r3 = df[df.Tumor_Sample_Barcode == "s3"].iloc[0]
    assert r3["n"] == 1 and r3["O"] == 1 and abs(r3["E"] - 0.25) < 1e-9


def test_parse_position():
    assert parse_position("175/393") == 175
    assert parse_position(204) == 204
    assert parse_position("p.?") is None
