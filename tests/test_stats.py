"""Statistics tests on synthetic data only (no local data required)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from run_confirmatory import fit_nb2, verdict
from refit_signatures import opportunity_adjust


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
    b, se, crit, note = fit_nb2(df, "pred_std")
    assert np.isfinite(b)


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
