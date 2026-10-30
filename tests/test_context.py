import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from context_utils import context_for_position, normalize_chrom


def test_context_for_position_basic():
    seq = "AAACCCGGGTTT"
    result = context_for_position(seq, 5, "C", "T")
    assert result is not None
    raw, norm, label = result
    assert raw == "CCC"
    assert label == "C[C>T]C"


def test_context_for_position_out_of_range():
    seq = "ACGT"
    assert context_for_position(seq, 1, "A", "T") is None


def test_normalize_chrom():
    assert normalize_chrom("23") == "X"
    assert normalize_chrom("24") == "Y"
    assert normalize_chrom("1") == "1"
    assert normalize_chrom("X") == "X"
