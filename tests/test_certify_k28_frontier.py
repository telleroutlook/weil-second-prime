"""Tests for the focused k=28 frontier certificate helper."""
from __future__ import annotations

from fractions import Fraction
import shutil
from pathlib import Path

from scripts.certify_k28_frontier import (
    certify,
    rayleigh_numerator,
    required_m0_keys,
)


def test_required_columns_share_cross_entry() -> None:
    keys = required_m0_keys()
    assert len(keys) == 55
    assert (0, 26) in keys
    assert (0, 27) in keys
    assert (26, 27) in keys
    assert (27, 26) not in keys


def test_integer_rayleigh_numerator_is_outward() -> None:
    c_aa = (Fraction(1), Fraction(3))
    c_ab = (Fraction(-8), Fraction(-7))
    c_bb = (Fraction(22), Fraction(24))
    value = rayleigh_numerator(c_aa, c_ab, c_bb)
    assert value[0] == 9 * 1 + 6 * (-8) + 22
    assert value[1] == 9 * 3 + 6 * (-7) + 24


def test_completed_checkpoint_replays_nonnegative_witness(tmp_path: Path) -> None:
    """Reassemble the focused result from its stored interval witnesses."""
    source = Path("pilots/cert_fp_second_N28_frontier_eta1_p512.frontier.ckpt.json")
    out = tmp_path / "replayed.json"
    shutil.copy(source, tmp_path / "replayed.frontier.ckpt.json")

    result = certify(
        56, 100, depth_2d=4, depth_3d=3, prec=512,
        out_path=out, resume=True,
    )
    lower, upper = map(Fraction, result["rayleigh_numerator_interval"])
    assert lower > 0
    assert result["rayleigh_witness_status"] == "nonnegative"
    assert result["certified_not_positive_definite"] is False
