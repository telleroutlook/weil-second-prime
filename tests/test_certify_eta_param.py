"""Tests for --eta parameter in certify_fp_second.

Uses N=3 odd (b_L < 0) so the certify returns immediately without any
matrix build. Verifies that eta appears in the result dict and that the
Fraction parsing of --eta is correct.
"""
from __future__ import annotations

from fractions import Fraction

from checker.fp_second.certify_fp_second import certify_sector


def test_eta_in_result_dict_default():
    result = certify_sector(5600, 10000, "odd", 3)  # N=3: b_L < 0
    assert "eta" in result
    assert "eta_float" in result
    assert result["eta"] == "1/2"
    assert abs(result["eta_float"] - 0.5) < 1e-12
    assert result["certified"] is False


def test_eta_custom_fraction():
    eta = Fraction(249, 100)
    result = certify_sector(5600, 10000, "odd", 3, eta=eta)
    assert result["eta"] == "249/100"
    assert abs(result["eta_float"] - 2.49) < 1e-10


def test_eta_unity():
    eta = Fraction(1, 1)
    result = certify_sector(5600, 10000, "odd", 3, eta=eta)
    assert result["eta"] == "1"
    assert abs(result["eta_float"] - 1.0) < 1e-12
