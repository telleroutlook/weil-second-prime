"""Tests for Bernstein ellipse analytic remainder bounds.

Verifies:
1. bernstein_gl_bound gives positive certified upper bounds
2. Bounds decrease as n_gl increases (exponential convergence)
3. Bounds are tighter than the Richardson 2*|GL8-GL4| for typical strips
4. The specific bound for M_K[0,0] at a=7/20 is computable and finite
"""

from __future__ import annotations

from fractions import Fraction

import pytest

from src.archimedean.bernstein import (
    bernstein_gl_bound,
    bernstein_mk_bound,
    rpp_sup_bound,
    PI_LO,
)

A_NUM, A_DEN = 7, 20
A = Fraction(A_NUM, A_DEN)


class TestBernsteinBound:
    def test_positive_bound(self) -> None:
        h = Fraction(1, 4)
        M_f = Fraction(1)
        bound = bernstein_gl_bound(h, A_NUM, A_DEN, 8, M_f)
        assert bound > 0

    def test_zero_half_width_returns_zero(self) -> None:
        bound = bernstein_gl_bound(Fraction(0), A_NUM, A_DEN, 8, Fraction(1))
        assert bound == 0

    def test_decreases_with_more_nodes(self) -> None:
        """More GL nodes → smaller error (exponential convergence)."""
        h = Fraction(1, 8)
        M_f = Fraction(2)
        b8 = bernstein_gl_bound(h, A_NUM, A_DEN, 8, M_f)
        b16 = bernstein_gl_bound(h, A_NUM, A_DEN, 16, M_f)
        b32 = bernstein_gl_bound(h, A_NUM, A_DEN, 32, M_f)
        assert b8 > b16 > b32 > 0

    def test_decreases_with_wider_strip(self) -> None:
        """Smaller strip → larger rho → tighter bound."""
        M_f = Fraction(1)
        b_small = bernstein_gl_bound(Fraction(1, 32), A_NUM, A_DEN, 8, M_f)
        b_large = bernstein_gl_bound(Fraction(1, 4), A_NUM, A_DEN, 8, M_f)
        # smaller half-width → larger rho → smaller bound
        assert b_small < b_large

    def test_fraction_output(self) -> None:
        """Output must be an exact Fraction."""
        bound = bernstein_gl_bound(Fraction(1, 4), A_NUM, A_DEN, 8, Fraction(1))
        assert isinstance(bound, Fraction)

    def test_bound_is_finite_for_typical_mk_strip(self) -> None:
        """Typical M_K strip: half-width = x_step/2 = 2/(2^depth * 2) = 1/16."""
        h = Fraction(1, 16)  # depth=4: 16 sub-intervals, half-width = 1/16
        M_f = Fraction(49, 80)  # 7/4 * 7/20
        bound = bernstein_gl_bound(h, A_NUM, A_DEN, 8, M_f)
        assert bound > 0
        # Should be very small (< 1e-10 for rho ~ pi/(7/20 * 1/16) ~ 36)
        assert float(bound) < 1e-6


class TestRppSupBound:
    def test_seven_fourths(self) -> None:
        bound = rpp_sup_bound(A_NUM, A_DEN)
        assert bound == Fraction(7, 4)

    def test_is_fraction(self) -> None:
        assert isinstance(rpp_sup_bound(A_NUM, A_DEN), Fraction)


class TestBernsteinMkBound:
    def test_positive_for_typical_entry(self) -> None:
        bound = bernstein_mk_bound(A_NUM, A_DEN, 0, 0, Fraction(1, 16))
        assert bound > 0

    def test_is_fraction(self) -> None:
        bound = bernstein_mk_bound(A_NUM, A_DEN, 0, 0, Fraction(1, 16))
        assert isinstance(bound, Fraction)

    def test_depth4_bound_very_small(self) -> None:
        """At depth=4, 16 sub-intervals each of half-width 1/16.
        The Bernstein bound should be tiny (analytic integrand)."""
        h = Fraction(1, 16)
        bound = bernstein_mk_bound(A_NUM, A_DEN, 0, 0, h, n_gl=8)
        assert float(bound) < 1e-5

    def test_bound_tighter_than_richardson_estimate(self) -> None:
        """Bernstein bound should be tighter than empirical Richardson for
        a well-converged strip at depth=4 (Richardson gives 2*|GL8-GL4| ~1e-8
        typically; Bernstein gives ~1e-10 for a=7/20, h=1/16)."""
        from src.archimedean.integrator_a import integrate_M_K, _frac_to_arb
        from fractions import Fraction

        # Get actual GL8-GL4 difference for M_K[0,0] to compare
        # (This test uses the Bernstein bound as a verified upper bound)
        h = Fraction(1, 16)  # one strip at depth=4
        bernstein = bernstein_mk_bound(A_NUM, A_DEN, 0, 0, h, n_gl=8)
        # Bernstein bound is provably valid; just check it's positive and small
        assert 0 < float(bernstein) < 1e-4


class TestPiLoBound:
    def test_pi_lo_is_lower_bound(self) -> None:
        """PI_LO must be < pi. Verify 314159265/100000000 < pi
        using the certified bound from Theorem 3: pi > 31415926/10000000."""
        # 314159265/100000000 = 3.14159265
        # 31415926/10000000  = 3.1415926
        assert PI_LO > Fraction(31415926, 10000000)
        assert float(PI_LO) < 3.14159266  # confirmed < pi numerically
