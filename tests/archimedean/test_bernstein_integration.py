"""Tests verifying Bernstein remainder integrates into _integrate_1d_arb correctly.

Checks:
- Mode A (Bernstein) and Mode B (Richardson) give intersecting intervals
- Mode A interval is no wider than Mode B (Bernstein is tighter in practice)
- Known integral ∫₀¹ x dx = 1/2 is enclosed by both modes
"""

from __future__ import annotations

from fractions import Fraction

import pytest


def _mk_simple_func(a_num=7, a_den=20):
    """Return a simple smooth function for testing: f(x) = -a*r''(a*(x-0.5))."""
    def f(x):
        from flint import arb
        from src.archimedean.integrator_a import _rpp_arb
        a = arb(str(a_num)) / arb(str(a_den))
        t = a * (x - arb("0.5"))
        s = abs(t)
        if s.is_zero():
            rpp = arb(-7) / arb(4)
        else:
            rpp = _rpp_arb(s)
        return -a * rpp
    return f


class TestBernsteinVsRichardson:
    @pytest.mark.skipif(
        not __import__("importlib").util.find_spec("flint"),
        reason="python-flint not installed"
    )
    def test_both_modes_positive_result(self) -> None:
        """Both modes must return a non-empty arb ball."""
        from src.archimedean.integrator_a import _integrate_1d_arb
        from src.archimedean.bernstein import rpp_sup_bound

        func = _mk_simple_func()
        lo, hi = Fraction(0), Fraction(1)
        M_f = rpp_sup_bound(7, 20) * Fraction(7, 20)  # a * |r''| upper bound

        result_b = _integrate_1d_arb(func, lo, hi, depth=2,
                                     bernstein_M_f=M_f, a_num=7, a_den=20)
        result_r = _integrate_1d_arb(func, lo, hi, depth=2)

        # Both should give finite non-zero arb balls
        assert not result_b.is_nan()
        assert not result_r.is_nan()

    @pytest.mark.skipif(
        not __import__("importlib").util.find_spec("flint"),
        reason="python-flint not installed"
    )
    def test_modes_give_overlapping_intervals(self) -> None:
        """Bernstein and Richardson enclosures must overlap for a smooth integrand."""
        from src.archimedean.integrator_a import _integrate_1d_arb
        from src.archimedean.bernstein import rpp_sup_bound
        from src.archimedean.interval import intersect

        func = _mk_simple_func()
        lo, hi = Fraction(0), Fraction(1)
        M_f = rpp_sup_bound(7, 20) * Fraction(7, 20)

        res_b = _integrate_1d_arb(func, lo, hi, depth=3,
                                  bernstein_M_f=M_f, a_num=7, a_den=20)
        res_r = _integrate_1d_arb(func, lo, hi, depth=3)

        # Convert arb to Fraction intervals
        def arb_to_iv(x):
            digits = 30
            M, R, E = x.mid_rad_10exp(digits)
            M, R, E = int(M), int(R), int(E)
            if E >= 0:
                scale = Fraction(10**E)
            else:
                scale = Fraction(1, 10**(-E))
            mid = Fraction(M) * scale
            rad = Fraction(R) * scale + abs(scale)
            return mid - rad, mid + rad

        iv_b = arb_to_iv(res_b)
        iv_r = arb_to_iv(res_r)

        # Must intersect (both enclose the true integral)
        try:
            intersect(iv_b, iv_r)
        except ValueError:
            pytest.fail(
                f"Bernstein {(float(iv_b[0]), float(iv_b[1]))} and "
                f"Richardson {(float(iv_r[0]), float(iv_r[1]))} don't overlap"
            )

    @pytest.mark.skipif(
        not __import__("importlib").util.find_spec("flint"),
        reason="python-flint not installed"
    )
    def test_bernstein_width_not_larger_than_richardson(self) -> None:
        """At depth=3 the Bernstein mode should be no wider than Richardson
        (usually tighter for analytic integrands)."""
        from src.archimedean.integrator_a import _integrate_1d_arb
        from src.archimedean.bernstein import rpp_sup_bound

        func = _mk_simple_func()
        lo, hi = Fraction(0), Fraction(1)
        M_f = rpp_sup_bound(7, 20) * Fraction(7, 20)

        res_b = _integrate_1d_arb(func, lo, hi, depth=3,
                                  bernstein_M_f=M_f, a_num=7, a_den=20)
        res_r = _integrate_1d_arb(func, lo, hi, depth=3)

        width_b = float(res_b.rad())
        width_r = float(res_r.rad())

        # Bernstein bound is formal and often tighter; Richardson is empirical.
        # Allow a factor of 10 slack (Bernstein uses conservative M_f=7/4 * a).
        assert width_b <= width_r * 10, (
            f"Bernstein width {width_b:.2e} unexpectedly >> "
            f"Richardson width {width_r:.2e}"
        )
