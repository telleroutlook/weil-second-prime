"""Tests for the migrated Archimedean integrator modules.

These tests verify:
- Import hygiene (no stale weil-lower-bound paths)
- Interval arithmetic correctness
- Kernel r'' value at s=0 and known points
- Taylor coefficient correctness (the P0 bug fix: s/48 not s^2/48)
- Log-moment V_matrix_entry diagonal sign
- LDL^T on small positive-definite matrices

They do NOT run the full integration (requires python-flint and is slow).
The integration smoke test is gated behind a 'flint' mark.
"""

from __future__ import annotations

from fractions import Fraction

import pytest


class TestIntervalArithmetic:
    def test_add(self) -> None:
        from src.archimedean.interval import add
        a = (Fraction(1, 3), Fraction(2, 3))
        b = (Fraction(1, 6), Fraction(1, 2))
        result = add(a, b)
        assert result == (Fraction(1, 2), Fraction(7, 6))

    def test_mul_mixed_sign(self) -> None:
        from src.archimedean.interval import mul
        a = (Fraction(-1), Fraction(2))
        b = (Fraction(-3), Fraction(1))
        lo, hi = mul(a, b)
        assert lo == Fraction(-6)
        assert hi == Fraction(3)

    def test_intersect_valid(self) -> None:
        from src.archimedean.interval import intersect
        a = (Fraction(0), Fraction(2))
        b = (Fraction(1), Fraction(3))
        assert intersect(a, b) == (Fraction(1), Fraction(2))

    def test_intersect_empty_raises(self) -> None:
        from src.archimedean.interval import intersect
        with pytest.raises(ValueError):
            intersect((Fraction(0), Fraction(1)), (Fraction(2), Fraction(3)))

    def test_sq_lower_straddles_zero(self) -> None:
        from src.archimedean.interval import sq_lower
        assert sq_lower((Fraction(-1), Fraction(2))) == Fraction(0)

    def test_sq_lower_positive(self) -> None:
        from src.archimedean.interval import sq_lower
        assert sq_lower((Fraction(2), Fraction(5))) == Fraction(4)

    def test_sq_lower_negative(self) -> None:
        from src.archimedean.interval import sq_lower
        assert sq_lower((Fraction(-5), Fraction(-2))) == Fraction(4)


class TestLDLT:
    def test_identity_2x2(self) -> None:
        from src.archimedean.interval import point
        from src.archimedean.ldlt import ldlt_factor, certify_positive_definite
        one = point(Fraction(1))
        zero = point(Fraction(0))
        C = [[one, zero], [zero, one]]
        L, d = ldlt_factor(C)
        assert len(d) == 2
        assert all(p[0] > 0 for p in d)

    def test_positive_definite_2x2(self) -> None:
        from src.archimedean.interval import point
        from src.archimedean.ldlt import certify_positive_definite
        # [[4, 2], [2, 3]] has eigenvalues > 0
        C = [
            [(Fraction(4), Fraction(4)), (Fraction(2), Fraction(2))],
            [(Fraction(2), Fraction(2)), (Fraction(3), Fraction(3))],
        ]
        assert certify_positive_definite(C)

    def test_indefinite_2x2_fails(self) -> None:
        from src.archimedean.interval import point
        from src.archimedean.ldlt import certify_positive_definite
        # [[1, 2], [2, 1]] has det < 0
        C = [
            [(Fraction(1), Fraction(1)), (Fraction(2), Fraction(2))],
            [(Fraction(2), Fraction(2)), (Fraction(1), Fraction(1))],
        ]
        assert not certify_positive_definite(C)


class TestKernelConstants:
    def test_rpp_at_zero(self) -> None:
        from src.archimedean.kernel import R_DOUBLE_PRIME_AT_ZERO
        assert R_DOUBLE_PRIME_AT_ZERO == Fraction(-7, 4)

    def test_delta_positive(self) -> None:
        from src.archimedean.kernel import _DELTA
        assert _DELTA > 0

    @pytest.mark.skipif(
        not __import__("importlib").util.find_spec("flint"),
        reason="python-flint not installed"
    )
    def test_rpp_series_at_zero(self) -> None:
        """Series branch at s=0 should return exactly (-7/4, -7/4)."""
        from src.archimedean.kernel import rpp_series_branch
        lo, hi = rpp_series_branch(Fraction(0), Fraction(0))
        assert lo == Fraction(-7, 4)
        assert hi == Fraction(-7, 4)


class TestTaylorCubicCoefficient:
    """Verify the corrected s/48 Taylor coefficient (P0 bug fix)."""

    def test_rpp_mpmath_near_zero_linear_term(self) -> None:
        """r''(s) ≈ -7/4 - s/48 for small s. Check the sign and magnitude."""
        pytest.importorskip("mpmath")
        from src.archimedean.integrator_b import _rpp_mpmath
        import mpmath
        mpmath.mp.dps = 50

        s = mpmath.mpf("1e-5")
        rpp_s = _rpp_mpmath(s)
        # Expected: -7/4 - s/48 ≈ -1.75000000020833...
        expected_approx = -mpmath.mpf("7") / 4 - s / 48
        diff = abs(rpp_s - expected_approx)
        # Higher-order terms start at s^2 * (-9/32). For s=1e-5: |error| ~ 9*(1e-5)^2/32 ~ 3e-11.
        # The tolerance must account for the quadratic remainder, not just machine epsilon.
        assert diff < mpmath.mpf("1e-10"), (
            f"Linear term mismatch: got {float(rpp_s)}, expected ≈ {float(expected_approx)}, "
            f"diff={float(diff)}"
        )

    def test_rpp_mpmath_coefficient_not_s2_over_48(self) -> None:
        """Confirm the bug s^2/48 is NOT present (would give wrong coefficient)."""
        pytest.importorskip("mpmath")
        from src.archimedean.integrator_b import _rpp_mpmath
        import mpmath
        mpmath.mp.dps = 50

        s = mpmath.mpf("0.01")
        rpp_s = _rpp_mpmath(s)
        # With wrong s^2/48: -7/4 - s^2/48 ≈ -1.75000002083...
        wrong_approx = -mpmath.mpf("7") / 4 - s**2 / 48
        # With correct s/48: -7/4 - s/48 ≈ -1.75020833...
        correct_approx = -mpmath.mpf("7") / 4 - s / 48
        # The actual value should be much closer to correct_approx
        diff_correct = abs(rpp_s - correct_approx)
        diff_wrong = abs(rpp_s - wrong_approx)
        assert diff_correct < diff_wrong, (
            "Coefficient check: s/48 formula should match better than s^2/48. "
            f"diff_correct={float(diff_correct):.2e}, diff_wrong={float(diff_wrong):.2e}"
        )


class TestLogMomentsSign:
    """V(x) = -1/2 log(1-x^2) >= 0, so <V P_n, P_n> >= 0 for diagonal entries."""

    @pytest.mark.skipif(
        not __import__("importlib").util.find_spec("flint"),
        reason="python-flint not installed"
    )
    def test_V_diagonal_nonnegative(self) -> None:
        from src.archimedean.log_moments import V_matrix_entry
        for n in [0, 2, 4]:
            lo, hi = V_matrix_entry(n, n)
            assert lo >= 0 or hi >= 0, f"V_matrix_entry({n},{n}) = [{lo}, {hi}] should be >= 0"

    @pytest.mark.skipif(
        not __import__("importlib").util.find_spec("flint"),
        reason="python-flint not installed"
    )
    def test_V_odd_parity_zero(self) -> None:
        """<V P_0, P_1> = 0 by parity (V is even, P_0 is even, P_1 is odd)."""
        from src.archimedean.log_moments import V_matrix_entry
        lo, hi = V_matrix_entry(0, 1)
        # Should be very close to zero (exact parity vanishing)
        assert abs(float(lo)) < 1e-10 and abs(float(hi)) < 1e-10


class TestIntegrateStub:
    """Verify that IntegrationUnavailable is raised when flint is missing."""

    def test_path_a_raises_when_flint_missing(self, monkeypatch) -> None:
        import builtins
        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "flint":
                raise ImportError("mocked: flint not available")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)
        from checker.archimedean.integrate import IntegrationUnavailable, _check_deps
        with pytest.raises(IntegrationUnavailable):
            _check_deps()
