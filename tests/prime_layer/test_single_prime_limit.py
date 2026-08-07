"""
S2 acceptance gate as pytest: single-prime-limit self-check (prime-layer level).

Fast (no archimedean integrals): verifies the two-prime prime layer reduces
EXACTLY to weil-first's single-prime layer when c3 = 0, element-wise < 1e-10.
Also pins the C11 discipline: S2 with c3 != 0 and no F provider must RAISE, never
silently drop the cross term.

The slow assembled-C check (full Schur, one archimedean build ~1 min) lives in
scripts/single_prime_limit_check.py --full.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest

from scripts.single_prime_limit_check import (
    TOL,
    check_prime_layer,
    single_prime_layer,
)
from src.prime_layer.legendre_shift_2prime import (
    M2_two_prime,
    S2_two_prime,
    tau2_at,
    tau3_at,
    window_check,
    C2,
    C3,
)

# Test points inside the second-prime window (1/2 log3, log2) ~ (0.549, 0.693).
WINDOW_POINTS = [0.55, 0.60, 0.65, 0.69]


class TestWindow:
    def test_window_bounds(self) -> None:
        assert not window_check(0.549)      # below 1/2 log3
        assert window_check(0.60)           # inside
        assert not window_check(0.70)       # above log2

    def test_tau_in_single_hop_regime(self) -> None:
        # Both tau_2, tau_3 must be in (0,2): single-hop regime for the window.
        for L in WINDOW_POINTS:
            assert 0 < float(tau2_at(L)) < 2
            assert 0 < float(tau3_at(L)) < 2

    def test_c3_larger_than_c2(self) -> None:
        # c_p = log p / sqrt p; c3 = log3/sqrt3 > c2 = log2/sqrt2.
        assert C3 > C2 > 0


class TestSinglePrimeLimit:
    """c3 = 0 must reproduce the weil-first single-prime layer element-wise."""

    @pytest.mark.parametrize("L", WINDOW_POINTS)
    @pytest.mark.parametrize("sector", ["even", "odd"])
    def test_prime_layer_matches(self, L: float, sector: str) -> None:
        r = check_prime_layer(L, sector, N=4)
        assert r["max_dM2"] < TOL, f"M2 mismatch {r['max_dM2']:.2e} at L={L} {sector}"
        assert r["max_dS2"] < TOL, f"S2 mismatch {r['max_dS2']:.2e} at L={L} {sector}"
        assert r["pass"]

    def test_M2_reduces_to_minus_c2_J(self) -> None:
        # Explicit: M2_two_prime(c3=0)[a,b] == -C2 * J_{ij}(tau2).
        indices = [0, 2, 4]
        L = 0.60
        tau2 = tau2_at(L)
        from src.prime_layer.legendre_shift import compute_J
        M = M2_two_prime(indices, L, c3=0.0)
        for a, i in enumerate(indices):
            for b, j in enumerate(indices):
                expected = -C2 * float(compute_J(i, j, tau2))
                assert abs(M[a][b] - expected) < TOL

    def test_S2_reduces_to_c2sq_E(self) -> None:
        indices = [0, 2, 4]
        L = 0.60
        tau2 = tau2_at(L)
        from src.prime_layer.legendre_shift import compute_E
        S = S2_two_prime(indices, L, c3=0.0)
        for a, i in enumerate(indices):
            for b, j in enumerate(indices):
                expected = (C2 ** 2) * float(compute_E(i, j, tau2))
                assert abs(S[a][b] - expected) < TOL


class TestC11Discipline:
    """The cross term must never be silently omitted (C11 bug class)."""

    def test_S2_raises_without_F_when_c3_nonzero(self) -> None:
        indices = [0, 2, 4]
        with pytest.raises(NotImplementedError):
            S2_two_prime(indices, 0.60, c3=C3, F_provider=None)

    def test_S2_uses_F_when_provided(self) -> None:
        # A stub F provider is accepted and contributes (proves the wire, not
        # the physics — the real F is S3). Sentinel value makes the term visible.
        indices = [0, 2]
        L = 0.60

        def F_stub(i: int, j: int, t2: Fraction, t3: Fraction) -> float:
            return 1.0  # sentinel

        S_with = S2_two_prime(indices, L, c3=C3, F_provider=F_stub)
        # Compare to c3=0 baseline plus the c3^2 E term and 2 c2 c3 * 1.
        from src.prime_layer.legendre_shift import compute_E
        tau2, tau3 = tau2_at(L), tau3_at(L)
        for a, i in enumerate(indices):
            for b, j in enumerate(indices):
                expected = (C2 ** 2) * float(compute_E(i, j, tau2))
                expected += (C3 ** 2) * float(compute_E(i, j, tau3))
                expected += C2 * C3 * (1.0 + 1.0)  # F_ij + F_ji, both stubbed 1
                assert abs(S_with[a][b] - expected) < TOL

    def test_M2_includes_both_shifts_when_c3_nonzero(self) -> None:
        # M2 has no cross term but must include BOTH -c2 J(tau2) and -c3 J(tau3).
        indices = [0, 2]
        L = 0.60
        from src.prime_layer.legendre_shift import compute_J
        tau2, tau3 = tau2_at(L), tau3_at(L)
        M = M2_two_prime(indices, L, c3=C3)
        for a, i in enumerate(indices):
            for b, j in enumerate(indices):
                expected = -C2 * float(compute_J(i, j, tau2)) - C3 * float(compute_J(i, j, tau3))
                assert abs(M[a][b] - expected) < TOL
