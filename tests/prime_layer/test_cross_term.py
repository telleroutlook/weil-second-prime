"""
Tests for the cross-prime term F_{ij}(tau_2, tau_3) (S3).

F is the exchange Gram between the two prime shifts,
    F_{ij}(tau_2,tau_3) = <C_{tau_3,1} P_j, C_{tau_2,1} P_i>,
implemented exactly in Q[tau_2,tau_3] as a sum of four shift-strip integrals.

Ground-truth invariants (each an independent check that no strip was dropped —
omitting one is the C11 bug class):
  1. F_{ij}(tau,tau) == E_{ij}(tau)   (both operators are C_{tau,1}).
  2. F_{ij} == 0 for i+j odd           (reflection x->-x gives F=(-1)^{i+j}F).
  3. F_{ij}(t2,t3) == F_{ji}(t3,t2)     (swapping the two operators).
  4. In-window (tau_2+tau_3 = log6/L > 2) the B,C strips are empty; F = A+D.
"""

from __future__ import annotations

from fractions import Fraction

import pytest

from src.prime_layer.legendre_shift import compute_E
from src.prime_layer.legendre_shift_2prime import (
    compute_F,
    default_F_provider,
    tau2_at,
    tau3_at,
)

INDICES = list(range(0, 8))
WINDOW_POINTS = [0.55, 0.60, 0.65, 0.69]


class TestCrossTermInvariants:
    @pytest.mark.parametrize("tau_num,tau_den", [(7, 5), (3, 2), (9, 5)])
    def test_F_tau_tau_equals_E(self, tau_num: int, tau_den: int) -> None:
        """F_{ij}(tau,tau) must equal E_{ij}(tau) exactly (same operator)."""
        tau = Fraction(tau_num, tau_den)
        for i in INDICES:
            for j in INDICES:
                assert compute_F(i, j, tau, tau) == compute_E(i, j, tau), (
                    f"F_{{{i},{j}}}(tau,tau) != E_{{{i},{j}}}(tau) at tau={tau}"
                )

    def test_F_parity_vanishing(self) -> None:
        t2, t3 = Fraction(7, 5), Fraction(9, 5)
        for i in INDICES:
            for j in INDICES:
                if (i + j) % 2 == 1:
                    assert compute_F(i, j, t2, t3) == Fraction(0)

    def test_F_operator_swap_symmetry(self) -> None:
        t2, t3 = Fraction(7, 5), Fraction(9, 5)
        for i in INDICES:
            for j in INDICES:
                assert compute_F(i, j, t2, t3) == compute_F(j, i, t3, t2)

    def test_F_is_exact_fraction(self) -> None:
        F = compute_F(0, 2, Fraction(7, 5), Fraction(9, 5))
        assert isinstance(F, Fraction)


class TestWindowBehaviour:
    @pytest.mark.parametrize("L", WINDOW_POINTS)
    def test_BC_strips_empty_in_window(self, L: float) -> None:
        # tau_2 + tau_3 = log 6 / L > 2 for all L < log 2, so the mixed strips
        # B (tau3-1, 1-tau2) and C (tau2-1, 1-tau3) are empty.
        t2, t3 = tau2_at(L), tau3_at(L)
        assert t2 + t3 > 2
        assert t3 - 1 >= 1 - t2   # B empty
        assert t2 - 1 >= 1 - t3   # C empty

    @pytest.mark.parametrize("L", WINDOW_POINTS)
    def test_F00_closed_form(self, L: float) -> None:
        # With P_0 = 1, F_00 = |A| + |D| = 2*(2 - max(tau2,tau3)).
        t2, t3 = tau2_at(L), tau3_at(L)
        tmax = max(t2, t3)
        assert compute_F(0, 0, t2, t3) == 2 * (2 - tmax)

    @pytest.mark.parametrize("L", WINDOW_POINTS)
    def test_F_cauchy_schwarz(self, L: float) -> None:
        # F_{ij} = <C_{tau3}P_j, C_{tau2}P_i> is an inner product between two
        # DIFFERENT functions, so it is NOT sign-definite. The real bound is
        # cross-Cauchy-Schwarz: F_{ij}^2 <= E_{ii}(tau2) * E_{jj}(tau3).
        t2, t3 = tau2_at(L), tau3_at(L)
        for i in INDICES:
            for j in INDICES:
                if (i + j) % 2 == 1:
                    continue
                F = compute_F(i, j, t2, t3)
                bound = compute_E(i, i, t2) * compute_E(j, j, t3)
                assert F * F <= bound, f"CS violated F_{{{i},{j}}} at L={L}"


class TestProvider:
    def test_default_F_provider_matches_compute_F(self) -> None:
        t2, t3 = tau2_at(0.60), tau3_at(0.60)
        for i in [0, 2, 4]:
            for j in [0, 2, 4]:
                assert default_F_provider(i, j, t2, t3) == float(compute_F(i, j, t2, t3))
