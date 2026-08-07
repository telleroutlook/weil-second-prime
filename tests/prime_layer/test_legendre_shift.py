"""Tests for Legendre shift algebra (Theorem 4 cross-checks)."""

from __future__ import annotations

from fractions import Fraction

import pytest

from src.prime_layer.legendre_shift import (
    compute_E,
    compute_J,
    legendre_poly,
    prime_legendre_matrices,
)

# tau = log2 / L, L = 7/20; rational midpoint for tests
TAU = Fraction(20 * 842, 7 * 1215)  # ≈ 1.9806..., uses LOG2_LO


class TestLegendrePolynomials:
    def test_p0(self) -> None:
        assert legendre_poly(0) == (Fraction(1),)

    def test_p1(self) -> None:
        assert legendre_poly(1) == (Fraction(0), Fraction(1))

    def test_p2(self) -> None:
        # P_2(x) = (3x^2 - 1)/2
        p2 = legendre_poly(2)
        assert p2[0] == Fraction(-1, 2)
        assert p2[1] == Fraction(0)
        assert p2[2] == Fraction(3, 2)

    def test_p3(self) -> None:
        # P_3(x) = (5x^3 - 3x)/2
        p3 = legendre_poly(3)
        assert p3[0] == Fraction(0)
        assert p3[1] == Fraction(-3, 2)
        assert p3[2] == Fraction(0)
        assert p3[3] == Fraction(5, 2)


class TestParityVanishing:
    """J_{ij} and E_{ij} vanish when i+j is odd (Theorem 4 parity invariant)."""

    @pytest.mark.parametrize("i,j", [(0, 1), (1, 0), (1, 2), (2, 3), (0, 3)])
    def test_J_odd_parity_zero(self, i: int, j: int) -> None:
        assert compute_J(i, j, TAU) == Fraction(0)

    @pytest.mark.parametrize("i,j", [(0, 1), (1, 0), (1, 2), (2, 3), (0, 3)])
    def test_E_odd_parity_zero(self, i: int, j: int) -> None:
        assert compute_E(i, j, TAU) == Fraction(0)


class TestSampleValues:
    """Theorem 4 explicit sample values."""

    def test_J00(self) -> None:
        tau = Fraction(3, 2)  # generic tau in (1,2)
        result = compute_J(0, 0, tau)
        expected = 4 - 2 * tau
        assert result == expected

    def test_J11(self) -> None:
        tau = Fraction(3, 2)
        result = compute_J(1, 1, tau)
        expected = tau**3 / 3 - 2 * tau + Fraction(4, 3)
        assert result == expected

    def test_J02(self) -> None:
        tau = Fraction(3, 2)
        result = compute_J(0, 2, tau)
        expected = -(tau**3) + 3 * tau**2 - 2 * tau
        assert result == expected

    def test_J_symmetric_for_even_indices(self) -> None:
        tau = Fraction(3, 2)
        # J is not necessarily symmetric, but E is (it's a Gram matrix)
        e02 = compute_E(0, 2, tau)
        e20 = compute_E(2, 0, tau)
        assert e02 == e20


class TestMatrixStructure:
    """Even and odd sector matrix shapes and sparsity."""

    def test_even_sector_8x8(self) -> None:
        indices = list(range(0, 16, 2))  # [0,2,4,6,8,10,12,14]
        assert len(indices) == 8
        tau = Fraction(3, 2)
        J, E = prime_legendre_matrices(indices, tau)
        assert len(J) == 8
        assert all(len(row) == 8 for row in J)
        # All diagonal elements should be non-zero for even sector
        for k in range(8):
            assert J[k][k] != Fraction(0) or E[k][k] != Fraction(0)

    def test_odd_sector_6x6(self) -> None:
        indices = list(range(1, 12, 2))  # [1,3,5,7,9,11]
        assert len(indices) == 6
        tau = Fraction(3, 2)
        J, E = prime_legendre_matrices(indices, tau)
        assert len(J) == 6
        assert all(len(row) == 6 for row in J)

    def test_E_positive_diagonal(self) -> None:
        """E_{ii} = 2 * integral_{-1}^{1-tau} P_i^2 >= 0."""
        tau = Fraction(3, 2)
        for i in range(0, 8, 2):
            e = compute_E(i, i, tau)
            assert e >= 0, f"E_{{{i},{i}}} = {e} is negative"

    def test_no_cross_parity_coupling(self) -> None:
        """Even indices should not couple to odd indices."""
        tau = Fraction(3, 2)
        for i in range(0, 6, 2):
            for j in range(1, 6, 2):
                assert compute_J(i, j, tau) == Fraction(0)
                assert compute_E(i, j, tau) == Fraction(0)


class TestRationalOutput:
    """All outputs must be exact rationals."""

    def test_J_returns_fraction(self) -> None:
        result = compute_J(0, 0, TAU)
        assert isinstance(result, Fraction)

    def test_E_returns_fraction(self) -> None:
        result = compute_E(0, 0, TAU)
        assert isinstance(result, Fraction)

    def test_all_matrix_entries_are_fractions(self) -> None:
        tau = Fraction(3, 2)
        indices = [0, 2, 4]
        J, E = prime_legendre_matrices(indices, tau)
        for row in J:
            for val in row:
                assert isinstance(val, Fraction)
        for row in E:
            for val in row:
                assert isinstance(val, Fraction)
