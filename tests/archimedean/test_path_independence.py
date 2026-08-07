"""O2 path independence verification: Path A ∩ Path B non-empty for all M_K entries.

These tests verify the Archimedean dual-path requirement (INV-16 equivalent):
  - Path A: Arb GL-8/GL-4 with Richardson remainder
  - Path B: mpmath independent quadrature with Bernstein ellipse remainder
  - Both paths must give intersecting intervals for every M_K[n_i, n_j]

All 100 basis-pair entries (64 even + 36 odd) were verified to intersect
at depth=4 prec=256 on 2026-08-05. Tests here check the same property at
depth=2 for speed, plus spot-checks at depth=4 for high-degree entries.

The depth=2 test catches gross errors; the depth=4 spot-checks catch the
convergence issue where low-depth Path A diverges from Path B.
"""

from __future__ import annotations

from fractions import Fraction

import pytest

from src.archimedean.interval import intersect

A_NUM, A_DEN = 7, 20


def _mk_a(ni: int, nj: int, depth: int = 2, prec: int = 128):
    from src.archimedean.integrator_a import integrate_M_K
    return integrate_M_K(ni, nj, A_NUM, A_DEN, depth=depth, prec=prec).to_interval()


def _mk_b(ni: int, nj: int):
    from src.archimedean.integrator_b import integrate_M_K_path_b
    return integrate_M_K_path_b(ni, nj, A_NUM, A_DEN).to_interval()


def _intersects(ni: int, nj: int, depth: int = 2) -> bool:
    iv_a = _mk_a(ni, nj, depth=depth)
    iv_b = _mk_b(ni, nj)
    try:
        intersect(iv_a, iv_b)
        return True
    except ValueError:
        return False


class TestPathIntersectionEvenDiagonal:
    """Even sector diagonal entries: low-degree at depth=2, high-degree at depth=4."""

    @pytest.mark.parametrize("n", [0, 2, 4, 6])
    def test_low_degree_depth2(self, n: int) -> None:
        assert _intersects(n, n, depth=2), f"M_K[{n},{n}] A∩B empty at depth=2"

    @pytest.mark.parametrize("n", [8, 10, 12, 14])
    def test_high_degree_depth4(self, n: int) -> None:
        """High-degree entries require depth=4 for Path A to converge."""
        assert _intersects(n, n, depth=4), f"M_K[{n},{n}] A∩B empty at depth=4"


class TestPathIntersectionOddDiagonal:
    """Odd sector diagonal entries: low-degree at depth=2, high-degree at depth=4."""

    @pytest.mark.parametrize("n", [1, 3, 5, 7])
    def test_low_degree_depth2(self, n: int) -> None:
        assert _intersects(n, n, depth=2), f"M_K[{n},{n}] A∩B empty at depth=2"

    @pytest.mark.parametrize("n", [9, 11])
    def test_high_degree_depth4(self, n: int) -> None:
        assert _intersects(n, n, depth=4), f"M_K[{n},{n}] A∩B empty at depth=4"


class TestPathIntersectionHighDegree:
    """High-degree off-diagonal entries require depth=4 for convergence.

    M_K[14,14] at depth=2 gives Path A value -8.55e-7 which does NOT
    intersect Path B's [-8.43e-7, -8.39e-7]. At depth=4, Path A converges
    to -8.412e-7 which IS inside Path B's interval. This tests that the
    certify tier uses adequate depth.
    """

    def test_14_14_depth4(self) -> None:
        assert _intersects(14, 14, depth=4), "M_K[14,14] A∩B empty at depth=4"

    def test_12_14_depth4(self) -> None:
        assert _intersects(12, 14, depth=4), "M_K[12,14] A∩B empty at depth=4"

    def test_14_0_depth4(self) -> None:
        assert _intersects(14, 0, depth=4), "M_K[14,0] A∩B empty at depth=4"

    def test_0_14_depth4(self) -> None:
        assert _intersects(0, 14, depth=4), "M_K[0,14] A∩B empty at depth=4"


class TestPathConsistency:
    """Cross-check that Path A and Path B midpoints agree to 4 significant figures."""

    @pytest.mark.parametrize("ni,nj", [(0, 0), (2, 2), (4, 4), (1, 1), (3, 3)])
    def test_midpoints_agree(self, ni: int, nj: int) -> None:
        iv_a = _mk_a(ni, nj, depth=4)
        iv_b = _mk_b(ni, nj)
        mid_a = (iv_a[0] + iv_a[1]) / 2
        mid_b = (iv_b[0] + iv_b[1]) / 2
        if abs(float(mid_a)) < 1e-30:
            # Both essentially zero — just check intersection
            assert _intersects(ni, nj, depth=4)
            return
        rel_diff = abs(float(mid_a - mid_b)) / abs(float(mid_a))
        assert rel_diff < 1e-3, (
            f"M_K[{ni},{nj}] relative midpoint difference {rel_diff:.2e} > 1e-3: "
            f"A={float(mid_a):.6e} B={float(mid_b):.6e}"
        )
