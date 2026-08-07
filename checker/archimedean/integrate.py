"""Archimedean primitive integrators for L = 7/20.

Implements Path A (Arb GL-8/GL-4 with Bernstein ellipse remainder) and
Path B (mpmath independent path) for all five primitive matrix blocks:
M_V, M_K, S_VV, S_VK, S_KK.

P0 bugs from weil-lower-bound are NOT present here:
  - integrate_M_K calls _integrate_1d_arb with GL-8/GL-4 Richardson remainder
  - _rpp_mpmath uses correct s/48 linear coefficient (not s^2/48)
  - All remainders use Bernstein ellipse analytic bounds

Both paths are independent: different arithmetic engines (Arb vs mpmath),
different splitting strategies, no shared approximation code.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

from src.archimedean.interval import Interval, add, intersect, point
from src.archimedean.integrator_a import (
    integrate_M_K as _a_mk,
    integrate_S_VK,
    integrate_S_KK,
    PathAResult,
)
from src.archimedean.log_moments import V_matrix_entry, V2_matrix_entry

# Frozen parameters for L = 7/20
A_NUM = 7
A_DEN = 20


class IntegrationUnavailable(Exception):
    """Raised when python-flint or mpmath is not installed."""


def _check_deps() -> None:
    try:
        import flint  # noqa: F401  # type: ignore[import]
        import mpmath  # noqa: F401  # type: ignore[import]
    except ImportError as exc:
        raise IntegrationUnavailable(
            "python-flint and mpmath are required for Archimedean integration"
        ) from exc


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_all_primitives_path_a(
    contract: dict[str, Any],
    precision: int = 256,
) -> dict[str, Any]:
    """Compute M_V, M_K, S_VV, S_VK, S_KK via Path A (Arb + GL-8/GL-4 remainder).

    Returns a dict with keys:
      M_V, M_K, S_VV, S_VK, S_KV, S_KK  ->  {(i,j): Interval}
      witnesses -> list of LeafWitnessA
    """
    _check_deps()

    sector = contract.get("sector", "even")
    index_set = contract.get("index_set", [0, 2, 4, 6, 8, 10, 12, 14])
    N = len(index_set)

    # M_V: exact formula via log-moment Beta derivatives
    M_V: dict[tuple[int, int], Interval] = {}
    for i, ni in enumerate(index_set):
        for j, nj in enumerate(index_set):
            M_V[(i, j)] = V_matrix_entry(ni, nj, precision)

    # M_K: Duffy 2D quadrature with GL-8/GL-4 Richardson remainder
    M_K: dict[tuple[int, int], Interval] = {}
    mk_witnesses = []
    for i, ni in enumerate(index_set):
        for j, nj in enumerate(index_set):
            result = _a_mk(ni, nj, A_NUM, A_DEN, depth=4, prec=precision)
            M_K[(i, j)] = result.to_interval()
            mk_witnesses.extend(result.leaves)

    # S_VV: exact formula
    S_VV: dict[tuple[int, int], Interval] = {}
    for i, ni in enumerate(index_set):
        for j, nj in enumerate(index_set):
            S_VV[(i, j)] = V2_matrix_entry(ni, nj, precision)

    # S_VK, S_KV, S_KK: compute directly per entry (Path A, depth=4/3)
    S_VK: dict[tuple[int, int], Interval] = {}
    S_KV: dict[tuple[int, int], Interval] = {}
    S_KK: dict[tuple[int, int], Interval] = {}
    for i, ni in enumerate(index_set):
        for j, nj in enumerate(index_set):
            svk = integrate_S_VK(ni, nj, A_NUM, A_DEN, depth=4, prec=precision)
            skv = integrate_S_VK(nj, ni, A_NUM, A_DEN, depth=4, prec=precision)
            skk = integrate_S_KK(ni, nj, A_NUM, A_DEN, depth=3, prec=precision)
            S_VK[(i, j)] = svk.to_interval()
            S_KV[(i, j)] = skv.to_interval()
            S_KK[(i, j)] = skk.to_interval()

    return {
        "M_V": M_V,
        "M_K": M_K,
        "S_VV": S_VV,
        "S_VK": S_VK,
        "S_KV": S_KV,
        "S_KK": S_KK,
        "witnesses": mk_witnesses,
    }


def compute_all_primitives_path_b(
    contract: dict[str, Any],
    precision: int = 256,
) -> dict[str, Any]:
    """Compute M_V, M_K, S_VV, S_VK, S_KV, S_KK via Path B (mpmath independent).

    Path B uses mpmath.quad with a different splitting strategy and the
    corrected _rpp_mpmath (s/48 linear term, not s^2/48) with Bernstein
    ellipse discretization remainder.
    """
    _check_deps()

    sector = contract.get("sector", "even")
    index_set = contract.get("index_set", [0, 2, 4, 6, 8, 10, 12, 14])

    from src.archimedean.integrator_b import (
        integrate_M_K_path_b,
        integrate_S_VK_path_b,
        integrate_S_KK_path_b,
        PathBResult,
    )

    # M_V: same exact formula (no quadrature needed, path-independent)
    M_V: dict[tuple[int, int], Interval] = {}
    for i, ni in enumerate(index_set):
        for j, nj in enumerate(index_set):
            M_V[(i, j)] = V_matrix_entry(ni, nj, precision)

    M_K: dict[tuple[int, int], Interval] = {}
    for i, ni in enumerate(index_set):
        for j, nj in enumerate(index_set):
            result = integrate_M_K_path_b(ni, nj, A_NUM, A_DEN)
            M_K[(i, j)] = result.to_interval()

    S_VV: dict[tuple[int, int], Interval] = {}
    for i, ni in enumerate(index_set):
        for j, nj in enumerate(index_set):
            S_VV[(i, j)] = V2_matrix_entry(ni, nj, precision)

    S_VK: dict[tuple[int, int], Interval] = {}
    S_KV: dict[tuple[int, int], Interval] = {}
    S_KK: dict[tuple[int, int], Interval] = {}
    for i, ni in enumerate(index_set):
        for j, nj in enumerate(index_set):
            S_VK[(i, j)] = integrate_S_VK_path_b(ni, nj, A_NUM, A_DEN).to_interval()
            S_KV[(i, j)] = integrate_S_VK_path_b(nj, ni, A_NUM, A_DEN).to_interval()
            S_KK[(i, j)] = integrate_S_KK_path_b(ni, nj, A_NUM, A_DEN).to_interval()

    return {
        "M_V": M_V,
        "M_K": M_K,
        "S_VV": S_VV,
        "S_VK": S_VK,
        "S_KV": S_KV,
        "S_KK": S_KK,
        "witnesses": [],
    }


def verify_intersection(
    primitives_a: dict[str, Any],
    primitives_b: dict[str, Any],
) -> dict[str, Any]:
    """Verify Path A ∩ Path B non-empty for every matrix entry.

    Returns:
      all_pass: bool
      checks: {key: bool}  — True if interval intersection is non-empty
      primitives: {key: Interval}  — intersection intervals (for CERTIFIED entries)
    """
    checks: dict[str, bool] = {}
    merged: dict[str, Any] = {}

    for block in ("M_V", "M_K", "S_VV", "S_VK", "S_KV", "S_KK"):
        a_block = primitives_a.get(block, {})
        b_block = primitives_b.get(block, {})
        for key in a_block:
            check_key = f"{block}_{key}"
            iv_a = a_block[key]
            iv_b = b_block.get(key)
            if iv_b is None:
                checks[check_key] = False
                continue
            try:
                merged_iv = intersect(iv_a, iv_b)
                checks[check_key] = True
                merged[check_key] = merged_iv
            except ValueError:
                checks[check_key] = False

    all_pass = all(checks.values())
    return {
        "all_pass": all_pass,
        "checks": checks,
        "primitives": merged,
    }
