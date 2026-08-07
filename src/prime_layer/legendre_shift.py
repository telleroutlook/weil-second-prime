"""Legendre shift algebra for the first prime layer.

Computes J_{ij}(tau) and E_{ij}(tau) exactly in Q[tau] using Legendre
recurrence and Python's Fraction type.  No quadrature, no sympy, no mpmath.

Key invariants (Theorem 4):
  J_{ij}(tau) = 2 * integral_{-1}^{1-tau} P_i(x) P_j(x+tau) dx
  E_{ij}(tau) = 2 * integral_{-1}^{1-tau} P_i(x) P_j(x)    dx
  Both vanish when i+j is odd.
  Both lie in Q[tau].

Sample values (used in unit tests):
  J_{00}(tau) = 4 - 2*tau
  J_{11}(tau) = tau^3/3 - 2*tau + 4/3
  J_{02}(tau) = -tau^3 + 3*tau^2 - 2*tau
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache


# ---------------------------------------------------------------------------
# Polynomial arithmetic over Q[x] (coefficients in ascending degree order)
# ---------------------------------------------------------------------------

Poly = list[Fraction]


def _zero(n: int = 1) -> Poly:
    return [Fraction(0)] * n


def _poly_add(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    result = _zero(n)
    for i, c in enumerate(a):
        result[i] += c
    for i, c in enumerate(b):
        result[i] += c
    return result


def _poly_scale(a: Poly, s: Fraction) -> Poly:
    return [c * s for c in a]


def _poly_mul(a: Poly, b: Poly) -> Poly:
    n = len(a) + len(b) - 1
    result = _zero(n)
    for i, ca in enumerate(a):
        for j, cb in enumerate(b):
            result[i + j] += ca * cb
    return result


def _poly_shift(poly: Poly, s: Fraction) -> Poly:
    """Return poly(x + s) as a polynomial in x."""
    # Horner + binomial: substitute x -> x + s coefficient by coefficient
    result: Poly = [Fraction(0)]
    for k in range(len(poly) - 1, -1, -1):
        # result = result * (x + s) + poly[k], but we build from highest degree
        # Alternative: direct composition via substitution
        pass

    # Direct substitution: poly(x+s) = sum_k c_k (x+s)^k
    result = _zero(1)
    x_plus_s_pow: Poly = [Fraction(1)]  # (x+s)^0 = 1
    x_monomial: Poly = [s, Fraction(1)]  # x + s
    for k, c in enumerate(poly):
        if k > 0:
            x_plus_s_pow = _poly_mul(x_plus_s_pow, x_monomial)
        if c != 0:
            result = _poly_add(result, _poly_scale(x_plus_s_pow, c))
    return result


def _poly_definite_integral(poly: Poly, lo: Fraction, hi: Fraction) -> Fraction:
    """Compute integral_{lo}^{hi} poly(x) dx exactly in Q."""
    total = Fraction(0)
    for k, c in enumerate(poly):
        total += c * (hi ** (k + 1) - lo ** (k + 1)) / (k + 1)
    return total


# ---------------------------------------------------------------------------
# Legendre polynomials
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def legendre_poly(n: int) -> tuple[Fraction, ...]:
    """Return P_n(x) coefficients in Q (ascending degree), cached."""
    if n == 0:
        return (Fraction(1),)
    if n == 1:
        return (Fraction(0), Fraction(1))
    p_prev = list(legendre_poly(n - 2))
    p_curr = list(legendre_poly(n - 1))
    # (n) P_n = (2n-1) x P_{n-1} - (n-1) P_{n-2}
    coeff_x = Fraction(2 * n - 1, n)
    coeff_prev = Fraction(n - 1, n)
    deg = n
    p_next: Poly = _zero(deg + 1)
    for i, c in enumerate(p_curr):
        if c != 0:
            p_next[i + 1] = _poly_scale([c], coeff_x)[0]
    # x * p_curr
    p_x_curr: Poly = _zero(len(p_curr) + 1)
    for i, c in enumerate(p_curr):
        p_x_curr[i + 1] += coeff_x * c
    # subtract coeff_prev * p_prev
    p_prev_scaled = _poly_scale(p_prev, coeff_prev)
    p_next = _poly_add(p_x_curr, [-c for c in p_prev_scaled])
    return tuple(p_next)


# ---------------------------------------------------------------------------
# J and E matrices
# ---------------------------------------------------------------------------

def compute_J(i: int, j: int, tau: Fraction) -> Fraction:
    """J_{ij}(tau) = <C_{tau,1} P_j, P_i> = 2 * int_{-1}^{1-tau} P_i(x) P_j(x+tau) dx.

    Returns 0 when i+j is odd (parity invariant, Theorem 4).
    """
    if (i + j) % 2 != 0:
        return Fraction(0)
    pi: Poly = list(legendre_poly(i))
    pj: Poly = list(legendre_poly(j))
    pj_shifted = _poly_shift(pj, tau)
    product = _poly_mul(pi, pj_shifted)
    return 2 * _poly_definite_integral(product, Fraction(-1), Fraction(1) - tau)


def compute_E(i: int, j: int, tau: Fraction) -> Fraction:
    """E_{ij}(tau) = <C_{tau,1} P_j, C_{tau,1} P_i> = 2 * int_{-1}^{1-tau} P_i(x) P_j(x) dx.

    Returns 0 when i+j is odd (parity invariant, Theorem 4).
    Uses C_{tau,1}^2 = 1_{E_- union E_+} from Theorem 1 exchange decomposition.
    """
    if (i + j) % 2 != 0:
        return Fraction(0)
    pi: Poly = list(legendre_poly(i))
    pj: Poly = list(legendre_poly(j))
    product = _poly_mul(pi, pj)
    return 2 * _poly_definite_integral(product, Fraction(-1), Fraction(1) - tau)


def prime_legendre_matrices(
    indices: list[int], tau: Fraction
) -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    """Return (J_matrix, E_matrix) for the given index set and tau.

    J[row][col] = J_{indices[row], indices[col]}(tau)
    E[row][col] = E_{indices[row], indices[col]}(tau)
    """
    N = len(indices)
    J_mat = [[compute_J(indices[r], indices[c], tau) for c in range(N)] for r in range(N)]
    E_mat = [[compute_E(indices[r], indices[c], tau) for c in range(N)] for r in range(N)]
    return J_mat, E_mat


# ---------------------------------------------------------------------------
# Verification of sample values (Theorem 4 cross-checks)
# ---------------------------------------------------------------------------

def _verify_sample_values(tau: Fraction) -> None:
    """Assert the three sample values from Theorem 4 hold exactly."""
    j00 = compute_J(0, 0, tau)
    expected_j00 = 4 - 2 * tau
    assert j00 == expected_j00, f"J_00 mismatch: {j00} != {expected_j00}"

    j11 = compute_J(1, 1, tau)
    expected_j11 = tau**3 / 3 - 2 * tau + Fraction(4, 3)
    assert j11 == expected_j11, f"J_11 mismatch: {j11} != {expected_j11}"

    j02 = compute_J(0, 2, tau)
    expected_j02 = -(tau**3) + 3 * tau**2 - 2 * tau
    assert j02 == expected_j02, f"J_02 mismatch: {j02} != {expected_j02}"
