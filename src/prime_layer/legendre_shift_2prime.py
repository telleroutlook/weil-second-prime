"""
legendre_shift_2prime.py — Legendre matrix algebra for two prime shifts.

Second-prime window L in (1/2 log 3, log 2): both n=2 (tau_2 = log2/L) and
n=3 (tau_3 = log3/L) are in the single-hop regime, so the prime layer is a
genuine TWO-SHIFT object.

The prime operator is  C = c_2 C_{tau_2,1} + c_3 C_{tau_3,1}  with
  c_p = Lambda(p)/sqrt(p):   c_2 = log2/sqrt2,   c_3 = log3/sqrt3.

Expanding the prime second moment gives (every term must be present — omitting
any one is the C11 false-positive bug class, see CLAUDE.md):

  M^{(2)}_{ij} = <(V+K)P_j, C P_i>-analogue
              = -( c_2 J_{ij}(tau_2) + c_3 J_{ij}(tau_3) )                (linear, no cross term)

  S^{(2)}_{ij} = <C P_j, C P_i>
              = c_2^2 E_{ij}(tau_2)
              + c_3^2 E_{ij}(tau_3)
              + c_2 c_3 ( F_{ij}(tau_2,tau_3) + F_{ji}(tau_2,tau_3) )      (cross-prime)

where
  J_{ij}(tau) = <C_{tau,1} P_j, P_i>              (single-shift, in legendre_shift.py)
  E_{ij}(tau) = <C_{tau,1} P_j, C_{tau,1} P_i>    (single-shift, in legendre_shift.py)
  F_{ij}(tau_2,tau_3) = <C_{tau_3,1} P_j, C_{tau_2,1} P_i>   (cross-prime, NEW — S3).

DISCIPLINE (C11): the cross term F is genuinely new to the two-prime window.
It is now implemented exactly (compute_F / default_F_provider). This module still
NEVER silently sets F=0: S2_two_prime requires an explicit F_provider when
c_3 != 0 (pass default_F_provider for the real term), else it RAISES. The
single-prime limit (c_3 = 0) needs no F and is used by the S2 acceptance
self-check.

Reuses the exact Fraction/Legendre arithmetic from legendre_shift.py.
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Callable, List, Optional

from .legendre_shift import (
    compute_J,
    compute_E,
    legendre_poly,
    _poly_mul,
    _poly_shift,
    _poly_definite_integral,
)

# ── Prime constants ─────────────────────────────────────────────────────────

LOG2 = math.log(2)
LOG3 = math.log(3)
C2 = LOG2 / math.sqrt(2)   # Lambda(2)/sqrt(2)
C3 = LOG3 / math.sqrt(3)   # Lambda(3)/sqrt(3)

# Rational-approximation denominator cap; matches recompute_schur convention.
_TAU_DEN = 10_000

# A cross-prime F provider maps (i, j, tau_2, tau_3) -> F_{ij}(tau_2, tau_3).
FProvider = Callable[[int, int, Fraction, Fraction], float]


def tau_at(L: float, prime: int) -> Fraction:
    """tau_p = log(p)/L as a rational approximation (denominator <= 10000)."""
    return Fraction(math.log(prime) / L).limit_denominator(_TAU_DEN)


def tau2_at(L: float) -> Fraction:
    return tau_at(L, 2)


def tau3_at(L: float) -> Fraction:
    return tau_at(L, 3)


def window_check(L: float) -> bool:
    """True iff L lies in the second-prime window (1/2 log3, log2)."""
    return LOG3 / 2 < L < LOG2


# ── Cross-prime term F (S3, the mathematical heart of the two-prime window) ──
#
# The scaled exchange operator (Theorem 1, b -> tau, L -> 1) acts on
# f in L^2(-1,1) as
#
#     (C_{tau,1} f)(x) = 1_{(-1, 1-tau)}(x) f(x+tau)   [forward strip]
#                      + 1_{(tau-1, 1)}(x) f(x-tau)     [backward strip].
#
# The cross-prime Gram is
#
#     F_{ij}(tau_2, tau_3) = <C_{tau_3,1} P_j, C_{tau_2,1} P_i>
#                          = integral_{-1}^{1} (C_{tau_2,1} P_i)(x) (C_{tau_3,1} P_j)(x) dx.
#
# Each operator contributes a forward + a backward strip, so the product is a
# sum of FOUR shift-combinations, each integrated over the intersection of the
# two strips involved:
#
#   A  P_i(x+tau_2) P_j(x+tau_3)  on (-1,           1 - max(tau_2,tau_3))
#   B  P_i(x+tau_2) P_j(x-tau_3)  on (tau_3 - 1,     1 - tau_2)
#   C  P_i(x-tau_2) P_j(x+tau_3)  on (tau_2 - 1,     1 - tau_3)
#   D  P_i(x-tau_2) P_j(x-tau_3)  on (max(tau_2,tau_3) - 1, 1)
#
# Every term is in Q[x] (Legendre shift by a rational tau), so F in Q[tau2,tau3].
# Ground-truth invariants (unit-tested): F_{ij}(tau,tau) == E_{ij}(tau) (same
# operator both sides), and F_{ij} == 0 for i+j odd (reflection x -> -x gives
# F = (-1)^{i+j} F). NO term is dropped — omitting one is the C11 bug class.


def _strip_integral(
    pi_shift: Fraction, pj_shift: Fraction, i: int, j: int, lo: Fraction, hi: Fraction
) -> Fraction:
    """integral_{lo}^{hi} P_i(x + pi_shift) P_j(x + pj_shift) dx, exact in Q.

    Returns 0 if the domain is empty (lo >= hi)."""
    if lo >= hi:
        return Fraction(0)
    pi = _poly_shift(list(legendre_poly(i)), pi_shift)
    pj = _poly_shift(list(legendre_poly(j)), pj_shift)
    product = _poly_mul(pi, pj)
    return _poly_definite_integral(product, lo, hi)


def compute_F(i: int, j: int, tau2: Fraction, tau3: Fraction) -> Fraction:
    """Cross-prime Gram F_{ij}(tau_2, tau_3) = <C_{tau_3,1} P_j, C_{tau_2,1} P_i>.

    Exact in Q[tau_2, tau_3]. Sums the four shift-combinations A,B,C,D over their
    (possibly empty) intersection domains. Returns 0 when i+j is odd (parity).
    """
    if (i + j) % 2 != 0:
        return Fraction(0)
    one = Fraction(1)
    tmax = max(tau2, tau3)
    # A: P_i(x+tau2) P_j(x+tau3) on (-1, 1 - max(tau2,tau3))
    A = _strip_integral(tau2, tau3, i, j, -one, one - tmax)
    # B: P_i(x+tau2) P_j(x-tau3) on (tau3 - 1, 1 - tau2)
    B = _strip_integral(tau2, -tau3, i, j, tau3 - one, one - tau2)
    # C: P_i(x-tau2) P_j(x+tau3) on (tau2 - 1, 1 - tau3)
    C = _strip_integral(-tau2, tau3, i, j, tau2 - one, one - tau3)
    # D: P_i(x-tau2) P_j(x-tau3) on (max(tau2,tau3) - 1, 1)
    D = _strip_integral(-tau2, -tau3, i, j, tmax - one, one)
    return A + B + C + D


def default_F_provider(i: int, j: int, tau2: Fraction, tau3: Fraction) -> float:
    """Float wrapper around compute_F, matching the FProvider signature."""
    return float(compute_F(i, j, tau2, tau3))


def M2_two_prime(
    indices: List[int],
    L: float,
    c2: float = C2,
    c3: float = C3,
) -> List[List[float]]:
    """Prime-layer first-moment matrix M^{(2)} for two primes.

        M^{(2)}_{ij} = -( c2 * J_{ij}(tau_2) + c3 * J_{ij}(tau_3) )

    Both shifts are always present. Setting c3 = 0 gives the single-prime limit
    (reduces to weil-first's -c2 * J_{ij}(tau_2)).
    """
    tau2 = tau2_at(L)
    tau3 = tau3_at(L)
    n = len(indices)
    M = [[0.0] * n for _ in range(n)]
    for a, i in enumerate(indices):
        for b, j in enumerate(indices):
            m = -c2 * float(compute_J(i, j, tau2))
            if c3 != 0.0:
                m += -c3 * float(compute_J(i, j, tau3))
            M[a][b] = m
    return M


def S2_two_prime(
    indices: List[int],
    L: float,
    c2: float = C2,
    c3: float = C3,
    F_provider: Optional[FProvider] = None,
) -> List[List[float]]:
    """Prime-layer second-moment matrix S^{(2)} for two primes.

        S^{(2)}_{ij} = c2^2 E_{ij}(tau_2) + c3^2 E_{ij}(tau_3)
                     + c2 c3 ( F_{ij} + F_{ji} )

    The cross-prime term F is NEW to the two-prime window and is implemented in
    compute_F. This function does NOT fabricate F=0 (the C11 omitted-cross-term
    bug). When c3 != 0 a real F_provider MUST be supplied (pass
    default_F_provider), else NotImplementedError is raised. The single-prime
    limit (c3 = 0) is exact without F and drives the S2 self-check.
    """
    tau2 = tau2_at(L)
    tau3 = tau3_at(L)
    n = len(indices)
    S = [[0.0] * n for _ in range(n)]
    for a, i in enumerate(indices):
        for b, j in enumerate(indices):
            s = c2 * c2 * float(compute_E(i, j, tau2))
            if c3 != 0.0:
                s += c3 * c3 * float(compute_E(i, j, tau3))
                if F_provider is None:
                    raise NotImplementedError(
                        "cross-prime F_{ij}(tau_2, tau_3) is required for c3 != 0 "
                        "(S3 deliverable); refusing to silently omit it (C11 bug class)."
                    )
                F_ij = F_provider(i, j, tau2, tau3)
                F_ji = F_provider(j, i, tau2, tau3)
                s += c2 * c3 * (F_ij + F_ji)
            S[a][b] = s
    return S
