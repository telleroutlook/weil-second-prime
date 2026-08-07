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

DISCIPLINE (C11): the cross term F is genuinely new to the two-prime window and
is NOT yet implemented (S3). This module DOES NOT silently set F=0. When c_3 != 0
and no F provider is supplied, S2_two_prime RAISES. The single-prime limit
(c_3 = 0) needs no F and is used by the S2 acceptance self-check.

Reuses the exact Fraction/Legendre arithmetic from legendre_shift.py.
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Callable, List, Optional

from .legendre_shift import compute_J, compute_E

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

    The cross-prime term F is NEW to the two-prime window. This function does
    NOT fabricate F=0 (that is the C11 omitted-cross-term bug). When c3 != 0 a
    real F_provider MUST be supplied, else NotImplementedError is raised. The
    single-prime limit (c3 = 0) is exact without F and drives the S2 self-check.
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
