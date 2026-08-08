"""Two-prime split-residual Schur recompute for the second window.

Reuses the trusted archimedean assembly of weil-first (checker.fp035.recompute_schur):
S0 = S_VV + S_VK + S_KV + S_KK (four terms), M0 = V + K, T, Gd, and the generic
R0/R2/F/C assembly in pivot_from_matrices. The ONLY thing that changes for the
second window is the prime layer, which becomes a genuine TWO-shift object:

    M^{(2)}_{ij} = -( c2 J_{ij}(tau2) + c3 J_{ij}(tau3) )
    S^{(2)}_{ij} = c2^2 E_{ij}(tau2) + c3^2 E_{ij}(tau3)
                 + c2 c3 ( F_{ij} + F_{ji} )     <- cross-prime term (S3)

with tau_p = log(p)/L, c_p = log(p)/sqrt(p).

Every prime shift and the cross term are always present (omitting any is the C11
false-positive bug class, CLAUDE.md). The profiling knobs (include_tau2 /
include_tau3 / include_cross) exist ONLY for the S4 influence probe and default
to all-on; a certificate run must use all-on.

NOTE: build_matrices below is float-CENTER (pilot grade). Verdicts require the
Arb-interval certify path (S5). This module is the pilot/profile engine.
"""
from __future__ import annotations

import math
from fractions import Fraction
from typing import Optional

import numpy as np

from checker.fp035.recompute_schur import build_matrices, pivot_from_matrices, _c_L
from src.prime_layer.legendre_shift_2prime import (
    C2,
    C3,
    M2_two_prime,
    S2_two_prime,
    default_F_provider,
    window_check,
)


def build_two_prime_matrices(
    L_num: int,
    L_den: int,
    sector: str,
    N: int,
    include_tau2: bool = True,
    include_tau3: bool = True,
    include_cross: bool = True,
    swap_c2_c3: bool = False,
) -> dict:
    """Assemble the two-prime matrix dict.

    Starts from the trusted archimedean build (S0 four-term, M0, T, Gd) and
    replaces the prime layer (M2, S2) with the two-shift version.

    Profiling knobs (S4 only; a certificate uses the defaults = all True):
      include_tau2 / include_tau3 : keep/zero the c2 (prime 2) / c3 (prime 3) shift.
      include_cross               : keep/zero the cross-prime F term in S2.
      swap_c2_c3                  : swap c2<->c3 (a mutation-style probe; a
                                    correct checker must reject a swapped cert).
    """
    L = L_num / L_den
    if not window_check(L):
        raise ValueError(f"L={L} outside second-prime window (1/2 log3, log2)")

    mats = build_matrices(L_num, L_den, sector, N)  # archimedean + throwaway 1-prime layer
    indices = mats["indices"]

    c2, c3 = (C3, C2) if swap_c2_c3 else (C2, C3)
    c2_eff = c2 if include_tau2 else 0.0
    c3_eff = c3 if include_tau3 else 0.0

    M2 = np.array(M2_two_prime(indices, L, c2=c2_eff, c3=c3_eff))

    if c3_eff != 0.0 and include_cross:
        F_provider = default_F_provider
    elif c3_eff != 0.0 and not include_cross:
        # cross term dropped ON PURPOSE (S4 probe). Not a silent omission: it is
        # an explicit, labelled experiment. A certificate never takes this path.
        F_provider = lambda i, j, t2, t3: 0.0  # noqa: E731
    else:
        F_provider = None  # c3 == 0: no cross term exists

    S2 = np.array(S2_two_prime(indices, L, c2=c2_eff, c3=c3_eff, F_provider=F_provider))

    mats["M2"] = M2
    mats["S2"] = S2
    mats["c2"] = c2_eff
    mats["c3"] = c3_eff
    mats["profile_flags"] = {
        "include_tau2": include_tau2,
        "include_tau3": include_tau3,
        "include_cross": include_cross,
        "swap_c2_c3": swap_c2_c3,
    }
    return mats


def verify_sector_2prime(
    L_num: int,
    L_den: int,
    sector: str,
    N: int,
    d: int,
    eta: float = 0.5,
    **flags,
) -> tuple[float, float, dict]:
    """Compute (min_pivot, b_L, info) for the two-prime second window (pilot grade)."""
    mats = build_two_prime_matrices(L_num, L_den, sector, N, **flags)
    piv, b_L = pivot_from_matrices(mats, d, eta, judge="pivot")
    info = {
        "L": mats["L"], "sector": sector, "N": N, "d": d, "eta": eta,
        "b_L": b_L, "c_L": _c_L(mats["L"]), "min_pivot": piv,
        "c2": mats["c2"], "c3": mats["c3"],
        "S0_definition": "S_VV+S_VK+S_KV+S_KK",
        "S2_definition": "c2^2 E(tau2) + c3^2 E(tau3) + c2 c3 (F_ij + F_ji)",
        "judge": "min LDL^T pivot",
        "profile_flags": mats["profile_flags"],
        "grade": "pilot (float center)",
    }
    return piv, b_L, info
