"""
S2 acceptance gate — single-prime limit self-check.

The two-prime prime layer (legendre_shift_2prime) must reduce EXACTLY to
weil-first's single-prime layer when prime 3 is switched off (c3 = 0). This is
the S2 acceptance criterion from HANDOFF/PLAN: catch a bad port or a formula
error (sign, transposition, wrong tau, spurious factor) BEFORE any two-prime
computation is trusted (do not carry a mistake into S3+).

Two levels of check:

  1. prime-layer element-wise (FAST): M2_two_prime(c3=0) and S2_two_prime(c3=0)
     vs the weil-first single-prime formula recomputed inline here. No
     archimedean integrals -> milliseconds. Runs in pytest.

  2. assembled-Schur element-wise (SLOW): build the full C = b_L F - R_eta once
     via the ported (byte-identical) recompute_schur.build_matrices, then swap
     only the prime layer (M2, S2) for the two-prime c3=0 version and re-assemble.
     max|C_first - C_second| must be < 1e-10. Run via run_and_wait (~1 min/build).

The archimedean part (S0 four-term, M0, T, Gd) is identical in both paths, so a
difference in C can only come from the prime layer — exactly what we isolate.
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction
from typing import List

import numpy as np

from checker.fp035.recompute_schur import (
    build_matrices,
    _c_L,
    _H,
    C2_FLOAT,
    KAPPA_FLOAT,
    L0,
    ETA,
)
from src.prime_layer.legendre_shift import compute_J, compute_E
from src.prime_layer.legendre_shift_2prime import M2_two_prime, S2_two_prime

TOL = 1e-10


def _max_abs_diff(A: List[List[float]] | np.ndarray, B: List[List[float]] | np.ndarray) -> float:
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    return float(np.max(np.abs(A - B)))


def single_prime_layer(indices: List[int], L: float) -> tuple[np.ndarray, np.ndarray]:
    """weil-first single-prime layer, recomputed inline (ground truth).

    Matches recompute_schur.build_matrices exactly:
      M2[a,b] = -C2 * J_{ij}(tau),  S2[a,b] = C2^2 * E_{ij}(tau),  tau = log2/L.
    """
    tau = Fraction(math.log(2) / L).limit_denominator(10000)
    n = len(indices)
    M2 = np.zeros((n, n))
    S2 = np.zeros((n, n))
    for a, i in enumerate(indices):
        for b, j in enumerate(indices):
            M2[a, b] = -C2_FLOAT * float(compute_J(i, j, tau))
            S2[a, b] = (C2_FLOAT ** 2) * float(compute_E(i, j, tau))
    return M2, S2


def check_prime_layer(L: float, sector: str, N: int) -> dict:
    """Level 1: prime-layer element-wise, single-prime limit (c3=0). Fast."""
    parity = 0 if sector == "even" else 1
    indices = list(range(parity, parity + 2 * N, 2))
    M2_first, S2_first = single_prime_layer(indices, L)
    # two-prime with prime 3 switched off (c3 = 0) -> no cross term needed
    M2_second = np.array(M2_two_prime(indices, L, c3=0.0))
    S2_second = np.array(S2_two_prime(indices, L, c3=0.0))
    dM = _max_abs_diff(M2_first, M2_second)
    dS = _max_abs_diff(S2_first, S2_second)
    return {
        "L": L, "sector": sector, "N": N,
        "max_dM2": dM, "max_dS2": dS,
        "pass": dM < TOL and dS < TOL,
    }


def _reassemble_C(mats: dict, M2: np.ndarray, S2: np.ndarray, d: int, eta: float) -> np.ndarray:
    """Re-run pivot_from_matrices' assembly with a swapped prime layer."""
    Gd = mats["Gd"]; T = mats["T"]; M0 = mats["M0"]; S0 = mats["S0"]; L = mats["L"]
    Ginv = np.diag([1.0 / g for g in Gd])
    R0 = S0 - M0.T @ Ginv @ M0
    R2 = S2 - M2.T @ Ginv @ M2
    R_eta = (1 + eta) * R0 + (1 + 1.0 / eta) * R2
    c_L = _c_L(L)
    b_L = _H(d) - c_L - L0 - KAPPA_FLOAT
    F = T + M0 + M2 - c_L * np.diag(Gd)
    return b_L * F - R_eta


def check_assembled_C(L: float, sector: str, N: int, d: int, eta: float = ETA) -> dict:
    """Level 2: assembled-Schur element-wise. Slow (one archimedean build)."""
    parity = 0 if sector == "even" else 1
    indices = list(range(parity, parity + 2 * N, 2))
    mats = build_matrices_for_L(L, sector, N)
    # Path A: weil-first single-prime layer (as build_matrices produced it)
    C_first = _reassemble_C(mats, mats["M2"], mats["S2"], d, eta)
    # Path B: two-prime layer with c3 = 0
    M2_second = np.array(M2_two_prime(indices, L, c3=0.0))
    S2_second = np.array(S2_two_prime(indices, L, c3=0.0))
    C_second = _reassemble_C(mats, M2_second, S2_second, d, eta)
    dC = _max_abs_diff(C_first, C_second)
    return {
        "L": L, "sector": sector, "N": N, "d": d,
        "max_dM2": _max_abs_diff(mats["M2"], M2_second),
        "max_dS2": _max_abs_diff(mats["S2"], S2_second),
        "max_dC": dC,
        "pass": dC < TOL,
    }


def build_matrices_for_L(L: float, sector: str, N: int) -> dict:
    """build_matrices takes L as a rational L_num/L_den; derive from float L."""
    f = Fraction(L).limit_denominator(10000)
    return build_matrices(f.numerator, f.denominator, sector, N)


def main() -> int:
    ap = argparse.ArgumentParser(description="S2 single-prime-limit self-check")
    ap.add_argument("--L", type=float, default=0.6, help="test point in second window")
    ap.add_argument("--N", type=int, default=3)
    ap.add_argument("--sector", choices=["even", "odd", "both"], default="both")
    ap.add_argument("--full", action="store_true",
                    help="also run the slow assembled-C check (archimedean build)")
    ap.add_argument("--d", type=int, default=1)
    args = ap.parse_args()

    sectors = ["even", "odd"] if args.sector == "both" else [args.sector]
    all_pass = True
    for sector in sectors:
        r = check_prime_layer(args.L, sector, args.N)
        status = "PASS" if r["pass"] else "FAIL"
        print(f"[prime-layer] L={r['L']} {sector} N={r['N']}: "
              f"max|dM2|={r['max_dM2']:.2e} max|dS2|={r['max_dS2']:.2e}  {status}",
              flush=True)
        all_pass &= r["pass"]

    if args.full:
        for sector in sectors:
            r = check_assembled_C(args.L, sector, args.N, args.d)
            status = "PASS" if r["pass"] else "FAIL"
            print(f"[assembled-C] L={r['L']} {sector} N={r['N']} d={args.d}: "
                  f"max|dM2|={r['max_dM2']:.2e} max|dS2|={r['max_dS2']:.2e} "
                  f"max|dC|={r['max_dC']:.2e}  {status}", flush=True)
            all_pass &= r["pass"]

    print(f"\nS2 acceptance gate: {'PASS' if all_pass else 'FAIL'} (tol={TOL:.0e})",
          flush=True)
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
