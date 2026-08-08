"""
Independent reproduction harness for CROSS_TERM_DOMINANCE_PROOF Layer-1/Layer-2.

Supervisor task (2026-08-08): the DOMINANCE report's Layer-2 (inf Delta-lambda >=
0.158 over 142 subintervals, N=6,8) is "method rigorous, self-reported, awaiting
independent reproduction". This harness recomputes the report's OWN judge
(min EIGENVALUE, not min-pivot) from the repo's trusted assembly, independently.

Efficiency: the archimedean block (S0,M0,T,Gd) depends only on (L,sector,N), NOT
on the prime flags. full and off differ ONLY in S2's cross term
    S2_full - S2_off = c2 c3 (F_ij + F_ji).
So we build the archimedean block ONCE per (L,sector,N) and derive both variants,
halving the cost versus a naive full+off rebuild.

Two judges reported per point:
  * min_pivot  : the repo's canonical positivity judge (LDL^T). Used ONLY to
                 CROSS-VALIDATE this assembly against the repo's certified S4
                 number (full pivot ~ -0.0197 at L=3/5,N=8,even). If this matches,
                 the assembly is trustworthy and the eigenvalue numbers are too.
  * min_eig    : the DOMINANCE report's judge. Delta-lambda = eig(full)-eig(off).

Grade here = float center (pilot). A certify-grade (Arb verified-eigenvalue)
confirmation of the anchor point is a separate step; this establishes WHICH set of
Layer-1 numbers is correct before spending the expensive Arb build.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from checker.fp035.recompute_schur import build_matrices, _c_L, _H, KAPPA_FLOAT, L0
from src.prime_layer.legendre_shift_2prime import (
    C2, C3, compute_J, compute_E, compute_F, tau2_at, tau3_at, window_check,
)

ETA = 0.5


def _min_pivot(C: np.ndarray) -> float:
    """Min LDL^T pivot of the symmetrized matrix (repo's positivity judge)."""
    A = 0.5 * (C + C.T)
    n = A.shape[0]
    Lm = np.eye(n)
    d = np.zeros(n)
    for j in range(n):
        s = A[j, j] - sum(Lm[j, k] * Lm[j, k] * d[k] for k in range(j))
        d[j] = s
        if abs(s) < 1e-300:
            return float(s)
        for i in range(j + 1, n):
            Lm[i, j] = (A[i, j] - sum(Lm[i, k] * Lm[j, k] * d[k] for k in range(j))) / s
    return float(np.min(d))


def build_prime_layers(indices, L, swap=False):
    """M2 (identical for full/off) and the cross-term matrix X = c2 c3 (F_ij+F_ji).

    S2_full = c2^2 E(tau2) + c3^2 E(tau3) + X ;  S2_off = S2_full - X.
    """
    c2, c3 = (C3, C2) if swap else (C2, C3)
    tau2 = tau2_at(L)
    tau3 = tau3_at(L)
    n = len(indices)
    M2 = np.zeros((n, n))
    S2_full = np.zeros((n, n))
    X = np.zeros((n, n))
    for a, i in enumerate(indices):
        for b, j in enumerate(indices):
            J2 = float(compute_J(i, j, tau2)); J3 = float(compute_J(i, j, tau3))
            E2 = float(compute_E(i, j, tau2)); E3 = float(compute_E(i, j, tau3))
            M2[a, b] = -(c2 * J2 + c3 * J3)
            xij = c2 * c3 * (float(compute_F(i, j, tau2, tau3)) + float(compute_F(j, i, tau2, tau3)))
            X[a, b] = xij
            S2_full[a, b] = c2 * c2 * E2 + c3 * c3 * E3 + xij
    return M2, S2_full, X


def schur_C(arch, M2, S2, d, L, eta=ETA):
    Gd = arch["Gd"]; T = arch["T"]; M0 = arch["M0"]; S0 = arch["S0"]
    Ginv = np.diag([1.0 / g for g in Gd])
    R0 = S0 - M0.T @ Ginv @ M0
    R2 = S2 - M2.T @ Ginv @ M2
    R_eta = (1 + eta) * R0 + (1 + 1.0 / eta) * R2
    cL = _c_L(L)
    b_L = _H(d) - cL - L0 - KAPPA_FLOAT
    F = T + M0 + M2 - cL * np.diag(Gd)
    C = b_L * F - R_eta
    return 0.5 * (C + C.T), b_L


def point_report(L_num, L_den, sector, N, swap=False):
    L = L_num / L_den
    parity = 0 if sector == "even" else 1
    d = 2 * N + parity
    arch = build_matrices(L_num, L_den, sector, N)
    indices = arch["indices"]
    M2, S2_full, X = build_prime_layers(indices, L, swap=swap)
    S2_off = S2_full - X

    Cf, b_L = schur_C(arch, M2, S2_full, d, L)
    Co, _ = schur_C(arch, M2, S2_off, d, L)

    ef = np.linalg.eigvalsh(Cf)
    eo = np.linalg.eigvalsh(Co)
    pf = _min_pivot(Cf)
    po = _min_pivot(Co)
    return {
        "L": f"{L_num}/{L_den}", "L_float": L, "sector": sector, "N": N, "d": d,
        "b_L": b_L, "swap_c2_c3": swap,
        "eig_full_min": float(ef[0]), "eig_full_2nd": float(ef[1]),
        "eig_off_min": float(eo[0]), "eig_off_2nd": float(eo[1]),
        "gap_full": float(ef[1] - ef[0]), "gap_off": float(eo[1] - eo[0]),
        "Delta_lambda": float(ef[0] - eo[0]),
        "pivot_full": pf, "pivot_off": po,
        "grade": "pilot (float center)",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", default="3/5")
    ap.add_argument("--sectors", default="even")
    ap.add_argument("--Ns", default="6,8")
    ap.add_argument("--out", default="pilots/independent_dominance.json")
    args = ap.parse_args()
    num, den = (int(x) for x in args.L.split("/"))
    L = num / den
    if not window_check(L):
        print(f"L={L} outside window"); return 2
    out = Path(args.out)
    results = []
    for sector in args.sectors.split(","):
        for N in (int(x) for x in args.Ns.split(",")):
            t0 = time.time()
            print(f"[indep {sector} N={N}] L={num}/{den} building...", flush=True)
            r = point_report(num, den, sector, N)
            r["elapsed_s"] = round(time.time() - t0, 1)
            results.append(r)
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps(results, indent=2))
            print(f"[indep {sector} N={N}] eig_full={r['eig_full_min']:+.6f} "
                  f"eig_off={r['eig_off_min']:+.6f} Dlam={r['Delta_lambda']:+.6f} "
                  f"| pivot_full={r['pivot_full']:+.6f} (repo S4 ref -0.0197) "
                  f"({r['elapsed_s']}s)", flush=True)
    print(f"\nWritten to {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
