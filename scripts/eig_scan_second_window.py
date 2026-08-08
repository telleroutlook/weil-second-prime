"""
Unified second-window eigenvalue L-scan serving BOTH outsourced verdicts.

At each rational L in the window ((1/2)log3, log2), do ONE archimedean build
(build_matrices: S0 four-term, M0, T, Gd) and derive THREE prime-layer variants,
sharing the identical archimedean block (the assembly's dominant cost). This makes
the three eigenvalue curves mutually consistent and cheap after one build:

  full : M2 = -(c2 J2 + c3 J3),  S2 = c2^2 E2 + c3^2 E3 + c2 c3 (F_ij+F_ji)
  off  : same M2,                 S2 = c2^2 E2 + c3^2 E3            (cross removed)
  p1   : M2 = -(c2 J2),           S2 = c2^2 E2                      (prime-2 ONLY)

Report 1 (CROSS_TERM_DOMINANCE Layer-2): Delta_lambda(L) = eig(full) - eig(off);
  inf over the window is the dominance margin. Report claims inf >= 0.158.
Report 2 (CROSS_TERM_MARGIN_COLLAPSE section 4): V2 = d eig(full)/dL and
  V1 = d eig(p1)/dL via symmetric difference; DeltaV = V2 - V1. p1 is EXACTLY the
  section-4 baseline C^(1) = b_L F0 - R_eta - c2 F2 (prime-2 only). The firewall
  holds: p1 is built without any prime-3 or cross data, so its slope cannot borrow
  from the dominance experiment.

min-swap guard: eig_2nd - eig_min is recorded per curve so a level crossing
(gap -> 0) that would break Hellmann-Feynman/finite-difference is visible.

Grade: float center (pilot) for the eig curves. A certify verified-eigenvalue
recheck at the anchor is done separately (authoritative_eig_check.py). Incremental-
durable: each L is written to disk as it completes (CLAUDE.md long-task rule),
--resume skips finished L.
"""
from __future__ import annotations

import argparse
import json
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

from checker.fp035.recompute_schur import build_matrices, _c_L, _H, KAPPA_FLOAT, L0
from src.prime_layer.legendre_shift_2prime import (
    C2, C3, compute_J, compute_E, compute_F, tau2_at, tau3_at, window_check,
)

ETA = 0.5


def _schur_eig(arch, M2, S2, d, L):
    Gd = arch["Gd"]; T = arch["T"]; M0 = arch["M0"]; S0 = arch["S0"]
    Ginv = np.diag([1.0 / g for g in Gd])
    R0 = S0 - M0.T @ Ginv @ M0
    R2 = S2 - M2.T @ Ginv @ M2
    R_eta = (1 + ETA) * R0 + (1 + 1.0 / ETA) * R2
    cL = _c_L(L)
    b_L = _H(d) - cL - L0 - KAPPA_FLOAT
    F = T + M0 + M2 - cL * np.diag(Gd)
    C = 0.5 * ((b_L * F - R_eta) + (b_L * F - R_eta).T)
    e = np.linalg.eigvalsh(C)
    return float(e[0]), float(e[1])


def scan_L(L_num: int, L_den: int, sector: str, N: int) -> dict:
    L = L_num / L_den
    parity = 0 if sector == "even" else 1
    d = 2 * N + parity
    arch = build_matrices(L_num, L_den, sector, N)
    indices = arch["indices"]
    tau2 = tau2_at(L); tau3 = tau3_at(L)
    n = len(indices)

    M2_full = np.zeros((n, n)); M2_p1 = np.zeros((n, n))
    S2_full = np.zeros((n, n)); S2_off = np.zeros((n, n)); S2_p1 = np.zeros((n, n))
    for a, i in enumerate(indices):
        for b, j in enumerate(indices):
            J2 = float(compute_J(i, j, tau2)); J3 = float(compute_J(i, j, tau3))
            E2 = float(compute_E(i, j, tau2)); E3 = float(compute_E(i, j, tau3))
            cross = C2 * C3 * (float(compute_F(i, j, tau2, tau3)) + float(compute_F(j, i, tau2, tau3)))
            M2_full[a, b] = -(C2 * J2 + C3 * J3)
            M2_p1[a, b] = -(C2 * J2)
            base = C2 * C2 * E2 + C3 * C3 * E3
            S2_full[a, b] = base + cross
            S2_off[a, b] = base
            S2_p1[a, b] = C2 * C2 * E2

    ef0, ef1 = _schur_eig(arch, M2_full, S2_full, d, L)
    eo0, eo1 = _schur_eig(arch, M2_full, S2_off, d, L)
    ep0, ep1 = _schur_eig(arch, M2_p1, S2_p1, d, L)
    return {
        "L": f"{L_num}/{L_den}", "L_float": L, "sector": sector, "N": N, "d": d,
        "eig_full_min": ef0, "eig_full_gap": ef1 - ef0,
        "eig_off_min": eo0, "eig_off_gap": eo1 - eo0,
        "eig_p1_min": ep0, "eig_p1_gap": ep1 - ep0,
        "Delta_lambda": ef0 - eo0,
        "grade": "pilot (float center)",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--Ls", required=True, help="comma-separated rationals a/b, e.g. 5500/10000,5600/10000")
    ap.add_argument("--sector", default="even")
    ap.add_argument("--N", type=int, default=8)
    ap.add_argument("--out", default="pilots/eig_scan.json")
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()

    out = Path(args.out)
    results = json.loads(out.read_text()) if (args.resume and out.exists()) else []
    done = {(r["L"], r["sector"], r["N"]) for r in results}

    for tok in args.Ls.split(","):
        num, den = (int(x) for x in tok.split("/"))
        L = num / den
        if not window_check(L):
            print(f"[skip] L={tok}={L:.6f} outside window", flush=True)
            continue
        if (f"{num}/{den}", args.sector, args.N) in done:
            print(f"[skip] {tok} {args.sector} N={args.N} done", flush=True)
            continue
        t0 = time.time()
        print(f"[scan {args.sector} N={args.N}] L={tok}={L:.6f} building...", flush=True)
        r = scan_L(num, den, args.sector, args.N)
        r["elapsed_s"] = round(time.time() - t0, 1)
        results.append(r)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(results, indent=2))
        print(f"[scan {args.sector} N={args.N}] L={L:.6f} "
              f"eig_full={r['eig_full_min']:+.6f}(gap {r['eig_full_gap']:.3f}) "
              f"eig_off={r['eig_off_min']:+.6f} eig_p1={r['eig_p1_min']:+.6f} "
              f"Dlam={r['Delta_lambda']:+.6f} ({r['elapsed_s']}s)", flush=True)
    print(f"\nWritten to {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
