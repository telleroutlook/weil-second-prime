"""
Authoritative eigenvalue reproduction from the TRUSTED certify assembly.

This does NOT use a re-implemented assembly. It calls the repo's own
scripts.certify_cross_influence.build_C_interval — the exact function whose min
LDL^T pivot the S4 certificate reports and the project trusts — for both the
full (real cross term) and cross_off variants, takes the outward-rounded interval
MIDPOINT matrix, and computes:

  * interval min-pivot via the repo's own min_pivot_signed  (must reproduce the
    S4 certified pivot ~ -0.0197 at L=3/5,N=8,even -> proves we are on the
    trusted assembly);
  * min EIGENVALUE via numpy on the symmetrized midpoint (the DOMINANCE report's
    judge) for full and off, and Delta-lambda = eig(full) - eig(off).

If eig_full here is NOT the report's +0.077 (and instead ~ -0.23), and if the
pivot reproduces -0.0197, then the report's Layer-1 eigenvalue numbers are
refuted on the project's OWN trusted assembly, independent of any reimplementation.

NOTE ON GRADE: the eigenvalue midpoint number is float-grade (numpy eigvalsh on
the interval center). A fully certified verified-eigenvalue (Kato-Temple +
Weyl-inflation over the interval radius) would only WIDEN bounds around this
center; it cannot move a -0.23 center to a +0.077 PD conclusion. The midpoint
already settles the sign dispute. The interval pivot IS certify-grade and anchors
trust in the assembly.
"""
from __future__ import annotations

import argparse
import json
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

from scripts.certify_cross_influence import build_C_interval, min_pivot_signed
from src.prime_layer.legendre_shift_2prime import window_check


def _mid_matrix(C):
    n = len(C)
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = 0.5 * (float(C[i][j][0]) + float(C[i][j][1]))
    return A


def _max_radius(C):
    n = len(C)
    r = 0.0
    for i in range(n):
        for j in range(n):
            r = max(r, 0.5 * (float(C[i][j][1]) - float(C[i][j][0])))
    return r


def check_point(L: Fraction, sector: str, N: int) -> dict:
    parity = 0 if sector == "even" else 1
    indices = list(range(parity, parity + 2 * N, 2))
    d = 2 * N + parity
    t0 = time.time()
    print(f"[auth {sector} N={N}] L={L}: build_C_interval FULL (trusted assembly)...", flush=True)
    Cf = build_C_interval(indices, L, d, include_cross=True)
    lo_f, hi_f, pos_f, note_f = min_pivot_signed(Cf)
    Af = _mid_matrix(Cf)
    ef = np.linalg.eigvalsh(0.5 * (Af + Af.T))
    radf = _max_radius(Cf)
    print(f"[auth {sector} N={N}] FULL: pivot in [{float(lo_f):.6e},{float(hi_f):.6e}] "
          f"(S4 ref -0.0197) | eig_min={ef[0]:+.6f} eig_2nd={ef[1]:+.6f} maxrad={radf:.2e} "
          f"({time.time()-t0:.0f}s)", flush=True)

    print(f"[auth {sector} N={N}] build_C_interval OFF...", flush=True)
    Co = build_C_interval(indices, L, d, include_cross=False)
    lo_o, hi_o, pos_o, note_o = min_pivot_signed(Co)
    Ao = _mid_matrix(Co)
    eo = np.linalg.eigvalsh(0.5 * (Ao + Ao.T))
    print(f"[auth {sector} N={N}] OFF: pivot in [{float(lo_o):.6e},{float(hi_o):.6e}] "
          f"| eig_min={eo[0]:+.6f}", flush=True)

    return {
        "L": str(L), "L_float": float(L), "sector": sector, "N": N, "d": d,
        "grade": "trusted-assembly midpoint (eig float, pivot certify)",
        "pivot_full": [float(lo_f), float(hi_f)], "pivot_full_note": note_f,
        "pivot_off": [float(lo_o), float(hi_o)], "pivot_off_note": note_o,
        "eig_full_min": float(ef[0]), "eig_full_2nd": float(ef[1]),
        "eig_off_min": float(eo[0]), "eig_off_2nd": float(eo[1]),
        "gap_full": float(ef[1] - ef[0]), "gap_off": float(eo[1] - eo[0]),
        "Delta_lambda_midpoint": float(ef[0] - eo[0]),
        "max_interval_radius_full": radf,
        "report_claim_eig_full": {"N6": 0.08341, "N8": 0.07689}.get(f"N{N}"),
        "report_claim_Delta_lambda": {"N6": 0.19526, "N8": 0.19521}.get(f"N{N}"),
        "elapsed_s": round(time.time() - t0, 1),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", default="3/5")
    ap.add_argument("--sector", default="even")
    ap.add_argument("--N", type=int, default=8)
    ap.add_argument("--out", default="pilots/authoritative_eig.json")
    args = ap.parse_args()
    num, den = (int(x) for x in args.L.split("/"))
    L = Fraction(num, den)
    if not window_check(float(L)):
        print(f"L={L} outside window"); return 2
    r = check_point(L, args.sector, args.N)
    out = Path(args.out)
    prev = json.loads(out.read_text()) if out.exists() else []
    prev.append(r)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(prev, indent=2))
    print(f"\n[auth] eig_full_min={r['eig_full_min']:+.6f} vs report {r['report_claim_eig_full']:+.6f}; "
          f"Delta_lambda={r['Delta_lambda_midpoint']:+.6f} vs report {r['report_claim_Delta_lambda']:+.6f}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
