"""
S5 — certify-grade positivity scan of the second window, from the LEFT END.

Supervisor steer (2026-08-08): the window is ((1/2)log3, log2) ~ (0.5493, 0.6931).
Positivity (if it survives anywhere) is most likely at the LEFT END L -> (1/2)log3+,
where c_L is smallest and prime 3 just entered the single-hop regime (tau3 -> 2-).
first-window experience: positivity lived only near the left edge and died as L
grew. So scan from the left, not from the middle (L=0.6 is already certified
NOT positive-definite, see docs/S4).

VERDICT GRADE: certify only. This script builds C = b_L F - R_eta as
outward-rounded Arb intervals and computes the interval min LDL^T pivot. A point
is declared positive-definite ONLY if every pivot's lower endpoint is > 0
(is_strictly_positive on the whole diagonal, no straddle). A pilot float number
NEVER decides positivity here (PROOF_CONSTITUTION A3).

Honest expected outcome (supervisor): the second window's certifiable range may
be SHORTER than the first window's (larger c_L + prime-3 negative contribution).
If even the left end is not positive-definite, that is a real E2-type boundary
result — "the method's reach in the second window is shorter than the first" —
not a failure.

Reuses the certified interval assembly of certify_cross_influence (full variant,
real cross term F). Long run: ~8 min per (sector) interval build. Incremental
per-point JSON (CLAUDE.md long-task rule).
"""
from __future__ import annotations

import argparse
import json
import math
import time
from fractions import Fraction
from pathlib import Path

from src.prime_layer.legendre_shift_2prime import window_check
from scripts.certify_cross_influence import build_C_interval, min_pivot_signed


def certify_positivity(L: Fraction, sector: str, N: int) -> dict:
    parity = 0 if sector == "even" else 1
    indices = list(range(parity, parity + 2 * N, 2))
    d = 2 * N + parity
    t0 = time.time()
    tau3 = math.log(3) / float(L)
    print(f"[S5 {sector}] L={L}={float(L):.6f} N={N} d={d} tau3={tau3:.5f} "
          f"(single-hop needs tau3<2): building interval C (full, real cross)...",
          flush=True)
    C = build_C_interval(indices, L, d, include_cross=True)
    lo, hi, all_pos, note = min_pivot_signed(C)
    verdict = "POSITIVE-DEFINITE" if all_pos else "NOT positive-definite"
    print(f"[S5 {sector}] min-pivot in [{float(lo):.6e}, {float(hi):.6e}] "
          f"all_pivots_pos={all_pos} note={note} -> {verdict} ({time.time()-t0:.0f}s)",
          flush=True)
    return {
        "L": str(L), "L_float": float(L), "sector": sector, "N": N, "d": d,
        "grade": "certify (Arb interval)",
        "tau3": tau3, "single_hop_ok": tau3 < 2,
        "min_pivot_lower": float(lo), "min_pivot_upper": float(hi),
        "min_pivot_note": note,
        "positive_definite": all_pos,
        "elapsed_s": round(time.time() - t0, 1),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="S5 certify-grade positivity scan (left-end first)")
    ap.add_argument("--points", default="11/20",
                    help="comma-separated rational L values, left-end first, e.g. 11/20,14/25")
    ap.add_argument("--sectors", default="even,odd")
    ap.add_argument("--even-N", type=int, default=8)
    ap.add_argument("--odd-N", type=int, default=7)
    ap.add_argument("--out", default="pilots/s5_positivity_scan.json")
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()

    out = Path(args.out)
    results = json.loads(out.read_text()) if (args.resume and out.exists()) else []
    done = {(r["L"], r["sector"]) for r in results}

    for pt in args.points.split(","):
        num, den = pt.split("/")
        L = Fraction(int(num), int(den))
        if not window_check(float(L)):
            print(f"[skip] L={L} outside window ((1/2)log3, log2)", flush=True)
            continue
        for sector in args.sectors.split(","):
            if (str(L), sector) in done:
                print(f"[skip] {sector}@{L} already done", flush=True)
                continue
            N = args.even_N if sector == "even" else args.odd_N
            r = certify_positivity(L, sector, N)
            results.append(r)
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps(results, indent=2))
            print(f"[saved] {sector}@{L} -> {args.out}", flush=True)

    print(f"\nS5 positivity scan written to {args.out}", flush=True)
    any_pd = False
    for r in results:
        print(f"  L={r['L']} {r['sector']}: min_pivot_lower={r['min_pivot_lower']:.4e} "
              f"PD={r['positive_definite']}", flush=True)
        any_pd |= r["positive_definite"]
    print(f"\nAny sector/point certified positive-definite: {any_pd}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
