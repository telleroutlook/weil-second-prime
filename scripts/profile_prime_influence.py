"""
S4 — per-sector, per-prime, per-cross-term influence profile on the Schur pivot.

The steer (weil-first): the first window's even-sector prime term was nearly
inert (zeroing it moved the pivot ~0.003). Do NOT budget compute symmetrically
for the second window; first MEASURE where the pivot actually moves, then spend
hard compute there.

Efficiency: the archimedean block (S0 four-term, M0, T, Gd) is IDENTICAL across
all prime-layer variants and is the expensive part (~tens of minutes at large
N/depth). The prime layer (M2, S2) is cheap exact Fraction arithmetic. So we
build the archimedean block ONCE per (sector, L) and swap the prime layer
instantly for each variant. The archimedean truncation is common-moded across
variants, so influence deltas isolate the prime-layer effect.

Variants profiled (each an explicit, labelled experiment — NOT a silent omission):
  full        : both shifts + cross term (the real certificate assembly)
  tau2_off    : c2 = 0 (prime 2 removed)
  tau3_off    : c3 = 0 (prime 3 removed)      -> reduces to single-prime
  cross_off   : both shifts, cross term F = 0  (isolates cross-prime influence)
  scale_cross : cross term multiplied by a large factor (E1 probe: a term with
                tiny true influence is exposed by a large-factor scaling mutant,
                not by demanding a sign flip it cannot cause).

PRECISION DISCIPLINE (PROOF_CONSTITUTION A3): this profiler is PILOT GRADE
(float centers, depth=4). Its ranking is a screening tool for compute
allocation. Any DIRECTIONAL VERDICT ("prime p / the cross term dominates sector
s") must be re-confirmed at certify grade (Arb interval) before it steers real
compute. Output is labelled accordingly.

Long-run discipline (CLAUDE.md): observable (per-step flush prints), pausable
(KeyboardInterrupt -> checkpoint), resumable (--resume skips done (sector,L)).
"""
from __future__ import annotations

import argparse
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from typing import Optional

import numpy as np

from checker.fp035.recompute_schur import build_matrices, pivot_from_matrices
from src.prime_layer.legendre_shift_2prime import (
    C2,
    C3,
    M2_two_prime,
    S2_two_prime,
    default_F_provider,
    window_check,
)

CKPT = Path("pilots/s4_profile_checkpoint.json")


def _prime_layer(indices, L, c2, c3, cross: bool, cross_scale: float = 1.0):
    """Return (M2, S2) for a prime-layer variant. cross=False sets F=0 (labelled)."""
    M2 = np.array(M2_two_prime(indices, L, c2=c2, c3=c3))
    if c3 == 0.0:
        F_provider = None
    elif not cross:
        F_provider = lambda i, j, t2, t3: 0.0  # noqa: E731  (explicit probe)
    elif cross_scale != 1.0:
        F_provider = lambda i, j, t2, t3: cross_scale * default_F_provider(i, j, t2, t3)  # noqa: E731
    else:
        F_provider = default_F_provider
    S2 = np.array(S2_two_prime(indices, L, c2=c2, c3=c3, F_provider=F_provider))
    return M2, S2


def profile_point(L_num: int, L_den: int, sector: str, N: int, d: int,
                  eta: float = 0.5, cross_scale: float = 100.0) -> dict:
    """Profile one (sector, L) point. One archimedean build, many prime swaps."""
    L = L_num / L_den
    if not window_check(L):
        raise ValueError(f"L={L} outside second-prime window")
    t0 = time.time()
    print(f"[{sector}] L={L_num}/{L_den} N={N} d={d}: building archimedean block "
          f"(one-time, expensive)...", flush=True)
    mats = build_matrices(L_num, L_den, sector, N)
    indices = mats["indices"]
    print(f"[{sector}] archimedean built in {time.time()-t0:.0f}s; swapping prime layers",
          flush=True)

    variants = {
        "full":        dict(c2=C2, c3=C3, cross=True),
        "tau2_off":    dict(c2=0.0, c3=C3, cross=True),
        "tau3_off":    dict(c2=C2, c3=0.0, cross=True),
        "cross_off":   dict(c2=C2, c3=C3, cross=False),
        "scale_cross": dict(c2=C2, c3=C3, cross=True, cross_scale=cross_scale),
    }
    results = {}
    for name, kw in variants.items():
        M2, S2 = _prime_layer(indices, L, kw["c2"], kw["c3"], kw["cross"],
                              kw.get("cross_scale", 1.0))
        m = dict(mats); m["M2"] = M2; m["S2"] = S2
        piv, b_L = pivot_from_matrices(m, d, eta, judge="pivot")
        results[name] = {"min_pivot": piv, "b_L": b_L}
        print(f"    {name:12s} min_pivot={piv:+.6f}  b_L={b_L:+.5f}", flush=True)

    base = results["full"]["min_pivot"]
    influence = {
        "d_tau2":  base - results["tau2_off"]["min_pivot"],
        "d_tau3":  base - results["tau3_off"]["min_pivot"],
        "d_cross": base - results["cross_off"]["min_pivot"],
        "d_scale_cross": results["scale_cross"]["min_pivot"] - base,
    }
    print(f"    influence: d_tau2={influence['d_tau2']:+.6f} "
          f"d_tau3={influence['d_tau3']:+.6f} d_cross={influence['d_cross']:+.6f} "
          f"d_scale_cross(x{cross_scale:g})={influence['d_scale_cross']:+.6f}",
          flush=True)
    return {
        "L": L, "L_num": L_num, "L_den": L_den, "sector": sector, "N": N, "d": d,
        "eta": eta, "cross_scale": cross_scale,
        "variants": results, "influence": influence,
        "grade": "PILOT (float center, depth=4) — directional verdicts need certify",
        "elapsed_s": round(time.time() - t0, 1),
    }


def _load_ckpt() -> dict:
    if CKPT.exists():
        return json.loads(CKPT.read_text())
    return {"done": {}, "results": []}


def _save_ckpt(state: dict) -> None:
    CKPT.parent.mkdir(parents=True, exist_ok=True)
    CKPT.write_text(json.dumps(state, indent=2))


def main() -> int:
    ap = argparse.ArgumentParser(description="S4 prime-influence profile (pilot grade)")
    ap.add_argument("--points", default="0.60",
                    help="comma-separated L values, e.g. 0.55,0.60,0.69")
    ap.add_argument("--sectors", default="even,odd")
    ap.add_argument("--even-N", type=int, default=8, help="even sector N (d=2N)")
    ap.add_argument("--odd-N", type=int, default=7, help="odd sector N (d=2N+1)")
    ap.add_argument("--cross-scale", type=float, default=100.0)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--out", default="pilots/s4_profile.json")
    args = ap.parse_args()

    state = _load_ckpt() if args.resume else {"done": {}, "results": []}
    Ls = [float(x) for x in args.points.split(",")]
    sectors = args.sectors.split(",")

    for L in Ls:
        f = Fraction(L).limit_denominator(10000)
        for sector in sectors:
            key = f"{sector}@{f.numerator}/{f.denominator}"
            if key in state["done"]:
                print(f"[skip] {key} (already done)", flush=True)
                continue
            N = args.even_N if sector == "even" else args.odd_N
            parity = 0 if sector == "even" else 1
            d = 2 * N + parity  # first complement degree = 2N (even) / 2N+1 (odd)
            try:
                r = profile_point(f.numerator, f.denominator, sector, N, d,
                                  cross_scale=args.cross_scale)
            except KeyboardInterrupt:
                print(f"\n[interrupted] checkpoint saved; --resume to continue", flush=True)
                _save_ckpt(state)
                return 130
            state["results"].append(r)
            state["done"][key] = True
            _save_ckpt(state)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(state["results"], indent=2))
    print(f"\nProfile written to {args.out} ({len(state['results'])} points). "
          f"GRADE: PILOT — confirm any compute-allocation verdict at certify grade.",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
