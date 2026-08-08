"""
Analyze eig_scan_second_window output into both outsourced verdicts.

Report 1 (DOMINANCE Layer-2): tabulate Delta_lambda(L) at each center L; report
  min over centers. Claim under test: inf Delta_lambda >= 0.158 (N=8). A single
  center below 0.158 refutes the UNIFORM bound (the report claims a window-wide
  infimum, so one interior counterexample suffices).

Report 2 (MARGIN_COLLAPSE section 4): at each center L with flanks L +/- delta,
  symmetric difference
     V2(L) = (eig_full(L+d) - eig_full(L-d)) / (2 d)
     V1(L) = (eig_p1(L+d)   - eig_p1(L-d))   / (2 d)
     DeltaV = V2 - V1.
  Conjecture "cross term delays collapse" predicts DeltaV > 0 (prime-2-only
  collapses FASTER than the full two-prime operator). Section-4 one-vote-veto:
  DeltaV <= 0 outside the error bar FALSIFIES the conjecture. Min-swap guard:
  skip / flag a center whose gap is smaller than a level-crossing threshold.

Error-bar model for DeltaV: the eig curve is a float-center pilot. The finite-
difference truncation error for a symmetric difference is O(d^2 * f''); with
d=1e-4 this is ~1e-8 * |f''|, negligible vs the signal. The dominant uncertainty
is the eig float rounding (~1e-12) amplified by 1/(2d)=5e3 -> ~5e-9 per slope.
So an |DeltaV| >> 1e-6 is comfortably outside the error bar. We report DeltaV and
the per-term slopes; a robust sign is the verdict.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def _find(rows, center, off, tol=1):
    """Locate the row whose L numerator is center+off (denominator assumed shared)."""
    target = center + off
    for r in rows:
        num = int(r["L"].split("/")[0])
        if abs(num - target) <= tol:
            return r
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", default="pilots/eig_scan_even_N8.json")
    ap.add_argument("--centers", default="55200,58000,60000,64000,68800")
    ap.add_argument("--delta-num", type=int, default=10, help="flank offset in L-numerator units")
    ap.add_argument("--den", type=int, default=100000)
    ap.add_argument("--dom-threshold", type=float, default=0.158)
    ap.add_argument("--gap-min", type=float, default=1e-3, help="min-swap guard threshold on eig gap")
    args = ap.parse_args()

    rows = json.loads(Path(args.scan).read_text())
    centers = [int(x) for x in args.centers.split(",")]
    delta = args.delta_num / args.den

    print("=" * 78)
    print("REPORT 1 — DOMINANCE Layer-2: Delta_lambda(L) coverage")
    print(f"  claim under test: inf Delta_lambda >= {args.dom_threshold}")
    print("=" * 78)
    print(f"{'L':>10} {'eig_full':>12} {'eig_off':>12} {'Delta_lambda':>14} {'gap_full':>10}")
    dom_vals = []
    for c in centers:
        r = _find(rows, c, 0)
        if r is None:
            print(f"{c/args.den:>10.5f}  (missing)")
            continue
        dom_vals.append((r["L_float"], r["Delta_lambda"]))
        print(f"{r['L_float']:>10.5f} {r['eig_full_min']:>+12.6f} {r['eig_off_min']:>+12.6f} "
              f"{r['Delta_lambda']:>+14.6f} {r['eig_full_gap']:>10.4f}")
    if dom_vals:
        Lmin, dmin = min(dom_vals, key=lambda t: t[1])
        verdict = "UPHELD" if dmin >= args.dom_threshold else "REFUTED"
        print(f"\n  min Delta_lambda over centers = {dmin:+.6f} at L={Lmin:.5f}")
        print(f"  vs claimed uniform bound {args.dom_threshold}: {verdict}")

    print("\n" + "=" * 78)
    print("REPORT 2 — MARGIN_COLLAPSE section 4: DeltaV = V2 - V1 sign")
    print(f"  symmetric difference delta={delta:.1e}; conjecture predicts DeltaV > 0")
    print("=" * 78)
    print(f"{'L':>10} {'V1(p1)':>12} {'V2(full)':>12} {'DeltaV':>12} {'gap_full':>9} {'guard':>8}")
    signs = []
    for c in centers:
        lo = _find(rows, c, -args.delta_num)
        mid = _find(rows, c, 0)
        hi = _find(rows, c, +args.delta_num)
        if not (lo and hi and mid):
            print(f"{c/args.den:>10.5f}  (missing flank)")
            continue
        V2 = (hi["eig_full_min"] - lo["eig_full_min"]) / (2 * delta)
        V1 = (hi["eig_p1_min"] - lo["eig_p1_min"]) / (2 * delta)
        dV = V2 - V1
        guard_ok = mid["eig_full_gap"] >= args.gap_min and mid["eig_p1_gap"] >= args.gap_min
        signs.append(dV)
        print(f"{mid['L_float']:>10.5f} {V1:>+12.5f} {V2:>+12.5f} {dV:>+12.5f} "
              f"{mid['eig_full_gap']:>9.4f} {'ok' if guard_ok else 'SWAP?':>8}")
    if signs:
        all_pos = all(s > 1e-6 for s in signs)
        all_nonpos = all(s <= 1e-6 for s in signs)
        if all_pos:
            v = "conjecture DIRECTION consistent (DeltaV>0 at all sampled L)"
        elif all_nonpos:
            v = "conjecture FALSIFIED (DeltaV<=0 outside error bar at all sampled L)"
        else:
            v = "MIXED sign across window (no uniform sign; conjecture not uniformly true)"
        print(f"\n  DeltaV signs: {[round(s,5) for s in signs]}")
        print(f"  verdict: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
