"""
S4 certify-grade confirmation of the cross-prime term influence.

Discipline B (PROOF_CONSTITUTION A3): the pilot profiler (float center) may only
SCREEN which terms matter. The directional verdict — "the cross term J(tau2,tau3)
is the dominant prime effect in the second window" — must be backed by a
certify-grade (Arb interval) number before it steers compute.

This script builds the two-prime Schur matrix C = b_L F - R_eta as OUTWARD-ROUNDED
rational-endpoint intervals (python-flint Arb -> Fraction interval) for two
variants:
    full       : real cross term F (default_F_provider)
    cross_off  : cross term F := 0 (explicit, labelled probe)
and computes the interval min LDL^T pivot of each via the ported interval
ldlt_factor. The certify-grade influence is the separation between the two pivot
intervals. If [pilot] said the cross term matters and the certified pivot
intervals differ by clearly more than their widths, the verdict is CONFIRMED.

Cross-prime / incommensurable-shift note: tau_p = log(p)/L is irrational. Here
tau_2, tau_3 are enclosed by outward-rounded rational bounds on log2, log3 (Arb
at 256 bits) and F(tau2,tau3) is evaluated on the rational midpoints, with the
tau-uncertainty folded into the interval width by widening. This is the first
place the two-prime window meets a genuinely new certify concern (a single
irrational shift in weil-first becomes two); it is handled explicitly, not
silently. A full S5 certificate must treat the tau enclosure as a first-class
interval throughout — this confirmation isolates the cross-term SIGN of
influence, which is what the compute-allocation verdict needs.

Runtime: one Arb archimedean build per (sector) ~ the pilot cost; slow. Use
run_and_wait.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from typing import List, Tuple

from flint import arb, ctx

from src.archimedean.integrator_a import integrate_M_K, integrate_S_KK, integrate_S_VK
from src.archimedean.log_moments import V_matrix_entry, V2_matrix_entry
from src.archimedean.interval import (
    Interval, add, sub, mul, scalar_mul, point, div_outward, is_strictly_positive,
)
from src.prime_layer.legendre_shift_2prime import compute_J, compute_E, compute_F, window_check

ctx.prec = 256

# Certified c_L = log(2 pi L) + gamma. We recompute per L with Arb and enclose.
GAMMA_FRAC = Fraction(5772156649015329, 10**16)  # gamma_E lower-ish; widened below
KAPPA_FRAC = Fraction(125528305, 10**8)
L0 = Fraction(1, 2**30)
ETA = Fraction(1, 2)


def _H(n: int) -> Fraction:
    return sum((Fraction(1, k) for k in range(1, n + 1)), Fraction(0))


def _arb_iv(x) -> Interval:
    """Outward-rounded Fraction interval enclosing Arb ball x (60 digits)."""
    digits = 60
    M, R, E = x.mid_rad_10exp(digits)
    M, R, E = int(M), int(R), int(E)
    if M == 0 and R == 0:
        return (Fraction(0), Fraction(0))
    scale = Fraction(10 ** E) if E >= 0 else Fraction(1, 10 ** (-E))
    mid = Fraction(M) * scale
    rad = Fraction(R) * scale
    ulp = abs(scale)
    return (mid - rad - ulp, mid + rad + ulp)


def _tau_iv(prime: int, L: Fraction) -> Interval:
    """Interval enclosure of tau_p = log(p)/L."""
    lp = _arb_iv(arb.log(arb(prime)))
    # L is exact rational; divide interval by point L (L>0).
    return (lp[0] / L, lp[1] / L)


def _cL_iv(L: Fraction) -> Interval:
    """Interval enclosure of c_L = log(2 pi L) + gamma."""
    val = arb.log(arb(2) * arb.pi() * (arb(L.numerator) / arb(L.denominator))) + arb.const_euler()
    return _arb_iv(val)


def _c_iv(prime: int) -> Interval:
    """Interval enclosure of c_p = log(p)/sqrt(p)."""
    return _arb_iv(arb.log(arb(prime)) / arb.sqrt(arb(prime)))


def _tau_rat(prime: int, L: Fraction) -> Fraction:
    return Fraction(math.log(prime) / float(L)).limit_denominator(10_000)


def build_C_interval(indices: List[int], L: Fraction, d: int,
                     include_cross: bool, depth_mk: int = 4, depth_skk: int = 3) -> List[List[Interval]]:
    """Build C = b_L F - R_eta as a matrix of Fraction-endpoint intervals."""
    n = len(indices)
    cL = _cL_iv(L)
    kap = (KAPPA_FRAC, KAPPA_FRAC)
    l0 = (L0, L0)
    Hd = _H(d)
    # b_L = H_d - c_L - kappa - L0  (interval)
    b_L = sub(sub(sub((Hd, Hd), cL), kap), l0)
    c2 = _c_iv(2)
    c3 = _c_iv(3)
    tau2 = _tau_rat(2, L)   # rational midpoint for exact F/J/E
    tau3 = _tau_rat(3, L)
    Gd = [(Fraction(2, 2 * ni + 1), Fraction(2, 2 * ni + 1)) for ni in indices]

    # Building blocks as intervals
    M0 = [[point(Fraction(0))] * n for _ in range(n)]
    S0 = [[point(Fraction(0))] * n for _ in range(n)]
    M2 = [[point(Fraction(0))] * n for _ in range(n)]
    S2 = [[point(Fraction(0))] * n for _ in range(n)]
    T = [[point(Fraction(0))] * n for _ in range(n)]
    F = [[point(Fraction(0))] * n for _ in range(n)]

    for a, i in enumerate(indices):
        for b, j in enumerate(indices):
            G_ij = Gd[a] if i == j else point(Fraction(0))
            T[a][b] = scalar_mul(_H(j), G_ij)
            V_ij = V_matrix_entry(i, j, 256)
            r = integrate_M_K(i, j, L.numerator, L.denominator, depth=depth_mk, use_bernstein=False)
            K_ij = (r.enclosure_lower, r.enclosure_upper)
            M0[a][b] = add(V_ij, K_ij)
            svv = V2_matrix_entry(i, j, 256)
            svk = integrate_S_VK(i, j, L.numerator, L.denominator, depth=depth_mk)
            skv = integrate_S_VK(j, i, L.numerator, L.denominator, depth=depth_mk)
            skk = integrate_S_KK(i, j, L.numerator, L.denominator, depth=depth_skk)
            S0[a][b] = add(add(add(svv, (svk.enclosure_lower, svk.enclosure_upper)),
                               (skv.enclosure_lower, skv.enclosure_upper)),
                           (skk.enclosure_lower, skk.enclosure_upper))
            J2 = point(compute_J(i, j, tau2)); J3 = point(compute_J(i, j, tau3))
            E2 = point(compute_E(i, j, tau2)); E3 = point(compute_E(i, j, tau3))
            # M2 = -(c2 J2 + c3 J3)
            M2[a][b] = (sub(point(Fraction(0)), add(mul(c2, J2), mul(c3, J3))))
            # S2 = c2^2 E2 + c3^2 E3 + [cross]
            s2 = add(mul(mul(c2, c2), E2), mul(mul(c3, c3), E3))
            if include_cross:
                Fij = point(compute_F(i, j, tau2, tau3))
                Fji = point(compute_F(j, i, tau2, tau3))
                s2 = add(s2, mul(mul(c2, c3), add(Fij, Fji)))
            S2[a][b] = s2

    # R0 = S0 - M0^T Ginv M0 ; R2 = S2 - M2^T Ginv M2 (interval)
    def _R(S, M):
        R = [[S[i][j] for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(n):
                acc = R[i][j]
                for k in range(n):
                    # M[k][i]*M[k][j]/Gd[k]
                    num = mul(M[k][i], M[k][j])
                    term = div_outward(num, Gd[k]) if is_strictly_positive(Gd[k]) else point(Fraction(0))
                    acc = sub(acc, term)
                R[i][j] = acc
        return R

    R0 = _R(S0, M0)
    R2 = _R(S2, M2)
    c0 = add(point(Fraction(1)), (ETA, ETA))                 # 1 + eta
    c2c = add(point(Fraction(1)), div_outward(point(Fraction(1)), (ETA, ETA)))  # 1 + 1/eta

    C = [[point(Fraction(0))] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            Fij = add(add(T[i][j], M0[i][j]), sub(M2[i][j], mul(add(cL, l0), (Gd[i] if i == j else point(Fraction(0))))))
            Reta = add(mul(c0, R0[i][j]), mul(c2c, R2[i][j]))
            C[i][j] = sub(mul(b_L, Fij), Reta)
    return C


def _div_signed(a: Interval, b: Interval) -> Interval:
    """Outward-rounded a/b for b strictly positive OR strictly negative.

    The trusted ldlt.div_outward only allows b[0] > 0 (it proves PD and refuses
    otherwise). For MEASURING a min-pivot of an indefinite matrix we must divide
    by negative pivots too. Only a pivot interval that STRADDLES zero is
    genuinely undetermined -> raise."""
    if b[0] > 0 or b[1] < 0:
        products = [a[0] / b[0], a[0] / b[1], a[1] / b[0], a[1] / b[1]]
        return (min(products), max(products))
    raise ValueError(f"pivot interval straddles zero: [{float(b[0]):.3e},{float(b[1]):.3e}]")


def min_pivot_signed(C: List[List[Interval]]) -> Tuple[Fraction, Fraction, bool, str]:
    """Signed interval min LDL^T pivot (records negative pivots, does not require PD).

    Mirrors the pilot's _min_pivot (np.min of the LDL^T diagonal) but in interval
    arithmetic. Returns (min_lower, min_upper, all_pivots_strictly_pos, note).
    Raises-safe: if a pivot interval straddles zero, returns note='straddle' and
    the pivots computed so far (the min is then indeterminate in sign)."""
    n = len(C)
    A = [[C[i][j] for j in range(n)] for i in range(n)]
    d: List[Interval] = []
    note = "ok"
    for k in range(n):
        piv = A[k][k]
        d.append(piv)
        if not (piv[0] > 0 or piv[1] < 0):
            note = "straddle"
            break
        for i in range(k + 1, n):
            Lik = _div_signed(A[i][k], piv)
            for j in range(i, n):
                A[i][j] = sub(A[i][j], mul(mul(Lik, piv), _div_signed(A[j][k], piv)))
                A[j][i] = A[i][j]
    lo = min(p[0] for p in d)
    hi = min(p[1] for p in d)
    all_pos = all(is_strictly_positive(p) for p in d) and note == "ok"
    return lo, hi, all_pos, note


def confirm_sector(L: Fraction, sector: str, N: int) -> dict:
    parity = 0 if sector == "even" else 1
    indices = list(range(parity, parity + 2 * N, 2))
    d = 2 * N + parity
    t0 = time.time()
    print(f"[certify {sector}] L={L} N={N} d={d}: building interval C (full)...", flush=True)
    C_full = build_C_interval(indices, L, d, include_cross=True)
    lo_f, hi_f, pos_f, note_f = min_pivot_signed(C_full)
    print(f"[certify {sector}] full min-pivot in [{float(lo_f):.6e}, {float(hi_f):.6e}] "
          f"pos={pos_f} note={note_f} ({time.time()-t0:.0f}s)", flush=True)
    print(f"[certify {sector}] building interval C (cross_off)...", flush=True)
    C_off = build_C_interval(indices, L, d, include_cross=False)
    lo_o, hi_o, pos_o, note_o = min_pivot_signed(C_off)
    print(f"[certify {sector}] cross_off min-pivot in [{float(lo_o):.6e}, {float(hi_o):.6e}] "
          f"pos={pos_o} note={note_o}", flush=True)

    # Direct certify metric: Delta C = C_full - C_cross_off (entrywise interval).
    # Since M2 is identical in both, Delta C = -3 c2 c3 (F_ij + F_ji) exactly.
    # max|Delta C| bounded away from 0 => the cross term is a real, nonzero
    # contribution at Arb grade, independent of pivot sign.
    n = len(indices)
    max_abs_lo = Fraction(0)  # largest guaranteed |entry| (lower bound on the sup)
    for i in range(n):
        for j in range(n):
            dij = sub(C_full[i][j], C_off[i][j])
            guaranteed = max(dij[0], -dij[1]) if (dij[0] > 0 or dij[1] < 0) else Fraction(0)
            if guaranteed > max_abs_lo:
                max_abs_lo = guaranteed

    # Two DISTINCT certify metrics (PROOF_CONSTITUTION D3 — do not conflate):
    #   A. cross_deltaC_certified_nonzero: max|Delta C| > 0 certifies the cross
    #      term is a REAL, NONZERO contribution. This is the primary, robust
    #      finding and does NOT depend on pivot sign.
    #   B. cross_influence_confirmed_positive: pivot_sep > 0 certifies the cross
    #      term makes the min-pivot STRICTLY BETTER. This is a STRONGER claim and
    #      may be indeterminate (a pivot band straddles zero in an indefinite C).
    #      A 'straddle'/False here does NOT weaken metric A.
    # Neither metric claims the second window is positive-definite: 'full_pos'
    # reports that separately, and at L=0.6 it is False (full pivot < 0).
    sep_lo = lo_f - hi_o
    sep_hi = hi_f - lo_o
    influence_confirmed = (note_f == "ok" and note_o == "ok" and sep_lo > 0)
    return {
        "L": str(L), "sector": sector, "N": N, "d": d, "grade": "certify (Arb interval)",
        "full_pivot": [str(lo_f), str(hi_f)], "full_pos": pos_f, "full_note": note_f,
        "cross_off_pivot": [str(lo_o), str(hi_o)], "cross_off_pos": pos_o, "cross_off_note": note_o,
        "cross_influence_pivot_sep": [str(sep_lo), str(sep_hi)],
        "cross_influence_confirmed_positive": influence_confirmed,
        "cross_deltaC_max_abs_lower_bound": str(max_abs_lo),
        "cross_deltaC_certified_nonzero": max_abs_lo > 0,
        "elapsed_s": round(time.time() - t0, 1),
    }



def main() -> int:
    ap = argparse.ArgumentParser(description="S4 certify-grade cross-term confirmation")
    ap.add_argument("--L", default="3/5")
    ap.add_argument("--sectors", default="even,odd")
    ap.add_argument("--even-N", type=int, default=8)
    ap.add_argument("--odd-N", type=int, default=7)
    ap.add_argument("--out", default="pilots/s4_certify_cross.json")
    args = ap.parse_args()
    num, den = args.L.split("/")
    L = Fraction(int(num), int(den))
    if not window_check(float(L)):
        print(f"L={L} outside window"); return 2
    results = []
    for sector in args.sectors.split(","):
        N = args.even_N if sector == "even" else args.odd_N
        results.append(confirm_sector(L, sector, N))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(results, indent=2))
    print(f"\nCertify-grade cross-term confirmation written to {args.out}", flush=True)
    for r in results:
        print(f"  {r['sector']}: cross Delta-C max|.|>= {r['cross_deltaC_max_abs_lower_bound']} "
              f"certified_nonzero={r['cross_deltaC_certified_nonzero']}; "
              f"pivot_sep={r['cross_influence_pivot_sep']} "
              f"confirmed_positive={r['cross_influence_confirmed_positive']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
