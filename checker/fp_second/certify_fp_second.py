"""
certify_fp_second.py — 2-prime Arb-interval certify path for weil-second-prime.

All arithmetic uses Interval = Tuple[Fraction, Fraction] with outward rounding.
No float() in the proof chain. Final judge: min LDL^T pivot (ldlt.ldlt_factor).

Math:
    C = b_L * F - R_eta
    F   = T + M0 + M2 - c_L * diag(Gd)
    R_eta = (1+η)*R0 + (1+1/η)*R2
    R0  = S0 - M0^T * Ginv * M0
    R2  = S2 - M2^T * Ginv * M2
    b_L = H(d) - c_L - L0 - κ        (d = first free dimension: 2N+1 odd, 2N even)
    M2[i,j] = -(c2·J(i,j,τ2) + c3·J(i,j,τ3))
    S2[i,j] = c2²·E(i,j,τ2) + c3²·E(i,j,τ3) + c2·c3·(F_ij + F_ji)

Long-task requirements (CLAUDE.md):
    observable   — prints "[sector] step N/total (elapsed Xs)" per (a,b) pair
    pausable     — KeyboardInterrupt saves checkpoint and exits cleanly
    resumable    — --resume loads checkpoint and skips already-computed (a,b) pairs
    incremental-durable — checkpoint written after EACH (a,b) pair

Usage:
    python3 -m checker.fp_second.certify_fp_second \\
        --L 5600 10000 --sector odd --N 13 \\
        --out pilots/cert_fp_second_N13.json

    # Resume after interruption:
    python3 -m checker.fp_second.certify_fp_second \\
        --L 5600 10000 --sector odd --N 13 \\
        --out pilots/cert_fp_second_N13.json --resume
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from src.archimedean.interval import (
    Interval,
    add, sub, mul, scalar_mul, neg, point, hull,
    is_strictly_positive,
)
from src.archimedean.ldlt import ldlt_factor
from src.archimedean.integrator_a import (
    integrate_M_K, integrate_S_VK, integrate_S_KK, _arb_to_interval,
)
from src.archimedean.log_moments import V_matrix_entry, V2_matrix_entry
from src.archimedean.kernel import kappa as compute_kappa
from src.prime_layer.legendre_shift import compute_J, compute_E
from src.prime_layer.legendre_shift_2prime import compute_F

Matrix = List[List[Interval]]

ETA = Fraction(1, 2)
L0 = Fraction(1, 2 ** 30)


# ── Certified irrational constants via Arb ────────────────────────────────────

def _arb_iv(arb_val) -> Interval:
    return _arb_to_interval(arb_val)


def c2_iv(prec: int = 256) -> Interval:
    """c2 = log(2)/sqrt(2), certified."""
    from flint import arb, ctx
    ctx.prec = prec
    return _arb_iv(arb(2).log() / arb(2).sqrt())


def c3_iv(prec: int = 256) -> Interval:
    """c3 = log(3)/sqrt(3), certified."""
    from flint import arb, ctx
    ctx.prec = prec
    return _arb_iv(arb(3).log() / arb(3).sqrt())


def cL_iv(L_num: int, L_den: int, prec: int = 256) -> Interval:
    """c_L = log(2·π·L) + γ (Euler-Mascheroni), certified."""
    from flint import arb, ctx
    ctx.prec = prec
    L_arb = arb(L_num) / arb(L_den)
    return _arb_iv((arb(2) * arb.pi() * L_arb).log() + arb.const_euler())


def b_L_iv(d: int, L_num: int, L_den: int, prec: int = 256) -> Interval:
    """b_L = H(d) - c_L - L0 - κ as Interval."""
    H_d = sum((Fraction(1, k) for k in range(1, d + 1)), Fraction(0))
    kappa = compute_kappa(L_num, L_den, prec)
    cL = cL_iv(L_num, L_den, prec)
    b_lo = H_d - cL[1] - L0 - kappa
    b_hi = H_d - cL[0] - L0 - kappa
    return (b_lo, b_hi)


def tau_frac(L_num: int, L_den: int, prime: int) -> Fraction:
    return Fraction(math.log(prime) / (L_num / L_den)).limit_denominator(10_000)


def _H(n: int) -> Fraction:
    return sum((Fraction(1, k) for k in range(1, n + 1)), Fraction(0))


# ── Checkpoint serialization ──────────────────────────────────────────────────

def _iv_to_list(iv: Interval) -> list:
    return [str(iv[0]), str(iv[1])]


def _list_to_iv(lst: list) -> Interval:
    return (Fraction(lst[0]), Fraction(lst[1]))


def _save_checkpoint(path: Path, state: dict) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(state, indent=2))
    tmp.replace(path)


def _load_checkpoint(path: Path) -> dict:
    return json.loads(path.read_text())


# ── Interval matrix helpers ───────────────────────────────────────────────────

def _zeros(n: int) -> Matrix:
    return [[point(Fraction(0))] * n for _ in range(n)]


def _add(A: Matrix, B: Matrix) -> Matrix:
    n = len(A)
    return [[add(A[i][j], B[i][j]) for j in range(n)] for i in range(n)]


def _sub(A: Matrix, B: Matrix) -> Matrix:
    n = len(A)
    return [[sub(A[i][j], B[i][j]) for j in range(n)] for i in range(n)]


def _mul(A: Matrix, B: Matrix) -> Matrix:
    n = len(A)
    C = _zeros(n)
    for i in range(n):
        for j in range(n):
            s = point(Fraction(0))
            for k in range(n):
                s = add(s, mul(A[i][k], B[k][j]))
            C[i][j] = s
    return C


def _scale_iv(c: Interval, A: Matrix) -> Matrix:
    n = len(A)
    return [[mul(c, A[i][j]) for j in range(n)] for i in range(n)]


def _scale_frac(c: Fraction, A: Matrix) -> Matrix:
    n = len(A)
    return [[scalar_mul(c, A[i][j]) for j in range(n)] for i in range(n)]


def _T(A: Matrix) -> Matrix:
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


def _symmetrize(A: Matrix) -> Matrix:
    n = len(A)
    S = _zeros(n)
    for i in range(n):
        for j in range(i + 1):
            h = hull(A[i][j], A[j][i])
            S[i][j] = h
            S[j][i] = h
    return S


# ── Building block assembly (pausable/resumable) ──────────────────────────────

def build_matrices_iv(
    L_num: int, L_den: int, sector: str, N: int,
    depth_2d: int = 4, depth_3d: int = 3, prec: int = 256,
    checkpoint_path: Optional[Path] = None,
    resume: bool = False,
    use_bernstein: bool = True,
) -> dict:
    """Assemble M0, S0, M2, S2, T, Gd as Interval matrices.

    Pausable: KeyboardInterrupt saves checkpoint and re-raises.
    Resumable: if resume=True, loads existing checkpoint and skips done (a,b) pairs.
    Incremental-durable: checkpoint written after EACH (a,b) pair.

    use_bernstein: if True (default), use Bernstein ellipse remainder for M_K
        (valid for n_row < ~38, i.e. N ≤ 19 second window).
        If False, use Richardson GL-8/GL-4 remainder (required for N ≥ 20).
    """
    parity = 0 if sector == "even" else 1
    indices = list(range(parity, parity + 2 * N, 2))
    n = len(indices)

    tau2 = tau_frac(L_num, L_den, 2)
    tau3 = tau_frac(L_num, L_den, 3)
    c2 = c2_iv(prec)
    c3 = c3_iv(prec)

    Gd = [Fraction(2, 2 * ni + 1) for ni in indices]

    T: Matrix = _zeros(n)
    for a, ni in enumerate(indices):
        T[a][a] = point(_H(ni) * Gd[a])

    # Load or initialise checkpoint state
    ckpt: Dict[str, dict] = {"M0": {}, "S0": {}}
    if resume and checkpoint_path is not None and checkpoint_path.exists():
        ckpt = _load_checkpoint(checkpoint_path)
        loaded = sum(1 for k in ckpt.get("M0", {}))
        print(f"  [resume] loaded {loaded}/{n*n} pairs from {checkpoint_path}", flush=True)

    M0: Matrix = _zeros(n)
    S0: Matrix = _zeros(n)

    # Restore already-computed entries from checkpoint
    for key, iv_list in ckpt.get("M0", {}).items():
        a, b = map(int, key.split(","))
        M0[a][b] = _list_to_iv(iv_list)
    for key, iv_list in ckpt.get("S0", {}).items():
        a, b = map(int, key.split(","))
        S0[a][b] = _list_to_iv(iv_list)

    done_keys = set(ckpt.get("M0", {}).keys())
    total = n * n
    t0 = time.time()

    try:
        for a, ni in enumerate(indices):
            for b, nj in enumerate(indices):
                key = f"{a},{b}"
                if key in done_keys:
                    step = a * n + b + 1
                    print(
                        f"  [{sector}/N={N}] ({a},{b})=(P{ni},P{nj}) "
                        f"{step}/{total} [skip/cached]",
                        flush=True,
                    )
                    continue

                step = a * n + b + 1
                print(
                    f"  [{sector}/N={N}] archimedean ({a},{b})=(P{ni},P{nj}) "
                    f"{step}/{total} ({time.time()-t0:.0f}s)",
                    flush=True,
                )

                M0[a][b] = add(
                    V_matrix_entry(ni, nj, prec),
                    integrate_M_K(ni, nj, L_num, L_den,
                                  depth=depth_2d, prec=prec,
                                  use_bernstein=use_bernstein).to_interval(),
                )
                svv = V2_matrix_entry(ni, nj, prec)
                svk = integrate_S_VK(ni, nj, L_num, L_den,
                                     depth=depth_2d, prec=prec).to_interval()
                skv = integrate_S_VK(nj, ni, L_num, L_den,
                                     depth=depth_2d, prec=prec).to_interval()
                skk = integrate_S_KK(ni, nj, L_num, L_den,
                                     depth=depth_3d, prec=prec).to_interval()
                S0[a][b] = add(add(add(svv, svk), skv), skk)

                # Incremental-durable: write after EACH completed pair
                ckpt["M0"][key] = _iv_to_list(M0[a][b])
                ckpt["S0"][key] = _iv_to_list(S0[a][b])
                ckpt["meta"] = {
                    "L": f"{L_num}/{L_den}", "sector": sector, "N": N,
                    "depth_2d": depth_2d, "depth_3d": depth_3d,
                }
                if checkpoint_path is not None:
                    _save_checkpoint(checkpoint_path, ckpt)

    except KeyboardInterrupt:
        if checkpoint_path is not None:
            _save_checkpoint(checkpoint_path, ckpt)
            done = sum(1 for k in ckpt.get("M0", {}))
            print(
                f"\n  [interrupt] checkpoint saved ({done}/{total} pairs) → {checkpoint_path}",
                flush=True,
            )
        raise

    # Prime layer (fast: Fraction arithmetic, no checkpoint needed)
    M2: Matrix = _zeros(n)
    S2: Matrix = _zeros(n)
    for a, ni in enumerate(indices):
        for b, nj in enumerate(indices):
            J2 = point(compute_J(ni, nj, tau2))
            J3 = point(compute_J(ni, nj, tau3))
            M2[a][b] = neg(add(mul(c2, J2), mul(c3, J3)))
            E2 = point(compute_E(ni, nj, tau2))
            E3 = point(compute_E(ni, nj, tau3))
            F_ij = point(compute_F(ni, nj, tau2, tau3))
            F_ji = point(compute_F(nj, ni, tau2, tau3))
            c2_sq = mul(c2, c2)
            c3_sq = mul(c3, c3)
            c2c3  = mul(c2, c3)
            S2[a][b] = add(
                add(mul(c2_sq, E2), mul(c3_sq, E3)),
                mul(c2c3, add(F_ij, F_ji)),
            )

    return {
        "indices": indices, "n": n,
        "M0": M0, "S0": S0, "M2": M2, "S2": S2,
        "Gd": Gd, "T": T,
        "tau2": tau2, "tau3": tau3,
        "c2": c2, "c3": c3,
    }


def assemble_schur_iv(
    mats: dict, d: int, L_num: int, L_den: int,
    eta: Fraction = ETA, prec: int = 256,
) -> Matrix:
    """Assemble C = b_L·F - R_eta as an Interval matrix and symmetrize."""
    n    = mats["n"]
    M0   = mats["M0"];  S0 = mats["S0"]
    M2   = mats["M2"];  S2 = mats["S2"]
    T    = mats["T"];   Gd = mats["Gd"]

    Ginv: Matrix = _zeros(n)
    for i in range(n):
        Ginv[i][i] = point(Fraction(1) / Gd[i])

    R0 = _sub(S0, _mul(_mul(_T(M0), Ginv), M0))
    R2 = _sub(S2, _mul(_mul(_T(M2), Ginv), M2))

    one_plus_eta     = Fraction(1) + eta
    one_plus_inv_eta = Fraction(1) + Fraction(1) / eta
    R_eta = _add(
        _scale_frac(one_plus_eta,     R0),
        _scale_frac(one_plus_inv_eta, R2),
    )

    cL  = cL_iv(L_num, L_den, prec)
    b_L = b_L_iv(d, L_num, L_den, prec)

    cL_diag: Matrix = _zeros(n)
    for i in range(n):
        cL_diag[i][i] = mul(cL, point(Gd[i]))

    F_mat = _sub(_add(_add(T, M0), M2), cL_diag)
    C     = _sub(_scale_iv(b_L, F_mat), R_eta)
    return _symmetrize(C)


# ── Top-level certify call ────────────────────────────────────────────────────

def certify_sector(
    L_num: int, L_den: int, sector: str, N: int,
    depth_2d: int = 4, depth_3d: int = 3, prec: int = 256,
    out_path: Optional[Path] = None,
    resume: bool = False,
    use_bernstein: bool = True,
    eta: Fraction = ETA,
) -> dict:
    parity = 0 if sector == "even" else 1
    d = 2 * N + parity

    t0 = time.time()
    print(
        f"[certify_fp_second] L={L_num}/{L_den}={L_num/L_den:.4f} "
        f"sector={sector} N={N} d={d}",
        flush=True,
    )

    b_L = b_L_iv(d, L_num, L_den, prec)
    if not is_strictly_positive(b_L):
        msg = f"b_L=[{float(b_L[0]):.4e}, {float(b_L[1]):.4e}] not strictly positive"
        print(f"  ✗ {msg}", flush=True)
        return {
            "L": f"{L_num}/{L_den}", "L_float": L_num / L_den,
            "sector": sector, "N": N, "d": d,
            "eta": str(eta), "eta_float": float(eta),
            "certified": False, "reason": msg,
        }
    print(
        f"  b_L ∈ [{float(b_L[0]):.6f}, {float(b_L[1]):.6f}]  (strictly positive ✓)",
        flush=True,
    )

    # Checkpoint lives next to the output file (or in cwd)
    ckpt_path: Optional[Path] = None
    if out_path is not None:
        ckpt_path = out_path.with_suffix(".ckpt.json")
    else:
        ckpt_path = Path(f"certify_fp_second_{sector}_N{N}.ckpt.json")

    mats = build_matrices_iv(
        L_num, L_den, sector, N,
        depth_2d=depth_2d, depth_3d=depth_3d, prec=prec,
        checkpoint_path=ckpt_path, resume=resume,
        use_bernstein=use_bernstein,
    )

    print(f"  assembling Schur C = b_L·F − R_eta (eta={float(eta):.4f}) ...", flush=True)
    C = assemble_schur_iv(mats, d, L_num, L_den, eta=eta, prec=prec)

    print(f"  LDL^T factorization (n={mats['n']}) ...", flush=True)
    certified = False
    min_pivot_lb: Optional[Fraction] = None
    error_msg: Optional[str] = None

    try:
        _, pivots = ldlt_factor(C)
        certified = True
        min_pivot_lb = min(p[0] for p in pivots)
    except ValueError as exc:
        error_msg = str(exc)

    elapsed = time.time() - t0

    result: dict = {
        "L": f"{L_num}/{L_den}",
        "L_float": L_num / L_den,
        "sector": sector,
        "N": N,
        "d": d,
        "eta": str(eta),
        "eta_float": float(eta),
        "certified": certified,
        "elapsed_s": round(elapsed, 1),
    }
    if certified and min_pivot_lb is not None:
        result["min_pivot_lower"] = str(min_pivot_lb)
        result["min_pivot_lower_float"] = float(min_pivot_lb)
        print(f"  ✓ CERTIFIED: min_pivot_lower = {float(min_pivot_lb):.6e}", flush=True)
        # Remove checkpoint on success
        if ckpt_path is not None and ckpt_path.exists():
            ckpt_path.unlink()
    else:
        result["error"] = error_msg or "unknown"
        print(f"  ✗ NOT CERTIFIED: {error_msg}", flush=True)

    # Incremental-durable: write partial result to out_path as soon as done
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(result, indent=2))
        print(f"[certify_fp_second] result → {out_path}", flush=True)

    return result


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(
        description="2-prime Schur positivity certifier (Arb intervals)"
    )
    p.add_argument("--L", nargs=2, type=int, default=[5600, 10000],
                   metavar=("NUM", "DEN"),
                   help="L as rational NUM/DEN, default 5600/10000 = 0.56")
    p.add_argument("--sector", choices=["even", "odd"], default="odd")
    p.add_argument("--N", type=int, default=13,
                   help="basis truncation size")
    p.add_argument("--depth2", type=int, default=4,
                   help="depth for 2-D integrals (M_K, S_VK)")
    p.add_argument("--depth3", type=int, default=3,
                   help="depth for S_KK (via Legendre expansion)")
    p.add_argument("--prec", type=int, default=256,
                   help="Arb working precision in bits")
    p.add_argument("--out", default=None,
                   help="write JSON result to this file")
    p.add_argument("--resume", action="store_true",
                   help="load checkpoint and skip already-computed (a,b) pairs")
    p.add_argument("--eta", type=str, default="1/2",
                   help="residual weight eta as fraction NUM/DEN or float, default 1/2")
    p.add_argument("--no-bernstein", action="store_true",
                   help="use Richardson GL-8/GL-4 remainder instead of Bernstein "
                        "(required for N >= 20, second window)")
    args = p.parse_args()

    out_path = Path(args.out) if args.out else None

    # Parse eta: accept "NUM/DEN" or float string
    eta_str = args.eta
    if "/" in eta_str:
        num_s, den_s = eta_str.split("/", 1)
        eta_frac = Fraction(int(num_s), int(den_s))
    else:
        eta_frac = Fraction(eta_str).limit_denominator(100000)

    try:
        result = certify_sector(
            args.L[0], args.L[1], args.sector, args.N,
            depth_2d=args.depth2, depth_3d=args.depth3, prec=args.prec,
            out_path=out_path, resume=args.resume,
            use_bernstein=not args.no_bernstein,
            eta=eta_frac,
        )
    except KeyboardInterrupt:
        print("\n[certify_fp_second] interrupted — checkpoint saved, re-run with --resume",
              flush=True)
        sys.exit(130)

    if args.out is None:
        print(json.dumps(result, indent=2))

    sys.exit(0 if result["certified"] else 1)


if __name__ == "__main__":
    main()
