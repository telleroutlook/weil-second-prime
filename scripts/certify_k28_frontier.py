"""Focused Arb certificate for the k=28 frontier mode.

The float k=28 chain found a negative eigenvector concentrated on P53 and P55.
This script verifies a strictly rational version of that diagnosis without
building the full 28x28 interval matrix.  It computes only the two M0 columns
needed by the P53/P55 Schur-complement entries, then evaluates the integer
Rayleigh witness v = 3*e(P53) + e(P55).  If the outward interval for

    v^T C v = 9*C[53,53] + 6*C[53,55] + C[55,55]

has strictly negative upper endpoint, C is certifiably not positive definite.

Long-task discipline: every computed integral is printed and checkpointed
individually; KeyboardInterrupt exits with a resumable checkpoint.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from checker.fp_second.certify_fp_second import (
    b_L_iv,
    c2_iv,
    c3_iv,
    cL_iv,
    tau_frac,
)
from src.archimedean.integrator_a import integrate_M_K, integrate_S_KK, integrate_S_VK
from src.archimedean.interval import (
    Interval,
    add,
    mul,
    point,
    scalar_mul,
    sub,
)
from src.archimedean.log_moments import V2_matrix_entry, V_matrix_entry
from src.prime_layer.legendre_shift import compute_E, compute_J
from src.prime_layer.legendre_shift_2prime import compute_F

Key = Tuple[int, int]

N = 28
TARGET_A = 26  # zero-based position of P53
TARGET_B = 27  # zero-based position of P55
WITNESS = (3, 1)
ETA = Fraction(1, 1)


def _indices() -> List[int]:
    return list(range(1, 1 + 2 * N, 2))


def _harmonic(n: int) -> Fraction:
    return sum((Fraction(1, k) for k in range(1, n + 1)), Fraction(0))


def required_m0_keys() -> List[Key]:
    """Return unique symmetric M0 entries needed for columns A and B."""
    keys = set()
    for row in range(N):
        for target in (TARGET_A, TARGET_B):
            keys.add((min(row, target), max(row, target)))
    return sorted(keys)


def _iv_to_list(iv: Interval) -> List[str]:
    return [str(iv[0]), str(iv[1])]


def _list_to_iv(value: List[str]) -> Interval:
    if len(value) != 2:
        raise ValueError(f"invalid interval list: {value!r}")
    return (Fraction(value[0]), Fraction(value[1]))


def _save_json(path: Path, value: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2))
    tmp.replace(path)


def rayleigh_numerator(
    c_aa: Interval, c_ab: Interval, c_bb: Interval,
    witness: Tuple[int, int] = WITNESS,
) -> Interval:
    """Return v^T C v for v=(witness[0], witness[1])."""
    wa, wb = witness
    aa = scalar_mul(Fraction(wa * wa), c_aa)
    ab = scalar_mul(Fraction(2 * wa * wb), c_ab)
    bb = scalar_mul(Fraction(wb * wb), c_bb)
    return add(add(aa, ab), bb)


def _m0_entry(
    ni: int, nj: int, L_num: int, L_den: int,
    depth_2d: int, prec: int,
) -> Interval:
    v = V_matrix_entry(ni, nj, prec)
    k = integrate_M_K(
        ni, nj, L_num, L_den,
        depth=depth_2d, prec=prec, use_bernstein=False,
    ).to_interval()
    return add(v, k)


def _s0_entry(
    ni: int, nj: int, L_num: int, L_den: int,
    depth_2d: int, depth_3d: int, prec: int,
) -> Interval:
    svv = V2_matrix_entry(ni, nj, prec)
    svk = integrate_S_VK(
        ni, nj, L_num, L_den, depth=depth_2d, prec=prec,
    ).to_interval()
    skv = integrate_S_VK(
        nj, ni, L_num, L_den, depth=depth_2d, prec=prec,
    ).to_interval()
    skk = integrate_S_KK(
        ni, nj, L_num, L_den, depth=depth_3d, prec=prec,
    ).to_interval()
    return add(add(add(svv, svk), skv), skk)


def _m2_s2_entries(
    L_num: int, L_den: int, prec: int,
) -> tuple[dict[Key, Interval], dict[Key, Interval]]:
    indices = _indices()
    tau2 = tau_frac(L_num, L_den, 2)
    tau3 = tau_frac(L_num, L_den, 3)
    c2 = c2_iv(prec)
    c3 = c3_iv(prec)

    required_m2 = set(required_m0_keys())
    for a, b in ((TARGET_A, TARGET_A), (TARGET_A, TARGET_B), (TARGET_B, TARGET_B)):
        required_m2.add((min(a, b), max(a, b)))
    required_s2 = {
        (TARGET_A, TARGET_A), (TARGET_A, TARGET_B), (TARGET_B, TARGET_B),
    }

    m2: dict[Key, Interval] = {}
    s2: dict[Key, Interval] = {}
    for key in sorted(required_m2 | required_s2):
        a, b = key
        ni, nj = indices[a], indices[b]
        j2 = point(compute_J(ni, nj, tau2))
        j3 = point(compute_J(ni, nj, tau3))
        m2[key] = add(
            scalar_mul(Fraction(-1), mul(c2, j2)),
            scalar_mul(Fraction(-1), mul(c3, j3)),
        )

        e2 = point(compute_E(ni, nj, tau2))
        e3 = point(compute_E(ni, nj, tau3))
        fij = point(compute_F(ni, nj, tau2, tau3))
        fji = point(compute_F(nj, ni, tau2, tau3))
        c2sq = mul(c2, c2)
        c3sq = mul(c3, c3)
        c2c3 = mul(c2, c3)
        s2[key] = add(
            add(mul(c2sq, e2), mul(c3sq, e3)),
            mul(c2c3, add(fij, fji)),
        )
    return m2, s2


def _residual_entry(
    s_entry: Interval,
    col_x: List[Interval],
    col_y: List[Interval],
    gram_diagonal: List[Fraction],
) -> Interval:
    product_sum = point(Fraction(0))
    for x, y, gd in zip(col_x, col_y, gram_diagonal):
        term = mul(mul(x, point(Fraction(1) / gd)), y)
        product_sum = add(product_sum, term)
    return sub(s_entry, product_sum)


def _columns(
    entries: dict[Key, Interval], target: int,
) -> List[Interval]:
    return [
        entries[(min(row, target), max(row, target))]
        for row in range(N)
    ]


def certify(
    L_num: int,
    L_den: int,
    depth_2d: int,
    depth_3d: int,
    prec: int,
    out_path: Path,
    resume: bool,
) -> dict:
    started = time.time()
    indices = _indices()
    gram_diagonal = [Fraction(2, 2 * ni + 1) for ni in indices]
    m0_keys = required_m0_keys()
    s0_keys = [
        (TARGET_A, TARGET_A),
        (TARGET_A, TARGET_B),
        (TARGET_B, TARGET_B),
    ]
    total_units = len(m0_keys) + len(s0_keys)
    ckpt_path = out_path.with_suffix(".frontier.ckpt.json")
    meta = {
        "L": f"{L_num}/{L_den}",
        "sector": "odd",
        "N": N,
        "target_positions": [TARGET_A, TARGET_B],
        "target_degrees": [indices[TARGET_A], indices[TARGET_B]],
        "witness": list(WITNESS),
        "depth_2d": depth_2d,
        "depth_3d": depth_3d,
        "prec": prec,
        "use_bernstein": False,
    }

    state: dict = {"meta": meta, "M0": {}, "S0": {}}
    if resume and ckpt_path.exists():
        loaded = json.loads(ckpt_path.read_text())
        if loaded.get("meta") != meta:
            raise ValueError(
                f"checkpoint metadata mismatch for {ckpt_path}; use a new output path"
            )
        state = loaded
        print(
            f"[resume] M0={len(state.get('M0', {}))}/{len(m0_keys)}, "
            f"S0={len(state.get('S0', {}))}/{len(s0_keys)}",
            flush=True,
        )

    completed = len(state.get("M0", {})) + len(state.get("S0", {}))
    try:
        for a, b in m0_keys:
            key = f"{a},{b}"
            if key in state["M0"]:
                continue
            completed += 1
            print(
                f"[k28-frontier] M0 ({a},{b})=(P{indices[a]},P{indices[b]}) "
                f"{completed}/{total_units} ({time.time()-started:.0f}s)",
                flush=True,
            )
            state["M0"][key] = _iv_to_list(
                _m0_entry(
                    indices[a], indices[b], L_num, L_den,
                    depth_2d, prec,
                )
            )
            _save_json(ckpt_path, state)

        for a, b in s0_keys:
            key = f"{a},{b}"
            if key in state["S0"]:
                continue
            completed += 1
            print(
                f"[k28-frontier] S0 ({a},{b})=(P{indices[a]},P{indices[b]}) "
                f"{completed}/{total_units} ({time.time()-started:.0f}s)",
                flush=True,
            )
            state["S0"][key] = _iv_to_list(
                _s0_entry(
                    indices[a], indices[b], L_num, L_den,
                    depth_2d, depth_3d, prec,
                )
            )
            _save_json(ckpt_path, state)
    except KeyboardInterrupt:
        _save_json(ckpt_path, state)
        done = len(state.get("M0", {})) + len(state.get("S0", {}))
        print(
            f"\n[interrupt] checkpoint saved ({done}/{total_units}) -> {ckpt_path}",
            flush=True,
        )
        raise

    m0: dict[Key, Interval] = {
        tuple(map(int, key.split(","))): _list_to_iv(value)
        for key, value in state["M0"].items()
    }
    s0: dict[Key, Interval] = {
        tuple(map(int, key.split(","))): _list_to_iv(value)
        for key, value in state["S0"].items()
    }
    m2, s2 = _m2_s2_entries(L_num, L_den, prec)

    col_a_m0 = _columns(m0, TARGET_A)
    col_b_m0 = _columns(m0, TARGET_B)
    col_a_m2 = _columns(m2, TARGET_A)
    col_b_m2 = _columns(m2, TARGET_B)

    r0_aa = _residual_entry(
        s0[TARGET_A, TARGET_A], col_a_m0, col_a_m0, gram_diagonal,
    )
    r0_ab = _residual_entry(
        s0[TARGET_A, TARGET_B], col_a_m0, col_b_m0, gram_diagonal,
    )
    r0_bb = _residual_entry(
        s0[TARGET_B, TARGET_B], col_b_m0, col_b_m0, gram_diagonal,
    )
    r2_aa = _residual_entry(
        s2[TARGET_A, TARGET_A], col_a_m2, col_a_m2, gram_diagonal,
    )
    r2_ab = _residual_entry(
        s2[TARGET_A, TARGET_B], col_a_m2, col_b_m2, gram_diagonal,
    )
    r2_bb = _residual_entry(
        s2[TARGET_B, TARGET_B], col_b_m2, col_b_m2, gram_diagonal,
    )

    c_l = cL_iv(L_num, L_den, prec)
    b_l = b_L_iv(2 * N + 1, L_num, L_den, prec)

    def c_entry(a: int, b: int, r0: Interval, r2: Interval) -> Interval:
        diagonal = a == b
        f_value = add(m0[a, b], m2[a, b])
        if diagonal:
            t_value = point(_harmonic(indices[a]) * gram_diagonal[a])
            f_value = add(add(t_value, f_value), scalar_mul(Fraction(-1), mul(c_l, point(gram_diagonal[a]))))
        residual = add(
            scalar_mul(Fraction(1) + ETA, r0),
            scalar_mul(Fraction(1) + Fraction(1) / ETA, r2),
        )
        return sub(mul(b_l, f_value), residual)

    c_aa = c_entry(TARGET_A, TARGET_A, r0_aa, r2_aa)
    c_ab = c_entry(TARGET_A, TARGET_B, r0_ab, r2_ab)
    c_bb = c_entry(TARGET_B, TARGET_B, r0_bb, r2_bb)
    numerator = rayleigh_numerator(c_aa, c_ab, c_bb)
    certified_negative = numerator[1] < 0
    if certified_negative:
        witness_status = "negative"
    elif numerator[0] >= 0:
        witness_status = "nonnegative"
    else:
        witness_status = "indeterminate"

    result = {
        "format_version": "frontier-rayleigh-pilot-1.0",
        "method": "two_column_schur_rayleigh_v1",
        "L": f"{L_num}/{L_den}",
        "sector": "odd",
        "N": N,
        "d": 2 * N + 1,
        "eta": str(ETA),
        "target_degrees": [indices[TARGET_A], indices[TARGET_B]],
        "witness": list(WITNESS),
        "remainder_mode": "richardson_gl8_gl4",
        "rayleigh_witness_status": witness_status,
        "rayleigh_numerator_interval": [
            str(numerator[0]), str(numerator[1]),
        ],
        "checkpoint_m0_entries": len(state["M0"]),
        "checkpoint_s0_entries": len(state["S0"]),
        "certified_not_positive_definite": certified_negative,
        "resume_assembly_elapsed_s": round(time.time() - started, 1),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    _save_json(out_path, result)
    print(
        f"v^T C v ∈ [{float(numerator[0]):.9e}, {float(numerator[1]):.9e}]",
        flush=True,
    )
    if certified_negative:
        print("CERTIFIED NOT POSITIVE DEFINITE: interval upper endpoint < 0", flush=True)
        ckpt_path.unlink(missing_ok=True)
    else:
        print(
            f"NOT CERTIFIED NEGATIVE: Rayleigh witness interval is {witness_status}",
            flush=True,
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", nargs=2, type=int, default=[56, 100], metavar=("NUM", "DEN"))
    parser.add_argument("--depth2", type=int, default=4)
    parser.add_argument("--depth3", type=int, default=3)
    parser.add_argument("--prec", type=int, default=512)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    try:
        result = certify(
            args.L[0], args.L[1], args.depth2, args.depth3,
            args.prec, args.out, args.resume,
        )
    except KeyboardInterrupt:
        print("\n[interrupt] rerun with --resume", flush=True)
        sys.exit(130)
    if not result["certified_not_positive_definite"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
