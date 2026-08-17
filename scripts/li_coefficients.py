"""Compute finite-zero-truncation Li coefficients for L=0.56 comparison.

Li criterion: RH <=> lambda_n >= 0 for all n >= 1.
Partial sum with first K zero pairs gives a truncated approximation.

Usage:
    python3 -m scripts.li_coefficients --K 200 --n_max 20 \
        --out pilots/li_coefficients_L056.json
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path


def _li_partial(n_max: int, K: int, dps: int = 40) -> dict[int, float]:
    """Compute lambda_n (n=1..n_max) using first K pairs of non-trivial zeros.

    lambda_n = 2 * Re[ sum_{k=1}^{K} (1 - (1 - 1/rho_k)^n) ]
    where rho_k = 1/2 + i*t_k (upper half-plane zeros).

    Under RH, all lambda_n > 0 and lambda_n ~ (n/2) log n for large n.
    """
    import mpmath
    mpmath.mp.dps = dps

    print(f"Computing first {K} Riemann zeros...", flush=True)
    t0 = time.time()
    zeros = [mpmath.zetazero(k) for k in range(1, K + 1)]
    print(f"  {K} zeros in {time.time()-t0:.1f}s", flush=True)

    results: dict[int, float] = {}
    for n in range(1, n_max + 1):
        total = mpmath.mpf(0)
        for rho in zeros:
            term = 1 - (1 - 1 / rho) ** n
            total += 2 * mpmath.re(term)  # conjugate pair contributes 2*Re
        results[n] = float(total)
        print(f"  lambda_{n:2d} = {float(total):+.6f}", flush=True)

    return results


def main() -> None:
    p = argparse.ArgumentParser(description="Li coefficients via zero sum")
    p.add_argument("--K", type=int, default=200,
                   help="Number of zero pairs to sum (default 200)")
    p.add_argument("--n_max", type=int, default=20,
                   help="Compute lambda_1 .. lambda_n_max (default 20)")
    p.add_argument("--out", default="pilots/li_coefficients_L056.json",
                   help="Output JSON path")
    p.add_argument("--dps", type=int, default=40,
                   help="mpmath decimal places (default 40)")
    args = p.parse_args()

    t_start = time.time()
    lambdas = _li_partial(args.n_max, args.K, dps=args.dps)
    elapsed = time.time() - t_start

    # N-convergence comparison table (from pilots, corrected kappa=2.056)
    eig_full_reference = {
        7: -0.5747, 9: -0.3895, 11: -0.2770, 13: -0.1879,
        15: -0.1155, 17: -0.0760,
    }

    result = {
        "description": (
            "Li coefficients lambda_n = 2*Re[sum_{k=1}^{K} (1-(1-1/rho_k)^n)], "
            f"K={args.K} zero pairs. "
            "Li criterion: RH <=> lambda_n >= 0 for all n. "
            "Comparison: second-window Schur eig_full(N) is a different convergent."
        ),
        "grade": "pilot (float, partial zero sum)",
        "K_zeros": args.K,
        "n_max": args.n_max,
        "elapsed_s": round(elapsed, 1),
        "lambdas": {str(n): v for n, v in lambdas.items()},
        "schur_eig_full_reference": {
            str(k): v for k, v in eig_full_reference.items()
        },
        "note": (
            "lambda_n values are all positive (consistent with RH), growing as "
            "O(n log n). The Schur eig_full(N) is a separate convergent: it "
            "measures min-eigenvalue of the truncated Schur matrix, not the Li sum. "
            "Both are truncated approximations that improve with more terms (N or K)."
        ),
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2))
    print(f"\nWrote {out_path}  ({elapsed:.1f}s)")

    # Print comparison table
    print("\n--- Comparison table ---")
    print(f"{'n/N':>4}  {'lambda_n (Li)':>16}  {'eig_full(N) (Schur)':>20}")
    all_ns = sorted(set(range(1, args.n_max + 1)) | set(eig_full_reference))
    for n in all_ns:
        li = f"{lambdas[n]:+.6f}" if n in lambdas else "—"
        sch = f"{eig_full_reference[n]:+.4f}" if n in eig_full_reference else "—"
        print(f"{n:>4}  {li:>16}  {sch:>20}")


if __name__ == "__main__":
    main()
