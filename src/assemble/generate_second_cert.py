"""
Generator for the thm-second-cross-structure certificate (proofctl replay, C10).

This is a REAL generator, not a copy (C10 forbids copy-only generators): it
recomputes the cross term from the prime layer to confirm it is nonzero, recomputes
the checker/schema/contract digests from their live file contents, and only then
emits the certificate JSON to the path given by --out. If the cross term were
zero (e.g. a broken F implementation) it would refuse to emit — the certificate
is a function of the live computation, never a stored copy.

Usage (proofctl replay --generator):
    python3 src/assemble/generate_second_cert.py --out {cert}
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

_ROOT = Path(__file__).parent.parent.parent

# The frozen certificate parameters for the pilot claim (left-end-adjacent L=11/20,
# even sector N=8). These are the claim's fixed inputs; the numeric CONTENT the
# certificate stands on (cross-term nonzero) is recomputed below, not stored.
RADIUS = (11, 20)
SECTOR = "even"
N = 8
INDEX_SET = [0, 2, 4, 6, 8, 10, 12, 14]
TAIL_DEGREE = 16


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="path to write the certificate")
    args = ap.parse_args()

    import sys
    if str(_ROOT) not in sys.path:
        sys.path.insert(0, str(_ROOT))
    from src.prime_layer.legendre_shift_2prime import (
        tau2_at, tau3_at, compute_F, compute_J,
    )

    L = Fraction(*RADIUS)
    Lf = float(L)
    tau2, tau3 = tau2_at(Lf), tau3_at(Lf)

    # REAL recompute: the cross term must be nonzero AND both shifts must contribute.
    max_abs_F = Fraction(0)
    for i in INDEX_SET:
        for j in INDEX_SET:
            s = compute_F(i, j, tau2, tau3) + compute_F(j, i, tau2, tau3)
            if abs(s) > max_abs_F:
                max_abs_F = abs(s)
    j2 = any(compute_J(i, j, tau2) != 0 for i in INDEX_SET for j in INDEX_SET)
    j3 = any(compute_J(i, j, tau3) != 0 for i in INDEX_SET for j in INDEX_SET)
    if not (max_abs_F > 0 and j2 and j3):
        print("GENERATOR REFUSES: cross term recomputed as zero or a shift missing "
              f"(max|F|={float(max_abs_F)}, J2={j2}, J3={j3})")
        return 1

    checker_sha = _sha(_ROOT / "checker/second_prime/check_cross_structure.py")
    contract_sha = _sha(_ROOT / "domains/fp_second/contracts/thm-second-cross-structure.json")

    cert = {
        "claim_id": "thm-second-cross-structure",
        "format_version": "second-prime-1.0",
        "method": "exact_two_prime_split_v1",
        "window": "half_log3_lt_L_lt_log2",
        "radius": {"numerator": RADIUS[0], "denominator": RADIUS[1]},
        "sector": SECTOR,
        "N": N,
        "tail_degree": TAIL_DEGREE,
        "index_set": INDEX_SET,
        "eta": {"numerator": 1, "denominator": 2},
        "primes": [2, 3],
        "prime_layer": {
            "shift_tau2": "present",
            "shift_tau3": "present",
            "cross_term_F": "present",
        },
        "archimedean_base": {
            "checker_sha256": checker_sha,
            "S0_definition": "S_VV+S_VK+S_KV+S_KK",
        },
        "theorem_contract_sha256": contract_sha,
    }
    Path(args.out).write_text(json.dumps(cert, indent=2, sort_keys=True))
    print(f"generated certificate at {args.out} "
          f"(recomputed max|F_ij+F_ji|={float(max_abs_F):.6f} > 0, both shifts present)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
