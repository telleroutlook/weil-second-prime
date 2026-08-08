"""Checker for thm-second-cross-structure — second-prime-window cross term.

Certifies the HONEST, certify-true result of the second window (S4): the
cross-prime term F_{ij}(tau_2, tau_3) = <C_{tau_3,1} P_j, C_{tau_2,1} P_i> is a
REAL, NONZERO contribution to the two-prime Schur second moment — the genuine new
structure the first window does not contain. This is metric A of S4
(sign-independent): Delta C = C_full - C_cross_off = -3 c_2 c_3 (F_ij + F_ji)
exactly, because M^{(2)} carries no cross term. So the cross-term magnitude is
PURELY the prime layer and is recomputed here in exact Fraction arithmetic
(fast, no archimedean integrals).

This checker does NOT claim positivity. The second window is NOT certified
positive-definite: at the left end L=11/20 the odd sector is certifiably negative
and the even sector straddles zero (see pilots/s5_positivity_L055.json). Those are
separate facts; this certificate is bounded to the cross-term new-structure
finding and to finite-scale scope (L < log2, no RH).

Obligations:
  1. second.cross-term-present-and-nonzero  (recompute max|F_ij+F_ji| > 0)
  2. second.both-shifts-present             (tau_2 AND tau_3 both contribute)
  3. second.window-bounds-hold              (1/2 log3 < L < log2 by certified rationals)
  4. second.four-term-S0-declared           (S_VV+S_VK+S_KV+S_KK, not S_KK-only)
  5. second.positivity-not-claimed          (this cert does not assert PD)
  6. second.conclusion-bounded-and-no-rh    (finite-scale, L<log2, no RH)

Exit codes: 0 CERTIFIED, 1 uncertified, 2 malformed/resource, 3 blocked.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

_ROOT = Path(__file__).parent.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

OBLIGATION_IDS = [
    "second.cross-term-present-and-nonzero",
    "second.both-shifts-present",
    "second.window-bounds-hold",
    "second.four-term-S0-declared",
    "second.positivity-not-claimed",
    "second.conclusion-bounded-and-no-rh",
]

# Certified rational bounds on log 2 and log 3, computed from the arctanh series
# (stdlib Fraction, auditable, no float). arctanh(x) = sum_k x^(2k+1)/(2k+1).
#   log2 = 2 arctanh(1/3):  (1+1/3)/(1-1/3) = 2.
#   log3 = 2 arctanh(1/2):  (1+1/2)/(1-1/2) = 3.
# Partial sum (K terms, all positive) is a LOWER bound; add a geometric tail
# bound for the UPPER bound. The tail after K terms of 2*sum x^(2k+1)/(2k+1) is
#   < 2/(2K+1) * x^(2K+1) / (1 - x^2)   (bounding 1/(2k+1) <= 1/(2K+1)).

def _log_bounds(inv_x_den: int, K: int = 24) -> tuple[Fraction, Fraction]:
    """Certified (lower, upper) for 2*arctanh(1/inv_x_den) via K-term series."""
    x = Fraction(1, inv_x_den)
    S = Fraction(0)
    for k in range(K):
        S += x ** (2 * k + 1) / (2 * k + 1)
    lo = 2 * S
    tail = Fraction(2, 2 * K + 1) * x ** (2 * K + 1) / (1 - x * x)
    hi = 2 * S + 2 * tail  # doubling the tail bound for safety margin
    return lo, hi


_LOG2_LO, _LOG2_HI = _log_bounds(3)   # log2 = 2 arctanh(1/3)
_LOG3_LO, _LOG3_HI = _log_bounds(2)   # log3 = 2 arctanh(1/2)


def _window_ok(L: Fraction) -> tuple[bool, str]:
    """Verify 1/2 log3 < L < log2 by certified rational bounds (not enumeration)."""
    # left:  1/2 log3 < L  <=>  log3 < 2L. Certified if UPPER bound on log3 < 2L.
    left = _LOG3_HI < 2 * L
    # right: L < log2. Certified if L < LOWER bound on log2.
    right = L < _LOG2_LO
    ok = left and right
    msg = (f"window (1/2 log3, log2): L={L}; "
           f"left(log3<2L via log3<{float(_LOG3_HI):.8f}<{float(2*L):.6f})={left}, "
           f"right(L<log2 via {float(L):.6f}<{float(_LOG2_LO):.8f}<log2)={right}")
    return ok, msg


def _recompute_cross(index_set: list[int], L: Fraction) -> tuple[Fraction, bool]:
    """Recompute max|F_ij + F_ji| over the index set (exact Fraction). Returns
    (max_abs, all_shifts_contribute)."""
    from src.prime_layer.legendre_shift_2prime import tau2_at, tau3_at, compute_F, compute_J
    Lf = float(L)
    tau2 = tau2_at(Lf)
    tau3 = tau3_at(Lf)
    max_abs = Fraction(0)
    for i in index_set:
        for j in index_set:
            s = compute_F(i, j, tau2, tau3) + compute_F(j, i, tau2, tau3)
            if abs(s) > max_abs:
                max_abs = abs(s)
    # both shifts present: some J(tau2) and some J(tau3) nonzero on the set
    j2 = any(compute_J(i, j, tau2) != 0 for i in index_set for j in index_set)
    j3 = any(compute_J(i, j, tau3) != 0 for i in index_set for j in index_set)
    return max_abs, (j2 and j3)


def verify(cert: dict) -> tuple[bool, list[bool], str]:
    msgs = []
    radius = cert.get("radius", {})
    L = Fraction(int(radius["numerator"]), int(radius["denominator"]))
    sector = cert.get("sector")
    index_set = cert.get("index_set", [])

    # 1 + 2: recompute the cross term (exact) and both-shift presence.
    max_abs, both_shifts = _recompute_cross(index_set, L)
    o1 = max_abs > 0
    msgs.append(f"cross-term recompute: max|F_ij+F_ji|={float(max_abs):.6f} (>0: {o1})")
    o2 = both_shifts and cert.get("prime_layer", {}).get("shift_tau2") == "present" \
        and cert.get("prime_layer", {}).get("shift_tau3") == "present"
    msgs.append(f"both shifts present (recomputed J(tau2),J(tau3) nonzero AND declared): {o2}")

    # 3: window bounds by certified rationals.
    o3, wmsg = _window_ok(L)
    msgs.append(wmsg)

    # 4: four-term S0 declared (not S_KK-only).
    o4 = cert.get("archimedean_base", {}).get("S0_definition") == "S_VV+S_VK+S_KV+S_KK" \
        and cert.get("prime_layer", {}).get("cross_term_F") == "present"
    msgs.append(f"four-term S0 + cross_term_F declared present: {o4}")

    # 5: positivity NOT claimed by this cert (no positivity/pivot/conclusion fields).
    forbidden = {"min_pivot", "min_eig", "positive_definite", "conclusion",
                 "lambda", "pivots", "eigenvalues"}
    o5 = not (forbidden & set(cert.keys()))
    msgs.append(f"positivity not self-claimed (no forbidden result fields): {o5}")

    # 6: scope bounded — window const is L<log2, method const, primes=[2,3].
    o6 = (cert.get("window") == "half_log3_lt_L_lt_log2"
          and cert.get("method") == "exact_two_prime_split_v1"
          and cert.get("primes") == [2, 3])
    msgs.append(f"scope bounded (L<log2 window const, two-prime method, primes=[2,3]): {o6}; "
                f"no RH implied (finite-scale second window only)")

    results = [o1, o2, o3, o4, o5, o6]
    return all(results), results, "; ".join(msgs)


def _mutation_metadata() -> dict:
    """Expose C11 fields from the committed mutation-catalog artifact.

    References a pre-run, auditable artifact (not re-running the catalog on every
    check). Absent/invalid artifact -> empty dict, so C11 correctly BLOCKS rather
    than silently claiming coverage. All values are STRINGS: proofctl replay
    unmarshals attestation metadata into map[string]string, and a non-string
    value makes the WHOLE metadata map fail to unmarshal (silently dropping every
    key). Learned 2026-08-08 from a release gate that saw metadata:null."""
    art = _ROOT / "pilots" / "mutation_catalog_second_cross.json"
    try:
        d = json.loads(art.read_text())
        if not d.get("baseline_certifies") or d.get("kill_rate") != 1.0:
            return {}
        return {
            "mutation_kill_rate": str(d.get("kill_rate_pct", "")),
            "mutation_catalog_digest": str(d.get("catalog_digest", "")),
        }
    except Exception:
        return {}


def main() -> int:
    # cert path from argv[1] (CLI) or stdin bridge
    cert_path = None
    claim_id = "thm-second-cross-structure"
    if len(sys.argv) >= 2:
        cert_path = Path(sys.argv[1])
    else:
        try:
            inp = json.load(sys.stdin)
            claim_id = inp.get("claim_id", claim_id)
            for ev in inp.get("evidence", []):
                hint = ev.get("local_path", "") or ev.get("path_hint", "")
                if hint and Path(hint).exists():
                    cert_path = Path(hint)
                    break
        except json.JSONDecodeError as exc:
            print(f"CHECKER PROTOCOL ERROR: {exc}", file=sys.stderr)
            return 2

    if cert_path is None or not cert_path.exists():
        print("CHECKER ERROR: certificate not found", file=sys.stderr)
        return 2
    try:
        cert = json.loads(cert_path.read_text())
    except Exception as exc:
        print(f"CHECKER ERROR: malformed certificate: {exc}", file=sys.stderr)
        return 2
    claim_id = cert.get("claim_id", claim_id)

    passed, results, explanation = verify(cert)
    out = {
        "protocol_version": 2,
        "claim_id": claim_id,
        "obligation_results": [
            {"id": oid, "verdict": "pass" if r else "fail"}
            for oid, r in zip(OBLIGATION_IDS, results)
        ],
        "status": "CERTIFIED" if passed else "UNCERTIFIED",
        "explanation": explanation,
        "metadata": {
            "format_version": "second-prime-1.0",
            "method": "exact_two_prime_split_v1",
            "grade": "certify (exact Fraction cross-term recompute)",
            "positivity_claimed": "false",
            "scope": "finite-scale second-window cross-term structure; L<log2; no RH",
            **_mutation_metadata(),
        },
    }
    print(json.dumps(out, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
