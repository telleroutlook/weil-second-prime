"""Archimedean base primitive checker.

Exit codes:
    0  CERTIFIED — all primitive integrals verified by both paths
    1  uncertified — some interval check failed
    2  malformed certificate or resource error
    3  O2_BLOCKED — required dependencies unavailable

This checker independently recomputes M_V, M_K, S_VV, S_VK, S_KK from first
principles using two independent integration paths and verifies their intersection.

P0 fixes applied (from weil-lower-bound audit):
  - integrate_M_K now calls _integrate_1d_arb with GL-8/GL-4 remainder (not raw GL-8)
  - Near-zero Taylor cubic coefficient is 7s^3/11520 (not s^3/2880)
  - All remainders use Bernstein ellipse analytic bounds (not empirical convergence)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_ROOT = Path(__file__).parent.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

SCHEMA = _ROOT / "schemas" / "certificate-archimedean-v1.schema.json"
MAX_CONTRACT_BYTES = 10 * 1024 * 1024


def _load_and_validate(path: Path) -> dict:
    raw = path.read_bytes()
    if len(raw) > MAX_CONTRACT_BYTES:
        raise ValueError("archimedean contract exceeds 1 MiB limit")
    try:
        value = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid JSON: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError("contract root must be an object")

    try:
        import jsonschema
    except ImportError as exc:
        raise RuntimeError("jsonschema is required for strict validation") from exc

    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    validator = jsonschema.Draft7Validator(schema)
    errors = sorted(validator.iter_errors(value), key=lambda e: list(e.path))
    if errors:
        msgs = [
            f"{'/'.join(str(p) for p in e.path) or '<root>'}: {e.message}"
            for e in errors
        ]
        raise ValueError("; ".join(msgs))
    return value


def _check_leaf_witnesses(contract: dict) -> list[str]:
    """Verify Path A leaf witnesses from mk_entries when present.

    For each leaf, checks that enclosure[lower] <= enclosure[upper] and
    that the remainder bound is non-negative. Returns list of failure keys.
    """
    from fractions import Fraction

    failures = []
    for entry in contract.get("mk_entries", []):
        for i, leaf in enumerate(entry.get("leaf_witnesses", [])):
            key = f"M_K[{entry['n_row']},{entry['n_col']}].leaf[{i}]"
            try:
                enc = leaf.get("enclosure", {})
                lo = Fraction(enc.get("lower", "0"))
                hi = Fraction(enc.get("upper", "0"))
                rem = Fraction(leaf.get("remainder", {}).get("bound", "0"))
                if lo > hi:
                    failures.append(f"{key}: enclosure inverted lo={lo} > hi={hi}")
                if rem < 0:
                    failures.append(f"{key}: negative remainder {rem}")
            except (ValueError, ZeroDivisionError) as exc:
                failures.append(f"{key}: parse error {exc}")
    return failures


def check(args: argparse.Namespace) -> int:
    try:
        contract = _load_and_validate(args.contract)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ARCHIMEDEAN CHECKER REJECT: {exc}", file=sys.stderr)
        return 2

    try:
        from checker.archimedean.integrate import (
            compute_all_primitives_path_a,
            compute_all_primitives_path_b,
            verify_intersection,
        )
    except ImportError as exc:
        output = {
            "status": "O2_BLOCKED",
            "obligation": "archimedean_primitives_o2_v1",
            "reason": f"integration module unavailable: {exc}",
        }
        print(json.dumps(output, sort_keys=True))
        return 3

    try:
        primitives_a = compute_all_primitives_path_a(contract, precision=256)
        primitives_b = compute_all_primitives_path_b(contract, precision=256)
        verified = verify_intersection(primitives_a, primitives_b)
    except Exception as exc:  # noqa: BLE001
        print(f"ARCHIMEDEAN CHECKER ERROR: {exc}", file=sys.stderr)
        return 1

    if not verified["all_pass"]:
        failing = [k for k, v in verified["checks"].items() if not v]
        print(
            f"ARCHIMEDEAN CHECKER UNCERTIFIED: failing checks: {failing}",
            file=sys.stderr,
        )
        return 1

    # Verify leaf witnesses when present in the certificate
    leaf_failures = _check_leaf_witnesses(contract)
    if leaf_failures:
        print(
            f"ARCHIMEDEAN CHECKER LEAF FAIL: {len(leaf_failures)} leaf witness failures: "
            f"{leaf_failures[:3]}",
            file=sys.stderr,
        )
        return 1

    # Serialise Fraction-valued intervals as [num/den, num/den] strings
    def _ser_iv(iv: object) -> list[str]:
        lo, hi = iv  # type: ignore[misc]
        return [f"{lo.numerator}/{lo.denominator}", f"{hi.numerator}/{hi.denominator}"]

    output = {
        "protocol_version": 2,
        "claim_id": "",
        "obligation_results": [
            {"id": "archimedean_primitives_o2_v1", "verdict": "pass"}
        ],
        "status": "CERTIFIED",
        "obligation": "archimedean_primitives_o2_v1",
        "primitives": {k: _ser_iv(v) for k, v in verified["primitives"].items()},
        "checks": verified["checks"],
    }
    print(json.dumps(output, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Archimedean base primitive checker for L=7/20."
    )
    parser.add_argument("contract", type=Path, help="Archimedean certificate JSON")
    parser.add_argument("--schema", type=Path, default=None,
                        help="Override schema path (default: schemas/certificate-archimedean-v1.schema.json)")
    return parser


def main() -> int:
    return check(build_parser().parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
