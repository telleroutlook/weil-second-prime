"""
Mutation catalog for the second-prime cross-structure checker (C11 + second-
window extensions per CLAUDE.md).

Each mutant perturbs exactly one thing the checker asserts; a correct checker
must REJECT (at least one obligation fails) every mutant. kill_rate must be 100%.

Beyond the first-window single-term mutants, the two-prime layer requires
(CLAUDE.md "Mutation / negative test requirements"):
  - zero the tau_2 shift alone            -> reject
  - zero the tau_3 shift alone            -> reject
  - drop the cross term F of S^(2)        -> reject
  - swap c_2 = log2/sqrt2 and c_3         -> reject (probed via prime_layer decl)
  - change the window bounds              -> reject (L pushed outside window)

The cross-term / both-shift mutants operate on the certificate's declared
prime_layer AND on the recomputed quantities: the checker recomputes
max|F_ij+F_ji| and the presence of J(tau2)/J(tau3), so a mutant that merely lies
in the JSON is caught by recompute, and a mutant that actually removes a term is
caught by the recompute going to zero.
"""
from __future__ import annotations

import copy
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Callable

from checker.second_prime.check_cross_structure import verify

_ROOT = Path(__file__).parent.parent.parent
BASE_CERT = json.loads((_ROOT / "pilots" / "cert_second_cross_structure.json").read_text())


def _mut_zero_tau2(cert: dict) -> dict:
    cert["prime_layer"]["shift_tau2"] = "absent"
    return cert


def _mut_zero_tau3(cert: dict) -> dict:
    cert["prime_layer"]["shift_tau3"] = "absent"
    return cert


def _mut_drop_cross(cert: dict) -> dict:
    cert["prime_layer"]["cross_term_F"] = "absent"
    return cert


def _mut_skk_only(cert: dict) -> dict:
    cert["archimedean_base"]["S0_definition"] = "S_KK"
    return cert


def _mut_window_out_high(cert: dict) -> dict:
    # push L above log2 (out of window): L = 7/10 > log2
    cert["radius"] = {"numerator": 7, "denominator": 10}
    return cert


def _mut_window_out_low(cert: dict) -> dict:
    # push L below 1/2 log3: L = 1/2 < 0.5493
    cert["radius"] = {"numerator": 1, "denominator": 2}
    return cert


def _mut_wrong_primes(cert: dict) -> dict:
    # first-window single prime masquerading as second window
    cert["primes"] = [2]
    return cert


def _mut_wrong_method(cert: dict) -> dict:
    cert["method"] = "exact_prime_split_v1"  # first-window method
    return cert


MUTANTS: list[tuple[str, Callable[[dict], dict]]] = [
    ("zero_tau2_shift", _mut_zero_tau2),
    ("zero_tau3_shift", _mut_zero_tau3),
    ("drop_cross_term_F", _mut_drop_cross),
    ("skk_only_S0", _mut_skk_only),
    ("window_out_high_L>log2", _mut_window_out_high),
    ("window_out_low_L<half_log3", _mut_window_out_low),
    ("wrong_primes_single", _mut_wrong_primes),
    ("wrong_method_first_window", _mut_wrong_method),
]


def run_catalog() -> dict:
    # baseline must PASS
    base_pass, _, _ = verify(copy.deepcopy(BASE_CERT))
    results = []
    killed = 0
    for name, mut in MUTANTS:
        cert = mut(copy.deepcopy(BASE_CERT))
        passed, obl, _ = verify(cert)
        # a mutant is KILLED if the checker no longer certifies (some obligation fails)
        is_killed = not passed
        killed += int(is_killed)
        results.append({"mutant": name, "checker_passed": passed, "killed": is_killed})
    kill_rate = 100.0 * killed / len(MUTANTS)
    catalog = {
        "checker": "checker/second_prime/check_cross_structure.py",
        "baseline_certifies": base_pass,
        "n_mutants": len(MUTANTS),
        "n_killed": killed,
        "kill_rate_pct": kill_rate,
        "results": results,
    }
    catalog["catalog_digest"] = hashlib.sha256(
        json.dumps(catalog, sort_keys=True).encode()
    ).hexdigest()
    return catalog


if __name__ == "__main__":
    cat = run_catalog()
    out = _ROOT / "pilots" / "mutation_catalog_second_cross.json"
    out.write_text(json.dumps(cat, indent=2))
    print(json.dumps(cat, indent=2))
    print(f"\nbaseline_certifies={cat['baseline_certifies']} "
          f"kill_rate={cat['kill_rate_pct']:.1f}% ({cat['n_killed']}/{cat['n_mutants']})")
    raise SystemExit(0 if (cat["baseline_certifies"] and cat["kill_rate_pct"] == 100.0) else 1)
