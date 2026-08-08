"""Pytest wrapper for the second-prime cross-structure mutation catalog (C11).

Ensures the checker (a) certifies the baseline cert and (b) kills 100% of the
mutation catalog (zeroing tau_2/tau_3, dropping the cross term, S_KK-only, both
window-bound violations, wrong primes/method). This is the C11 obligation for
the second window and the CLAUDE.md two-prime negative-test requirement.
"""
from __future__ import annotations

import pytest

from tests.mutation.mutation_catalog_second_cross import run_catalog


@pytest.fixture(scope="module")
def catalog():
    return run_catalog()


def test_baseline_certifies(catalog):
    assert catalog["baseline_certifies"] is True


def test_kill_rate_100(catalog):
    assert catalog["kill_rate"] == 1.0, catalog["results"]
    assert catalog["kill_rate_pct"] == "100%"


def test_all_mutants_killed(catalog):
    unkilled = [r["mutant"] for r in catalog["results"] if not r["killed"]]
    assert not unkilled, f"unkilled mutants: {unkilled}"


@pytest.mark.parametrize("required_mutant", [
    "zero_tau2_shift",
    "zero_tau3_shift",
    "drop_cross_term_F",
    "skk_only_S0",
    "window_out_high_L>log2",
    "window_out_low_L<half_log3",
])
def test_required_two_prime_mutants_present_and_killed(catalog, required_mutant):
    rec = next((r for r in catalog["results"] if r["mutant"] == required_mutant), None)
    assert rec is not None, f"missing required mutant {required_mutant}"
    assert rec["killed"], f"{required_mutant} not killed"
