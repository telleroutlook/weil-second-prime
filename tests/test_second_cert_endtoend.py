"""End-to-end test: the pilot second-window certificate is schema-valid and the
checker certifies it (exit 0 / all obligations pass). Also pins the honest scope:
the certificate must NOT carry positivity/pivot/conclusion fields."""
from __future__ import annotations

import json
from pathlib import Path

import jsonschema
import pytest

_ROOT = Path(__file__).parent.parent
CERT = _ROOT / "pilots" / "cert_second_cross_structure.json"
SCHEMA = _ROOT / "schemas" / "certificate-second-prime-v1.schema.json"


@pytest.fixture(scope="module")
def cert():
    return json.loads(CERT.read_text())


def test_cert_schema_valid(cert):
    jsonschema.validate(cert, json.loads(SCHEMA.read_text()))


def test_cert_carries_no_positivity_fields(cert):
    forbidden = {"min_pivot", "min_eig", "positive_definite", "conclusion",
                 "lambda", "pivots", "eigenvalues"}
    assert not (forbidden & set(cert.keys()))


def test_cert_declares_both_shifts_and_cross(cert):
    pl = cert["prime_layer"]
    assert pl["shift_tau2"] == "present"
    assert pl["shift_tau3"] == "present"
    assert pl["cross_term_F"] == "present"


def test_cert_four_term_s0(cert):
    assert cert["archimedean_base"]["S0_definition"] == "S_VV+S_VK+S_KV+S_KK"


def test_checker_certifies(cert):
    from checker.second_prime.check_cross_structure import verify
    passed, results, _ = verify(cert)
    assert passed, results
