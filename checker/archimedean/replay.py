"""Archimedean base primitive replay for the first-prime checker."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


class O2Blocked(Exception):
    """Raised when the O2 primitive replay is unavailable or incomplete."""


def sha256_file(path: Path) -> str:
    """Return lowercase hex SHA256 of a file's bytes."""
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return digest


def replay_archimedean_base(
    archimedean_base: dict[str, Any],
    base_certificate: Path,
    base_checker: Path,
    base_schema: Path,
) -> dict[str, Any]:
    """Replay the Archimedean base certificate and return verified primitive matrices.

    Verifies SHA256 digests for certificate, checker, and schema before invocation.
    Raises O2Blocked if the base checker is unavailable or fails.
    """
    cert_digest = sha256_file(base_certificate)
    if cert_digest != archimedean_base["certificate_sha256"]:
        raise ValueError(
            f"archimedean_base.certificate_sha256 mismatch: "
            f"expected {archimedean_base['certificate_sha256']}, got {cert_digest}"
        )

    checker_digest = sha256_file(base_checker)
    if checker_digest != archimedean_base["checker_sha256"]:
        raise ValueError(
            f"archimedean_base.checker_sha256 mismatch: "
            f"expected {archimedean_base['checker_sha256']}, got {checker_digest}"
        )

    schema_digest = sha256_file(base_schema)
    if schema_digest != archimedean_base["schema_sha256"]:
        raise ValueError(
            f"archimedean_base.schema_sha256 mismatch: "
            f"expected {archimedean_base['schema_sha256']}, got {schema_digest}"
        )

    if archimedean_base["obligation"] != "archimedean_primitives_o2_v1":
        raise ValueError(
            f"unexpected obligation: {archimedean_base['obligation']!r}"
        )

    # Invoke the base checker as a subprocess; it must exit 0 and emit JSON.
    try:
        result = subprocess.run(
            [
                sys.executable,
                str(base_checker),
                str(base_certificate),
                "--schema", str(base_schema),
            ],
            capture_output=True,
            text=True,
            timeout=3600,
        )
    except FileNotFoundError as exc:
        raise O2Blocked(f"base checker not found: {base_checker}") from exc
    except subprocess.TimeoutExpired as exc:
        raise O2Blocked("base checker timed out after 600 seconds") from exc

    if result.returncode == 3:
        raise O2Blocked(f"base checker O2_BLOCKED: {result.stdout.strip()}")
    if result.returncode != 0:
        raise O2Blocked(
            f"base checker failed (exit {result.returncode}): {result.stderr.strip()}"
        )

    try:
        output = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise O2Blocked(f"base checker produced invalid JSON: {exc}") from exc

    if output.get("status") != "CERTIFIED":
        raise O2Blocked(
            f"base checker did not return CERTIFIED: status={output.get('status')!r}"
        )

    # Return the verified primitive blocks for use in the first-prime checker
    primitives = output.get("primitives")
    if not primitives:
        raise O2Blocked("base checker output missing 'primitives' block")

    return primitives
