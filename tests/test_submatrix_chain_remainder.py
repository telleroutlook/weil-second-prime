"""Safety tests for the submatrix-chain pilot."""
from __future__ import annotations

import ast
from pathlib import Path


def test_submatrix_chain_does_not_skip_mk_remainder() -> None:
    """Raw GL8 values produced a false k=28 sign; do not re-enable that path."""
    tree = ast.parse(
        Path("scripts/submatrix_chain.py").read_text(), filename="submatrix_chain.py"
    )
    calls = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "integrate_M_K"
    ]
    assert calls
    for call in calls:
        skipped = next(
            (kw.value for kw in call.keywords if kw.arg == "skip_remainder"),
            None,
        )
        assert not (
            isinstance(skipped, ast.Constant) and skipped.value is True
        )
