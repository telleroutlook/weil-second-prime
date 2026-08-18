"""Regression tests for the pilot B0 branch fitter."""
from __future__ import annotations

import numpy as np
from pytest import approx

from scripts.fit_b0 import fit_mode, load_modes


def test_fit_mode_scales_ir_branch() -> None:
    """The optimizer must recover an IR-scale asymptote, not its initial guess."""
    ks = np.arange(25.0, 31.0)
    expected = np.array([-2.0e-6, 0.65, 1.5e-6])
    lams = expected[0] * expected[1] ** ks + expected[2]

    popt, rms, count = fit_mode(
        ks, lams, 25, -1, (0.0, 1.0e-5), "IR"
    )

    assert count == len(ks)
    assert rms < 1.0e-14
    assert float(popt[1]) == approx(0.65, rel=1.0e-3)
    assert float(popt[2]) == approx(1.5e-6, abs=1.0e-11)


def test_load_modes_separates_k28_frontier_mode() -> None:
    """k=28's negative λ₀ must not be merged into the continuing UV branch."""
    uv_k, uv_lam, ir_k, ir_lam, frontier_k, frontier_lam = load_modes()

    assert uv_k[-1] == 28
    assert ir_k[-1] == 28
    assert frontier_k.tolist() == [28]
    assert float(uv_lam[-1]) == approx(0.012972974484854764, abs=1e-12)
    assert float(ir_lam[-1]) == approx(1.8655016981625037e-6, rel=1e-10)
    assert float(frontier_lam[-1]) == approx(-0.1818156514988331, abs=1e-12)
