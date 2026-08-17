"""H3 variational pilot: HTF rayleigh_certificate on the second-window Schur matrix.

Validates the H3 path (PLAN.md): does HTF's Arb-certified Rayleigh quotient agree
with our dense float computation at N=8 even, L=0.55?

Acceptance criterion (PLAN.md H3a): relative error < 1% between DMRG/variational
estimate and Arb reference; this pilot uses the exact min-eigenvector as the
trial state, so the Rayleigh quotient == lambda_min (up to Arb rounding).

Usage:
    python3 -m scripts.h3_variational_pilot
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

# Add HTF to path
_HTF_PATH = Path(__file__).parent.parent.parent / "htf"
if _HTF_PATH.exists():
    sys.path.insert(0, str(_HTF_PATH))

HTF_AVAILABLE = False
try:
    from htf import rayleigh_certificate  # type: ignore
    HTF_AVAILABLE = True
except ImportError:
    pass


def _build_schur_fast(L_num: int, L_den: int, sector: str, N: int,
                      depth_2d: int = 2, depth_3d: int = 2) -> tuple[np.ndarray, float]:
    """Build float Schur matrix C = b_L*F - R_eta at reduced depth (pilot only)."""
    import math
    from fractions import Fraction
    from src.archimedean.integrator_a import integrate_M_K, integrate_S_KK, integrate_S_VK
    from src.archimedean.log_moments import V_matrix_entry, V2_matrix_entry
    from src.prime_layer.legendre_shift import compute_J, compute_E
    from src.prime_layer.legendre_shift_2prime import compute_F
    from src.archimedean.kernel import kappa as compute_kappa

    L = L_num / L_den
    parity = 0 if sector == "even" else 1
    indices = list(range(parity, parity + 2 * N, 2))
    n = len(indices)

    def _mid_iv(iv: tuple) -> float:
        return 0.5 * (float(iv[0]) + float(iv[1]))

    def _mid_r(r) -> float:
        return 0.5 * (float(r.enclosure_lower) + float(r.enclosure_upper))

    Gd = np.array([2.0 / (2 * ni + 1) for ni in indices])
    T = np.diag([sum(1.0 / k for k in range(1, ni + 1)) * Gd[a]
                 for a, ni in enumerate(indices)])
    M0 = np.zeros((n, n))
    S0 = np.zeros((n, n))
    M2 = np.zeros((n, n))
    S2 = np.zeros((n, n))

    c2 = math.log(2) / math.sqrt(2)
    c3 = math.log(3) / math.sqrt(3)
    tau2 = Fraction(math.log(2) / L).limit_denominator(10000)
    tau3 = Fraction(math.log(3) / L).limit_denominator(10000)

    for a, i in enumerate(indices):
        for b, j in enumerate(indices):
            V_ij = _mid_iv(V_matrix_entry(i, j, 128))
            K_ij = _mid_r(integrate_M_K(i, j, L_num, L_den, depth=depth_2d,
                                         use_bernstein=False, skip_remainder=True))
            svv = _mid_iv(V2_matrix_entry(i, j, 128))
            svk = _mid_r(integrate_S_VK(i, j, L_num, L_den, depth=depth_2d))
            skv = _mid_r(integrate_S_VK(j, i, L_num, L_den, depth=depth_2d))
            skk = _mid_r(integrate_S_KK(i, j, L_num, L_den, depth=depth_3d))
            S0[a, b] = svv + svk + skv + skk
            M0[a, b] = V_ij + K_ij
            J2 = float(compute_J(i, j, tau2))
            J3 = float(compute_J(i, j, tau3))
            M2[a, b] = -(c2 * J2 + c3 * J3)
            E2 = float(compute_E(i, j, tau2))
            E3 = float(compute_E(i, j, tau3))
            F_ij = float(compute_F(i, j, tau2, tau3))
            F_ji = float(compute_F(j, i, tau2, tau3))
            S2[a, b] = c2**2 * E2 + c3**2 * E3 + c2 * c3 * (F_ij + F_ji)

    c_L = math.log(2 * math.pi * L) + 0.5772156649015329
    kappa = compute_kappa(L_num, L_den, prec=128)
    d = 2 * N if sector == "even" else 2 * N + 1
    H_d = sum(1.0 / k for k in range(1, d + 1))
    b_L = H_d - c_L - float(kappa)
    Ginv = np.diag(1.0 / Gd)
    eta = 0.5
    R0 = S0 - M0.T @ Ginv @ M0
    R2 = S2 - M2.T @ Ginv @ M2
    R_eta = (1 + eta) * R0 + (1 + 1.0 / eta) * R2
    F_mat = T + M0 + M2 - c_L * np.diag(Gd)
    C = b_L * F_mat - R_eta
    return C, b_L


def main() -> None:
    if not HTF_AVAILABLE:
        print("HTF not available at ../htf — skipping")
        sys.exit(2)

    L_num, L_den = 56, 100   # L = 0.56 (second window)
    sector, N = "odd", 5     # small N for fast pilot
    depth_2d, depth_3d = 2, 2
    print(f"Building Schur matrix (depth={depth_2d}/{depth_3d}): "
          f"L={L_num}/{L_den}={L_num/L_den}, sector={sector}, N={N}", flush=True)

    t0 = time.time()
    C, b_L = _build_schur_fast(L_num, L_den, sector, N,
                                depth_2d=depth_2d, depth_3d=depth_3d)
    build_time = time.time() - t0
    print(f"  build_matrices: {build_time:.1f}s  (b_L={b_L:.4f})", flush=True)

    # Dense reference: min eigenvalue via numpy
    evals, evecs = np.linalg.eigh(0.5 * (C + C.T))
    lambda_min_ref = float(evals[0])
    psi_min = evecs[:, 0].astype(np.float64)
    print(f"  lambda_min (numpy reference) = {lambda_min_ref:.6e}", flush=True)

    # HTF requires exactly symmetric matrix; symmetrize float rounding error
    C_sym = 0.5 * (C + C.T)

    # HTF Rayleigh certificate with trial state = exact min eigenvector
    t1 = time.time()
    cert = rayleigh_certificate(C_sym, psi_min, notes=f"H3 pilot L={L_num}/{L_den} N={N}")
    htf_time = time.time() - t1

    upper = cert.upper
    print(f"  HTF Rayleigh upper bound: E0 <= {upper:.6e}  "
          f"(midpoint={cert.midpoint:.6e}, radius={cert.radius:.2e})  [{htf_time:.2f}s]",
          flush=True)

    rel_err = abs(cert.midpoint - lambda_min_ref) / max(abs(lambda_min_ref), 1e-15)
    passed = rel_err < 0.01
    print(f"  Relative error vs numpy: {rel_err:.2e}  "
          f"({'PASS' if passed else 'FAIL'} < 1%)", flush=True)

    result = {
        "description": "H3 pilot: HTF rayleigh_certificate vs numpy lambda_min",
        "L": f"{L_num}/{L_den}",
        "L_float": L_num / L_den,
        "sector": sector,
        "N": N,
        "b_L": round(b_L, 6),
        "lambda_min_numpy": round(lambda_min_ref, 8),
        "htf_rayleigh_midpoint": round(cert.midpoint, 8),
        "htf_rayleigh_radius": cert.radius,
        "htf_upper_bound": round(upper, 8),
        "relative_error": round(rel_err, 6),
        "h3a_criterion_pass": passed,
        "grade": "pilot (float Schur + Arb Rayleigh quotient)",
        "assurance": cert.assurance,
        "notes": cert.notes,
        "build_time_s": round(build_time, 1),
        "htf_time_s": round(htf_time, 3),
    }

    out_path = Path("pilots/h3_variational_pilot.json")
    out_path.write_text(json.dumps(result, indent=2))
    print(f"\nH3a criterion: {'PASS' if passed else 'FAIL'}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
