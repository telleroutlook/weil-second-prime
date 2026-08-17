"""
Extend sub-matrix analysis from k=17 to k=18 by computing row 17 (P35).

Loads the N=17 Arb checkpoint (rows 0-16, P1-P33), computes the new row
corresponding to P35 (index 17), assembles the 18×18 Galerkin matrix,
and saves the eigenvalue spectrum to pilots/submatrix_k18.json.

This gives the k=18 data point for convergence model discrimination
in ~70 minutes, much faster than waiting for the full N=21 build.
"""
import json, math, time, pathlib
import numpy as np
from fractions import Fraction

from src.archimedean.integrator_a import integrate_M_K, integrate_S_KK, integrate_S_VK
from src.archimedean.log_moments import V_matrix_entry, V2_matrix_entry
from src.prime_layer.legendre_shift import compute_J, compute_E
from src.prime_layer.legendre_shift_2prime import compute_F
from src.archimedean.kernel import kappa as compute_kappa

CKPT17 = pathlib.Path("pilots/cert_fp_second_N17.ckpt.json")
OUT    = pathlib.Path("pilots/submatrix_k18.json")
ROW_SAVE = pathlib.Path("pilots/submatrix_k18_row17.npz")

L = 0.56; L_NUM, L_DEN = 56, 100; PARITY = 1
N_OLD = 17  # existing checkpoint size
N_NEW = 18  # new size: add P35
INDICES = list(range(PARITY, PARITY + 2*N_NEW, 2))  # [1,3,...,35]
n = len(INDICES)
assert INDICES[-1] == 35, f"Expected P35 as last, got P{INDICES[-1]}"

c2 = math.log(2) / math.sqrt(2)
c3 = math.log(3) / math.sqrt(3)
TAU2 = Fraction(math.log(2) / L).limit_denominator(10000)
TAU3 = Fraction(math.log(3) / L).limit_denominator(10000)


def mid_iv(iv):
    return 0.5 * (float(Fraction(iv[0])) + float(Fraction(iv[1])))

def mid_r(r):
    return 0.5 * (float(r.enclosure_lower) + float(r.enclosure_upper))


def load_n17_checkpoint():
    with open(CKPT17) as f:
        ckpt = json.load(f)
    M0 = np.zeros((N_NEW, N_NEW))
    S0 = np.zeros((N_NEW, N_NEW))
    for k2, iv in ckpt['M0'].items():
        a, b = map(int, k2.split(','))
        if a < N_OLD and b < N_OLD:
            M0[a, b] = mid_iv(iv)
    for k2, iv in ckpt['S0'].items():
        a, b = map(int, k2.split(','))
        if a < N_OLD and b < N_OLD:
            S0[a, b] = mid_iv(iv)
    return M0, S0


def compute_row17():
    """Compute M0[17,:], S0[17,:], M2[17,:], S2[17,:] for i=P35."""
    i = INDICES[17]  # P35
    assert i == 35
    m0 = np.zeros(n); s0 = np.zeros(n)
    m2 = np.zeros(n); s2 = np.zeros(n)
    t0 = time.time()
    for b, j in enumerate(INDICES):
        tb = time.time()
        V   = mid_iv(V_matrix_entry(i, j, 128))
        K   = mid_r(integrate_M_K(i, j, L_NUM, L_DEN, depth=2,
                                   use_bernstein=False, skip_remainder=True))
        m0[b] = V + K
        svv = mid_iv(V2_matrix_entry(i, j, 128))
        svk = mid_r(integrate_S_VK(i, j, L_NUM, L_DEN, depth=2))
        skv = mid_r(integrate_S_VK(j, i, L_NUM, L_DEN, depth=2))
        skk = mid_r(integrate_S_KK(i, j, L_NUM, L_DEN, depth=3))
        s0[b] = svv + svk + skv + skk
        J2  = float(compute_J(i, j, TAU2)); J3 = float(compute_J(i, j, TAU3))
        m2[b] = -(c2 * J2 + c3 * J3)
        E2  = float(compute_E(i, j, TAU2)); E3 = float(compute_E(i, j, TAU3))
        Fij = float(compute_F(i, j, TAU2, TAU3)); Fji = float(compute_F(j, i, TAU2, TAU3))
        s2[b] = c2**2 * E2 + c3**2 * E3 + c2 * c3 * (Fij + Fji)
        print(f"  b={b:2d} P{j:2d} done ({time.time()-tb:.1f}s, total {time.time()-t0:.1f}s)",
              flush=True)
    return m0, s0, m2, s2


def compute_m2_s2_full(M0, S0):
    """Compute the full M2 and S2 matrices analytically (no slow integrals)."""
    M2 = np.zeros((n, n)); S2 = np.zeros((n, n))
    for a, i in enumerate(INDICES):
        for b, j in enumerate(INDICES):
            J2 = float(compute_J(i, j, TAU2)); J3 = float(compute_J(i, j, TAU3))
            M2[a, b] = -(c2 * J2 + c3 * J3)
            E2 = float(compute_E(i, j, TAU2)); E3 = float(compute_E(i, j, TAU3))
            Fij = float(compute_F(i, j, TAU2, TAU3)); Fji = float(compute_F(j, i, TAU2, TAU3))
            S2[a, b] = c2**2 * E2 + c3**2 * E3 + c2 * c3 * (Fij + Fji)
    return M2, S2


def do_spectrum(M0, S0, M2, S2):
    kappa_val = float(compute_kappa(L_NUM, L_DEN, prec=128))
    c_L = math.log(2 * math.pi * L) + 0.5772156649015329
    Gd = np.array([2.0 / (2 * ni + 1) for ni in INDICES])
    Ginv = np.diag(1.0 / Gd)
    R0 = S0 - M0.T @ Ginv @ M0
    R2 = S2 - M2.T @ Ginv @ M2
    d = 2 * n + 1
    b_L = sum(1 / k for k in range(1, d + 1)) - c_L - kappa_val
    T_mat = np.diag([sum(1.0 / k for k in range(1, ni + 1)) * Gd[aa]
                     for aa, ni in enumerate(INDICES)])
    F_mat = T_mat + M0 + M2 - c_L * np.diag(Gd)

    results = {}
    for eta in [0.5, 1.0, 1.22, 2.0, 2.49, 4.0]:
        R_eta = (1 + eta) * R0 + (1 + 1.0 / eta) * R2
        C = b_L * F_mat - R_eta; C = 0.5 * (C + C.T)
        evals = np.linalg.eigvalsh(C)
        n_neg = int(np.sum(evals < -1e-9))
        results[eta] = {"lambda_min": float(evals[0]), "n_neg": n_neg,
                        "top4": [float(e) for e in evals[:4]]}
        mark = " *** POSITIVE ***" if evals[0] > 0 else ""
        print(f"  eta={eta:.3f}  lambda_min={evals[0]:+.8f}  n_neg={n_neg}{mark}", flush=True)
    return b_L, results


def main():
    print(f"Computing k=18 sub-matrix (adding P35 to N=17 checkpoint)", flush=True)

    # Load N=17 checkpoint
    M0, S0 = load_n17_checkpoint()
    print(f"Loaded N=17 checkpoint ({N_OLD}×{N_OLD} blocks)", flush=True)

    # Load or compute row 17
    if ROW_SAVE.exists():
        print("Loading saved row 17 from checkpoint...", flush=True)
        data = np.load(str(ROW_SAVE))
        m0_row = data['m0']; s0_row = data['s0']
        m2_row = data['m2']; s2_row = data['s2']
    else:
        print(f"Computing row 17 (P35 × P1..P35, {n} entries)...", flush=True)
        m0_row, s0_row, m2_row, s2_row = compute_row17()
        np.savez(str(ROW_SAVE), m0=m0_row, s0=s0_row, m2=m2_row, s2=s2_row)
        print(f"Row 17 saved to {ROW_SAVE}", flush=True)

    # Fill in row/col 17 (symmetric for M0 and S0)
    M0[17, :] = m0_row; M0[:, 17] = m0_row  # M0 symmetric
    S0[17, :] = s0_row; S0[:, 17] = s0_row  # S0 symmetric

    # Compute full M2, S2 analytically (fast)
    print("Computing M2, S2 analytically...", flush=True)
    M2, S2 = compute_m2_s2_full(M0, S0)
    M2[17, :] = m2_row; M2[:, 17] = m2_row  # analytic, symmetric
    S2[17, :] = s2_row; S2[:, 17] = s2_row

    # Eigenvalue spectrum
    print(f"\nEigenvalue spectrum for k=18 (P1..P35), η scan:", flush=True)
    b_L, spectrum = do_spectrum(M0, S0, M2, S2)
    print(f"\nb_L = {b_L:.6f}", flush=True)

    result = {"k": N_NEW, "b_L": b_L, "indices": INDICES, "spectrum": {
        str(eta): v for eta, v in spectrum.items()
    }}
    OUT.write_text(json.dumps(result, indent=2))
    print(f"\nSaved to {OUT}")


if __name__ == "__main__":
    main()
