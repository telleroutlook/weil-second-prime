"""
Extend sub-matrix analysis by computing a single new row (P_{2*row_idx+1}).

Usage:
  python3 scripts/compute_submatrix_row.py --row ROW

  ROW is the row index (0-based). Row 17 = P35, Row 18 = P37, etc.
  Each row extends the sub-matrix by one dimension.

Loads N=17 Arb checkpoint (rows 0-16) plus any already-computed extension
rows from pilots/submatrix_row_*.npz, assembles the (ROW+1)×(ROW+1) Galerkin
matrix, computes eigenvalues, saves to pilots/submatrix_k{ROW+1}.json.
"""
import argparse, json, math, time, pathlib
import numpy as np
from fractions import Fraction
from fractions import Fraction

from src.archimedean.integrator_a import integrate_M_K, integrate_S_KK, integrate_S_VK
from src.archimedean.log_moments import V_matrix_entry, V2_matrix_entry
from src.prime_layer.legendre_shift import compute_J, compute_E
from src.prime_layer.legendre_shift_2prime import compute_F
from src.archimedean.kernel import kappa as compute_kappa

CKPT17 = pathlib.Path("pilots/cert_fp_second_N17.ckpt.json")
N_OLD = 17  # existing checkpoint size (P1..P33)

L = 0.56; L_NUM, L_DEN = 56, 100; PARITY = 1
c2 = math.log(2) / math.sqrt(2)
c3 = math.log(3) / math.sqrt(3)
TAU2 = Fraction(math.log(2) / L).limit_denominator(10000)
TAU3 = Fraction(math.log(3) / L).limit_denominator(10000)


def mid_iv(iv):
    return 0.5 * (float(Fraction(iv[0])) + float(Fraction(iv[1])))

def mid_r(r):
    return 0.5 * (float(r.enclosure_lower) + float(r.enclosure_upper))


def all_indices(n):
    return list(range(PARITY, PARITY + 2 * n, 2))


def load_n17(n):
    with open(CKPT17) as f:
        ckpt = json.load(f)
    M0 = np.zeros((n, n)); S0 = np.zeros((n, n))
    for k2, iv in ckpt['M0'].items():
        a, b = map(int, k2.split(','))
        if a < N_OLD and b < N_OLD and a < n and b < n:
            M0[a, b] = mid_iv(iv)
    for k2, iv in ckpt['S0'].items():
        a, b = map(int, k2.split(','))
        if a < N_OLD and b < N_OLD and a < n and b < n:
            S0[a, b] = mid_iv(iv)
    return M0, S0


def compute_row(row_idx, indices):
    """Compute M0[row_idx,:], S0[row_idx,:] for new row i = indices[row_idx]."""
    i = indices[row_idx]
    n = len(indices)
    m0 = np.zeros(n); s0 = np.zeros(n)
    t0 = time.time()
    for b, j in enumerate(indices):
        tb = time.time()
        V  = mid_iv(V_matrix_entry(i, j, 128))
        K  = mid_r(integrate_M_K(i, j, L_NUM, L_DEN, depth=2,
                                  use_bernstein=False, skip_remainder=True))
        m0[b] = V + K
        svv = mid_iv(V2_matrix_entry(i, j, 128))
        svk = mid_r(integrate_S_VK(i, j, L_NUM, L_DEN, depth=2))
        skv = mid_r(integrate_S_VK(j, i, L_NUM, L_DEN, depth=2))
        skk = mid_r(integrate_S_KK(i, j, L_NUM, L_DEN, depth=3))
        s0[b] = svv + svk + skv + skk
        print(f"  P{i}×P{j:2d} done ({time.time()-tb:.1f}s, total {time.time()-t0:.1f}s)",
              flush=True)
    return m0, s0


def compute_m2_s2_full(indices):
    n = len(indices)
    M2 = np.zeros((n, n)); S2 = np.zeros((n, n))
    for a, i in enumerate(indices):
        for b, j in enumerate(indices):
            J2 = float(compute_J(i, j, TAU2)); J3 = float(compute_J(i, j, TAU3))
            M2[a, b] = -(c2 * J2 + c3 * J3)
            E2 = float(compute_E(i, j, TAU2)); E3 = float(compute_E(i, j, TAU3))
            Fij = float(compute_F(i, j, TAU2, TAU3)); Fji = float(compute_F(j, i, TAU2, TAU3))
            S2[a, b] = c2**2 * E2 + c3**2 * E3 + c2 * c3 * (Fij + Fji)
    return M2, S2


def do_spectrum(M0, S0, M2, S2, indices):
    n = len(indices)
    kappa_val = float(compute_kappa(L_NUM, L_DEN, prec=128))
    c_L = math.log(2 * math.pi * L) + 0.5772156649015329
    Gd = np.array([2.0 / (2 * ni + 1) for ni in indices])
    Ginv = np.diag(1.0 / Gd)
    R0 = S0 - M0.T @ Ginv @ M0
    R2 = S2 - M2.T @ Ginv @ M2
    d = 2 * n + 1
    b_L = sum(1 / k for k in range(1, d + 1)) - c_L - kappa_val
    T_mat = np.diag([sum(1.0 / k for k in range(1, ni + 1)) * Gd[aa]
                     for aa, ni in enumerate(indices)])
    F_mat = T_mat + M0 + M2 - c_L * np.diag(Gd)
    results = {}
    for eta in [0.5, 1.0, 1.22, 2.0, 4.0]:
        R_eta = (1 + eta) * R0 + (1 + 1.0 / eta) * R2
        C = b_L * F_mat - R_eta; C = 0.5 * (C + C.T)
        evals = np.linalg.eigvalsh(C)
        n_neg = int(np.sum(evals < -1e-9))
        results[eta] = {"lambda_min": float(evals[0]), "n_neg": n_neg,
                        "top6": [float(e) for e in evals[:6]]}
        mark = " *** POSITIVE ***" if evals[0] > 0 else ""
        print(f"  eta={eta:.2f}  lambda_min={evals[0]:+.8f}  n_neg={n_neg}{mark}", flush=True)
    return b_L, results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--row", type=int, required=True,
                        help="Row index to compute (17=P35, 18=P37, ...)")
    args = parser.parse_args()
    row_idx = args.row
    n = row_idx + 1  # new matrix size
    indices = all_indices(n)
    new_degree = indices[row_idx]

    print(f"Computing k={n} sub-matrix (adding P{new_degree}, row {row_idx})", flush=True)

    row_file = pathlib.Path(f"pilots/submatrix_row{row_idx:02d}.npz")
    out_file = pathlib.Path(f"pilots/submatrix_k{n:02d}.json")

    # Load base N=17 checkpoint
    M0, S0 = load_n17(n)

    # Load previously-computed extension rows (rows N_OLD..row_idx-1)
    for r in range(N_OLD, row_idx):
        rf = pathlib.Path(f"pilots/submatrix_row{r:02d}.npz")
        if rf.exists():
            data = np.load(str(rf))
            m0_r = data['m0']; s0_r = data['s0']
            M0[r, :len(m0_r)] = m0_r; M0[:len(m0_r), r] = m0_r
            S0[r, :len(s0_r)] = s0_r; S0[:len(s0_r), r] = s0_r
            print(f"Loaded extension row {r} (P{all_indices(n)[r]})", flush=True)
        else:
            raise FileNotFoundError(
                f"Extension row {r} not found at {rf}. "
                f"Run --row {r} first."
            )

    # Compute or load the new row
    if row_file.exists():
        print(f"Loading saved row {row_idx} from {row_file}...", flush=True)
        data = np.load(str(row_file))
        m0_row = data['m0']; s0_row = data['s0']
    else:
        print(f"Computing row {row_idx} (P{new_degree} × P1..P{new_degree}, {n} entries)...",
              flush=True)
        m0_row, s0_row = compute_row(row_idx, indices)
        np.savez(str(row_file), m0=m0_row, s0=s0_row)
        print(f"Saved row {row_idx} to {row_file}", flush=True)

    # Fill row/col
    M0[row_idx, :n] = m0_row; M0[:n, row_idx] = m0_row
    S0[row_idx, :n] = s0_row; S0[:n, row_idx] = s0_row

    # M2, S2 analytical
    print("Computing M2, S2 analytically...", flush=True)
    M2, S2 = compute_m2_s2_full(indices)

    print(f"\nEigenvalue spectrum for k={n} (P1..P{new_degree}), η scan:", flush=True)
    b_L, spectrum = do_spectrum(M0, S0, M2, S2, indices)
    print(f"\nb_L = {b_L:.6f}  (k={n})", flush=True)

    result = {
        "k": n, "b_L": b_L,
        "new_degree": new_degree,
        "indices": indices,
        "spectrum": {str(eta): v for eta, v in spectrum.items()},
    }
    out_file.write_text(json.dumps(result, indent=2))
    print(f"Saved to {out_file}")


if __name__ == "__main__":
    main()
