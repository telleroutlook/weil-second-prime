"""
Sequential sub-matrix chain: compute k=18..25 eigenvalue spectra.

Extends the N=17 Arb checkpoint with new rows (P35..P49) one at a time.
Each row is saved to pilots/submatrix_row{r:02d}.npz for resumability.
Each spectrum is saved to pilots/submatrix_k{k:02d}.json.

Row numbering: row r corresponds to basis function P_{2r+1}.
  r=17 → P35 (k=18), r=18 → P37 (k=19), ..., r=24 → P49 (k=25)

Usage:
  python3 scripts/submatrix_chain.py [--start ROW] [--end ROW]

  Default: --start 17 --end 24  (k=18..25)

If pilots/submatrix_k18_row17.npz exists (from compute_submatrix_k18.py),
row 17 is loaded from it automatically (no recompute).

Checkpointing: each completed row saved immediately. Restart resumes from
the last saved row.
"""
import argparse, json, math, time, pathlib
import numpy as np
from fractions import Fraction

from src.archimedean.integrator_a import integrate_M_K, integrate_S_KK, integrate_S_VK
from src.archimedean.log_moments import V_matrix_entry, V2_matrix_entry
from src.prime_layer.legendre_shift import compute_J, compute_E
from src.prime_layer.legendre_shift_2prime import compute_F
from src.archimedean.kernel import kappa as compute_kappa

CKPT17 = pathlib.Path("pilots/cert_fp_second_N17.ckpt.json")
N_OLD = 17

L = 0.56; L_NUM, L_DEN = 56, 100; PARITY = 1
c2 = math.log(2) / math.sqrt(2)
c3 = math.log(3) / math.sqrt(3)
TAU2 = Fraction(math.log(2) / L).limit_denominator(10000)
TAU3 = Fraction(math.log(3) / L).limit_denominator(10000)
kappa_val = float(compute_kappa(L_NUM, L_DEN, prec=128))
c_L = math.log(2 * math.pi * L) + 0.5772156649015329


def mid_iv(iv):
    return 0.5 * (float(Fraction(iv[0])) + float(Fraction(iv[1])))


def mid_r(r):
    return 0.5 * (float(r.enclosure_lower) + float(r.enclosure_upper))


def all_indices(n):
    return list(range(PARITY, PARITY + 2 * n, 2))


def load_n17_block(n_max):
    with open(CKPT17) as f:
        ckpt = json.load(f)
    M0 = np.zeros((n_max, n_max)); S0 = np.zeros((n_max, n_max))
    for k2, iv in ckpt['M0'].items():
        a, b = map(int, k2.split(','))
        if a < N_OLD and b < N_OLD and a < n_max and b < n_max:
            M0[a, b] = mid_iv(iv)
    for k2, iv in ckpt['S0'].items():
        a, b = map(int, k2.split(','))
        if a < N_OLD and b < N_OLD and a < n_max and b < n_max:
            S0[a, b] = mid_iv(iv)
    return M0, S0


def compute_row(row_idx, indices):
    i = indices[row_idx]
    n = len(indices)
    m0 = np.zeros(n); s0 = np.zeros(n)
    t0 = time.time()
    for b, j in enumerate(indices):
        tb = time.time()
        V = mid_iv(V_matrix_entry(i, j, 128))
        K = mid_r(integrate_M_K(i, j, L_NUM, L_DEN, depth=2,
                                  use_bernstein=False, skip_remainder=True))
        m0[b] = V + K
        svv = mid_iv(V2_matrix_entry(i, j, 128))
        svk = mid_r(integrate_S_VK(i, j, L_NUM, L_DEN, depth=2))
        skv = mid_r(integrate_S_VK(j, i, L_NUM, L_DEN, depth=2))
        skk = mid_r(integrate_S_KK(i, j, L_NUM, L_DEN, depth=3))
        s0[b] = svv + svk + skv + skk
        print(f"  [row {row_idx}] b={b:2d} P{j:2d}: {time.time()-tb:.1f}s  total={time.time()-t0:.1f}s",
              flush=True)
    return m0, s0


def compute_m2_s2(indices):
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


def do_spectrum(M0, S0, M2, S2, indices, k):
    n = len(indices)
    Gd = np.array([2.0 / (2 * ni + 1) for ni in indices])
    Ginv = np.diag(1.0 / Gd)
    R0 = S0 - M0.T @ Ginv @ M0
    R2 = S2 - M2.T @ Ginv @ M2
    d = 2 * k + 1
    b_L = sum(1 / i for i in range(1, d + 1)) - c_L - kappa_val
    T_mat = np.diag([sum(1.0 / i for i in range(1, ni + 1)) * Gd[aa]
                     for aa, ni in enumerate(indices)])
    F_mat = T_mat + M0 + M2 - c_L * np.diag(Gd)
    results = {}
    best_eta = None; best_lmin = -1e9
    for eta in [0.5, 0.65, 0.80, 0.90, 1.0, 1.10, 1.22, 2.0]:
        R_eta = (1 + eta) * R0 + (1 + 1.0 / eta) * R2
        C = b_L * F_mat - R_eta; C = 0.5 * (C + C.T)
        evals = np.linalg.eigvalsh(C)
        n_neg = int(np.sum(evals < -1e-9))
        lmin = float(evals[0])
        results[eta] = {"lambda_min": lmin, "n_neg": n_neg,
                        "top4": [float(e) for e in evals[:4]]}
        mark = " *** POSITIVE ***" if evals[0] > 0 else ""
        print(f"  k={k} eta={eta:.2f}  λ_min={evals[0]:+.8f}  n_neg={n_neg}{mark}", flush=True)
        if lmin > best_lmin:
            best_lmin = lmin; best_eta = eta
    print(f"  Best: eta={best_eta:.2f}  λ_min_opt={best_lmin:+.8f}", flush=True)
    results["_best_eta"] = best_eta; results["_best_lmin"] = best_lmin
    return b_L, results


def row_save_path(r):
    # Try k18-specific file first (from compute_submatrix_k18.py)
    if r == 17:
        alt = pathlib.Path("pilots/submatrix_k18_row17.npz")
        if alt.exists():
            return alt
    return pathlib.Path(f"pilots/submatrix_row{r:02d}.npz")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=17,
                        help="First row to compute (17=P35, k=18)")
    parser.add_argument("--end", type=int, default=24,
                        help="Last row to compute (24=P49, k=25)")
    args = parser.parse_args()

    row_start = args.start
    row_end = args.end
    n_max = row_end + 1

    print(f"Sub-matrix chain: rows {row_start}..{row_end}  (k={row_start+1}..{row_end+1})",
          flush=True)
    print(f"Loading N=17 checkpoint...", flush=True)
    M0 = np.zeros((n_max, n_max)); S0 = np.zeros((n_max, n_max))
    m0_base, s0_base = load_n17_block(n_max)
    M0[:N_OLD, :N_OLD] = m0_base[:N_OLD, :N_OLD]
    S0[:N_OLD, :N_OLD] = s0_base[:N_OLD, :N_OLD]

    # Load already-completed extension rows
    for r in range(N_OLD, row_start):
        rp = row_save_path(r)
        if rp.exists():
            data = np.load(str(rp))
            m0r = data['m0']; s0r = data['s0']
            M0[r, :len(m0r)] = m0r; M0[:len(m0r), r] = m0r
            S0[r, :len(s0r)] = s0r; S0[:len(s0r), r] = s0r
            print(f"  Loaded extension row {r} from {rp}", flush=True)
        else:
            print(f"WARNING: extension row {r} not found at {rp}", flush=True)

    t_chain = time.time()
    for row_idx in range(row_start, row_end + 1):
        n = row_idx + 1  # current matrix size
        indices = all_indices(n)
        new_degree = indices[row_idx]
        out_json = pathlib.Path(f"pilots/submatrix_k{n:02d}.json")

        print(f"\n{'='*60}", flush=True)
        print(f"Row {row_idx} (P{new_degree}), k={n}  [chain elapsed {time.time()-t_chain:.0f}s]",
              flush=True)

        # Load or compute row
        rp = row_save_path(row_idx)
        if rp.exists():
            print(f"Loading row {row_idx} from {rp}", flush=True)
            data = np.load(str(rp))
            m0_row = data['m0']; s0_row = data['s0']
        else:
            print(f"Computing row {row_idx} (P{new_degree} × P1..P{new_degree}, {n} entries)...",
                  flush=True)
            t_row = time.time()
            m0_row, s0_row = compute_row(row_idx, indices)
            save_path = pathlib.Path(f"pilots/submatrix_row{row_idx:02d}.npz")
            np.savez(str(save_path), m0=m0_row, s0=s0_row)
            print(f"Row {row_idx} done in {time.time()-t_row:.0f}s, saved to {save_path}",
                  flush=True)

        # Fill symmetric matrix
        M0[row_idx, :n] = m0_row[:n]; M0[:n, row_idx] = m0_row[:n]
        S0[row_idx, :n] = s0_row[:n]; S0[:n, row_idx] = s0_row[:n]

        # M2, S2 analytical
        M2, S2 = compute_m2_s2(indices)

        # Spectrum
        b_L, spectrum = do_spectrum(M0[:n, :n], S0[:n, :n], M2, S2, indices, n)
        print(f"b_L={b_L:.6f}", flush=True)

        # Save
        result = {
            "k": n, "b_L": b_L, "new_degree": new_degree,
            "indices": indices,
            "spectrum": {str(eta): v for eta, v in spectrum.items()},
        }
        out_json.write_text(json.dumps(result, indent=2))
        print(f"Saved to {out_json}", flush=True)

    print(f"\nChain complete. Total time: {time.time()-t_chain:.0f}s", flush=True)


if __name__ == "__main__":
    main()
