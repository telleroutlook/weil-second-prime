"""
Process N=19 certify result when cert_fp_second_N19.ckpt.json is complete.

Loads the N=19 Arb checkpoint (float centers), assembles the Schur matrix
at optimal eta, computes eigenvalues, and reports the convergence model fit.

Usage:
  python3 scripts/process_N19_result.py
  python3 scripts/process_N19_result.py --ckpt pilots/cert_fp_second_N19.ckpt.json
"""
import argparse, json, math, pathlib
import numpy as np
from fractions import Fraction
from src.prime_layer.legendre_shift import compute_J, compute_E
from src.prime_layer.legendre_shift_2prime import compute_F
from src.archimedean.kernel import kappa as compute_kappa

L = 0.56; L_NUM, L_DEN = 56, 100; PARITY = 1; N = 19
INDICES = list(range(PARITY, PARITY + 2*N, 2))
n = len(INDICES)
c2 = math.log(2) / math.sqrt(2)
c3 = math.log(3) / math.sqrt(3)
TAU2 = Fraction(math.log(2) / L).limit_denominator(10000)
TAU3 = Fraction(math.log(3) / L).limit_denominator(10000)


def mid_iv(iv):
    return 0.5 * (float(Fraction(iv[0])) + float(Fraction(iv[1])))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ckpt", default="pilots/cert_fp_second_N19.ckpt.json")
    args = parser.parse_args()

    ckpt_path = pathlib.Path(args.ckpt)
    if not ckpt_path.exists():
        print(f"Checkpoint not found: {ckpt_path}")
        return

    print(f"Loading N=19 checkpoint from {ckpt_path}...")
    with open(ckpt_path) as f:
        ckpt = json.load(f)

    M0 = np.zeros((n, n)); S0 = np.zeros((n, n))
    n_loaded_M0 = 0; n_loaded_S0 = 0
    for k2, iv in ckpt.get('M0', {}).items():
        a, b = map(int, k2.split(','))
        if a < n and b < n:
            M0[a, b] = mid_iv(iv); n_loaded_M0 += 1
    for k2, iv in ckpt.get('S0', {}).items():
        a, b = map(int, k2.split(','))
        if a < n and b < n:
            S0[a, b] = mid_iv(iv); n_loaded_S0 += 1

    expected = n * n
    print(f"Loaded M0: {n_loaded_M0}/{expected}  S0: {n_loaded_S0}/{expected}")
    if n_loaded_M0 < expected or n_loaded_S0 < expected:
        print("Symmetrizing to fill missing upper-triangle entries...")
        M0 = np.maximum(M0, M0.T); S0 = np.maximum(S0, S0.T)
        M0 = 0.5 * (M0 + M0.T); S0 = 0.5 * (S0 + S0.T)
        # Fill any remaining zeros on diagonal from submatrix npz if available
        npz_row = pathlib.Path("pilots/submatrix_row18.npz")
        if npz_row.exists():
            r18 = np.load(str(npz_row))
            for col in range(n):
                if M0[n-1, col] == 0.0 and col < len(r18['m0']):
                    M0[n-1, col] = r18['m0'][col]; M0[col, n-1] = r18['m0'][col]
                if S0[n-1, col] == 0.0 and col < len(r18['s0']):
                    S0[n-1, col] = r18['s0'][col]; S0[col, n-1] = r18['s0'][col]
            print(f"Patched last row from submatrix_row18.npz")

    # Compute M2, S2 analytically
    M2 = np.zeros((n, n)); S2 = np.zeros((n, n))
    for a, i in enumerate(INDICES):
        for b, j in enumerate(INDICES):
            J2 = float(compute_J(i, j, TAU2)); J3 = float(compute_J(i, j, TAU3))
            M2[a, b] = -(c2 * J2 + c3 * J3)
            E2 = float(compute_E(i, j, TAU2)); E3 = float(compute_E(i, j, TAU3))
            Fij = float(compute_F(i, j, TAU2, TAU3)); Fji = float(compute_F(j, i, TAU2, TAU3))
            S2[a, b] = c2**2 * E2 + c3**2 * E3 + c2 * c3 * (Fij + Fji)

    kappa_val = float(compute_kappa(L_NUM, L_DEN, prec=128))
    c_L = math.log(2 * math.pi * L) + 0.5772156649015329
    d = 2 * N + 1
    b_L = sum(1/k for k in range(1, d+1)) - c_L - kappa_val
    Gd = np.array([2.0 / (2*ni+1) for ni in INDICES])
    Ginv = np.diag(1.0 / Gd)
    R0 = S0 - M0.T @ Ginv @ M0
    R2 = S2 - M2.T @ Ginv @ M2
    T_mat = np.diag([sum(1.0/k for k in range(1, ni+1)) * Gd[aa] for aa, ni in enumerate(INDICES)])
    F_mat = T_mat + M0 + M2 - c_L * np.diag(Gd)

    print(f"\nN=19 float analysis: b_L = {b_L:.6f}")
    print(f"||R0||_F = {np.linalg.norm(R0,'fro'):.4e}  ||R2||_F = {np.linalg.norm(R2,'fro'):.4e}")
    frob_eta = math.sqrt(np.linalg.norm(R2,'fro') / np.linalg.norm(R0,'fro'))
    print(f"Frobenius eta* = {frob_eta:.4f}")

    print("\nEta scan:")
    print(f"{'eta':>6} | {'lambda_min':>12} | {'n_neg':>5} | {'lambda_1':>12}")
    best_eta = 0.5; best_lmin = -1e9
    results = {}
    for eta in [0.5, 1.0, frob_eta, 1.22, 2.0, 2.49, 4.0]:
        R_eta = (1+eta)*R0 + (1+1.0/eta)*R2
        C = b_L*F_mat - R_eta; C = 0.5*(C+C.T)
        evals = np.linalg.eigvalsh(C)
        n_neg = int(np.sum(evals < -1e-9))
        lmin = float(evals[0]); l1 = float(evals[1])
        mark = " ***POS***" if lmin > 0 else ""
        print(f"{eta:6.3f} | {lmin:+12.8f} | {n_neg:5d} | {l1:+12.8f}{mark}")
        results[eta] = {"lambda_min": lmin, "lambda_1": l1, "n_neg": n_neg}
        if lmin > best_lmin:
            best_lmin = lmin; best_eta = eta

    print(f"\nBest: eta={best_eta:.4f}  lambda_min={best_lmin:+.8f}")

    # Convergence model comparison
    print("\n--- Convergence model check ---")
    # Sub-matrix data k=13..17
    k_prev = [13, 14, 15, 16, 17]
    lam_prev = [-0.16040557, -0.12815207, -0.09898596, -0.07925052, -0.06277780]
    lam19 = best_lmin
    print(f"lambda_0(19) actual = {lam19:+.6f}")
    print("Model predictions at N=19:")
    print(f"  exp+B (5-pt fit): -0.039253")
    print(f"  pow+B (5-pt fit): -0.037556")
    print(f"  exp(B=0):         -0.038968")
    print(f"  1/k+1/k^2 (B=0):  -0.037223")

    # 6-pt fit with N=19 included
    from scipy.optimize import curve_fit
    k_all = np.array([13, 14, 15, 16, 17, 19], dtype=float)
    lam_all = np.array([*lam_prev, lam19])
    def exp_model(k, A, r, B): return A * r**k + B
    def exp0_model(k, A, r): return A * r**k
    try:
        popt_exp, _ = curve_fit(exp_model, k_all, lam_all, p0=[-3.5, 0.79, -0.001], maxfev=10000)
        A, r, B = popt_exp
        print(f"\n6-pt exp+B fit (k=13..17,19): A={A:.4f} r={r:.5f} B={B:+.6f}")
        for Np in [21, 23, 25, 30]:
            print(f"  N={Np}: {exp_model(Np, *popt_exp):+.6f}")
    except Exception as e:
        print(f"6-pt exp+B fit failed: {e}")

    try:
        popt_exp0, _ = curve_fit(exp0_model, k_all, lam_all, p0=[-3.5, 0.79], maxfev=10000)
        A0, r0 = popt_exp0
        print(f"6-pt exp(B=0) fit: A={A0:.4f} r={r0:.5f}")
        for Np in [21, 23, 25, 30]:
            print(f"  N={Np}: {exp0_model(Np, *popt_exp0):+.6f}")
    except Exception as e:
        print(f"6-pt exp(B=0) fit failed: {e}")

    # Min eigenvector at N=19
    eta = best_eta
    R_eta = (1+eta)*R0 + (1+1.0/eta)*R2
    C = b_L*F_mat - R_eta; C = 0.5*(C+C.T)
    evals, evecs = np.linalg.eigh(C)
    v0 = evecs[:, 0]
    idx_sorted = np.argsort(np.abs(v0))[::-1]
    print(f"\nMinimum eigenvector (lambda_0={evals[0]:+.6f}) top-5 components:")
    for k in range(5):
        a = idx_sorted[k]
        print(f"  P{INDICES[a]:2d}: {v0[a]:+.4f}")
    print(f"  (UV mode: expect P{INDICES[-1]} = P{2*N-1} to dominate)")


if __name__ == "__main__":
    main()
