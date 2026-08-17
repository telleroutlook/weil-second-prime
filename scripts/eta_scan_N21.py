"""
Eta scan for N=21, L=0.56: build M0/S0/M2/S2 row-by-row with checkpointing.
Extends eta_scan_N15.py to N=21 for convergence model discrimination.

N=19 diagnostic window: 3-pt model (b_L>0) predicts N=19 ~ -0.054,
6-pt model (all N) predicts ~ -0.043. N=21 float result further constrains B.

Usage:
  python3 scripts/eta_scan_N21.py             # fresh run
  python3 scripts/eta_scan_N21.py --resume    # resume from checkpoint
  python3 scripts/eta_scan_N21.py --scan-only # skip build, just do eta scan from saved matrices
"""
import argparse, json, math, time, pathlib
import numpy as np
from fractions import Fraction

from src.archimedean.integrator_a import integrate_M_K, integrate_S_KK, integrate_S_VK
from src.archimedean.log_moments import V_matrix_entry, V2_matrix_entry
from src.prime_layer.legendre_shift import compute_J, compute_E
from src.prime_layer.legendre_shift_2prime import compute_F
from src.archimedean.kernel import kappa as compute_kappa

CKPT = pathlib.Path("pilots/eta_scan_N21.ckpt.npz")
OUT  = pathlib.Path("pilots/eta_scan_N21.json")

L_NUM, L_DEN = 56, 100
L = 0.56; N = 21; PARITY = 1
INDICES = list(range(PARITY, PARITY + 2*N, 2))
n = len(INDICES)
Gd = np.array([2.0 / (2*ni + 1) for ni in INDICES])
c2 = math.log(2) / math.sqrt(2)
c3 = math.log(3) / math.sqrt(3)
TAU2 = Fraction(math.log(2) / L).limit_denominator(10000)
TAU3 = Fraction(math.log(3) / L).limit_denominator(10000)


def mid_iv(iv):
    return 0.5 * (float(iv[0]) + float(iv[1]))

def mid_r(r):
    return 0.5 * (float(r.enclosure_lower) + float(r.enclosure_upper))

def build_row(a: int, i: int) -> tuple:
    """Build row a of M0, S0, M2, S2. Returns (m0_row, s0_row, m2_row, s2_row)."""
    m0 = np.zeros(n); s0 = np.zeros(n)
    m2 = np.zeros(n); s2 = np.zeros(n)
    for b, j in enumerate(INDICES):
        V   = mid_iv(V_matrix_entry(i, j, 128))
        K   = mid_r(integrate_M_K(i, j, L_NUM, L_DEN, depth=2, use_bernstein=False, skip_remainder=True))
        m0[b] = V + K
        svv = mid_iv(V2_matrix_entry(i, j, 128))
        svk = mid_r(integrate_S_VK(i, j, L_NUM, L_DEN, depth=2))
        skv = mid_r(integrate_S_VK(j, i, L_NUM, L_DEN, depth=2))
        skk = mid_r(integrate_S_KK(i, j, L_NUM, L_DEN, depth=3))
        s0[b] = svv + svk + skv + skk
        J2   = float(compute_J(i, j, TAU2)); J3 = float(compute_J(i, j, TAU3))
        m2[b] = -(c2*J2 + c3*J3)
        E2   = float(compute_E(i, j, TAU2)); E3 = float(compute_E(i, j, TAU3))
        Fij  = float(compute_F(i, j, TAU2, TAU3)); Fji = float(compute_F(j, i, TAU2, TAU3))
        s2[b] = c2**2*E2 + c3**2*E3 + c2*c3*(Fij + Fji)
    return m0, s0, m2, s2


def do_eta_scan(M0, S0, M2, S2):
    kappa_val = float(compute_kappa(L_NUM, L_DEN, prec=128))
    c_L = math.log(2 * math.pi * L) + 0.5772156649015329
    d = 2*N + 1
    b_L = sum(1/k for k in range(1, d+1)) - c_L - kappa_val
    Ginv = np.diag(1.0 / Gd)
    R0 = S0 - M0.T @ Ginv @ M0
    R2 = S2 - M2.T @ Ginv @ M2
    T_mat = np.diag([sum(1.0/k for k in range(1, ni+1))*Gd[aa] for aa, ni in enumerate(INDICES)])
    F_mat = T_mat + M0 + M2 - c_L * np.diag(Gd)

    print(f"b_L={b_L:.6f}  kappa={kappa_val:.6f}", flush=True)

    norm_R0 = np.linalg.norm(R0, 'fro')
    norm_R2 = np.linalg.norm(R2, 'fro')
    k_full = norm_R2 / norm_R0
    eta_star = math.sqrt(k_full)
    print(f"||R0||_F={norm_R0:.6e}  ||R2||_F={norm_R2:.6e}  k={k_full:.4f}  eta*={eta_star:.4f}", flush=True)

    etas = [0.5, 0.75, 0.887, 1.0, 1.25, 1.5, 2.0, 2.49, eta_star, 3.0, 4.0, 5.0, 6.97]
    results = []
    best_eta = 0.5; best_eig = -1e9
    print("\nEta scan: lambda_min(C(eta))")
    for eta in etas:
        R_eta = (1 + eta)*R0 + (1 + 1.0/eta)*R2
        C = b_L*F_mat - R_eta; C = 0.5*(C + C.T)
        eig = float(np.linalg.eigvalsh(C)[0])
        mark = " *** POSITIVE ***" if eig > 0 else ""
        print(f"  eta={eta:.3f}  lambda_min={eig:+.8f}{mark}", flush=True)
        results.append({"eta": eta, "lambda_min": eig})
        if eig > best_eig: best_eig = eig; best_eta = eta

    print(f"\nBest: eta={best_eta:.3f}  lambda_min={best_eig:+.8f}")
    print("\nConvergence model check:")
    print(f"  3-pt model (b_L>0) N=21 prediction: -0.042732")
    print(f"  6-pt model (all N) N=21 prediction: -0.020706")
    if best_eig > -0.020706:
        print("  => ABOVE 6-pt: B > 0, form POSITIVE DEFINITE in limit")
    elif best_eig < -0.042732:
        print("  => BELOW 3-pt: B < 0, form INDEFINITE")
    else:
        print("  => BETWEEN models: inconclusive")

    return {
        "N": N, "b_L": b_L, "kappa": kappa_val,
        "norm_R0_F": norm_R0, "norm_R2_F": norm_R2,
        "k_frob": k_full, "eta_star_frob": eta_star,
        "eta_scan": results,
        "best_eta": best_eta, "best_lambda_min": best_eig,
        "model_3pt_N21": -0.042732,
        "model_6pt_N21": -0.020706,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--scan-only", action="store_true")
    args = parser.parse_args()

    M0 = np.zeros((n, n)); S0 = np.zeros((n, n))
    M2 = np.zeros((n, n)); S2 = np.zeros((n, n))
    start_row = 0

    if (args.resume or args.scan_only) and CKPT.exists():
        data = np.load(CKPT)
        M0[:] = data["M0"]; S0[:] = data["S0"]
        M2[:] = data["M2"]; S2[:] = data["S2"]
        start_row = int(data.get("done_rows", 0))
        print(f"Resumed from checkpoint: {start_row}/{n} rows done", flush=True)
        if args.scan_only:
            if start_row < n:
                print(f"WARNING: only {start_row}/{n} rows available, eta scan on partial data", flush=True)
            result = do_eta_scan(M0, S0, M2, S2)
            OUT.write_text(json.dumps(result, indent=2))
            print(f"Saved to {OUT}")
            return

    t_total = time.time()
    for a in range(start_row, n):
        i = INDICES[a]
        ta = time.time()
        m0, s0, m2, s2 = build_row(a, i)
        M0[a] = m0; S0[a] = s0; M2[a] = m2; S2[a] = s2
        print(f"row a={a:2d} P{i:2d} done ({time.time()-ta:.1f}s, total {time.time()-t_total:.1f}s)", flush=True)
        np.savez(CKPT, M0=M0, S0=S0, M2=M2, S2=S2, done_rows=a+1)

    print(f"\nAll {n} rows complete. Running eta scan...", flush=True)
    result = do_eta_scan(M0, S0, M2, S2)
    OUT.write_text(json.dumps(result, indent=2))
    CKPT.unlink(missing_ok=True)
    print(f"Done. Results saved to {OUT}")


if __name__ == "__main__":
    main()
