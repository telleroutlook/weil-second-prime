"""
General sub-matrix analysis script for any k in the chain k=19..25.

Usage: PYTHONPATH=. python3 scripts/process_submatrix_k.py --k 20

Computes:
1. Total positivity checks (S0, M0, R0, R_eta)
2. η scan and η_opt
3. UV-cross decomposition at η=1.0
4. Rate r(k-1 → k) and model comparison
5. Updated B₀ estimate and predictions
"""
import argparse, json, math, pathlib, numpy as np, sys
from fractions import Fraction
from src.prime_layer.legendre_shift import compute_J, compute_E
from src.prime_layer.legendre_shift_2prime import compute_F
from src.archimedean.kernel import kappa as compute_kappa

parser = argparse.ArgumentParser()
parser.add_argument("--k", type=int, required=True, help="submatrix size (k=19..25)")
args_cli = parser.parse_args()

k = args_cli.k
n = k
CKPT17 = pathlib.Path("pilots/cert_fp_second_N17.ckpt.json")
K18_ROW17 = pathlib.Path("pilots/submatrix_k18_row17.npz")
K_JSON = pathlib.Path(f"pilots/submatrix_k{k:02d}.json")
PREV_JSON = pathlib.Path(f"pilots/submatrix_k{k-1:02d}_analysis.json")

if not K_JSON.exists():
    print(f"pilots/submatrix_k{k:02d}.json not found. Chain still running.")
    sys.exit(1)

L = 0.56; L_NUM, L_DEN = 56, 100; PARITY = 1
c2 = math.log(2) / math.sqrt(2); c3 = math.log(3) / math.sqrt(3)
TAU2 = Fraction(math.log(2) / L).limit_denominator(10000)
TAU3 = Fraction(math.log(3) / L).limit_denominator(10000)
kappa_val = float(compute_kappa(L_NUM, L_DEN, prec=128))
c_L = math.log(2 * math.pi * L) + 0.5772156649015329

def mid_iv(iv):
    return 0.5 * (float(Fraction(iv[0])) + float(Fraction(iv[1])))

# Load N17 checkpoint
with open(CKPT17) as f:
    ckpt = json.load(f)
N17 = 17
M0f = np.zeros((N17, N17)); S0f = np.zeros((N17, N17))
for k2, iv in ckpt['M0'].items():
    a, b = map(int, k2.split(','))
    if a < N17 and b < N17: M0f[a, b] = mid_iv(iv)
for k2, iv in ckpt['S0'].items():
    a, b = map(int, k2.split(','))
    if a < N17 and b < N17: S0f[a, b] = mid_iv(iv)

M0 = np.zeros((n, n)); S0 = np.zeros((n, n))
M0[:N17, :N17] = M0f; S0[:N17, :N17] = S0f

# Row 17 (computed in k=18 context, may have fewer than n entries)
r17 = np.load(str(K18_ROW17))
nr17 = len(r17['m0'])
M0[17, :nr17] = r17['m0'][:nr17]; M0[:nr17, 17] = r17['m0'][:nr17]
S0[17, :nr17] = r17['s0'][:nr17]; S0[:nr17, 17] = r17['s0'][:nr17]

# Rows 18..n-1: each row r was computed in k=r+1 context (r+1 entries available)
for r in range(18, n):
    rfile = pathlib.Path(f"pilots/submatrix_row{r:02d}.npz")
    if not rfile.exists():
        print(f"Missing {rfile}"); sys.exit(1)
    rd = np.load(str(rfile))
    nr = len(rd['m0'])
    M0[r, :nr] = rd['m0'][:nr]; M0[:nr, r] = rd['m0'][:nr]
    S0[r, :nr] = rd['s0'][:nr]; S0[:nr, r] = rd['s0'][:nr]

# Build prime-layer matrices
indices = list(range(PARITY, PARITY + 2 * k, 2))
M2 = np.zeros((n, n)); S2 = np.zeros((n, n))
for a, i in enumerate(indices):
    for b, j in enumerate(indices):
        J2 = float(compute_J(i, j, TAU2)); J3 = float(compute_J(i, j, TAU3))
        M2[a, b] = -(c2 * J2 + c3 * J3)
        E2 = float(compute_E(i, j, TAU2)); E3 = float(compute_E(i, j, TAU3))
        Fij = float(compute_F(i, j, TAU2, TAU3)); Fji = float(compute_F(j, i, TAU2, TAU3))
        S2[a, b] = c2**2 * E2 + c3**2 * E3 + c2 * c3 * (Fij + Fji)

d = 2 * k + 1
b_L = sum(1 / i for i in range(1, d + 1)) - c_L - kappa_val
Gd = np.array([2.0 / (2 * ni + 1) for ni in indices])
Ginv = np.diag(1.0 / Gd)
R0 = S0 - M0.T @ Ginv @ M0
R2 = S2 - M2.T @ Ginv @ M2
T_mat = np.diag([sum(1.0 / j for j in range(1, ni + 1)) * Gd[aa] for aa, ni in enumerate(indices)])
F_mat = T_mat + M0 + M2 - c_L * np.diag(Gd)

print(f"=== Total positivity checks (k={k}) ===")
R_eta1 = 2 * R0 + 2 * R2
print(f"S0 ALL pos:          {np.all(S0 >= -1e-12)}  (neg={np.sum(S0 < -1e-12)}/{n*n})")
print(f"M0 ALL pos:          {np.all(M0 >= -1e-12)}  (neg={np.sum(M0 < -1e-12)}/{n*n})")
print(f"R0 ALL nonneg:       {np.all(R0 >= -1e-12)}  (neg={np.sum(R0 < -1e-12)}/{n*n})")
print(f"R_eta(eta=1) nonneg: {np.all(R_eta1 >= -1e-12)}  (neg={np.sum(R_eta1 < -1e-12)}/{n*n})")

print(f"\n=== k={k} analysis (b_L={b_L:.6f}) ===\n")

print("η scan:")
best_lam = None; best_eta = None
for eta in [0.50, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 1.00, 1.10, 1.22, 2.00]:
    R_eta = (1 + eta) * R0 + (1 + 1.0 / eta) * R2
    C = b_L * F_mat - R_eta; C = 0.5 * (C + C.T)
    evals = np.linalg.eigvalsh(C)
    n_neg = int(np.sum(evals < -1e-9))
    print(f"  η={eta:.2f}  λ_min={evals[0]:+.8f}  n_neg={n_neg}")
    if best_lam is None or evals[0] > best_lam:
        best_lam = evals[0]; best_eta = eta

print(f"\n  Best: η={best_eta:.2f}  λ_min_opt={best_lam:+.8f}")

print(f"\n=== UV decomposition at η=1.0 ===")
R_eta = 2.0 * R0 + 2.0 * R2
C = b_L * F_mat - R_eta; C = 0.5 * (C + C.T)
evals, evecs = np.linalg.eigh(C)
lam_k = evals[0]; v0 = evecs[:, 0]; uv = k - 1
v_UV = v0[uv]; v_IR = v0[:uv]
C_UV_diag = C[uv, uv]
C_diag = C_UV_diag * v_UV**2
C_cross = 2.0 * np.dot(C[:uv, uv], v_IR) * v_UV
C_block = np.dot(v_IR, C[:uv, :uv] @ v_IR)
loc3 = sum(v0[np.argsort(np.abs(v0))[::-1][:3]]**2)
print(f"λ₀(k={k}, η=1.0)     = {lam_k:+.8f}")
print(f"UV-diag              = {C_diag:+.8f}  (C_UV={C_UV_diag:+.6f})")
print(f"UV-cross             = {C_cross:+.8f}")
print(f"IR-block             = {C_block:+.8f}")
print(f"|v_UV|²              = {v_UV**2:.4f}   top-3 loc = {loc3:.3f}")

col_UV = C[:uv, uv]
v0_check = -v0 if v0[uv] < 0 else v0
all_col_neg = all(x < 0 for x in col_UV)
all_v_pos = all(x > 0 for x in v0_check)
print(f"\n=== Sign-definiteness check ===")
print(f"C[j,UV] ALL<0: {all_col_neg}")
print(f"v₀ ALL same sign: {all_v_pos}")
print(f"UV-cross sign-definite negative: {all_col_neg and all_v_pos}")
if not all_col_neg:
    pos_entries = [(j, col_UV[j]) for j in range(len(col_UV)) if col_UV[j] >= 0]
    print(f"  VIOLATION: positive C[j,UV] entries: {pos_entries[:5]}")

print(f"\n=== Model comparison ===")
if PREV_JSON.exists():
    with open(PREV_JSON) as f:
        prev = json.load(f)
    lam_prev = prev.get('lambda0_eta1', None)
    uvcross_prev = prev.get('uvcross_eta1', None)
    irblock_prev = prev.get('irblock_eta1', None)
    if lam_prev:
        print(f"Rate r(k-1→k)         = {lam_k / lam_prev:.4f}")
        print(f"Expected (geometric)  = 0.7897")
    if uvcross_prev:
        print(f"UV-cross ratio        = {C_cross / uvcross_prev:.4f}")
    if irblock_prev:
        print(f"IR-block ratio        = {C_block / irblock_prev:.4f}")
else:
    print(f"No previous JSON found at {PREV_JSON}")

# Save results to JSON for next step
out = {
    'k': k, 'n': n, 'b_L': b_L,
    'lambda0_eta1': float(lam_k),
    'lambda0_opt': float(best_lam),
    'eta_opt': float(best_eta),
    'uvcross_eta1': float(C_cross),
    'irblock_eta1': float(C_block),
    'uvdiag_eta1': float(C_diag),
    'vUV_sq': float(v_UV**2),
    'loc3': float(loc3),
}
K_JSON_OUT = pathlib.Path(f"pilots/submatrix_k{k:02d}_analysis.json")
with open(K_JSON_OUT, 'w') as f:
    json.dump(out, f, indent=2)
print(f"\nSaved analysis to {K_JSON_OUT}")
