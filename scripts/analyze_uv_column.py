"""
Analyze the UV column of C[k=18] to understand the constant UV-cross term.

The UV-cross term: 2 × v_UV × Σ_{j<UV} C[j,UV] × v_j

If UV-cross → constant as k→∞, there must be an asymptotic balance between:
  - C[j,UV] decay as UV mode degree increases
  - v_j / v_UV ratio (eigenvector structure)

This script:
1. Reconstructs full C matrix for k=18
2. Analyzes C[:,UV] column by layer (d = UV - j)
3. Checks whether C[j,UV] ~ f(d) × g(j) (separable structure)
4. Estimates analytical UV-cross bound from layer structure
5. Probes whether UV-cross ≈ -integral of a kernel
"""
import json, math, pathlib, numpy as np
from fractions import Fraction
from src.prime_layer.legendre_shift import compute_J, compute_E
from src.prime_layer.legendre_shift_2prime import compute_F
from src.archimedean.kernel import kappa as compute_kappa

K18_ROW17 = pathlib.Path("pilots/submatrix_k18_row17.npz")
CKPT17 = pathlib.Path("pilots/cert_fp_second_N17.ckpt.json")

L=0.56; L_NUM,L_DEN=56,100; PARITY=1
c2=math.log(2)/math.sqrt(2); c3=math.log(3)/math.sqrt(3)
TAU2=Fraction(math.log(2)/L).limit_denominator(10000)
TAU3=Fraction(math.log(3)/L).limit_denominator(10000)
kappa_val=float(compute_kappa(L_NUM,L_DEN,prec=128))
c_L=math.log(2*math.pi*L)+0.5772156649015329

def mid_iv(iv): return 0.5*(float(Fraction(iv[0]))+float(Fraction(iv[1])))

k=18; n=18
indices=list(range(PARITY, PARITY+2*n, 2))  # P1,P3,...,P35

# Load matrices
with open(CKPT17) as f: ckpt=json.load(f)
M0=np.zeros((n,n)); S0=np.zeros((n,n))
for k2,iv in ckpt['M0'].items():
    a,b=map(int,k2.split(','))
    if a<17 and b<17: M0[a,b]=mid_iv(iv)
for k2,iv in ckpt['S0'].items():
    a,b=map(int,k2.split(','))
    if a<17 and b<17: S0[a,b]=mid_iv(iv)
r17=np.load(str(K18_ROW17))
M0[17,:n]=r17['m0'][:n]; M0[:n,17]=r17['m0'][:n]
S0[17,:n]=r17['s0'][:n]; S0[:n,17]=r17['s0'][:n]

M2=np.zeros((n,n)); S2=np.zeros((n,n))
for a,i in enumerate(indices):
    for b,j in enumerate(indices):
        J2=float(compute_J(i,j,TAU2)); J3=float(compute_J(i,j,TAU3))
        M2[a,b]=-(c2*J2+c3*J3)
        E2=float(compute_E(i,j,TAU2)); E3=float(compute_E(i,j,TAU3))
        Fij=float(compute_F(i,j,TAU2,TAU3)); Fji=float(compute_F(j,i,TAU2,TAU3))
        S2[a,b]=c2**2*E2+c3**2*E3+c2*c3*(Fij+Fji)

d_val=2*n+1
b_L=sum(1/i for i in range(1,d_val+1))-c_L-kappa_val
Gd=np.array([2.0/(2*ni+1) for ni in indices]); Ginv=np.diag(1.0/Gd)
R0=S0-M0.T@Ginv@M0; R2=S2-M2.T@Ginv@M2
T_mat=np.diag([sum(1.0/j for j in range(1,ni+1))*Gd[aa] for aa,ni in enumerate(indices)])
F_mat=T_mat+M0+M2-c_L*np.diag(Gd)

eta=1.0
R_eta=(1+eta)*R0+(1+1.0/eta)*R2
C=b_L*F_mat-R_eta; C=0.5*(C+C.T)
evals,evecs=np.linalg.eigh(C)
v0=evecs[:,0]; uv=n-1  # UV index = 17

print(f"=== UV column analysis k={n} (UV=P{indices[uv]}) ===")
print(f"b_L={b_L:.6f}, λ₀={evals[0]:.8f}")
print(f"|v_UV|={abs(v0[uv]):.6f},  |v_UV|²={v0[uv]**2:.6f}")
print()

# C[:,UV] column
col_UV = C[:,uv]
print("C[j,UV] column (d = UV-j, layer depth):")
print(f"{'j':>3} {'P_j':>5} {'d=UV-j':>7} {'C[j,UV]':>12} {'v_j':>10} {'contrib':>12}")
total_cross = 0.0
layer_sums = {}
for j in range(n-1):
    d_layer = (n-1) - j
    contrib = 2*col_UV[j]*v0[j]*v0[uv]
    total_cross += contrib
    if d_layer not in layer_sums:
        layer_sums[d_layer] = 0.0
    layer_sums[d_layer] += contrib
    print(f"{j:3d} {indices[j]:5d} {d_layer:7d} {col_UV[j]:12.6f} {v0[j]:10.6f} {contrib:12.8f}")

print(f"\nUV-diag: {C[uv,uv]*v0[uv]**2:+.8f}")
print(f"UV-cross total: {total_cross:+.8f}")
print(f"IR-block: {np.dot(v0[:uv], C[:uv,:uv]@v0[:uv]):+.8f}")
print(f"Sum (λ₀ check): {C[uv,uv]*v0[uv]**2 + total_cross + np.dot(v0[:uv], C[:uv,:uv]@v0[:uv]):+.8f}")

print("\nLayer sums (near-UV=d=1..5, far=d>=6):")
near_uv = 0.0; far = 0.0
for d_layer in sorted(layer_sums.keys()):
    s = layer_sums[d_layer]
    tag = "near-UV" if d_layer<=5 else "far"
    if d_layer<=5: near_uv += s
    else: far += s
    print(f"  d={d_layer:2d}: {s:+.8f}  ({tag})")
print(f"  near-UV total (d=1..5): {near_uv:+.8f}")
print(f"  far total (d>=6):       {far:+.8f}")

print("\n=== Separability check: C[j,UV] vs layer d ===")
print("Does C[j,UV] ≈ A × decay(d)^j for some A, decay?")
col_vals = col_UV[:n-1]
col_abs = np.abs(col_vals)
# Fit log-linear vs j
j_arr = np.arange(n-1)
valid = col_abs > 1e-10
if valid.sum() > 5:
    coeffs = np.polyfit(j_arr[valid], np.log(col_abs[valid]), 1)
    print(f"  log|C[j,UV]| ~ {coeffs[0]:.4f}*j + {coeffs[1]:.4f}")
    print(f"  → decay rate per j: {math.exp(coeffs[0]):.4f}")
    print(f"  → implied 'radius': {math.exp(-coeffs[0]):.4f}")

print("\n=== F vs R decomposition of UV column ===")
F_col = F_mat[:,uv]
R_col = R_eta[:,uv]
print(f"{'j':>3} {'b_L*F[j,UV]':>14} {'R_η[j,UV]':>12} {'net C[j,UV]':>12}")
for j in range(min(n-1, 18)):
    print(f"{j:3d} {b_L*F_col[j]:14.8f} {R_col[j]:12.8f} {col_UV[j]:12.8f}")

print("\n=== M0 vs M2 decomposition of F column ===")
print(f"{'j':>3} {'T[j,UV]':>10} {'M0[j,UV]':>10} {'M2[j,UV]':>10} {'F[j,UV]':>10}")
for j in range(min(n,18)):
    t_val = T_mat[j,uv]
    m0_val = M0[j,uv]
    m2_val = M2[j,uv]
    f_val = F_mat[j,uv]
    print(f"{j:3d} {t_val:10.6f} {m0_val:10.6f} {m2_val:10.6f} {f_val:10.6f}")

print("\n=== R0 vs R2 decomposition of R_η column ===")
print(f"{'j':>3} {'(1+η)R0[j,UV]':>15} {'(1+1/η)R2[j,UV]':>17} {'R_η[j,UV]':>12}")
for j in range(min(n,18)):
    r0_part = (1+eta)*R0[j,uv]
    r2_part = (1+1/eta)*R2[j,uv]
    print(f"{j:3d} {r0_part:15.8f} {r2_part:17.8f} {r0_part+r2_part:12.8f}")

print("\n=== UV-cross asymptotic estimate ===")
# Near-UV block (d=1..5) contributes ≈ -0.020 (stable)
# Far block (d>=6) contributes ≈ -0.009 at k=18, rate ≈ 0.92
# As k→∞: far→0, near-UV≈-0.020 → B₀ ≥ -0.020 (near-UV lower bound)
# But near-UV also evolves... check d=1 term specifically
d1_j = n-2  # j=n-2, d=1
print(f"d=1 (j={d1_j}=P{indices[d1_j]}): C[j,UV]={col_UV[d1_j]:.8f}, v_j={v0[d1_j]:.8f}")
print(f"  contrib = {2*col_UV[d1_j]*v0[d1_j]*v0[uv]:.8f}")
print(f"  This is the P_{{2k-3}} × P_{{2k-1}} entry — nearest-UV coupling")
print()
print(f"Fraction of UV-cross from d=1 alone: {layer_sums.get(1,0)/total_cross:.3f}")
print(f"Fraction from d=1..3:                {sum(layer_sums.get(d,0) for d in [1,2,3])/total_cross:.3f}")
