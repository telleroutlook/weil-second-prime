"""
Compare UV column structure across k=13,15,17,18 to study convergence of UV-cross.

Key question: Do C[UV-d, UV] and v_{UV-d}/v_UV converge as k→∞ for fixed d?
If yes → UV-cross → constant → B₀ < 0 is analytically provable.

Uses only the N=17 checkpoint (submatrix approach) + k=18 row17.
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
eta=1.0

def mid_iv(iv): return 0.5*(float(Fraction(iv[0]))+float(Fraction(iv[1])))

def build_matrix(n, M0_full, S0_full):
    """Build C matrix for given n using top-left n×n of M0_full, S0_full."""
    indices=list(range(PARITY,PARITY+2*n,2))
    M0=M0_full[:n,:n]; S0=S0_full[:n,:n]
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
    T_mat=np.diag([sum(1.0/i for i in range(1,ni+1))*Gd[aa] for aa,ni in enumerate(indices)])
    F_mat=T_mat+M0+M2-c_L*np.diag(Gd)
    R_eta=(1+eta)*R0+(1+1.0/eta)*R2
    C=b_L*F_mat-R_eta; C=0.5*(C+C.T)
    return C, b_L, indices, R0, R2, F_mat

# Load full N=18 matrices
N17=17; N18=18
with open(CKPT17) as f: ckpt=json.load(f)
M0_full=np.zeros((N18,N18)); S0_full=np.zeros((N18,N18))
for k2,iv in ckpt['M0'].items():
    a,b=map(int,k2.split(','))
    if a<N17 and b<N17: M0_full[a,b]=mid_iv(iv)
for k2,iv in ckpt['S0'].items():
    a,b=map(int,k2.split(','))
    if a<N17 and b<N17: S0_full[a,b]=mid_iv(iv)
r17=np.load(str(K18_ROW17))
M0_full[17,:N18]=r17['m0'][:N18]; M0_full[:N18,17]=r17['m0'][:N18]
S0_full[17,:N18]=r17['s0'][:N18]; S0_full[:N18,17]=r17['s0'][:N18]

print("=== UV column convergence across k=13,15,17,18 ===\n")

# Analyze each k
results = {}
for n in [13, 15, 17, 18]:
    C, b_L, indices, R0, R2, F_mat = build_matrix(n, M0_full, S0_full)
    evals,evecs=np.linalg.eigh(C)
    v0=evecs[:,0]; uv=n-1
    lam0=evals[0]
    # Normalize so v_UV > 0 for comparison
    if v0[uv] < 0: v0 = -v0
    col_UV = C[:,uv]
    # UV-cross contributions
    cross = 2.0*np.dot(col_UV[:uv], v0[:uv])*v0[uv]
    results[n] = {'C': C, 'v': v0, 'uv': uv, 'lam0': lam0, 'cross': cross,
                  'col_UV': col_UV, 'b_L': b_L, 'indices': indices,
                  'R0': R0, 'R2': R2, 'F_mat': F_mat}

print(f"{'k':>4} {'b_L':>8} {'λ₀':>10} {'UV-cross':>10} {'|v_UV|²':>9}")
for n in [13, 15, 17, 18]:
    r = results[n]
    print(f"{n:4d} {r['b_L']:8.6f} {r['lam0']:10.6f} {r['cross']:10.6f} {r['v'][r['uv']]**2:9.4f}")

print("\n=== C[UV-d, UV] convergence for d=1..5 ===")
print(f"{'d':>3}", end="")
for n in [13, 15, 17, 18]:
    print(f"  k={n:2d}:C[UV-d,UV]", end="")
print()
for d in range(1, min(6, 13)):  # d=1..5
    print(f"{d:3d}", end="")
    for n in [13, 15, 17, 18]:
        r = results[n]
        j = r['uv'] - d
        if j >= 0:
            print(f"      {r['col_UV'][j]:+.8f}", end="")
        else:
            print(f"             N/A", end="")
    print()

print("\n=== v[UV-d]/v[UV] ratio convergence for d=1..5 ===")
print("(This is the relative amplitude of near-UV eigenvector components)")
print(f"{'d':>3}", end="")
for n in [13, 15, 17, 18]:
    print(f"  k={n:2d}: ratio  ", end="")
print()
for d in range(1, min(6, 13)):
    print(f"{d:3d}", end="")
    for n in [13, 15, 17, 18]:
        r = results[n]
        j = r['uv'] - d
        if j >= 0:
            ratio = r['v'][j] / r['v'][r['uv']]
            print(f"      {ratio:+.6f}    ", end="")
        else:
            print(f"             N/A", end="")
    print()

print("\n=== Product C[UV-d,UV] × v[UV-d]/v[UV] × |v_UV|² convergence ===")
print("(This is the per-d contribution to UV-cross, normalized by v_UV² = |v_UV|²)")
print("If this converges for each d, UV-cross → B₀×something < 0")
print(f"{'d':>3}", end="")
for n in [13, 15, 17, 18]:
    print(f"  k={n:2d}: contrib ", end="")
print()
for d in range(1, min(14, 13)):
    print(f"{d:3d}", end="")
    for n in [13, 15, 17, 18]:
        r = results[n]
        j = r['uv'] - d
        if j >= 0:
            contrib = 2*r['col_UV'][j]*r['v'][j]*r['v'][r['uv']]
            print(f"      {contrib:+.8f}", end="")
        else:
            print(f"             N/A", end="")
    print()

print("\n=== R0[UV-d,UV] convergence (main driver of C[j,UV] < 0) ===")
print(f"{'d':>3}", end="")
for n in [13, 15, 17, 18]:
    print(f"  k={n:2d}: R0[j,UV]  ", end="")
print()
for d in range(1, min(6, 13)):
    print(f"{d:3d}", end="")
    for n in [13, 15, 17, 18]:
        r = results[n]
        j = r['uv'] - d
        if j >= 0:
            print(f"      {r['R0'][j,r['uv']]:+.8f}  ", end="")
        else:
            print(f"             N/A", end="")
    print()

print("\n=== Toy model: does v[UV-d]/v[UV] follow a geometric sequence? ===")
print("If v[UV-d] ≈ r^d × v[UV], then sum_d contrib ≈ 2|v_UV|² × (sum_d C[UV-d,UV]×r^d)")
print()
for n in [17, 18]:
    r = results[n]
    ratios = []
    for d in range(1, min(8, n)):
        j = r['uv'] - d
        if j >= 0:
            ratios.append(r['v'][j] / r['v'][r['uv']])
    # Check if ratios form geometric sequence
    if len(ratios) >= 3:
        geom_rates = [ratios[i+1]/ratios[i] for i in range(len(ratios)-1)]
        print(f"k={n}: v[UV-d]/v[UV] ratios: {[f'{x:.4f}' for x in ratios[:6]]}")
        print(f"  consecutive ratios (r_{d+1}/r_d): {[f'{x:.4f}' for x in geom_rates[:5]]}")
        print()

print("\n=== b_L × F[UV-d,UV] vs R_η[UV-d,UV] for d=1..5 across k ===")
print("(Check whether R_η dominates over b_L×F for near-UV entries)")
for n in [13, 15, 17, 18]:
    r = results[n]
    print(f"\nk={n} (b_L={r['b_L']:.4f}):")
    print(f"  {'d':>2}  {'b_L×F':>12}  {'R_η':>12}  {'ratio R_η/(b_L×F)':>20}")
    for d in range(1, min(6, n)):
        j = r['uv'] - d
        if j >= 0:
            f_val = r['b_L'] * r['F_mat'][j, r['uv']]
            R_eta_val = (1+eta)*r['R0'][j,r['uv']] + (1+1.0/eta)*r['R2'][j,r['uv']]
            if abs(f_val) > 1e-12:
                rat = R_eta_val / f_val
            else:
                rat = float('inf')
            print(f"  {d:2d}  {f_val:12.8f}  {R_eta_val:12.8f}  {rat:20.4f}")
