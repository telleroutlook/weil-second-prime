"""
Process k=19 sub-matrix result when pilots/submatrix_k19.json appears.

Computes:
1. λ₀ at all η values
2. UV-cross decomposition at η=1.0 and η_opt
3. Convergence rate r(18→19) — key model discriminator
4. Updated B₀ estimate and predictions for k=20..25

Run: python3 scripts/process_k19_result.py
"""
import json, math, pathlib, numpy as np, sys
from fractions import Fraction
from src.prime_layer.legendre_shift import compute_J, compute_E
from src.prime_layer.legendre_shift_2prime import compute_F
from src.archimedean.kernel import kappa as compute_kappa

K19_JSON = pathlib.Path("pilots/submatrix_k19.json")
K18_ROW17 = pathlib.Path("pilots/submatrix_k18_row17.npz")
K19_ROW18 = pathlib.Path("pilots/submatrix_row18.npz")
CKPT17 = pathlib.Path("pilots/cert_fp_second_N17.ckpt.json")

if not K19_JSON.exists():
    print("pilots/submatrix_k19.json not found. Chain still running.")
    sys.exit(1)

L=0.56; L_NUM,L_DEN=56,100; PARITY=1
c2=math.log(2)/math.sqrt(2); c3=math.log(3)/math.sqrt(3)
TAU2=Fraction(math.log(2)/L).limit_denominator(10000)
TAU3=Fraction(math.log(3)/L).limit_denominator(10000)
kappa_val=float(compute_kappa(L_NUM,L_DEN,prec=128))
c_L=math.log(2*math.pi*L)+0.5772156649015329

def mid_iv(iv): return 0.5*(float(Fraction(iv[0]))+float(Fraction(iv[1])))

with open(CKPT17) as f: ckpt=json.load(f)
N17=17
M0f=np.zeros((N17,N17)); S0f=np.zeros((N17,N17))
for k2,iv in ckpt['M0'].items():
    a,b=map(int,k2.split(','))
    if a<N17 and b<N17: M0f[a,b]=mid_iv(iv)
for k2,iv in ckpt['S0'].items():
    a,b=map(int,k2.split(','))
    if a<N17 and b<N17: S0f[a,b]=mid_iv(iv)

# Load all rows
k = 19; n = 19
M0 = np.zeros((n,n)); S0 = np.zeros((n,n))
M0[:N17,:N17] = M0f; S0[:N17,:N17] = S0f

# row 17 (k18) — only has 18 entries; entry [17,18] filled by row18 symmetry
r17 = np.load(str(K18_ROW17))
nr17 = len(r17['m0'])
M0[17,:nr17]=r17['m0'][:nr17]; M0[:nr17,17]=r17['m0'][:nr17]
S0[17,:nr17]=r17['s0'][:nr17]; S0[:nr17,17]=r17['s0'][:nr17]

# row 18 (k19)
r18 = np.load(str(K19_ROW18))
M0[18,:n]=r18['m0'][:n]; M0[:n,18]=r18['m0'][:n]
S0[18,:n]=r18['s0'][:n]; S0[:n,18]=r18['s0'][:n]

indices=list(range(PARITY,PARITY+2*k,2))
M2=np.zeros((n,n)); S2=np.zeros((n,n))
for a,i in enumerate(indices):
    for b,j in enumerate(indices):
        J2=float(compute_J(i,j,TAU2)); J3=float(compute_J(i,j,TAU3))
        M2[a,b]=-(c2*J2+c3*J3)
        E2=float(compute_E(i,j,TAU2)); E3=float(compute_E(i,j,TAU3))
        Fij=float(compute_F(i,j,TAU2,TAU3)); Fji=float(compute_F(j,i,TAU2,TAU3))
        S2[a,b]=c2**2*E2+c3**2*E3+c2*c3*(Fij+Fji)

d=2*k+1; b_L=sum(1/i for i in range(1,d+1))-c_L-kappa_val
Gd=np.array([2.0/(2*ni+1) for ni in indices]); Ginv=np.diag(1.0/Gd)
R0=S0-M0.T@Ginv@M0; R2=S2-M2.T@Ginv@M2
T_mat=np.diag([sum(1.0/j for j in range(1,ni+1))*Gd[aa] for aa,ni in enumerate(indices)])
F_mat=T_mat+M0+M2-c_L*np.diag(Gd)

# Total positivity checks (verified k=13..18: always True)
print("=== Total positivity checks ===")
R_eta1=2*R0+2*R2
print(f"S0 ALL pos: {np.all(S0>=-1e-12)}  (neg={np.sum(S0<-1e-12)}/{n*n})")
print(f"M0 ALL pos: {np.all(M0>=-1e-12)}  (neg={np.sum(M0<-1e-12)}/{n*n})")
print(f"R0 ALL nonneg: {np.all(R0>=-1e-12)}  (neg={np.sum(R0<-1e-12)}/{n*n})")
print(f"R_eta(eta=1) ALL nonneg: {np.all(R_eta1>=-1e-12)}  (neg={np.sum(R_eta1<-1e-12)}/{n*n})")


print(f"=== k=19 analysis (b_L={b_L:.6f}) ===\n")

print("η scan:")
lam18_eta1 = -0.048918
best_lam = 1e9; best_eta = None
for eta in [0.5, 0.65, 0.80, 0.90, 1.0, 1.10, 1.22, 2.0]:
    R_eta=(1+eta)*R0+(1+1.0/eta)*R2
    C=b_L*F_mat-R_eta; C=0.5*(C+C.T)
    evals=np.linalg.eigvalsh(C)
    n_neg=int(np.sum(evals<-1e-9))
    print(f"  η={eta:.2f}  λ_min={evals[0]:+.8f}  n_neg={n_neg}")
    if evals[0] > best_lam or best_eta is None:
        best_lam = evals[0]; best_eta = eta

print(f"\n  Best: η={best_eta:.2f}  λ_min_opt={best_lam:+.8f}")

# UV decomposition at η=1.0
print(f"\n=== UV decomposition at η=1.0 ===")
R_eta=(1+1.0)*R0+(1+1.0)*R2
C=b_L*F_mat-R_eta; C=0.5*(C+C.T)
evals,evecs=np.linalg.eigh(C)
lam19_eta1=evals[0]; v0=evecs[:,0]; uv=k-1
v_UV=v0[uv]; v_IR=v0[:uv]
C_UV_diag=C[uv,uv]
C_diag=C_UV_diag*v_UV**2
C_cross=2.0*np.dot(C[:uv,uv],v_IR)*v_UV
C_block=np.dot(v_IR,C[:uv,:uv]@v_IR)
loc3=sum(v0[np.argsort(np.abs(v0))[::-1][:3]]**2)
print(f"λ₀(k=19, η=1.0)  = {lam19_eta1:+.8f}")
print(f"UV-diag           = {C_diag:+.8f}  (C_UV={C_UV_diag:+.6f})")
print(f"UV-cross          = {C_cross:+.8f}")
print(f"IR-block          = {C_block:+.8f}")
print(f"|v_UV|²           = {v_UV**2:.4f}   top-3 loc = {loc3:.3f}")

# Sign-definiteness check (empirical lemma verification)
col_UV = C[:uv, uv]
if v0[uv] < 0: v0_check = -v0
else: v0_check = v0
all_col_neg = all(x < 0 for x in col_UV)
all_v_pos = all(x > 0 for x in v0_check)
print(f"\n=== Sign-definiteness check ===")
print(f"C[j,UV] ALL<0 (k=13..18: always True): {all_col_neg}")
print(f"v₀ ALL same sign (k=13..18: always True): {all_v_pos}")
print(f"UV-cross sign-definite negative: {all_col_neg and all_v_pos}")
if not all_col_neg:
    pos_entries = [(j, col_UV[j]) for j in range(len(col_UV)) if col_UV[j] >= 0]
    print(f"  VIOLATION: positive C[j,UV] entries: {pos_entries[:5]}")
if not all_v_pos:
    neg_entries = [(j, v0_check[j]) for j in range(len(v0_check)) if v0_check[j] <= 0]
    print(f"  VIOLATION: non-positive v₀ entries: {neg_entries[:5]}")

# Model comparison
print(f"\n=== Model comparison ===")
rate1819 = lam19_eta1 / lam18_eta1
print(f"Rate r(18→19)         = {rate1819:.4f}")
print(f"Expected (geometric)  = 0.7897")
print(f"Expected (component)  = 0.8885")
print(f"")
print(f"UV-cross at k=19      = {C_cross:.5f}")
print(f"UV-cross at k=18      = -0.02876")
print(f"UV-cross ratio 18→19  = {C_cross/-0.02876:.4f}")
print(f"")
print(f"IR-block at k=19      = {C_block:.5f}")
print(f"IR-block at k=18      = -0.018018")
print(f"IR-block ratio 18→19  = {C_block/-0.018018:.4f}")
print(f"")

# Updated B0 estimate
B0_new = np.mean([C_cross, -0.028927, -0.028759])  # k=17,18,19
print(f"Updated B₀ estimate (UV-cross mean k=17..19): {B0_new:.5f}")

# k=20..25 predictions
print(f"\n=== Updated k=20..25 predictions ===")
ir_rate = C_block / (-0.018018)  # post-k19 IR rate
B0 = B0_new
ud_rate = 0.879  # stable
lam_prev = lam19_eta1
print(f"Using B₀={B0:.5f}, IR_rate={ir_rate:.4f} (per step)")
print(f"{'k':>3}  {'pred_lam0':>10}  {'rate':>7}")
ir_val = C_block
ud_val = C_diag
for kk in range(20, 26):
    ir_val *= ir_rate
    ud_val *= ud_rate
    lam_pred = B0 + ir_val + ud_val
    rate = lam_pred / lam_prev
    print(f"{kk:3d}  {lam_pred:+10.5f}  {rate:.4f}")
    lam_prev = lam_pred
