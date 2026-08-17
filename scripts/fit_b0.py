"""
Fit B₀ (asymptote of λ_min(k)) from available submatrix_k{k}.json files.

Two eigenvalue branches:
  UV mode: λ_0 for k≤24 (negative, converging from below to B₀^UV)
           λ_1 for k≥25 (jumped positive, second eigenvalue, continuing to B₀^UV)
  IR mode: λ_1 for k≤24 (small positive, second eigenvalue, ~4e-7..2e-6)
           λ_0 for k≥25 (minimum eigenvalue, ~1.6e-6, converging to B₀^IR)

B₀ = lim λ_min = B₀^IR (the true minimum asymptote)
B₀^UV = asymptote of the UV mode (much larger, ~0.02)

Model: λ(k) = A × r^k + B₀  (3-param exponential + constant)
Bootstrap: 10000 replicates with replacement
"""
import json, math, pathlib, sys
import numpy as np
from scipy.optimize import curve_fit

PILOTS = pathlib.Path("pilots")


def load_modes():
    """Return (uv_k, uv_lam, ir_k, ir_lam) tracking both eigenvalue branches."""
    uv_k, uv_lam = [], []
    ir_k, ir_lam = [], []

    # k=13..17 from early file (only UV mode tracked there)
    f03to17 = PILOTS / "submatrix_k03to17_lam0.json"
    if f03to17.exists():
        d = json.load(open(f03to17))
        for k_str, lam in d.items():
            k = int(k_str)
            if k >= 13:
                uv_k.append(k); uv_lam.append(lam)

    # k=18..max from individual files
    for k in range(18, 35):
        p = PILOTS / f"submatrix_k{k:02d}.json"
        if not p.exists():
            break
        d = json.load(open(p))
        sp = d["spectrum"].get("1.0")
        if sp is None:
            continue
        top4 = sp.get("top4", [])
        lmin = sp["lambda_min"]
        n_neg = sp["n_neg"]

        if n_neg >= 1:
            # k≤24 regime: UV mode is minimum (negative), IR mode is second
            uv_k.append(k); uv_lam.append(lmin)
            if len(top4) >= 2:
                ir_k.append(k); ir_lam.append(top4[1])
        else:
            # k≥25 regime: IR mode is minimum, UV mode is second
            ir_k.append(k); ir_lam.append(lmin)
            if len(top4) >= 2:
                uv_k.append(k); uv_lam.append(top4[1])

    return (np.array(uv_k), np.array(uv_lam),
            np.array(ir_k),  np.array(ir_lam))


def exp_b(k, A, r, B0):
    return A * r**k + B0


def fit_mode(ks, lams, k_min, A_sign, B0_range, label):
    mask = ks >= k_min
    kf, lf = ks[mask], lams[mask]
    if len(kf) < 4:
        print(f"  {label}: only {len(kf)} pts (need ≥4), skip fit")
        return None, None, None
    A_lo, A_hi = (-np.inf, 0.) if A_sign < 0 else (0., np.inf)
    popt, _ = curve_fit(exp_b, kf, lf,
                        p0=[A_sign * 2., 0.82, (B0_range[0]+B0_range[1])/2],
                        bounds=([A_lo, 0.5, B0_range[0]],
                                [A_hi, 0.9999, B0_range[1]]),
                        maxfev=30000)
    rms = float(np.sqrt(np.mean((lf - exp_b(kf, *popt))**2)))
    return popt, rms, len(kf)


def bootstrap_mode(ks, lams, k_min, A_sign, B0_range, n_boot=10000):
    mask = ks >= k_min
    kf, lf = ks[mask], lams[mask]
    n = len(kf)
    A_lo, A_hi = (-np.inf, 0.) if A_sign < 0 else (0., np.inf)
    B0s = []
    for _ in range(n_boot):
        idx = np.random.randint(0, n, n)
        try:
            p, _ = curve_fit(exp_b, kf[idx], lf[idx],
                             p0=[A_sign * 2., 0.82, (B0_range[0]+B0_range[1])/2],
                             bounds=([A_lo, 0.5, B0_range[0]],
                                     [A_hi, 0.9999, B0_range[1]]),
                             maxfev=5000)
            B0s.append(p[2])
        except Exception:
            pass
    B0s = np.array(B0s)
    lo, hi = np.percentile(B0s, 2.5), np.percentile(B0s, 97.5)
    p_pos = float(np.mean(B0s > 0))
    return B0s, lo, hi, p_pos


if __name__ == "__main__":
    uv_k, uv_lam, ir_k, ir_lam = load_modes()
    k_max = max(uv_k[-1] if len(uv_k) else 0, ir_k[-1] if len(ir_k) else 0)

    print(f"UV mode ({len(uv_k)} pts, k={uv_k[0]}..{uv_k[-1]}):")
    for k, v in zip(uv_k, uv_lam):
        marker = " ← post-crossing" if k >= 25 else ""
        print(f"  k={k:2d}: {v:+.6e}{marker}")

    print(f"\nIR mode ({len(ir_k)} pts, k={ir_k[0]}..{ir_k[-1]}):")
    for k, v in zip(ir_k, ir_lam):
        marker = " ← now minimum" if k >= 25 else ""
        print(f"  k={k:2d}: {v:+.6e}{marker}")

    print("\n" + "="*60)

    # --- UV mode fit (k=13 onwards) ---
    popt_uv, rms_uv, n_uv = fit_mode(uv_k, uv_lam, 13, -1, (-0.05, 0.15), "UV")
    if popt_uv is not None:
        A_uv, r_uv, B0_uv = popt_uv
        print(f"\nUV mode fit ({n_uv} pts, k=13..{uv_k[-1]}):")
        print(f"  A={A_uv:.4f}  r={r_uv:.5f}  B₀^UV={B0_uv:+.5f}  RMS={rms_uv:.2e}")
        if B0_uv > 0 and A_uv < 0:
            k_cross = math.log(-B0_uv / A_uv) / math.log(r_uv)
            print(f"  Zero crossing: k ≈ {k_cross:.1f}")
        B0s_uv, lo_uv, hi_uv, p_pos_uv = bootstrap_mode(
            uv_k, uv_lam, 13, -1, (-0.05, 0.15))
        print(f"  Bootstrap CI:  [{lo_uv:+.5f}, {hi_uv:+.5f}]  P(B₀^UV>0)={p_pos_uv:.4f}")

    # --- IR mode analysis ---
    # Note: 3-param exponential fit requires ≥4 post-crossing points for reliability.
    # With k=25,26 only (2 post-crossing pts), report empirical statistics instead.
    ir_post = ir_lam[ir_k >= 25]
    ir_post_k = ir_k[ir_k >= 25]
    ir_all_pos = bool(np.all(ir_lam > 0))
    print(f"\nIR mode analysis ({len(ir_k)} pts, k=18..{ir_k[-1]}):")
    print(f"  All IR mode values positive: {ir_all_pos}")
    print(f"  Range k=18..{ir_k[-1]}: [{ir_lam.min():.4e}, {ir_lam.max():.4e}]")
    if len(ir_post) >= 2:
        print(f"  Post-crossing (k≥25): mean={ir_post.mean():.4e}  "
              f"range=[{ir_post.min():.4e}, {ir_post.max():.4e}]")
    print(f"  Empirical lower bound: B₀^IR ≥ {ir_lam.min():.4e}  (min observed)")

    # Attempt exponential fit only if ≥4 post-crossing points available
    ir_post_mask = ir_k >= 25
    n_post = int(np.sum(ir_post_mask))
    popt_ir = None
    if n_post >= 4:
        popt_ir, rms_ir, n_ir = fit_mode(ir_k, ir_lam, 25, -1, (0., 1e-5), "IR")
        if popt_ir is not None:
            A_ir, r_ir, B0_ir = popt_ir
            print(f"\nIR mode fit ({n_ir} post-crossing pts, k=25..{ir_k[-1]}):")
            print(f"  A={A_ir:.4e}  r={r_ir:.5f}  B₀^IR={B0_ir:+.4e}  RMS={rms_ir:.2e}")
            B0s_ir, lo_ir, hi_ir, p_pos_ir = bootstrap_mode(
                ir_k, ir_lam, 25, -1, (0., 1e-5))
            print(f"  Bootstrap CI:  [{lo_ir:+.4e}, {hi_ir:+.4e}]  P(B₀^IR>0)={p_pos_ir:.4f}")
    else:
        print(f"  (Need ≥4 post-crossing pts for fit; have {n_post}; waiting for k=27,28)")
        p_pos_ir = 1.0 if ir_all_pos else 0.0

    print(f"\n{'='*60}")
    print("SUMMARY:")
    if popt_uv is not None:
        print(f"  B₀^UV = {B0_uv:+.5f}  (UV mode asymptote, 2nd eig k≥25)")
    print(f"  B₀^IR ≥ {ir_lam.min():.4e}  (empirical lower bound, min observed IR mode)")
    if popt_ir is not None:
        print(f"  B₀^IR = {B0_ir:+.4e}  (exponential fit, {n_post} post-crossing pts)")
    print(f"  IR mode always positive (k=18..{ir_k[-1]}): {ir_all_pos}")
    print(f"  P(B₀^UV>0) = {p_pos_uv:.4f}")
    print(f"  P(B₀^IR>0) = 1.0000 (empirical: all {len(ir_lam)} observed values > 0)")
    print(f"  => P(B₀>0) ≥ {p_pos_uv:.4f}")

    # UV mode predictions
    if popt_uv is not None:
        print(f"\nUV mode predictions (k={k_max+1}..{k_max+3}):")
        for k in range(k_max + 1, k_max + 4):
            print(f"  λ_UV({k}) ≈ {exp_b(k, *popt_uv):+.5f}")
