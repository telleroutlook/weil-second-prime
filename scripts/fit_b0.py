"""
Fit B₀ (asymptote of λ_min(k)) from available submatrix_k{k}.json files.

Three mode classes after k=28:
  UV mode: λ_0 for k≤24 (negative, converging from below to B₀^UV)
           λ_1 for k≥25 (jumped positive, second eigenvalue, continuing to B₀^UV)
  IR mode: λ_1 for k≤24 (small positive, second eigenvalue, ~4e-7..2e-6)
           λ_0 for k≥25 (minimum eigenvalue, ~1.6e-6, converging to B₀^IR)
  Frontier mode: a distinct new negative mode can enter at λ_0 (observed at k=28);
           it is neither the continuing UV branch nor the IR branch.

B₀ = lim λ_min = B₀^IR (the true minimum asymptote)
B₀^UV = asymptote of the UV mode (much larger, ~0.02)

Model: λ(k) = A × r^k + B₀  (3-param exponential + constant)
Bootstrap: 10000 replicates with replacement
"""
import json, math, pathlib, sys
import numpy as np
from scipy.optimize import curve_fit

PILOTS = pathlib.Path("pilots")


def spectrum_prefix() -> str:
    """Return the corrected chain prefix, or the legacy raw-GL prefix."""
    for k in range(18, 35):
        if (PILOTS / f"submatrix_rich_k{k:02d}.json").exists():
            return "submatrix_rich_k"
    return "submatrix_k"


def load_modes():
    """Return tracked UV, IR, and new frontier-minimum modes.

    IR is identified by its O(1e-6) scale. UV is then matched to the remaining
    candidate nearest its previous value. This matters at k=28: the old UV branch
    is λ_2, while λ_0 is a new negative frontier mode. Treating λ_0 as UV would
    silently merge distinct spectral branches.
    """
    uv_k, uv_lam = [], []
    ir_k, ir_lam = [], []
    frontier_k, frontier_lam = [], []

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
        p = PILOTS / f"{spectrum_prefix()}{k:02d}.json"
        if not p.exists():
            break
        d = json.load(open(p))
        sp = d["spectrum"].get("1.0")
        if sp is None:
            continue
        top4 = sp.get("top4", [])
        lmin = sp["lambda_min"]
        candidates = [float(x) for x in top4]
        if not candidates:
            continue

        # Identify IR by scale, not by ordinal position in the spectrum.
        ir_candidates = [x for x in candidates if abs(x) < 1.0e-4]
        if not ir_candidates:
            raise ValueError(f"no IR-scale eigenvalue found at k={k}")
        ir_value = min(ir_candidates, key=abs)
        ir_k.append(k); ir_lam.append(ir_value)

        remaining = [x for x in candidates if x is not ir_value]
        if remaining:
            uv_value = (min(remaining, key=lambda x: abs(x - uv_lam[-1]))
                        if uv_lam else min(remaining))
            uv_k.append(k); uv_lam.append(uv_value)

            # A frontier minimum is λ_0 only when it belongs to neither tracked
            # branch. At k=28 this is the decisive negative mode.
            if lmin != ir_value and lmin != uv_value:
                frontier_k.append(k); frontier_lam.append(lmin)

    return (np.array(uv_k), np.array(uv_lam),
            np.array(ir_k),  np.array(ir_lam),
            np.array(frontier_k), np.array(frontier_lam))


def exp_b(k, A, r, B0):
    return A * r**k + B0


def _fit_exp_scaled(kf: np.ndarray, lf: np.ndarray,
                    A_sign: int, B0_range: tuple[float, float]):
    """Fit the exponential model with data and constants scaled to O(1).

    The IR branch is O(1e-6), while the original starting guess used A=2.
    Unscaled ``curve_fit`` finite differences can then fail to move away from
    the initial point and produce a misleading asymptote.  Scaling only A and
    B0 (not r) keeps the model unchanged while making the optimizer well-behaved
    for both the UV and IR branches.
    """
    y_scale = max(float(np.max(np.abs(lf))), np.finfo(float).tiny)
    y_scaled = lf / y_scale
    b_lo = B0_range[0] / y_scale
    b_hi = B0_range[1] / y_scale
    b_mid = (b_lo + b_hi) / 2.0

    A_lo, A_hi = (-np.inf, 0.) if A_sign < 0 else (0., np.inf)
    p0 = [A_sign * 2.0, 0.82, b_mid]
    popt_scaled, _ = curve_fit(
        exp_b, kf, y_scaled, p0=p0,
        bounds=([A_lo, 0.5, b_lo], [A_hi, 0.9999, b_hi]),
        maxfev=30000,
    )
    A_scaled, r, B0_scaled = popt_scaled
    return np.array([A_scaled * y_scale, r, B0_scaled * y_scale])


def fit_mode(ks, lams, k_min, A_sign, B0_range, label):
    mask = ks >= k_min
    kf, lf = ks[mask], lams[mask]
    if len(kf) < 4:
        print(f"  {label}: only {len(kf)} pts (need ≥4), skip fit")
        return None, None, None
    popt = _fit_exp_scaled(kf, lf, A_sign, B0_range)
    rms = float(np.sqrt(np.mean((lf - exp_b(kf, *popt))**2)))
    return popt, rms, len(kf)


def bootstrap_mode(ks, lams, k_min, A_sign, B0_range, n_boot=10000):
    mask = ks >= k_min
    kf, lf = ks[mask], lams[mask]
    n = len(kf)
    B0s = []
    for _ in range(n_boot):
        idx = np.random.randint(0, n, n)
        try:
            p = _fit_exp_scaled(
                kf[idx], lf[idx], A_sign, B0_range,
            )
            B0s.append(p[2])
        except Exception:
            pass
    B0s = np.array(B0s)
    lo, hi = np.percentile(B0s, 2.5), np.percentile(B0s, 97.5)
    p_pos = float(np.mean(B0s > 0))
    return B0s, lo, hi, p_pos


if __name__ == "__main__":
    # Bootstrap intervals are part of the pilot record; make them replayable.
    # Keep this seed explicit so concurrent k=28 watchers report identical CIs.
    np.random.seed(20260818)
    uv_k, uv_lam, ir_k, ir_lam, frontier_k, frontier_lam = load_modes()
    k_max = max(uv_k[-1] if len(uv_k) else 0,
                ir_k[-1] if len(ir_k) else 0,
                frontier_k[-1] if len(frontier_k) else 0)
    print("Grade: discovery pilot (float eigenvalues); not certificate evidence.")
    prefix = spectrum_prefix()
    if prefix == "submatrix_k":
        print(
            "WARNING: legacy submatrix_k*.json used raw-GL8 M_K centers "
            "(skip_remainder=True). Sign/mode fits from this input are "
            "quarantined until submatrix_rich_k*.json replaces them."
        )
    print("Bootstrap seed: 20260818")

    print(f"UV mode ({len(uv_k)} pts, k={uv_k[0]}..{uv_k[-1]}):")
    for k, v in zip(uv_k, uv_lam):
        marker = " ← post-crossing" if k >= 25 else ""
        print(f"  k={k:2d}: {v:+.6e}{marker}")

    print(f"\nIR mode ({len(ir_k)} pts, k={ir_k[0]}..{ir_k[-1]}):")
    for k, v in zip(ir_k, ir_lam):
        marker = " ← now minimum" if k >= 25 else ""
        print(f"  k={k:2d}: {v:+.6e}{marker}")

    if len(frontier_k):
        print(f"\nFrontier minima not assigned to UV or IR "
              f"({len(frontier_k)} pts):")
        for k, v in zip(frontier_k, frontier_lam):
            print(f"  k={k:2d}: {v:+.6e}  ← NEW NEGATIVE MODE"
                  if v < 0 else f"  k={k:2d}: {v:+.6e}")

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
    ir_b0_range = (0., 1e-5)
    if n_post >= 4:
        popt_ir, rms_ir, n_ir = fit_mode(ir_k, ir_lam, 25, -1, ir_b0_range, "IR")
        if popt_ir is not None:
            A_ir, r_ir, B0_ir = popt_ir
            print(f"\nIR mode fit ({n_ir} post-crossing pts, k=25..{ir_k[-1]}):")
            print(f"  A={A_ir:.4e}  r={r_ir:.5f}  B₀^IR={B0_ir:+.4e}  RMS={rms_ir:.2e}")
            B0s_ir, lo_ir, hi_ir, p_pos_ir = bootstrap_mode(
                ir_k, ir_lam, 25, -1, ir_b0_range)
            print(f"  Bootstrap CI ({len(B0s_ir)} fits):  "
                  f"[{lo_ir:+.4e}, {hi_ir:+.4e}]")
            print("  P(B₀^IR>0) is not reported from this bootstrap: "
                  "the fit bound already imposes B₀^IR ≥ 0.")
            if B0_ir >= 0.99 * ir_b0_range[1]:
                print("  WARNING: B₀^IR is at/near the upper fit bound; "
                      "the 4-point asymptote is under-identified.")
            if hi_ir >= 0.99 * ir_b0_range[1]:
                print("  WARNING: the bootstrap upper endpoint is at the fit bound; "
                      "treat the 4-point IR CI as upper-censored, not a sharp interval.")
    else:
        print(f"  (Need ≥4 post-crossing pts for fit; have {n_post}; waiting for k=27,28)")
        p_pos_ir = 1.0 if ir_all_pos else 0.0

    print(f"\n{'='*60}")
    print("SUMMARY:")
    if popt_uv is not None:
        print(f"  B₀^UV = {B0_uv:+.5f}  (UV mode asymptote, 2nd eig k≥25)")
    print(f"  B₀^IR ≥ {ir_lam.min():.4e}  (empirical lower bound, min observed IR mode)")
    if popt_ir is not None:
        print(f"  B₀^IR = {B0_ir:+.4e}  (bounded 4-pt model fit; upper-censored pilot estimate)")
    print(f"  IR mode always positive (k=18..{ir_k[-1]}): {ir_all_pos}")
    print(f"  P(B₀^UV>0) = {p_pos_uv:.4f}")
    print(f"  P(B₀^IR>0) = 1.0000 (empirical: all {len(ir_lam)} observed values > 0)")
    if len(frontier_k):
        print("  Global B₀>0 is NOT supported: a new negative frontier minimum exists.")
    else:
        print(f"  => P(B₀>0) ≥ {p_pos_uv:.4f}")

    # UV mode predictions
    if popt_uv is not None:
        print(f"\nUV mode predictions (k={k_max+1}..{k_max+3}):")
        for k in range(k_max + 1, k_max + 4):
            print(f"  λ_UV({k}) ≈ {exp_b(k, *popt_uv):+.5f}")
