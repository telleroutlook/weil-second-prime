# Method Boundary Analysis — Second-Window Split-Residual Schur Method

**Date:** 2026-08-18  
**Grade:** analysis (float pilot data + certified anchors)  
**Status:** Working document; **RAW-GL PILOT SIGN QUARANTINE (2026-08-18)**. The earlier k=18..28 chain used `integrate_M_K(..., skip_remainder=True)`. Its high-degree M0 centers are not sign-safe: raw GL8 gave `M0(P55,P55)=0.422735`, while the focused Richardson-Arb enclosure is `[0.012416,0.012561]`. A focused audit of the integer witness `v=3P53+P55` gives `v^T C v ∈ [+0.2292033466,+0.3953780542]`, entirely positive. Thus the reported k=28 negative frontier mode and all legacy branch/B0 fits are quarantined; they are not evidence about positivity. The focused audit does not prove the full matrix positive definite. N=25 completed indeterminate at pivot `(1,1)` with interval `[-2.353866e-2,+5.044900e-2]`; N=27 is running.

---

## 1. Setup

The Weil quadratic form positivity criterion at parameter $L$ reduces to showing
that the Schur matrix

$$C = b_L \cdot F - R_\eta, \qquad b_L = H_d - c_L - \kappa(L)$$

is positive definite, where:
- $d = 2N+1$ (odd sector) or $2N$ (even sector); $N$ = basis truncation
- $H_d = \sum_{k=1}^d 1/k$ (harmonic number)
- $c_L = \log(2\pi L) + \gamma_E$ (logarithmic constant, $\gamma_E \approx 0.5772$)
- $\kappa(L) = $ `compute_kappa(L_num, L_den, prec=128)` (Weil kernel maximum)
- $R_\eta = (1+\eta)R_0 + (1+1/\eta)R_2$ (split residual, $\eta = 1/2$ default)
- $R_0 = S^{(0)} - M_0^T G^{-1} M_0$, $R_2 = S^{(2)} - M_2^T G^{-1} M_2$

**Critical structure**: $C \succ 0$ requires (a) $b_L > 0$ (criterion applies) AND
(b) the remaining contribution $b_L \cdot F - R_\eta \succ 0$.

---

## 2. First vs. Second Window: The Kappa Gap

The fundamental difference between the windows is in $\kappa(L)$:

| Window | L example | $\kappa(L)$ | $c_L$ | Threshold $c_L + \kappa$ | $N_{\min}$ (odd, $b_L > 0$) |
|--------|-----------|-------------|-------|---------------------------|------------------------------|
| First  | 0.35      | 1.255       | 1.365 | 2.620                     | 5                            |
| Second | 0.56      | 2.056       | 1.835 | 3.891                     | 13                           |
| Second | 0.60      | 2.215       | 1.904 | 4.119                     | 17                           |
| Second | 0.62      | ~2.27       | ~1.94 | ~4.21                     | ~21                          |

The second window requires $H_d > 3.89$ at minimum, meaning $d \geq 27$ ($N \geq 13$
for odd sector). In contrast, the first window ($L = 7/20$) reaches $b_L > 0$ at
$N = 5$ ($d = 11$, $H(11) = 3.020 > 2.620$).

**Root cause of the gap**: In the second window, both primes $p=2,3$ lie in the
single-hop regime, so $\kappa(L) = \kappa_2(L) + \kappa_3(L)$ accumulates
contributions from two prime layers. In the first window, only $p=2$ contributes.
The 64% increase in $\kappa$ (1.255 → 2.056 at $L=0.56$) directly raises the
harmonic threshold by 1.271 units, pushing the effective $N_{\min}$ from 5 to 13.

---

## 3. Certify Data: L = 0.56, Odd Sector

## 3. Certify Data: L = 0.56, Odd Sector

### 3.1 Legacy float chain (raw GL8 M_K; forensic only)

**Quality warning:** rows k=18..28 below were extended with `skip_remainder=True`.
The kappa constant was corrected, but the raw GL8 quadrature center omits the
truncation remainder. Do not use these eigenvalues, branch labels, sign changes,
or B0 fits as evidence. They are retained only to document the source of the
false k=28 signal.

| N | d | b_L | eig_full (η=0.5) | eig_full (η=1.0) | Δ (η=0.5) | r² |
|---|---|-----|---------|---------|---|-----------------|
|  7 | 15 | −0.573 | −0.5747 | — | — | — |
|  9 | 19 | −0.343 | −0.3895 | — | +0.185 | — |
| 11 | 23 | −0.157 | −0.2770 | — | +0.113 | 0.607 |
| 13 | 27 | +0.000 | −0.1879 | — | +0.089 | 0.792 |
| 15 | 31 | +0.136 | −0.1155 | **−0.09899** | +0.072 | 0.812 |
| 17 | 35 | +0.256 | −0.0760 | **−0.0628** | +0.040 | 0.658 |
| 19 | 39 | +0.362 | **−0.04280** | **−0.03764** | +0.025 | — |
| **20** | **41** | **+0.412** | **−0.03144** | **−0.02741** | — | — |
| **21** | **43** | **+0.459** | **−0.02150** | **−0.01883** | — | — |
| **22** | **45** | **+0.504** | **−0.01209** | **−0.01066** | — | — |
| **23** | **47** | **+0.547** | **−0.00580** | **−0.00579** | — | — |
| **24** | **49** | **+0.588** | **−0.00121** | **−0.00068** | — | — |
| **25** | **51** | **+0.628** | **+1.47e-6** | **+1.54e-6 *** | — | — |
| **26** | **53** | **+0.666** | **+1.58e-6** | **+1.64e-6** | — | — |
| **27** | **55** | **+0.702** | **+1.70e-6** | **+1.76e-6** | — | — |
| **28** | **57** | **+0.738** | **INVALID** | **INVALID** | — | — |

*k=25 (2026-08-17): **ZERO CROSSING** — n_neg=0 at ALL η∈{0.5,0.65,0.8,0.9,1.0,1.1,1.22}. Only η=2.0 still n_neg=1.*
*k=26 (2026-08-17): All η∈[0.5,1.22] n_neg=0. λ_min is IR mode ~1.6e-6; UV mode (λ₁) = +7.43e-3.*
*k=27 (2026-08-18): **ALL η n_neg=0** (including η=2.0 for first time). λ_IR=+1.763e-6; λ_UV=+12.29e-3.*
*k=28 (2026-08-18, corrected): the apparent spectrum λ₀=-0.181816, λ₁=+1.865502e-6, λ₂=+0.012973 was a raw-GL8 artifact. A focused Richardson-Arb audit contradicts the negative witness; see §3.4.*

**N=15 full-float validation** (2026-08-16, bwgqjkw5z, 10838s):
Full 15×15 matrix from scratch (not sub-matrix): λ_min(η=1.0)=**−0.098986**, Frobenius η*=1.1715 gives
λ_min=−0.100778. Sub-matrix k=15 gives −0.098986. **Perfect agreement** — validates sub-matrix
extraction approach. Best η=1.0 (not Frobenius η*), same pattern as N=17.

### 3.1b Sub-matrix sweep and N-convergence analysis (2026-08-16, updated)

**Method**: Extracting the top-left k×k block from the N=17 Arb checkpoint gives the exact
k-dimensional Galerkin system (entries $M_0[a,b]$, $S_0[a,b]$ are independent of basis size N).
This provides a full sweep from k=3 to k=17 with a single 5-second computation.

**Full eigenvalue spectrum vs k** (η=1.0, L=0.56, odd sector):

| k | b_L | n_neg | λ₀ | λ₁ | λ₂ | λ₃ |
|---|-----|-------|-----|-----|-----|-----|
| 12 | −0.075 | 12 | −0.2022 | −0.0862 | −0.0487 | −0.0383 |
| 13 | +0.000 | 8 | −0.1604 | −0.0633 | −0.0351 | −0.0254 |
| 14 | +0.070 | 6 | −0.1282 | −0.0437 | −0.0181 | −0.0153 |
| 15 | +0.136 | 5 | −0.0990 | −0.0244 | −0.0065 | −0.0046 |
| 16 | +0.198 | 3 | −0.0793 | −0.0114 | −3e−7 | +0.0012 |
| 17 | +0.256 | 2 | −0.0628 | −0.0009 | +6e−7 | +0.0034 |
| 18 | +0.310 | **1** | **−0.0489** | +0.004 | +0.009 | +0.016 |
| **19** | **+0.362** | **1** | **−0.0376** | +7e−7 | +0.006 | +0.016 |

Key observations:
- At k=12 (b_L<0): **all 12 eigenvalues negative**. The diagonal b_L·F term is negative, driving the entire form negative.
- At k=13 (b_L≈0): **n_neg drops from 12 to 8** — a phase transition as b_L crosses zero.
- k=17→18: n_neg drops 2→**1**; **λ₁ crosses zero**, confirming B₁=+0.045 (positive). Only λ₀ remains negative.
- λ₂ passes through zero at k≈16, consistent with the LDL^T pivot(2,2) straddling zero.

**Convergence rates (λ₀ and λ₁ at k=13..19):**

| quantity | k=13 | k=16 | k=18 | k=19 | k=20 | k=21 | **k=22** | **k=23** | **k=24** | **k=25** | trend |
|----------|------|------|------|------|------|------|---------|---------|---------|---------|-------|
| λ₀ | −0.1604 | −0.0793 | −0.0489 | −0.0376 | −0.02741 | −0.01883 | **−0.01066** | **−0.00579** | **−0.00068** | **+1.5e-6** | → 0↗ |
| λ₀ ratio | — | — | 0.779 | 0.770 | 0.728 | 0.687 | **0.566** | **0.544** | **0.117** | **−0.002** | → 0 ↘ |
| k×λ₀ | −2.09 | −1.27 | −0.88 | −0.72 | −0.548 | −0.395 | **−0.234** | **−0.133** | **−0.016** | **≈0** | → 0↗ |

### 3.2 Withdrawn raw-GL convergence and B0 fits

The acceleration diagnostic, UV/IR branch assignment, zero-crossing claim, and
all B0/bootstrap fits derived from k=18..28 legacy rows are **withdrawn**. The
input M0 rows omitted the Richardson truncation remainder, so neither signs nor
mode ordinals are reliable. In particular, the former statements
"zero crossing confirmed at k=25", "B0_UV=+0.02361", and
"B0_IR=+8.4655e-6" are no longer evidence claims.

`fit_b0.py` now prints a quarantine warning for legacy `submatrix_k*.json`;
once corrected rows exist it will prefer `submatrix_rich_k*.json` without mixing
the two sequences.

### 3.3 Focused Richardson-Arb audit of the false k=28 signal

The false raw-GL k=28 minimum suggested the integer Rayleigh witness
$v=3P_{53}+P_{55}$. A focused Arb run computed only the required two M0 columns
and the three S0 entries, with GL-8/GL-4 Richardson remainder:

```text
v^T C v ∈ [+0.2292033466, +0.3953780542]
rayleigh_witness_status = nonnegative
certified_not_positive_definite = false
```

This does **not** prove the full N=28 matrix positive definite. It does refute
the proposed negative witness and exposes the magnitude of the raw-GL error:

| M0 entry | legacy raw-GL center | focused Richardson-Arb center/radius |
|---|---:|---:|
| `P53,P53` | `0.050151` | `0.012955 ± 0.000144` |
| `P53,P55` | `-0.138160` | `0.009174 ± 0.000127` |
| `P55,P55` | `0.422735` | `0.012488 ± 0.000073` |

The high-order quadrature error, not a new spectral branch, produced the apparent
k=28 negative mode.

<!-- WITHDRAWN RAW-GL BLOCK BEGIN (2026-08-18):
     Everything through the "Zero crossing ... confirmed" sentence below used
     skip_remainder=True M0 rows. It is forensic history only. -->
**Withdrawn historical B₀ sensitivity analysis** (raw GL8 input):

| Fit window | n pts | B₀ (exp+B) | r | AIC | k=25 pred |
|------------|-------|-----------|---|-----|-----------|
| k≥10 | 8 | **+0.031** | 0.833 | −97.1 | +0.010 |
| k≥11 | 7 | +0.024 | 0.824 | −84.5 | +0.005 |
| k≥12 | 6 | +0.003 | 0.796 | −80.1 | −0.007 |
| k≥13 | 5 | −0.001 | 0.788 | −65.3 | −0.010 |
| k≥14 | 4 | −0.019 | 0.738 | −55.2 | −0.023 |

**Key insight**: Including k=10..12 (where b_L<0, qualitatively different regime) pushes
B₀ to +0.031. The k≥13 window (b_L>0, physically consistent) gives B₀≈−0.001.
The b_L<0 regime should NOT be mixed with the b_L>0 convergence model.

**Bootstrap 95% CI for B₀^UV** (tracked UV mode, k=13..28, 10000 replicates, seed 20260818):
**B₀^UV ∈ [+0.02028, +0.02917] (95%) — entirely positive.** P(B₀^UV>0) = **1.0000**.
16-pt best fit: A=−2.031, r=0.831, **B₀^UV=+0.02361**, RMS=1.33e-3.
Historical 15-pt fit (k=13..27): B₀^UV=+0.02319, CI=[+0.01916,+0.03009], RMS=1.35e-3.
Earlier (k=13..26, 14-pt, **UV mode only**): P(B₀^UV>0) = 1.0000, CI=[+0.01790,+0.02751].
Earlier (k=13..25, 13-pt, **UV mode only**): P(B₀^UV>0) = 1.0000, CI=[+0.0179,+0.0276].
Earlier (k=13..22, 10-pt): P(B₀>0) = 0.9957, CI=[+0.006,+0.034].
Trend: 0.44 → 0.780 → 0.932 → 0.9957 → **1.0000**.

**B₀^UV > 0 confirmed at P=1.000.** Zero-crossing: UV mode crosses zero between k=24 and k=25, confirmed empirically (λ_UV(24)=−0.00068, λ_UV(25)=+3.25e-3).

**k=28 frontier-mode decomposition (discovery only).** In the η=1 minimum eigenvector, the largest weights are on $P_{53}$ ($|v|^2=0.8896$) and the new $P_{55}$ ($|v|^2=0.0956$). The Rayleigh split is:

$$\lambda_0 = (+2.07336)_{\text{block }P_1..P_{53}}
          + (-4.47559)_{P_{53}\text{--}P_{55}\text{ cross}}
          + (+2.22041)_{P_{55}\text{ diagonal}} = -0.181816.
$$

All diagonal entries remain positive in this float reconstruction; the negative value is caused by the new high-degree coupling. A four-point IR-only exponential fit gives B₀^IR=+8.4655e-6 with bootstrap interval [+2.4882e-6,+1.0000e-5], but the upper endpoint equals the imposed fit bound. Treat it as an upper-censored discovery estimate, not a sharp asymptote. It also does not address the new negative frontier mode.

**Zero-crossing prediction** (geometric model, k=13..19):
- λ₀(k) ≈ 3.68 × 0.786^k → 0 geometrically (B₀=0 forced)
- λ₀(25) ≈ −0.009, λ₀(30) ≈ −0.003, λ₀(35) ≈ −0.001
- Stays negative for all finite k; approaches 0 exponentially
- If B₀ = +0.004 (3-param best fit): λ₀ crosses zero at k≈22

**N=19 discriminating test — RESULT (2026-08-16 16:48):**

| N | Geometric r=0.790 (predicted) | UV-cross constant B₀=−0.029 (predicted) | **Actual** | Verdict |
|---|-------------------------------|------------------------------------------|------------|---------|
| 19 | −0.039 | −0.047 | **−0.03764** | **Geometric wins** (3.4% off vs 25% off) |

- **Rate r(18→19) = 0.7695** — geometric predicted 0.7897 ✓; component predicted 0.8885 ✗
- **η_opt = 0.80**, λ₀(η_opt) = −0.03640; n_neg = 1 (still 1 negative eigenvalue)
- **Geometric model (B₀≈0) is the favored interpretation** at k=19.

**UV-cross at k=19 = −0.02569** — dropped 10.7% from −0.02876 at k=18. **At k=20 = −0.02395** (further 6.7% drop). UV-cross monotone decay confirmed through k=20; "UV-cross constant at −0.029" hypothesis fully rejected.

Updated k=21..28 predictions using **13-pt exp+B model** (A=−2.347, r=0.820, B₀=+0.0177):

| k | Pred λ₀(13-pt) | **Actual** | UV-cross pred |
|---|----------------|------------|----------------------|
| 20 | — | **−0.02741** | **−0.02395 (actual)** |
| 21 | — | **−0.01883 (actual)** | **−0.02169 (actual)** |
| 22 | — | **−0.01066 (actual)** | **−0.01947 (actual)** |
| 23 | **−0.007** | **−0.00579 (actual)** | **−0.01933 (actual)** |
| 24 | **−0.002** | **−0.00068 (actual)** | **−0.01701 (actual)** |
| **25** | **+0.001** | **+1.54e-6 (ZERO CROSSING ✓)** | n/a (n_neg=0) |
| **26** | **+0.004** | **+1.64e-6 (IR mode; UV mode=+7.43e-3)** | — |
| 27 | +0.007 | pending | — |
| 28 | +0.009 | pending | — |

Zero crossing: k ≈ 24.6 (legacy fit_b0.py). This former confirmation is **withdrawn**.
<!-- WITHDRAWN RAW-GL BLOCK END -->

### 3.1d UV-mode decomposition and rate transition (2026-08-16)

**Method**: Decompose $\lambda_0 = v^T C v$ into three components using the min-eigenvector $v$:
- **UV-diag**: $C_{kk}|v_{UV}|^2$ (diagonal of UV basis function P_{2k-1})
- **UV-cross**: $2v_{UV}\sum_{j<k} C_{jk}v_j$ (coupling between UV and all interior modes)
- **IR-block**: $\sum_{j,l<k} C_{jl}v_jv_l$ (interior–interior quadratic form)

**Results (η=1.0, k=13..19):**

| k | λ₀ | UV-diag | UV-cross | IR-block | |v_UV|² | loc-3 |
|---|-----|---------|---------|---------|--------|-------|
| 13 | −0.1604 | −0.0056 | −0.0279 | −0.1269 | 0.122 | 0.299 |
| 14 | −0.1282 | −0.0057 | −0.0296 | −0.0928 | 0.160 | 0.326 |
| 15 | −0.0990 | −0.0057 | −0.0300 | −0.0632 | 0.210 | 0.429 |
| 16 | −0.0793 | −0.0050 | −0.0306 | −0.0436 | 0.256 | 0.479 |
| 17 | −0.0628 | −0.0038 | −0.0289 | −0.0301 | 0.291 | 0.538 |
| 18 | −0.0489 | −0.0021 | −0.0288 | −0.0180 | 0.338 | 0.605 |
| **19** | −0.0376 | −0.0008 | −0.0257 | −0.0111 | 0.363 | 0.633 |
| **20** | **−0.02741** | **+0.00165** | **−0.02395** | **−0.00511** | **0.377** | **0.663** |
| **21** | **−0.01883** | **+0.00251** | **−0.02169** | **+0.000349** | — | — |
| **22** | **−0.01066** | **+0.00474** | **−0.01947** | **+0.00408** | — | — |
| **23** | **−0.00579** | **+0.00658** | **−0.01933** | **+0.00695** | — | — |
| **24** | **−0.00068** | **+0.00814** | **−0.01701** | **+0.00819** | — | — |
| **25** | **+1.5e-6** | — | — | — | — | — |

**IR-block sign flip at k=21**: IR-block crossed zero (k=20: −0.00511 → k=21: +0.000349 → k=22: +0.00408 → k=23: +0.00695 → k=24: +0.00819). IR-block growing steadily positive.

**UV-diag continuing positive growth**: k=20: +0.00165, k=21: +0.00251, k=22: +0.00474, k=23: +0.00658, k=24: +0.00814. UV-diag is a positive and growing contribution.

**UV-cross monotone decay confirmed through k=24**: −0.02395(k=20), −0.02169(k=21), −0.01947(k=22), −0.01933(k=23), −0.01701(k=24). Decay ratio ≈ 0.90–0.89/step.

**Zero crossing mechanism (k=24→25)**: At k=24: UV-diag+IR = +0.00814+0.00819 = +0.01633 vs |UV-cross| = 0.01701. Net: −0.00068 (just barely negative). At k=25: n_neg=0 for ALL η∈[0.5,1.22] — the sum tips positive. λ₀(k=25)=+1.54e-6.

**k=25 UV decomp**: returns zeros because n_neg=0 (no negative eigenvector to decompose).

**UV-diag**: **SIGN FLIP at k=20!** At k=19: −0.0008 (C_UV≈−0.002), at k=20: +0.00165 (C_UV=+0.00437). The UV diagonal element b_L·F[UV,UV]−R_η[UV,UV] crossed zero between k=19 and k=20. UV-diag now contributes POSITIVELY to λ₀, slightly counteracting the UV-cross term.

**UV-cross series k=13..20**: −0.0279, −0.0296, −0.0300, −0.0306, −0.0289, −0.0288, −0.0257, **−0.0239**. Per-step ratios: 1.061, 1.014, 1.019, 0.945, 0.994, 0.893, **0.932**. Both k=19 and k=20 confirm steady decline; UV-cross is NOT constant — it is decaying. IR-block ratio 19→20 = 0.459.

**k=20 resolves the discriminator**: UV-cross(20) = −0.0239 ≈ −0.023 (B₀=0 prediction was −0.023). The B₀<0 scenario required UV-cross ≥ −0.025; actual is −0.024. **B₀=0 scenario confirmed by k=20.**

**Rate transition (updated with k=20)**: r(19→20) = 0.728, FASTER than geometric (0.790). Component breakdown at k=20: UV-diag −6% (positive, subtracts), UV-cross 87%, IR-block 19%. IR-block decaying rapidly.

- **Confirmed**: UV-cross → 0 at rate ≈ 0.93/step (k=19→20); B₀ ≥ 0
- **UV-diag sign flip**: C[UV,UV] crossed zero between k=19 and k=20; this UV-diagonal contribution now partially compensates UV-cross, but the net λ₀ keeps decreasing
- **Convergence accelerating**: r(17→18)=0.779, r(18→19)=0.728, r(19→20)=0.728 — stabilized near 0.73 (faster than initial 0.79)

**Localization trend**: top-3 squared norm grows 0.299 → **0.663** across k=13..20, suggesting ~0.75 by k=25.

### 3.1e UV column structure and sign mechanism (2026-08-16)

**Key finding**: The UV column C[:,UV] is **entirely negative** for all j < UV at every k tested.
This forces UV-cross < 0 independently of the eigenvector signs.

**Why C[j,UV] < 0 for all j**: Decompose $C_{j,UV} = b_L F_{j,UV} - R_\eta(j,UV)$.
The term $R_\eta(j,UV) = (1+\eta)R_0(j,UV) + (1+1/\eta)R_2(j,UV)$ is positive for all near-UV entries
and **dominates** $b_L F_{j,UV}$ by a factor of 4–20× at k=18,19. So $C_{j,UV} \approx -R_\eta(j,UV) < 0$.

| k | b_L | ratio $R_\eta / (b_L F)$ at d=1 | C[UV-1,UV] |
|---|-----|----------------------------------|------------|
| 13 | 0.0002 | ≫ 1000 | −0.0296 |
| 15 | 0.136 | 10.6 | −0.0232 |
| 17 | 0.256 | 5.0 | −0.0180 |
| 18 | 0.310 | 5.8 | −0.0179 |

As b_L increases with k, the ratio shrinks from ≫1000 toward ~5, but R_η still dominates at k=18,19.

**Layer contribution breakdown** (k=18, η=1.0):

| d=UV−j | C[UV-d,UV] | v ratio | contribution | note |
|--------|-----------|---------|-------------|------|
| 1 | −0.01789 | 0.701 | −0.00847 | nearest-UV coupling |
| 2 | −0.01072 | 0.549 | −0.00397 | |
| 3 | −0.00716 | 0.446 | −0.00216 | |
| 4 | −0.00902 | 0.422 | −0.00257 | |
| 5 | −0.01041 | 0.394 | −0.00277 | near-UV block |
| 6..17 | varies | decaying | −0.00882 | far (shrinking) |
| **total** | | | **−0.02876** | |

**Near-UV block evolution** across k (d=1..5 sum):

| k | near-UV sum | far sum | UV-cross |
|---|------------|---------|---------|
| 13 | −0.0162 | −0.0117 | −0.0279 |
| 15 | −0.0198 | −0.0102 | −0.0300 |
| 17 | −0.0196 | −0.0093 | −0.0289 |
| 18 | −0.0199 | −0.0088 | −0.0288 |

**Asymptotic behaviour (revised after k=19)**: Near-UV sum appeared to stabilize at −0.020 for k=15..18, but k=19 shows it dropped to −0.0177 (−11%). The far sum also shrinks. Therefore:

$$\text{UV-cross}(k) \xrightarrow{k\to\infty} 0 \quad (\text{geometric decay, both near-UV and far sums decaying})$$

**The earlier bound B₀ ≤ −0.020 is retracted.** k=19 breaks the plateau assumption; the near-UV block itself is decaying. B₀ is consistent with 0 (PSD limit).

**Eigenvector sign argument**: The minimum eigenvector v₀ of C(η) has all components of the same
sign (verified k=13..19). Since C[j,UV] < 0 for all j, every term 2C[j,UV]v_j v_UV is negative.
UV-cross < 0 is not an accident of oscillating signs — it is sign-definite.

**Analytical foundation — total positivity and V-dominance (2026-08-16)**: The sign structure of M₀ separates cleanly into archimedean (V) and prime-shift (K) parts:

| Component | neg entries (k=13) | property |
|-----------|-------------------|---------|
| V_part = ⟨VP_a, VP_b⟩ | 0/169 | ALL positive (V(x) > 0 everywhere) |
| K_part = ⟨KP_a, P_b⟩ | 67/169 | has negatives (prime-shift oscillations) |
| M₀ = V_part + K_part | 0/169 | ALL positive (V dominates K) |
| SVV = ⟨VP_a, VP_b⟩ | 0/169 | ALL positive |
| SVK = ⟨VP_a, KP_b⟩ | 79/169 | has negatives |
| S₀ (all terms) | 0/169 | ALL positive (SVV dominates SVK) |

At the (0,0) entry: V_part = 0.4268, K_part = −0.0476 (K is 11% of V). The archimedean contribution dominates by 9× at the strongest entry, and K_part decays rapidly off-diagonal. This gives the proof sketch: **M₀ > 0 entrywise because the archimedean V integral is positive everywhere and overwhelms the prime-shift K corrections.**

Direct inspection of the component matrices also confirms:

| Matrix | k=13 | k=15 | k=17 | k=18 | k=19 | **k=20** | Property |
|--------|------|------|------|------|------|---------|---------|
| S₀ | ALL pos | ALL pos | ALL pos | ALL pos | ALL pos | **ALL pos** | Totally positive Gram matrix |
| M₀ | ALL pos | ALL pos | ALL pos | ALL pos | ALL pos | **ALL pos** | Totally positive |
| R₀ = S₀ − M₀ᵀG⁻¹M₀ | ALL ≥ 0 | ALL ≥ 0 | ALL ≥ 0 | ALL ≥ 0 | ALL ≥ 0 | **ALL ≥ 0** | Totally non-negative |
| R_η(η=1) = 2R₀+2R₂ | ALL ≥ 0 | ALL ≥ 0 | ALL ≥ 0 | ALL ≥ 0 | ALL ≥ 0 | **ALL ≥ 0** | Totally non-negative |

Since R_η has all non-negative entries, R_η[j,UV] ≥ 0 for every (j,UV) pair. The UV column negativity C[j,UV] < 0 then follows from R_η[j,UV] > b_L F[j,UV], observed to hold at k=18,19,20. Note:
- R₂ has some negative entries (40/324 at k=18), but they are dominated by 2R₀ so R_η = 2R₀ + 2R₂ ≥ 0 entrywise.
- C is NOT a Z-matrix (3-12% positive off-diagonal entries in the IR block); the single-sign eigenvector property follows from a different mechanism (investigated but not yet analytically proved).

**Open question (updated 2026-08-16)**: Why is S₀ totally positive? M₀ totally positive? Both confirmed numerically through k=20; M₀ analytically via V-dominance (V_part > |K_part|). **Revised question**: is B₀ exactly 0 (geometric decay) or slightly positive (~+0.008, zero-crossing at k≈27)?

**Long-range sign question (updated)**: C[UV-1,UV] will flip positive near k≈84 (when b_L≈1.81). But by k=84, UV-cross will be negligible (~A×0.786^84 ≈ A×10^{-16} under geometric model). The sign flip at k=84 is irrelevant for B₀ — both λ₀ and UV-cross are exponentially small by then. The long-range sign question is resolved: it doesn't affect B₀.

**Compensation mechanism (updated)**: k×R₀[UV-d,UV] at k=19 compared to k=18:
- d=1: k=19: 0.1267 (k=18: 0.1268) — constant ✓
- d=2: k=19: 0.0967 (k=18: 0.0969) — constant ✓
- d=3: k=19: 0.0804 (k=18: 0.0806) — constant ✓

The c_d/k law holds. However, the eigenvector product v[UV-d]×v[UV] is NOT growing fast enough to maintain UV-cross constant — hence UV-cross is decaying (B₀→0).

### 3.1f UV-cross convergence mechanism: R₀ decay vs. eigenvector growth (2026-08-16)

**Observation**: R₀[UV-d, UV] decays precisely as c_d/k:

| d | k×R₀[UV-d,UV] (k=13) | k×R₀[UV-d,UV] (k=18) | Δ per step |
|---|---|---|---|
| 1 | 0.12745 | 0.12676 | −0.0001 |
| 2 | 0.09819 | 0.09691 | −0.0002 |
| 3 | 0.08237 | 0.08060 | −0.0003 |
| 4 | 0.07242 | 0.07020 | −0.0004 |
| 5 | 0.06563 | 0.06297 | −0.0004 |

The product k×R₀[UV-d,UV] is essentially constant (drift <0.2% per k-step), confirming **R₀[UV-d,UV] ~ c_d/k** to high precision.

**Compensation mechanism**: As k grows, R₀[UV-d,UV] ~ c_d/k shrinks, but the eigenvector product v[UV-d]×v[UV] grows (UV localization increasing: |v[UV]|² = 0.122→0.338 from k=13→18). The contribution 2C[UV-d,UV]×v[UV-d]×v[UV] stays approximately constant for each d because:
1. C[UV-d,UV] ≈ −R_η[UV-d,UV] (since b_L F is sub-dominant by 4–20×)
2. R_η[UV-d,UV] ~ c_d/k (same rate as R₀)
3. v[UV-d]×v[UV] grows, but C[UV-d,UV] shrinks, and their product per-d is near-constant

**Per-d contribution stability (k=13..19 with k×UV-cross diagnostic)**:

| k | UV-cross | cumsum d=1..5 | tail d>5 | **k×UV-cross** |
|---|---------|--------------|---------|---------------|
| 13 | −0.0279 | −0.0162 | −0.0117 | −0.363 |
| 14 | −0.0296 | −0.0174 | −0.0123 | −0.415 |
| 15 | −0.0300 | −0.0198 | −0.0102 | −0.451 |
| 16 | −0.0306 | −0.0197 | −0.0109 | −0.490 |
| 17 | −0.0289 | −0.0196 | −0.0093 | −0.492 |
| 18 | −0.0288 | −0.0199 | −0.0088 | **−0.518** |
| **19** | **−0.0257** | **−0.0177** | **−0.0080** | **−0.488 ↓** |

**k×UV-cross diagnostic**: If B₀ < 0 (UV-cross stabilizes at constant), k×UV-cross → ∞. If B₀ = 0 (UV-cross decays geometrically), k×UV-cross → 0. The series peaked at k=18 (−0.518) and dropped at k=19 (−0.488), suggesting k×UV-cross is bounded and possibly → 0. This strongly disfavors a stable B₀ < 0.

**Per-d at k=19 (full decomposition)**:
d=1: C=−0.01403, v=+0.413, contrib=−0.00698 (cumsum=−0.00698)
d=2: C=−0.00875, v=+0.316, contrib=−0.00334 (cumsum=−0.01032)
d=3: C=−0.00844, v=+0.272, contrib=−0.00277 (cumsum=−0.01308)
d=4: C=−0.00899, v=+0.240, contrib=−0.00261 (cumsum=−0.01569)
d=5: C=−0.00749, v=+0.222, contrib=−0.00200 (cumsum=−0.01769)
d=6..18: −0.00800 (tail, cumsum=−0.02569)

Near-UV block (d=1..5) = −0.01769 at k=19, DOWN from −0.0199 at k=18 (−11% drop).
The near-UV block is NOT stabilized — it is also decaying.

Key observations:
- **n_neg decreasing**: 8→6→5→3→2→1. At k=18, only λ₀ remains negative.
- **UV-cross**: oscillates in [−0.031, −0.028], never approaching 0.
- **cumsum d=1..5**: stabilized at −0.0196 to −0.0199 for k=15..18 (converged near-UV block).
- **tail d>5**: monotonically decreasing from −0.012 to −0.009, consistent with tail→0.

**Sign argument for B₀ < 0**: All individual per-d contributions are negative:
- 2C[UV-d,UV]v[UV-d]v[UV] < 0 for **every d** (verified k=13..19, all d up to k-1)
- Therefore UV-cross = Σ_d (negative) < 0 is a sum of same-sign terms
- No cancellation is possible; UV-cross is sign-definite negative for all finite k
- **HOWEVER**: sign-definiteness only prevents UV-cross from going positive — it does NOT prevent UV-cross → 0 from below
- Each per-d term is negative AND shrinking: UV-cross_d = 2×C[UV-d,UV]×v_{UV-d}×v_UV where C[UV-d,UV] ~ −c_d/k → 0
- The sum of shrinking negative terms can still converge to 0
- k=19 data shows UV-cross dropped to −0.02569 (10.7% below k=18), consistent with decay-to-0

**Revised conclusion**: the sign argument establishes UV-cross(k) < 0 for all finite k, but **does NOT establish B₀ < 0**. The near-UV cumsum (d=1..5) stabilized at −0.020 for k=15..18, but k=19 suggests it is also decaying. B₀ may be 0 (geometric model) or a small negative constant. The earlier bound "B₀ ≤ −0.020" was premature — it assumed the near-UV block had converged, but k=19 breaks that assumption.

**Status after k=20**: B₀ = +0.008 (8-pt fit), P(B₀>0)=0.932. Both B₀=0 and B₀≈+0.008 consistent with data; k=21..25 will discriminate. UV-cross monotone decay confirmed.

**k=19→20 actual values** (2026-08-16):
- UV-cross: k=19: **−0.02569**, k=20: **−0.02395** (ratio 0.932, monotone decline confirmed)
- UV-diag: k=20 **flipped positive** (+0.00165; C_UV=+0.00437 > 0)
- n_neg: 1 for all η at k=20 ✓ (sign-definiteness lemma holds ✓)
- Total positivity: S0, M0, R0, R_eta ALL positive/non-negative at k=20 ✓

**Implication**: UV-cross monotone decay is confirmed through k=20. The "B₀ ≤ −0.020" bound is definitively retracted; B₀ ≥ −0.004 at 95% CI. The B₀=0 vs B₀=+0.008 question requires k=25+.

### 3.2 Certify results (Arb intervals, odd, L=0.56)

| N | b_L | Certify verdict | Pivot interval |
|---|-----|-----------------|----------------|
| 15 | +0.136 | **CERTIFIABLY NOT POSITIVE** | C[0,0] ∈ [−5.39e−3, −5.38e−3] |
| 17 | +0.256 | **INDETERMINATE** | C[2,2] ∈ [−9.03e−4, +3.27e−3] |
| 19 | +0.362 | **CERTIFIABLY NOT PD** (η=0.1) | pivot (0,0) ∈ [−3.804e−02, −5.152e−04]; upper bound < 0 ✓ (6.5 s, cached ckpt) |
| 19 | +0.362 | *INDET* (η=0.5) | pivot (1,1) ∈ [−1.56e+02, +7.42e-02]; Schur blow-up: C[0,0]≈+0.02 tiny → 1/C[0,0] explodes |
| 19 | +0.362 | *INDET* (η=1.0) | pivot (0,0) ∈ [−1.54e−03, +6.67e−02]; straddles zero |

**N=19 float analysis (complete 361/361 entries, 2026-08-16 18:03):**
- η=0.5: λ_min=−0.04280, **λ₁=−0.00190** (n_neg=2); C[0,0]≈+0.02 → pivot (1,1) blows up
- η=1.0: λ_min=−0.03764, n_neg=1; λ₁=+7.3e-7 (near-null); C[0,0]≈+0.030 → pivot (0,0) straddles zero
- **η=0.1: n_neg=7; C[0,0]≈−0.021 (no near-cancellation) → pivot (0,0) certified negative (6.5 s)**
- η_opt=0.80: λ_min=−0.03640 (best from η scan); Frob η*=1.269
- 6-pt exp+B (k=13..17,19): A=−3.168, r=0.797, B₀=+0.0043 → N=25: −0.006, N=30: +0.001 (zero crossing ~N=29)
- Eigenvector top-5: P37(−0.60), P35(−0.41), P33(−0.32), P31(−0.27), P29(−0.24)

**η=0.1 certify mechanism (2026-08-16)**: R₀[0,0]=6.75e-4, R₂[0,0]=5.716e-3, b_L·F[0,0]=0.04311.
At η=0.1: R_η[0,0]=(1.1)×R₀+(11)×R₂=0.06362, so C[0,0]=0.04311−0.06362=−0.02051. No near-cancellation
(ratio 0.043/0.064=0.68, far from 1). Arb interval at prec=256: [−0.03804, −5.15e-4], upper bound < 0 → CNPD.
The (1+1/η)=11 factor amplifies R₂ interval errors (interval width=0.038 vs float=0.021), but upper bound stays
negative. Contrast: η∈[0.5, 1.0] has b_L·F[0,0]≈R_η[0,0] (ratio ~0.95–2.4) → near-cancellation → INDET.

**Near-cancellation structure**: C[0,0]=b_L·F[0,0]−R_η[0,0]. C[0,0]<0 iff R_η[0,0]>b_L·F[0,0]=0.04311.
At η=0.1: R_η[0,0]=0.0636 (>>0.0431) → well negative. Zero crossing at η≈0.14 (between float -0.001 at η=0.15 and +0.008 at η=0.20). η≥0.15 brings R_η[0,0] below b_L·F[0,0], causing near-cancellation or positive C[0,0].

### 3.1c Eigenvector structure at N=17 (2026-08-16)

At N=17, η=1.0, three eigenvectors near zero were extracted.

**Full spectrum**: only **2 negative eigenvalues** out of 17 (n_neg=2 at k=17):

| index | eigenvalue | character |
|-------|-----------|-----------|
| λ₀ | −0.0628 | UV boundary mode (dominant) |
| λ₁ | −0.0009 | mixed UV/mid mode (near-zero) |
| λ₂ | +1×10⁻⁶ | near-null IR mode |
| λ₃..λ₁₆ | +0.003 to +0.093 | positive, bulk spectrum |

**v₀ (UV boundary mode, λ₀=−0.0628)**:
All components negative, near-monotone decay from P33→P1:

| P33 | P31 | P29 | P27 | P25 | P23 | ... | P3 | P1 |
|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| −0.539 | −0.392 | −0.306 | −0.271 | −0.252 | −0.237 | ... | −0.072 | −0.096 |

Top-3 components carry 54% of squared norm (0.291+0.154+0.094=0.539).
This is the signature of a **UV truncation artifact**: the mode is peaked at the highest basis function P33 and decays monotonically toward lower degree. As N increases, the dominant component shifts to P_{2N−1}, and λ₀ decreases in magnitude.

**v₁ (mixed UV/mid mode, λ₁=−0.0009)**:
Sign change in the middle range — negative for P1 through P23, then mixed for P25-P33:

| P33 | P31 | P29 | ... | P17 | P15 | P13 | ... | P1 |
|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| +0.462 | +0.294 | +0.160 | ... | −0.332 | −0.345 | −0.368 | ... | −0.290 |

This mode has large contributions from both the UV boundary (P33,+0.462) and the mid-range (P13−P17,≈−0.34). The near-zero eigenvalue (−0.0009) and rapid convergence toward positive (λ₁→+0.045 from the 5-point fit) suggest this mode is already resolved by N=17.

**v₂ (near-null IR mode, λ₂=+1e−6)**:
Highly concentrated at low degree, alternating signs: P3(+0.784), P1(−0.442), P5(−0.427).
This is the algebraic near-null mode of the IR sector; its near-zero eigenvalue is related to the P1/P3/P5 interaction structure.

**Interpretation summary**:
- v₀ is a UV truncation artifact; UV-cross analysis (§3.1d) shows the rate is transitioning
  from r≈0.79 (crossover regime) toward a slower asymptotic rate as UV-cross dominates
- v₁ **confirmed crossed zero at k=18** (n_neg drops 2→1); B₁=+0.045 clearly positive
- v₂ is essentially zero (≈numerical precision) at k=16,17,18

**k=19 result (2026-08-16 16:48)**: λ₀(k=19) = −0.03764 — geometric model validated; UV-cross dropped to −0.02569. See §3.1b, §3.1d for full analysis.
- N=17: pivot (2,2) center ≈ +1.2e-3 at all tested η values, but interval straddles zero.
  **Arb eta scan of N=17 (from checkpoint, 6s each):**
  | η | Pivot (2,2) interval | Center | Lower bound |
  |---|---|---|---|
  | 0.5 | [-9.03e-4, +3.27e-3] | +1.18e-3 | -9.03e-4 |
  | 1.0 | [-7.88e-4, +3.59e-3] | +1.40e-3 | **-7.88e-4** (tightest) |
  | 2.49 | [-2.33e-3, +4.90e-3] | +1.28e-3 | -2.33e-3 |
  | 4.0 | [-4.24e-3, +6.22e-3] | +0.99e-3 | -4.24e-3 |
  
  The pivot center is positive at all η; the Arb issue is **precision**, not sign.

  **Float λ_min(C(η)) at N=17 — full matrix scan (2026-08-16):**
  From Arb checkpoint centers, recomputing M2/S2 in float:
  b_L=0.2555, ‖R₀‖_F=3.43e-2, ‖R₂‖_F=5.11e-2, k=1.491, Frobenius η*=1.221

  | η | λ_min(C) | note |
  |---|---|---|
  | 0.5 | -0.07372 | default |
  | 0.75 | -0.06292 | |
  | **1.0** | **-0.06278** | **global minimum (best η)** |
  | 1.221 (Frob η*) | -0.06580 | |
  | 2.49 (entry η*) | -0.09861 | entry-wise optimal → worse globally |
  | 4.0 | -0.14591 | |

  **Critical diagnosis**: The full matrix λ_min = −0.063 at best. Although pivot(2,2)
  center is ≈+1.2e-3, the LDL^T decomposition stops there — the remaining
  Schur complements carry the deeply negative eigenvalues of the full matrix.
  N=17 is clearly NOT close to positive definite; the positive pivot(2,2) center
  reflects only the leading 3×3 block behavior.

  The entry-wise η* = 2.49 is counterproductive globally (λ_min worsens). The
  global minimizer is η ≈ 1.0, close to the Frobenius η* = 1.221.

  **η_opt shift with k** (sub-matrix sweep k=13..20, 2026-08-16):
  True optimal η (maximizing λ_min) is DECREASING with k:

  | k | η_opt | λ_min(η_opt) | λ_min(η=1) | improve | Frob η* |
  |---|-------|-------------|------------|---------|---------|
  | 13 | 1.056 | −0.16027 | −0.16041 | +0.09% | 1.163 |
  | 14 | 1.018 | −0.12814 | −0.12815 | +0.01% | 1.193 |
  | 15 | 0.934 | −0.09880 | −0.09899 | +0.18% | 1.172 |
  | 16 | 0.904 | −0.07888 | −0.07925 | +0.46% | 1.206 |
  | 17 | 0.871 | −0.06216 | −0.06278 | +0.98% | 1.221 |
  | 18 | 0.800 | −0.04788 | −0.04892 | +2.1% | ~1.22 |
  | 19 | 0.800 | −0.03640 | −0.03764 | +3.3% | ~1.22 |
  | **20** | **0.800** | **−0.02609** | −0.02741 | **+4.8%** | ~1.22 |
  | **21** | **0.750** | **−0.01725** | −0.01883 | **+8.4%** | ~1.22 |
  | **22** | **0.750** | **−0.00889** | −0.01066 | **+16.6%** | ~1.22 |
  | **23** | **0.700** | **−0.00389** | −0.00579 | **+32.8%** | ~1.22 |
  | **24** | **0.700** | **+0.000001** | −0.00068 | — | ~1.22 |
  | **25** | **1.000** | **+1.54e-6** | +1.54e-6 | 0% | ~1.22 |

  η_opt shifted from 0.800→0.750 at k=21, then 0.750→0.700 at k=23. At k=24 the η_opt path reached positive. At k=25 all η∈[0.5,1.22] give positive λ_min. Improvement growing: 4.8%→8.4%→16.6%→32.8% — as λ₀→0, η-optimization increasingly effective.

- N=19: λ_min(η=0.80) = −0.03640, λ_min(η=1.0) = −0.03764, n_neg=1. The
  N-convergence ratio is ~0.617 per step (geometric B₀=0 model); zero crossing
  requires many more increments beyond N=19.

---

## 4. Method Boundary Estimate

### 4.1 Derivation

At N=15, L=0.56, odd sector (full first-column computation, 2026-08-16):

$$F[0,0] = 0.11899, \quad b_L = 0.1360, \quad b_L \cdot F[0,0] = 1.618 \times 10^{-2}$$
$$R_\eta[0,0]_{(\eta=0.5)} = 2.157 \times 10^{-2}, \quad C[0,0] = -5.39 \times 10^{-3} \checkmark$$

To make $C[0,0] = 0$ at $\eta = 0.5$:
$$b_L^* \cdot F[0,0] = R_{\eta=0.5}[0,0] = 2.157\times10^{-2} \Rightarrow b_L^* = 2.157\times10^{-2}/0.11899 = 0.1813$$

However, with optimal $\eta$:
$$b_L^* \cdot F[0,0] = R_{\eta^*}[0,0] = (\sqrt{R_0[0,0]}+\sqrt{R_2[0,0]})^2 = 1.308\times10^{-2} \Rightarrow b_L^* = 0.1099$$

The current $b_L = 0.136 > 0.1099$: **C[0,0] is already positive at N=15 with $\eta^*=2.49$!**

**Method boundary with $\eta=0.5$ (standard)**:
To make C[0,0] > 0 at $\eta=0.5$: need $b_L > 0.1813$.
$H_d^* = c_L + \kappa + b_L^* = 1.835 + 2.056 + 0.1813 = 4.072$
Solving $H_d \geq 4.072$: requires $d \geq 36$, i.e. **N ≥ 17** (odd sector) for (0,0) pivot only.
The bottleneck is the other pivots; full $C \succ 0$ requires larger N (certify in progress).

**Note**: The earlier estimate "N≈34" used incorrect $F[0,0] \approx 6.7 \times 10^{-3}$ (a factor-of-18 error).
The correct $F[0,0] = 0.119$ makes the analysis much more favorable.

### 4.2 Revised Compute Cost (2026-08-16 correction)

**With $\eta=0.5$**: (0,0) pivot requires $N \geq 17$; full $C \succ 0$ needs higher N
(empirical: N=17 certify in ~6.5h, N=19 in ~18h).

**With $\eta^* \approx 2.49$**: (0,0) pivot already positive at N=15. Full matrix
result (2026-08-16): λ_min(η=1.0) = −0.09899 at N=15; $C \not\succ 0$ at N=15.
The (0,0) pivot is positive but all other pivots are not yet.

**Earlier "N≈34 / 22-24h" estimate was wrong**: it used $F[0,0] \approx 6.7\times10^{-3}$
(a factor-of-18 error vs the correct $0.119$). All N_boundary estimates in this document
should be treated as pilot-grade lower bounds pending a full certify at larger N.

---

## 5. Why Larger L Does Not Help

A natural hypothesis: move to larger $L$ (closer to $\log 2 \approx 0.693$)
to exploit larger cross-prime coupling $J(\tau_2, \tau_3)$.

**Refutation — kappa profile across the second window:**

| $L$ | $\kappa(L)$ | $c_L$ | $c_L+\kappa$ | $N_{\min}$ (odd, $b_L>0$) | $N_{\text{boundary}}$ (est.) |
|-----|-------------|-------|--------------|--------------------------|------------------------------|
| 0.55 | 2.017 | 1.817 | 3.834 | 13 | **33** (optimal) |
| 0.56 | 2.056 | 1.835 | 3.891 | 13 | 35 |
| 0.57 | 2.096 | 1.853 | 3.949 | 14 | 37 |
| 0.58 | 2.135 | 1.870 | 4.006 | 15 | 39 |
| 0.60 | 2.215 | 1.904 | 4.120 | 17 | 44 |
| 0.62 | 2.296 | 1.937 | 4.233 | 19 | 49 |
| 0.65 | 2.418 | 1.984 | 4.402 | 23 | 58 |
| 0.68 | 2.542 | 2.029 | 4.571 | 27 | 69 |
| 0.69 | 2.584 | 2.044 | 4.628 | 28 | 73 |

$N_{\text{boundary}}$ = estimated N where $b_L \approx 0.94$ (calibrated from L=0.56 analysis
with the incorrect $F[0,0]\approx6.7\times10^{-3}$). **These estimates are provisional
and likely too pessimistic**: the corrected $F[0,0]=0.119$ lowers the required $b_L$ substantially.
With $\eta^*$-optimization, the (0,0) pivot is already positive at $N=15$ for $L=0.56$;
full-matrix N_boundary at L=0.56: sub-matrix chain shows λ_min→0 geometrically (B₀=0 model, r≈0.786). The monotone trend (larger L → harder) remains valid.

**Consequence:** Both $\kappa(L)$ and $c_L$ increase with $L$, raising the threshold
$c_L + \kappa$. The cross-prime coupling increment $\Delta\lambda$ is larger at bigger $L$
but is overwhelmed by the higher positivity threshold. The optimal L for this method
within the second window is near the **left endpoint** $L \approx \tfrac{1}{2}\log 3 \approx 0.549$.

The N_boundary column will be recomputed once the geometric decay rate is confirmed through k=25.

---

## 6. The Eta Optimization Direction — New Structural Finding (2026-08-16)

The residual $R_\eta = (1+\eta)R_0 + (1+1/\eta)R_2$ can be minimized over scalar $\eta > 0$.

### 6.1 (0,0) Near-Cancellation: A New Second-Window Structure

Computing the **full first column** of $M_0$ and $M_2$ at $N=15$, $L=0.56$ (using the
exact formula $R_0[0,0] = S_0[0,0] - \sum_k M_0[k,0]^2 / G_d[k]$):

| Quantity | Value |
|----------|-------|
| $S_0[0,0]$ | $4.084 \times 10^{-1}$ |
| $\sum_k M_0[k,0]^2/G_d[k]$ | $4.074 \times 10^{-1}$ |
| $R_0[0,0]$ | $\mathbf{1.076 \times 10^{-3}}$ (near-zero!) |
| $S_2[0,0]$ | $1.995 \times 10^{-1}$ |
| $\sum_k M_2[k,0]^2/G_d[k]$ | $1.929 \times 10^{-1}$ |
| $R_2[0,0]$ | $\mathbf{6.652 \times 10^{-3}}$ |
| $k[0,0] = R_2[0,0]/R_0[0,0]$ | **6.18** |

The archimedean Schur complement $R_0[0,0]$ is nearly zero because $S_0[0,0] \approx \sum_k M_0[k,0]^2/G_d[k]$
to 4 significant figures. This near-complete cancellation makes the prime-layer residual
$R_2[0,0]$ dominate by a factor of $\approx 6$ at the $P_1$–$P_1$ entry.

**Contrast**: At $a=7$ ($P_{15}$, diagonal-only approximation): $R_2/R_0 \approx 0.36$ —
the near-cancellation is specific to the lowest basis function $P_1$.

### 6.2 The Critical Eta Range for C[0,0] > 0

With $b_L \cdot F[0,0] = 1.618 \times 10^{-2}$, the condition $C[0,0] > 0$ requires:
$$\eta \in \bigl(0.887,\ 6.97\bigr)$$

(Solving the quadratic $(1+\eta) \cdot R_0[0,0] + (1+1/\eta) \cdot R_2[0,0] < b_L F[0,0]$.)

At $\eta = 0.5$ (current default): $C[0,0] = -5.39 \times 10^{-3} < 0$ — **outside the range**.

At $\eta = 1.0$: $R_\eta[0,0] = 2 \cdot R_0[0,0] + 2 \cdot R_2[0,0] = 1.546 \times 10^{-2} < 1.618 \times 10^{-2}$
→ $C[0,0] = +0.72 \times 10^{-3} > 0$. ✓

At $\eta^* = 2.49$: $R_\eta^*[0,0] = (\sqrt{R_0[0,0]}+\sqrt{R_2[0,0]})^2 = 1.308 \times 10^{-2}$
→ $C[0,0] = +3.10 \times 10^{-3} > 0$ (maximum at this entry). ✓

**Summary table for $C[0,0]$ at $N=15$, $L=0.56$:**

| $\eta$ | $R_\eta[0,0]$ | $C[0,0]$ | Sign |
|--------|----------------|----------|------|
| 0.5 (default) | $2.157 \times 10^{-2}$ | $-5.39 \times 10^{-3}$ | **negative** |
| 0.887 (critical) | $1.618 \times 10^{-2}$ | $0$ | zero |
| 1.0 | $1.545 \times 10^{-2}$ | $+7.3 \times 10^{-4}$ | **positive** |
| 2.49 ($\eta^*$) | $1.308 \times 10^{-2}$ | $+3.10 \times 10^{-3}$ | **positive** (max) |
| 6.97 (upper critical) | $1.618 \times 10^{-2}$ | $0$ | zero |

### 6.3 Full Matrix Eta Scan (2026-08-16 update)

**N=17 float full-matrix scan (from Arb checkpoint centers, 4 s):**

$\|R_0\|_F = 3.43\times10^{-2}$, $\|R_2\|_F = 5.11\times10^{-2}$,
Frobenius $\eta^*_F = \sqrt{\|R_2\|/\|R_0\|} = 1.221$

| $\eta$ | $\lambda_\min(C)$ |
|---|---|
| 0.5 | -0.0737 |
| 0.75 | -0.0629 |
| **1.0** | **-0.0628 (global min)** |
| 1.221 (Frob η*) | -0.0658 |
| 2.49 (entry η*) | -0.0986 |
| 4.0 | -0.1459 |
| 6.97 | -0.2439 |

**Critical finding**: The entry-wise $\eta^* = 2.49$ (which maximises $C[0,0]$, §6.2)
is counterproductive for the global minimum eigenvalue. At $\eta=2.49$, $\lambda_\min$
is **57% more negative** than at $\eta=1.0$. The global minimiser of $\|R_\eta\|$ is
the Frobenius ratio $\eta^*_F = \sqrt{\|R_2\|_F/\|R_0\|_F} \approx 1.22$, not the
scalar diagonal entry ratio $\sqrt{R_2[0,0]/R_0[0,0]} = \sqrt{6.18} \approx 2.49$.

The $(0,0)$ entry-wise optimisation is necessary but not sufficient:
$C[0,0](\eta^*) > 0$ at $N=15$ does NOT imply $C(\eta^*) \succ 0$.

**N=15 float scan (COMPLETED, bwgqjkw5z, 10838s)**: Full 15×15 matrix from scratch.
Result: $\lambda_\min(\eta=1.0) = -0.09899$, Frobenius $\eta^*=1.172$ gives $\lambda_\min=-0.10078$.
Best $\eta=1.0$ (not Frobenius $\eta^*$). **Prediction confirmed: $\lambda_\min < 0$ at all $\eta$ for $N=15$.**

The sub-matrix $k=15$ gives $-0.098986$ — perfect agreement, validating the sub-matrix approach.

### 6.4 Why the Near-Cancellation at P1?

Physically: $P_1$ is the constant function on $[-L,L]$, matching the zero mode of the
archimedean kernel. The archimedean Schur complement $R_0[0,0]$ measures the residual
orthogonality of $P_1$ after projecting onto the kernel eigenfunctions. Near-zero $R_0[0,0]$
means $P_1$ is nearly in the span of the projected archimedean kernel — a structural
feature of the Weil kernel at the scale $L \approx \frac{1}{2}\log 3$.

This near-cancellation is a **new numerical structure** specific to the second window
(two-prime regime), not observed in the first window ($L \approx 0.35$, single prime).

---

## 7. Comparison: First Window Method Boundary

For context, the first window ($L = 7/20 = 0.35$, single prime $p=2$):
- $\kappa(0.35) = 1.255$, $c_L(0.35) = 1.365$, threshold = 2.620
- N=5 already has $b_L > 0$
- Certify positive at N=8 (even) and N=6 (odd): **N ≤ 8**
- Method boundary: effectively N ≈ 8–10

The second window's method boundary (with $\eta=0.5$) is much lower than the earlier
"N≈34" estimate — that estimate used an incorrect $F[0,0]$. The corrected analysis
(§4.1) shows (0,0) pivot requires only $N \geq 17$ at $\eta=0.5$, and $N \geq 15$
already works at $\eta^* \approx 2.49$. The additional $\kappa_3$ contribution from $p=3$
still imposes a significant barrier compared to the first window, but the barrier is
lower than previously thought.

---

## 8. Conclusions (Updated 2026-08-16)

1. **Method boundary corrected**: Earlier "N≈34" estimate used $F[0,0] \approx 6.7\times10^{-3}$
   (18× error vs correct $0.119$). Corrected: (0,0) pivot requires $N \geq 17$ at $\eta=0.5$,
   or $N \geq 15$ at $\eta^* \approx 2.49$.

2. **Near-cancellation structure at (0,0)**: $R_0[0,0] \approx 1.08\times10^{-3}$ (nearly zero,
   archimedean Schur complement near-cancels). $R_2[0,0]/R_0[0,0] = 6.18$ — a new structural
   feature of the second window not present in the first window.

3. **Eta optimization: C[0,0] > 0 already at N=15, but full matrix remains negative:**
   For any $\eta \in (0.887, 6.97)$, the (0,0) pivot is positive at $N=15$. However,
   the float full-matrix $\eta$ scan at $N=17$ shows that the globally optimal $\eta$
   for $\lambda_\min$ is the Frobenius minimiser $\eta^*_F \approx 1.22$, NOT the
   entry-wise $\eta^*=2.49$. At best ($\eta=1.0$), $\lambda_\min(C) = -0.063$ at $N=17$
   — far from positive-definite. The N=15 full matrix eta scan **completed** (10838s):
   $\lambda_\min(\eta=1.0) = -0.09899$ — $\lambda_\min < 0$ at all $\eta$ (prediction confirmed).

4. **Larger L is strictly worse** (confirmed): moving toward $L = \log 2$ increases $\kappa$
   and raises the threshold. Optimal L for the current method is near $L = \tfrac{1}{2}\log 3$.

5. **T1 pilot extrapolations (L=0.62, 0.65) were invalid**: computed with $b_L < 0$;
   the "positive $\lambda_\infty$" predictions have no bearing on true convergence.

6. **Current certify results and B₀ status** (corrected 2026-08-18):
   - N=19 η=0.1: pivot(0,0)∈[−0.038,−5.15e-4] — CNPD certified.
   - N=25 η=1, prec=512, Richardson mode (completed 2026-08-20):
     pivot(1,1)∈[−2.353866e-2,+5.044900e-2] — **INDETERMINATE**. This is not
     a positivity certificate and not a certified negative pivot.
   - **Legacy dual-mode/B0 narrative withdrawn**: the k=18..28 chain omitted the
     M_K Richardson remainder. Its UV/IR/frontier labels, zero crossing, B0 fits,
     and bootstrap probabilities are not evidence.
   - **Focused Richardson-Arb audit** at N=28 for `v=3P53+P55`:
     `v^T C v ∈ [+0.2292033466,+0.3953780542]`, so this witness does not certify
     non-positive-definiteness. The full-matrix sign remains unknown from this audit.
   - `submatrix_chain.py` now uses the Richardson remainder and writes to a new
     `submatrix_rich_*` namespace; it will not silently resume legacy raw rows.

7. **Honest research narrative**: The second-window investigation documents the
   honest boundary of the split-residual Schur method when extended to two primes,
   including new structural phenomena (near-cancellation, eta-sensitive boundary).

---

## 9. Next Steps

| Priority | Action | Trigger |
|----------|--------|---------|
| **DONE** | N=19 certify η=0.5: INDET (pivot(1,1) blow-up) | `cert_fp_second_N19.json` ✓ |
| **DONE** | N=19 certify η=0.1: **CNPD** pivot(0,0)∈[−0.038,−5.15e-4] | `cert_fp_second_N19_eta01.json` ✓ |
| **DONE** | k=20..22: λ₀→0, UV structural transitions | `submatrix_k20..22.json` ✓ |
| **DONE** | k=23: λ₀=−0.00579, r=0.544 | `submatrix_k23.json` ✓ |
| **DONE** | k=24: λ₀=−0.00068, r=0.117 | `submatrix_k24.json` ✓ |
| **DONE** | **k=25: zero crossing confirmed** (n_neg=0, λ_IR=+1.54e-6, λ_UV=+3.25e-3) | `submatrix_k25.json` ✓ |
| **DONE** | **k=26: dual-mode confirmed** (λ_IR=+1.64e-6, λ_UV=+7.44e-3) | `submatrix_k26.json` ✓ |
| **DONE** | **k=27: ALL η n_neg=0** incl η=2.0; λ_IR=+1.763e-6, λ_UV=+12.29e-3 | `submatrix_k27.json` ✓ |
| **DONE** | **Raw-GL k=28 negative mode refuted for witness 3P53+P55**; legacy chain quarantined | `cert_fp_second_N28_frontier_eta1_p512.json` |
| **DONE** | `submatrix_chain.py` fixed to retain Richardson remainder and isolate corrected outputs | `tests/test_submatrix_chain_remainder.py` |
| **DONE** | Certify N=25, η=1/1, prec=512, --no-bernstein: **INDETERMINATE**, pivot(1,1) interval straddles zero | `cert_fp_second_N25_eta1_p512.json` |
| **Active** | Certify N=27, η=1/1, prec=512, --no-bernstein (checkpoint running) | `cert_fp_second_N27_eta1_p512.ckpt.json` |
| **High** | Treat N=25/N=27 as finite-scale min-pivot checks; no legacy chain narrative may override them | checker output |
| **High** | Rebuild the k=18..28 chain under `submatrix_rich_*` only if mode tracking is still needed | corrected chain |
| **High** | Add a pilot-sign firewall to proofctl candidate C13: raw-center outputs cannot enter sign narratives | upstream proposal |

---

*See also:* `PLAN.md §第四编 H1`, `docs/SECOND_WINDOW_PAPER_DRAFT.md §6.3`,
`paper/PAPER_LINT.md P-W2, P-W3, P-W8`.
