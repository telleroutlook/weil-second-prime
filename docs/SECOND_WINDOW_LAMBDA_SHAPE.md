# Second-window λ_min(L) shape characterization

**Date:** 2026-08-08
**Scope:** shape of the Weil-form min-eigenvalue inside the second prime window
$L\in(\tfrac12\log3,\ \log2)\approx(0.5493,\,0.6931)$, two primes $2,3$ + cross term.
**Not in scope:** RH, $L\to\infty$, extrapolation past $L=\log2$. This is an
*empirical shape description of a finite-scale object inside one window*, nothing more.

---

## 0. Grade discipline (read first)

- All curve points are **pilot** grade: `numpy.eigvalsh` on the **float center** of
  the assembled Schur matrix $C=b_LF-R_\eta$ (`eig_scan_second_window.scan_L`).
- Two points are **certify-anchored**: rebuilt on the trusted outward-rounded Arb
  interval assembly (`authoritative_eig_check` → `build_C_interval`), whose interval
  min-pivot reproduces the S4 certified anchor $-0.0197$ at $L=3/5$ exactly. At those
  two points the float-center eigenvalue reproduces the certify-assembly eigenvalue to
  all printed digits, so the pilot curve is trustworthy **as an eigenvalue center**.
- **λ_min (eigenvalue) ≠ min-pivot (positivity judge).** The task asks for λ_min shape,
  so every curve here is the **eigenvalue**. Note the two judges differ in magnitude
  (e.g. $L{=}0.6$ even: eig $=-0.231$, pivot $=-0.0197$) but agree in sign here (both
  negative), so no positivity verdict is affected. Positivity itself is judged by
  min-pivot elsewhere (S5); this document does **not** re-decide positivity.
- **No finite-sample → whole-window assertion.** Between sampled $L$ the curve is
  interpolated by eye only. "Infimum at the right end" is a statement about the
  *sampled grid + monotone trend*, not a certified window-wide extremum.

---

## 1. The two sector curves (pilot, eigenvalue)

**EVEN** $N=8,\ d=16$:

| L | 0.550 | 0.552 | 0.560 | 0.580 | 0.600 | 0.620 | 0.640 | 0.650 | 0.660 | 0.670 | 0.688 | 0.690 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| λ_min | −0.093 | −0.116 | −0.161 | −0.192 | −0.231 | −0.277 | −0.282 | −0.283 | −0.286 | −0.286 | −0.280 | −0.279 |

Shape: **steep descent** from the left end, **flattening** into a very shallow basin
past $L\approx0.64$, with a barely-perceptible up-turn at the extreme right
(≤ 0.01 over 0.64→0.69). Sample argmin $-0.286$ at $L\approx0.66$–$0.67$.

**ODD** $N=7,\ d=15$:

| L | 0.550 | 0.560 | 0.580 | 0.600 | 0.620 | 0.640 | 0.660 | 0.680 | 0.690 |
|---|---|---|---|---|---|---|---|---|---|
| λ_min | −0.118 | −0.167 | −0.234 | −0.247 | −0.281 | −0.303 | −0.324 | −0.353 | **−0.369** |

Shape: **monotone decreasing across the whole window**, no flattening, no turn-up.
Infimum at the **right end** $L\to\log2^-$ and still falling.

**WINDOW** $\lambda_{\min}^{\text{window}}(L)=\min(\text{even},\text{odd})$: the odd
sector is **more negative than even at every sampled $L$**, so odd dominates the
window λ_min throughout. Window infimum among samples: $-0.369$ at $L=0.69$ (right end).

---

## 2. d-refinement (anti-Gibbs) — a real correction

The apparent even-sector *interior minimum + right-end turn-up* was tested at three
first-corrector dimensions $d=2N\in\{12,16,20\}$ at $L\in\{0.66,0.67,0.69\}$:

| L | N6 / d12 | N8 / d16 | N10 / d20 |
|---|---|---|---|
| 0.66 | −0.4917 | −0.2856 | −0.1409 |
| 0.67 | −0.4879 | −0.2860 | −0.1358 |
| 0.69 | −0.4685 | −0.2788 | −0.1417 |

- The **absolute magnitude is strongly d-dependent** (shrinks as $d$ grows: $b_L=H_d-c_L-\dots$
  rises with $H_d$). Magnitude is a $b_L$-level effect, not shape — expected, harmless.
- **The even "interior minimum" is NOT d-robust.** argmin over $\{0.66,0.67,0.69\}$
  drifts with $d$: **0.66 (N6) → 0.67 (N8) → 0.69 (N10)**, and at $N{=}10$ the tiny
  turn-up *inverts* (0.67 becomes a local max). The wiggle amplitude ($<0.01$) is within
  d-drift. **Verdict: the even interior-min / right-end turn-up is a Gibbs-type truncation
  artifact — do NOT claim it.** (This is exactly the artifact the d-refinement was for;
  it caught one.)
- **The odd right-end infimum IS d-robust:** decreasing-to-the-right holds at both
  $N{=}6$ ($-0.439\to-0.461$, 0.68→0.69) and $N{=}7$ ($-0.353\to-0.369$).

What survives d-refinement:
1. Steep left-end descent (both sectors, all d).
2. Odd monotone-decreasing → **odd infimum at the right end** $L\to\log2^-$.
3. Odd dominates → **window λ_min infimum at the right end.**
The even shallow-basin fine structure does **not** survive and is dropped.

---

## 3. Comparison with the first window

First-window even λ_min (`weil-first-prime`, task2/task4, d=10/12/14): **monotone
decreasing in $L$** at every $d$; positive at the left end ($L{=}0.35$: $+0.00078$),
crosses zero near $L_c\approx0.42$, increasingly negative to the right
(d=12: $-0.021@0.38\to-0.083@0.45$). Infimum at the **right end** of the studied range;
positivity lives only at the **left edge**.

**Cross-window shape rule (empirical, NOT a theorem):**

| feature | first window | second window (dominant = odd) |
|---|---|---|
| trend in L | monotone decreasing | monotone decreasing |
| infimum location | right end | right end $(L\to\log2^-)$ |
| most-positive end | left end | left end |
| positivity at left end | **positive** (+0.0008…+0.05) | **negative** (odd $-0.118$; even indeterminate S5) |

The two windows share the **same qualitative shape**: λ_min decreases with $L$, infimum at
the right end, positivity (if any) concentrated at the left end. The difference is a
**downward level shift**: the second window's larger $c_L$ (≈1.8 vs 1.37) and the added
negative prime-3 contribution push the whole curve down, so the left end that was
*positive* in the first window is *already negative* (odd) or *knife-edge* (even, S5) in
the second. **Same shape, shorter reach** — consistent with the S5 boundary result.

---

## 4. Honest boundaries

- These are **pilot eigenvalue curves** (2 certify-anchored). "Right-end infimum" is a
  sampled-grid + monotone-trend statement, not a certified window-wide extremum.
- The cross-window shape agreement is an **empirical regularity across two windows**, not
  a proven law and **carries no implication for $L\to\log2^+$, larger primes, or RH.**
- No positivity verdict is issued or changed here; positivity is min-pivot (S5), and the
  second window is not certified positive-definite anywhere (odd negative at left end,
  even indeterminate).
- Nothing here approaches the wall. It is window-interior shape description only.

## 5. Artifacts

- `pilots/shape_even_N8_merged.json` — even N8 curve (22 pts)
- `pilots/shape_odd_N7.json` — odd N7 curve (9 pts, fills the S5 data gap)
- `pilots/dref_even_N6.json`, `dref_even_N10.json`, `dref_odd_N6.json` — d-refinement
- `pilots/shape_certify_anchor_even.json` (even L=0.67, FULL variant only — the OFF
  build hit walltime; the FULL certify result is complete and reproduces the pilot
  eigenvalue center), `shape_certify_anchor_odd.json` (odd L=0.69)
- `scripts/shape_scan_parallel.py` — parallel pilot driver (incremental-durable, --resume)
