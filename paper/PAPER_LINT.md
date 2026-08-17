# PAPER_LINT.md — pre-submission self-audit for weil-second-prime

> **Provenance.** This file is the second-window adaptation of
> `weil-first-prime/paper/PAPER_LINT.md` (ported 2026-08-16).
> Parts I–V (P1–P54, S1–S5) are domain-agnostic and maintained in the
> first-window file; run them from there using:
> ```bash
> TEX=paper/range_termination_paper.tex  # or the draft .tex when it exists
> # then execute P1–P54, S1–S5 from weil-first-prime/paper/PAPER_LINT.md
> ```
> **This file adds Part VI (P-W1…P-W13)** — the second-window-specific defect
> catalog. All items here are mandatory and encode bugs actually observed in this
> project. Run on every `paper/*.tex` and on `docs/SECOND_WINDOW_PAPER_DRAFT.md`
> (as a prose stand-in) before external review.

---

## PART VI — weil-second-prime domain checklist (P-W1…P-W13, mandatory)

Each item must be **run** (grep / script), not eyeballed. Precedents reference
events recorded in `PLAN.md`, `CLAUDE.md`, and the certify logs.

---

### P-W1 — Full S0 (four terms) AND full two-prime S2 (all cross terms)

The Schur matrix uses:
$$S^{(0)} = S_{VV}+S_{VK}+S_{KV}+S_{KK} \quad\text{(four terms)}$$
$$S^{(2)} = c_2^2 E_2 + c_3^2 E_3 + c_2 c_3 (F_{ij}+F_{ji}) \quad\text{(full cross terms)}$$

Both `F(i,j,τ₂,τ₃)` and `F(j,i,τ₂,τ₃)` must appear; swapping order is NOT the
same because `compute_F` is not symmetric in its first two arguments.

```bash
grep -n 'S\^{(0)}\|S_{VV}\|S_{KK}\|S\^{(2)}\|compute_F\|F_{ij}\|F_{ji}' "$TEX"
```

**Manual check:** confirm all four S0 blocks and the two compute_F directions are
present. Confirm $c_2^2 E_2 + c_3^2 E_3 + c_2 c_3 (F_{ij}+F_{ji})$, not just $E_2+E_3$.

**Precedent:** `legendre_shift_2prime` has four invariant tests for `compute_F`
(`F(τ,τ)=E`, `i+j` odd→0, symmetry, cross-Cauchy-Schwarz). If the paper's formula
differs from the code, one of them is wrong.

---

### P-W2 — kappa must be computed for the SECOND window, never imported from weil-first

The second-window κ(L) = `compute_kappa(L_num, L_den, prec=128)`. For L=0.56 this
gives 2.056, not 1.255 (the first-window constant `KAPPA_FLOAT`).

**This single error produced the N=13 false-positive signal** (b_L was inflated by
0.801, flipping C[0][0] from −5.4×10⁻³ to a spurious positive).

```bash
grep -n 'kappa\|KAPPA\|1\.255\|1\.25528\|compute_kappa' "$TEX"
grep -rn 'KAPPA_FLOAT\|kappa.*1\.255\|1\.25528' scripts/ checker/ src/
```

**Pass:** no occurrence of `1.255` or `1.25528` (first-window value) in a κ context;
every κ(L) value is annotated with its L and shown to come from `compute_kappa`.

**Precedent (2026-08-15):** `eig_scan_second_window.py` line 39 imported
`KAPPA_FLOAT=1.25528305` from weil-first. All N=7..17 "positive signal" results
were immediately withdrawn and re-scanned with correct kappa.

---

### P-W3 — INDETERMINATE ≠ NON-POSITIVE (Arb interval straddles zero)

When an Arb pivot interval straddles zero (e.g. N=17 pivot(2,2)∈[−9.03×10⁻⁴, +3.27×10⁻³]),
the result is **INDETERMINATE** — neither certified positive nor certified negative.
The paper must not label this "certifiably not positive" or "negative."

The two admissible verdicts are:
- **CERTIFIABLY NOT POSITIVE:** the entire interval is strictly negative (e.g. N=15 pivot ∈ [−5.39×10⁻³, −5.38×10⁻³]).
- **INDETERMINATE:** interval straddles zero — Arb precision insufficient, not a sign verdict.

```bash
grep -n 'INDETERMINATE\|straddles\|not positive\|non-positive\|certif' "$TEX"
```

**Manual check:** for each verdict displayed, confirm which case applies and that the
prose matches the interval data exactly.

**Precedent (2026-08-16):** initial write-up of N=17 used the phrase "non-positive";
corrected to "INDETERMINATE (precision insufficient)" after reviewing the interval
[−9.03×10⁻⁴, +3.27×10⁻³]. The interval is 4.2×10⁻³ wide with center near +1.2×10⁻³.

---

### P-W4 — Bernstein-blowup bug: S_KK/S_VK internal calls must use use_bernstein=False

For the second window (a=14/25), k≥43 makes the Bernstein bound (2R)^k ≫ 1.
Any `integrate_S_KK` or `integrate_S_VK` call that internally uses `use_bernstein=True`
(the default) will produce intervals of width ±10¹⁷ for n_row+n_col > 38.

The paper must not cite any S0 result computed with the blowup bug active.

```bash
grep -rn 'use_bernstein' src/archimedean/integrator_a.py checker/fp_second/
```

**Pass:** every internal `integrate_M_K` call inside `integrate_S_KK` and
`integrate_S_VK` passes `use_bernstein=False`. Regression tests:
`TestSkkSvkBernsteinBlowupRegression` in `tests/archimedean/test_integrators.py`.

**Precedent (2026-08-15):** discovered during N=17 certify. S0 entries with n_row≥43
had interval widths ≫1 before the fix.

---

### P-W5 — Richardson mode must be labeled epistemic grade "empirical coverage"

Results certified with `use_bernstein=False` (GL-8/GL-4 Richardson remainder) are
**not** formal Bernstein-ellipse certificates. The bound 2|I₈−I₄| is an empirical
coverage over-estimate, not an analytic theorem. Any displayed certify result must
carry this label when Richardson mode was used.

```bash
grep -n 'Richardson\|use_bernstein.*False\|GL.8\|GL-8\|empirical' "$TEX"
```

**Pass:** every Richardson-mode result is labeled "empirical truncation coverage"
or "[certify, Richardson-mode]". The phrase "formally certified" may not appear
without a Bernstein or analytic-bound qualifier.

**See:** `docs/richardson_remainder_analysis.md` §5 for the full epistemic table.

---

### P-W6 — Float-negative ≠ Arb-certified-negative

When a float pilot shows λ_min < 0 but the corresponding Arb interval is too wide
to certify the sign, the paper must write "float-negative, Arb-certification needed,"
not "negative (certified)." The second window has heavier integrands; Arb intervals
can be orders of magnitude wider than the float signal.

```bash
grep -n 'negative\|float.*neg\|Arb.*pend\|interval.*exp\|radius' "$TEX"
```

---

### P-W7 — No depth=2 quick-scan verdict promoted to a stated result

Reduced-depth pilot runs (`skip_remainder=True`, `depth_2d=2`, `depth_3d=2`) are
discovery-tier only. Certify runs use `depth=4` / `depth_2d=4` / `depth_3d=3` with
full remainder computation. A quick-scan number may not appear in a certified bound.

```bash
grep -n 'depth=2\|skip_remainder\|pilot\|quick.*scan' "$TEX"
```

---

### P-W8 — Cross-window numeric contamination (κ, c_L, τ, r, b_L)

Any constant or ratio measured in the first window (L≤7/20) may not be pasted into
second-window results. Known risky constants:
- κ: 1.255 (first) vs 2.056 (second at L=0.56) — 64% difference
- c_L: 1.365 (first) vs ~1.835 (second at L=0.56)
- r (increment ratio): first-window r≈0.79–0.80; second-window r≈0.55–0.81

```bash
grep -n '1\.365\|1\.36527\|1\.255\|0\.79\|0\.80\|0\.49\|7/20\|7.*20\|0\.35' "$TEX"
```

**Manual check:** every constant has a computation attached showing it was computed
at the target L in the **second** window. No first-window value is imported.

---

### P-W9 — b_L > 0 is necessary, not sufficient; state both conditions

$b_L > 0$ ensures $H_d > c_L + \kappa$ (the Schur criterion applies). But the form
$C = b_L F - R_\eta$ must also be positive-definite. Papers must not conflate
"b_L > 0 (criterion applies)" with "form is positive-definite."

```bash
grep -n 'b_L\|b_L.*>.*0\|criterion.*appli\|sufficient\|positive.*definite' "$TEX"
```

**Pass:** every displayed b_L value is paired with either the certify pivot result
(showing C ≻ 0 or C ⊁ 0) or an explicit statement "b_L > 0 is necessary but
the certify result is pending/negative."

---

### P-W10 — Conclusion boundary: L < log 2, not L ≤ 7/20

The second-window conclusion is bounded to $L < \log 2 \approx 0.693$. The
first-window bound $L \le 7/20 = 0.35$ must not appear as a second-window result.

```bash
grep -n 'log.2\|\\log2\|0\.693\|7/20\|7.*20\|0\.35\b' "$TEX" | grep -iv 'first.window\|weil.first\|cite\|bib'
```

**Manual check:** any displayed conclusion interval uses (log(3)/2, log(2)), not
the first-window interval. $L = \log 2$ is a hard boundary: p=2 exits the single-hop
regime, and Theorem 3 requires genuine extension. Do not extrapolate past log(2).

---

### P-W11 — No RH / near-RH / unbridged upgrade

The only admissible conclusion is finite-scale Weil positivity for $L$ in the
second prime window, $L < \log 2$. Two gaps separate this from RH: (1) finite-scale
→ full-interval positivity (no known path), (2) Weil positivity → RH (Weil equivalence
requires the full window $L \to \infty$, which this project does not claim).

```bash
grep -n 'RH\|Riemann Hypothesis\|near RH\|implies.*RH\|critical line\|all zeros\|Weil.*equiv' "$TEX" | grep -iv '%\|bib\|cite\|disclaimer'
```

**Manual check:** every hit is either a disclaimer or a literature citation, never
an asserted consequence of this paper's theorems.

---

### P-W12 — N-convergence narrative: extrapolation must be graded discovery-tier

The float N-convergence table (N=7..17) shows a positive increment trend but the
increment ratio r fluctuates (0.545–0.812). Any extrapolation to $\lambda_\infty$
from this table is **discovery-tier** and must be labeled as such. Specifically:
- A "geometric extrapolation" is valid only in the range where λ_N is already positive.
- The table shows λ_N < 0 throughout (N=7..17); a limit extrapolation that claims the
  limit is positive is unsupported by the certified data.

```bash
grep -n 'extrap\|geometric.*limit\|lambda.*infty\|N.*infty\|converge.*positive' "$TEX"
```

**Pass:** any extrapolated λ_∞ estimate is labeled "float pilot extrapolation, not
certified; see certify results in §6.3 for the authoritative sign verdict."

---

### P-W13 — Difficulty conservation (难度守恒) and no uniform lower bound over the full window

Second-window positivity at finite N does not imply full-window Weil positivity, which
does not imply RH. Two gaps remain: finite-scale → full window, and Weil positivity
→ RH. Any step that makes RH or "full-window positivity" appear to follow automatically
has evaporated difficulty — it is wrong (PROOF_CONSTITUTION Part B).

Additionally: RH ⟺ Λ=0 is a zero-margin critical system. Any derivation of a
uniform $\ge\varepsilon$ lower bound over the ENTIRE second window (all L ∈ (log3/2, log2))
would imply Λ < 0, contradicting Rodgers–Tao (Λ ≥ 0). Only finite-scale, specific-L
bounds are admissible.

```bash
grep -n '\\ge\\varepsilon\|uniform.*lower.*bound\|for all L\|whole window\|automatic\|Lambda\|Rodgers' "$TEX"
```

**Precedent:** PROOF_CONSTITUTION §B (no-放缩, C″); any uniform ≥ε floor over the
full window is the C″ tell.

---

## Running order for pre-submission (second window)

### Phase 0 — Universal checks
Run **all** items from `weil-first-prime/paper/PAPER_LINT.md` Parts I–V (P1–P54, S1–S5).
Use `TEX=paper/range_termination_paper.tex` (or the current draft `.tex`).

### Phase 1 — Second-window-specific checks (this file, in order)

| # | Item | What to run |
|---|---|---|
| P-W1 | S0 four terms + S2 full cross terms | `grep compute_F` + manual formula check |
| **P-W2** | kappa(L) must be compute_kappa, never 1.255 | `grep 1.255` + `grep KAPPA_FLOAT` |
| **P-W3** | INDETERMINATE ≠ NON-POSITIVE | `grep certif` + interval sign check |
| **P-W4** | S_KK/S_VK use_bernstein=False | `grep use_bernstein` in integrator_a.py |
| P-W5 | Richardson mode labeled empirical | `grep Richardson` |
| P-W6 | Float-negative ≠ Arb-certified | `grep negative` |
| P-W7 | No depth=2 result promoted | `grep depth=2` |
| **P-W8** | No cross-window constant contamination | `grep 1.255`, `grep 1.365` |
| P-W9 | b_L > 0 is necessary, not sufficient | `grep b_L` |
| **P-W10** | Conclusion boundary L < log2 | `grep 7/20`, `grep 0.35` |
| **P-W11** | No RH / near-RH | `grep RH` |
| P-W12 | N-convergence extrapolation graded discovery | `grep extrap` |
| **P-W13** | Difficulty conservation / no uniform lower bound | `grep uniform.*lower` |

Items in **bold** were triggered by actual bugs in this project (P-W2: kappa bug;
P-W3: N=17 INDETERMINATE misread; P-W4: Bernstein blowup; P-W8: cross-window
contamination; P-W10/W11/W13: conclusion scope).

---

## Notes for adaptation

- **P-W1 compute_F paths:** `src/prime_layer/legendre_shift_2prime.py`.
- **P-W2 compute_kappa path:** `src/archimedean/kernel.py`.
- **P-W4 regression test:** `tests/archimedean/test_integrators.py::TestSkkSvkBernsteinBlowupRegression`.
- **P-W5 epistemic table:** `docs/richardson_remainder_analysis.md §5`.
- **proofctl C12 candidate** (parameter-default blowup detection): proposed at
  `~/github/proofctl/PLAN.md` T37-13 — pending implementation upstream.
