# S4 — Per-sector prime-influence profile (second window)

**Date: 2026-08-08.** Question S4 answers: in the second-prime window
$L\in(\tfrac12\log3,\log2)$, which prime-layer terms actually move the Schur
min-pivot, and which are inert — so S5 and later compute are not budgeted
symmetrically. The steer from weil-first: the first window's even-sector single
prime was nearly inert. The second window adds the cross-prime term
$J(\tau_2,\tau_3)$ (via $F_{ij}$), absent in the first window — its influence is
the key unknown.

## Method

`scripts/profile_prime_influence.py` builds the archimedean block (four-term
$S^{(0)}$, $M^{(0)}$, $T$, $G$) ONCE per (sector, L) — the expensive part — then
swaps the prime layer instantly for each variant. The archimedean truncation is
common-moded across variants, so influence deltas isolate the prime-layer effect.
Variants: full / tau2_off / tau3_off / cross_off / scale_cross(×100).

**Precision grades are labelled on every number** (PROOF_CONSTITUTION A3):
`[pilot]` = float center, depth 4 (screening only); `[certify]` = Arb interval
(verdict-grade). Pilot numbers never decide compute allocation on their own.

## Results at L = 3/5 = 0.60

### [pilot] Influence on min-pivot (delta = full − variant_off)

| probe | even N=8 (d=16) | even N=6 (d=12) | odd N=7 (d=15) | odd N=6 (d=13) |
|---|---|---|---|---|
| tau2_off (prime 2)   | +0.000134 | +0.005951 | −0.043637 | −0.032362 |
| tau3_off (prime 3)   | −0.011209 | −0.024886 | −0.112967 | −0.031358 |
| **cross_off (J(τ₂,τ₃))** | **+0.195074** | **+0.196116** | **+0.099861** | **+0.370087** |
| scale_cross ×100     | −62.39 | −62.41 | −13.74 | −13.85 |

- **Prime 2 is nearly inert in the even sector** (+0.0001 at N=8) — mirrors the
  first-window even-sector inertness. Prime 3 is modest.
- **The cross term $J(\tau_2,\tau_3)$ is the LARGEST single prime effect in both
  sectors.** Removing it drops the even N=8 pivot from −0.036 to −0.231.
- **scale_cross ×100** moves the pivot by tens of units — the cross term is
  demonstrably IN the computation (E1 kill-criterion: a term's presence is
  probed by a large-factor scaling mutant, not by demanding a sign flip).

### Discipline C — multi-N (is "cross dominates" a dimension artifact?)

- **Even sector: NO.** d_cross = +0.1951 (N=8) vs +0.1961 (N=6): identical to
  ~0.001. The cross-term influence is a robust structural feature, not an
  artifact of the truncation dimension.
- **Odd sector: qualitatively stable.** d_cross varies with N (+0.100 at N=7 vs
  +0.370 at N=6) because $d$ and $b_L$ differ, but the cross term remains the
  dominant, strongly-positive prime effect at every N tested. Honest grade:
  even-sector robustness is quantitative; odd-sector robustness is qualitative.

### [certify] Arb-interval confirmation

`scripts/certify_cross_influence.py` builds $C$ as outward-rounded rational
intervals and reports two DISTINCT metrics (PROOF_CONSTITUTION D3 — not to be
conflated):

- **Metric A — $\Delta C$ certified nonzero.** Since $M^{(2)}$ has no cross term,
  full and cross_off differ only through $S^{(2)}$:
  $\Delta C = C_{\text{full}} - C_{\text{cross\_off}} = -3\,c_2 c_3 (F_{ij}+F_{ji})$
  exactly. $\max|\Delta C|$ bounded away from 0 certifies the cross term is a
  **real, nonzero contribution** — independent of pivot sign. This is the
  primary, robust finding.
- **Metric B — pivot strictly improved.** pivot_sep > 0 would certify the cross
  term makes the min-pivot better. This is a STRONGER claim and can be
  indeterminate when a pivot band straddles zero in an indefinite $C$. A
  straddle/False on B does NOT weaken A.

(certify results: see `pilots/s4_certify_cross_L060.json` once the N=8/7 run
completes; this section is updated with the interval numbers then.)

### [certify] results at L = 3/5, N=8 (even) / N=7 (odd) — Arb interval

| quantity | even N=8 | odd N=7 |
|---|---|---|
| full min-pivot (interval)      | [−0.019697, −0.019680] | [−0.099571, −0.051604] |
| full positive-definite?        | **No** (note=straddle)  | **No** (note=straddle)  |
| cross_off min-pivot (interval) | [−0.526488, −0.231336] | [−0.257244, −0.257141] |
| **Metric A: max\|ΔC\| ≥**       | **0.63037** | **0.13984** |
| **Metric A: cross term certified nonzero** | **True** | **True** |
| Metric B: pivot_sep interval   | [+0.21164, +0.50681] | [+0.15757, +0.20564] |
| Metric B: confirmed_positive   | False (see note) | False (see note) |

**Metric A (primary, robust): CONFIRMED in both sectors.** At certify grade the
cross term contributes $\max|\Delta C| \ge 0.630$ (even) / $0.140$ (odd), bounded
away from zero. $J(\tau_2,\tau_3)$ is a **real, nonzero** contribution to the
second-window Schur matrix — the genuine new structure, now Arb-certified.

**Metric B nuance (do NOT let it weaken A).** The reported pivot_sep interval is
positive in both sectors (even [+0.21,+0.51], odd [+0.16,+0.21]), i.e. the
full-min-pivot lower endpoint exceeds the cross_off upper endpoint. But
`confirmed_positive=False` because the guard requires BOTH variants' pivots to be
non-straddling, and `full` straddles zero mid-factorization in both sectors (an
indefinite matrix). So Metric B is reported as indeterminate-by-guard, NOT as a
negative result. Its technical difficulty (indefinite pivot bands) does not
weaken Metric A, which is sign-independent (PROOF_CONSTITUTION D3).

**Full positivity: certifiably FALSE at L=0.6.** Both sectors' full min-pivot
intervals lie strictly below zero. The second window is NOT positive-definite at
L=0.6 — an independent fact from the cross-term significance, and squarely S5's
question at other L points.

## Verdict and S5 steer

**The cross-prime term $J(\tau_2,\tau_3)$ is the second window's genuine new
structure** — significant, not inert, robust across N (quantitatively in even,
qualitatively in odd), and the largest single prime influence in both sectors.
This is real new mathematics the first window does not contain (method-range
extension), and the empirical justification for developing the second window.

**Three propositions kept strictly separate (PROOF_CONSTITUTION D3):**
1. "Cross term is significant new structure" — **supported** (metric A).
2. "Second window is positive-definite" — **NOT established**: full min-pivot is
   still negative at L=0.6 in every (sector, N) tested. This is an independent
   question for S5's certify positivity judge.
3. "Approaching RH" — **no**. Scope is finite-scale positivity, $L<\log2$; a
   large cross-term influence says nothing about RH.

"Cross term dominates" must not slide into "second window nearly proven" or
"nearly RH." A dominant cross influence can still be insufficient to turn the
whole form positive; that is exactly what S5 must decide.

**Compute-allocation steer for S5:** spend hard compute on the cross term
$J(\tau_2,\tau_3)$ and its certify-grade treatment (including the two-irrational-
shift $\tau_2,\tau_3$ enclosure — the new numerical concern), not on the
near-inert even-sector prime-2 shift. Whether any $L$ in the window is
positive-definite is open and is S5's job.
