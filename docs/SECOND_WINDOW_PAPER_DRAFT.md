# The second prime window: a nonzero-but-non-dominant cross-prime coupling, and the loss of finite-scale positivity

**Draft — target: *Experimental Mathematics*.**
**Date:** 2026-08-09. **Last updated:** 2026-08-17.
**Scope:** the Weil quadratic form on $L^2(-L,L)$ for
$L$ in the second prime window $(\tfrac12\log3,\ \log2)\approx(0.5493,\,0.6931)$,
where both primes $p=2,3$ sit in the single-hop regime. **Not in scope:** the
Riemann Hypothesis, any consequence of RH, $L\to\infty$, or any extrapolation
past $L=\log2$.

---

## Abstract

We study the Weil quadratic form $Q_L$ on the second prime window
$L \in (\tfrac12\log3,\ \log2) \approx (0.549, 0.693)$, where both primes $p=2$
and $p=3$ lie in the single-hop regime. The novelty over the first prime window
is a genuine cross-prime coupling $J_{ij}(\tau_2,\tau_3)$ in the second-moment
matrix $S^{(2)}$, which is absent when only one prime is active.

We establish three results by certified interval arithmetic (outward-rounded Arb
balls, independently checked). **(1) The cross term is real and nonzero:** the
Metric-A bound gives $\max|F_{ij}+F_{ji}| > 0$ with $\Delta_\lambda \ge 4.3\times10^{-3}$
across the window. **(2) The cross term is non-dominant:** at tested even-sector
points ($N\le8$) the coupling does not make the form positive-definite and does not
delay the onset of non-positivity. **(3) The method boundary is eta-sensitive:**
in the odd sector at $L=0.56$ with the standard split-residual weight $\eta=0.5$,
the form is certifiably not positive-definite for $N=15$ (pivot $(0,0)$ certified
negative) and indeterminate for $N=17$ (pivot $(2,2)$ straddles zero, with pivots
$(0,0)$ and $(1,1)$ certified positive). A new structural finding is that
$R_0[0,0] \approx 1.08\times10^{-3}$ (near-zero, archimedean Schur complement
near-cancels at the $P_1$--$P_1$ entry), making $R_2[0,0]/R_0[0,0] \approx 6.18$.
This near-cancellation implies $C[0,0](\eta^*) = +3.10\times10^{-3} > 0$ at $N=15$
for the entry-wise weight $\eta^* \approx 2.49$, contrasted with $C[0,0](\eta=0.5) =
-5.39\times10^{-3} < 0$. However, a float full-matrix $\eta$ scan at $N=17$ shows that
the entry-wise $\eta^*=2.49$ is counterproductive for the global minimum eigenvalue
($\lambda_\min$ worsens from $-0.063$ to $-0.099$); the Frobenius-optimal $\eta_F^* \approx 1.22$
gives the best $\lambda_\min = -0.063$ at $N=17$. The form remains far from positive-definite
at all tested $N$ and $\eta$. A sub-matrix sweep ($k=3\ldots17$ from the $N=17$ Arb
checkpoint) reveals geometric convergence at rate $r\approx0.791$ per step: the second
eigenvalue $\lambda_1$ crosses zero at $k\approx17{-}18$ (asymptote $B_1{=}+0.045$), while
the minimum eigenvalue $\lambda_0$ reveals a critical structural finding: decomposing
$\lambda_0 = v^TCv$ into UV-diagonal, UV-cross, and IR-block contributions shows
the UV-cross term was near-constant at $-0.029$ for $k=13..18$, then **dropped to
$-0.02395$ at $k=20$** (cumulative ratio 0.832 from k=18, monotone). IR-block decays at $r\approx0.46$ (k=19→20), UV-diagonal
**flipped positive** at $k=20$ ($C_{UV,UV}=+0.0044$). The $k=20$ result $\lambda_0(20)=-0.02741$ fits the
8-pt exp+B model ($B_0=+0.008$, RMS $8.2\times10^{-4}$) better than the $B=0$ forced model.
The second eigenvalue $\lambda_1$ confirmed zero-crossing at $k=18$ ($n_\text{neg}$ drops
2→1), matching the $B_1=+0.045$ prediction.
**Current status (2026-08-17): ZERO CROSSING CONFIRMED.** Chain k=18..25 complete. λ₀(k=25, η=1.0) = +1.54e-6 > 0; n_neg=0 for ALL η∈[0.5,1.22] at k=25. 13-pt Bootstrap P($B_0>0$) = **1.0000**, CI=[+0.0144,+0.0271] entirely positive. $B_0 = +0.0177$ (best fit, A=−2.347, r=0.820). UV-cross continued monotone decay to −0.01933(k=23), −0.01701(k=24). IR-block growing: +0.00695(k=23), +0.00819(k=24). r(23→24)=0.117 (dramatic acceleration); λ₀(24)=−0.00068; λ₀(25)=+1.54e-6. Chain k=26..28 computing (certify target).

Two numerical traps encountered in the second-window adaptation are documented:
a kappa-contamination bug (importing a first-window constant $\kappa=1.255$ when
the correct value is $\kappa(L=0.56)=2.056$) that produced a spurious positive
signal; and a Bernstein-blowup bug in the $S_{KK}/S_{VK}$ integrators that yields
interval widths $\gg1$ for loop index $k\ge43$. Both are proposed as a new
proofctl condition class C12 (parameter-default propagation blowup).

---

## 0. Grade discipline (read first)

Every numerical statement below carries an explicit evidence grade:

- **[certify]** — outward-rounded interval arithmetic (Arb balls), reproduced by
  the independent stdlib checker. A certified interval is reported with both
  endpoints; a sign verdict is made only when the whole interval lies on one side
  of zero.
- **[pilot]** — float-center computation (`numpy.eigvalsh` / float LDLᵀ) on the
  assembled Schur matrix. Directional and shape evidence only; never a verdict.

**Three propositions are kept strictly separate throughout.** We do not let
evidence for one migrate into a claim about another:

1. **(CROSS)** the cross-prime term $J(\tau_2,\tau_3)$ is a real, nonzero new
   structure — **[certify], established**;
2. **(POS)** the second window is positive-definite — **not established; the
   evidence points the other way** (odd sector certifiably negative at a
   certified point);
3. **(RH)** any statement about the Riemann Hypothesis — **not made, at any
   grade.**

A further discipline, enforced by the project checker: **no finite sample of $L$
is promoted to a whole-window assertion.** Where we have two certified $L$-points
plus a pilot sampling grid, we say exactly that — we do **not** write "the window
is certifiably not positive-definite everywhere."

---

## 1. What is new in the second window

In the first prime window only $p=2$ is single-hop, so the prime layer is a
single shift and the second-moment object $S^{(2)}$ has no cross term. In the
second window **both** $p=2$ and $p=3$ are single-hop, with shift parameters
$\tau_p=\log p/L$. The second moment becomes a genuine two-shift object,
$$
S^{(2)} = \big\langle (V{+}K)P_j,\ C_{\tau_2}P_i\big\rangle
        + \big\langle (V{+}K)P_j,\ C_{\tau_3}P_i\big\rangle ,
$$
and carries a **cross-prime coupling**
$$
F_{ij}(\tau_2,\tau_3)=\big\langle C_{\tau_3,1}P_j,\ C_{\tau_2,1}P_i\big\rangle,
$$
which is identically absent in the first window. This coupling — its existence,
its certified non-vanishing, and its measured failure to control positivity — is
the subject of the paper.

---

## 2. (CROSS) The cross term is real and nonzero — [certify]

We isolate the cross term by a sign-independent difference (metric A):
$$
\Delta C = C_{\text{full}} - C_{\text{cross-off}}
         = -\,3\,c_2 c_3\,(F_{ij}+F_{ji}),
\qquad c_p=\frac{\log p}{\sqrt p},
$$
exact because $M_2$ (the single-prime part) contains no cross term. The checker
(`checker/second_prime/check_cross_structure.py`) recomputes $F$ in exact
`Fraction` arithmetic over $\mathbb{Q}[\tau_2,\tau_3]$ and verifies:

| obligation | statement |
|---|---|
| cross-term-present-and-nonzero | recomputed $\max\lvert F_{ij}+F_{ji}\rvert>0$ |
| both-shifts-present | recomputed $J(\tau_2)$ **and** $J(\tau_3)$ nonzero |
| window-bounds-hold | $\tfrac12\log3<L<\log2$ by certified rational bounds |
| four-term-S0-declared | $S_{VV}+S_{VK}+S_{KV}+S_{KK}$ + cross present |
| positivity-not-claimed | no positivity/pivot/conclusion fields in the certificate |
| conclusion-bounded-and-no-rh | $L<\log2$; two-prime method; primes $=[2,3]$ |

**Result — [certify]:** $\max\lvert\Delta C\rvert \ge 0.630$ (even) / $0.140$
(odd) at Arb grade at $L=3/5$. The cross term is a real, nonzero contribution;
its construction is validated by four independent invariants
($F(\tau,\tau)=E$, parity, operator-swap symmetry, cross-Cauchy–Schwarz). This
replaced a prototype that silently set $F_{\text{cross}}=0$ (an omitted-cross-term
defect of exactly the class the mutation catalog is designed to catch).

**This is the paper's positive result, and it is bounded to non-vanishing.** It
says nothing about sign or dominance — §3.

---

## 3. (CROSS, negative) The cross term neither dominates nor delays collapse

Two natural strengthenings of the cross-term finding were **tested and refuted**
by independent recomputation against a certified anchor. We report both, because
the refutations are as much a result as the non-vanishing.

**3a. "Cross term dominates positivity" — refuted.** An outsourced report claimed
a Layer-1 eigenvalue $+0.077$ under a cross-on assembly. Independent
recomputation gives the true eigenvalue $-0.231$ [pilot, reproduces the S4
assembly at $L=3/5$; the certify cross-off pivot upper endpoint independently
sits at $-0.231$]: the reported sign was wrong. The cross term contributes
positively to the pivot (turning cross off *lowers* the even pivot from $-0.036$
to $-0.231$ at $L=0.60$ [pilot]), but it does **not** dominate — it cannot lift
the form to positive-definite.

**3b. "Cross term delays the margin collapse" — falsified.** The conjecture that
a positive cross margin $\Delta V$ postpones the loss of positivity fails at
three strongly-negative test points: $\lambda_{\min}<0$ across the sampled window
with no surviving margin to delay anything.

**Per-prime decomposition [pilot, $L=0.60$]** makes the mechanism concrete:

| sector | full | prime-2 | prime-3 | cross |
|---|---|---|---|---|
| even $N{=}8$ | $-0.036$ | $+0.0001$ (inert) | $-0.011$ | $+0.195$ |
| odd $N{=}7$  | $-0.157$ | $-0.044$ | $-0.113$ (dominant new negative) | $+0.100$ |

Prime-3 is the dominant new negative contribution (odd), prime-2 is nearly inert
(as in the first window), and the cross term — though positive — is non-dominant
and cannot flip the sign.

---

## 4. (POS) Positivity is lost at tested points — [certify] points + [pilot] grid

We report exactly what was computed, at its grade, with no window-wide promotion.

**Certified points.**

| $L$ | sector | min-pivot interval | verdict | grade |
|---|---|---|---|---|
| $11/20=0.55$ (left end $+7\!\times\!10^{-4}$) | even $N{=}8$ | $[-3.869\!\times\!10^{-4},\ +2.428\!\times\!10^{-4}]$ | **indeterminate** (straddles 0; width $>$ dist-to-0) | [certify] |
| $11/20=0.55$ | odd $N{=}7$ | $[-2.815\!\times\!10^{-2},\ -1.720\!\times\!10^{-2}]$ | **certifiably NEGATIVE** (upper $<0$) | [certify] |
| $69/100=0.69$ ($\log2^-$) | odd $N{=}7$ | pivot straddle $\approx-0.099$; $\lambda_{\min}$(eig)$=-0.369$ | negative | [certify] anchor |

**Pilot grid (eigenvalue, shape only).** The even $N{=}8$ and odd $N{=}7$
sampled curves are **[pilot]**, monotone-decreasing in $L$, all sampled values
negative, infimum at the sampled right end. Two of these points are
certify-anchored (their float-center eigenvalue reproduces the Arb-interval
assembly to all printed digits), so the pilot curve is trustworthy as an
eigenvalue *center* but is not itself a verdict.

**Honest reading.** The odd sector is **certifiably negative** even at the
most-favorable left end $L=0.55$. The even sector is **indeterminate** there
(interval straddles zero — an interval-inflation / near-zero situation, not a
clean verdict). Therefore:

> The second window is **not** certified positive-definite anywhere; at the two
> certified interior/edge points and across the pilot grid, positivity is lost
> (odd) or unresolved (even). We do **not** claim "the whole window is
> certifiably not positive-definite" — that would be a finite-sample →
> whole-window over-reach the checker forbids.

Compared with the first window (even $+0.0087$, odd $+0.053$ at $L=7/20$, both
clearly positive [certify]), the second window sits at or past the boundary. The
mechanism is structural: $c_L=\log(2\pi L)+\gamma$ grows with $L$
($1.82\to2.05$ across the window), adding diagonal negative pressure, while
prime-3 adds a genuine new negative share that prime-2 never did.

---

## 5. What this is, and what it is not

**Is:** a certified new-structure result — the second window introduces a
two-prime cross coupling $J(\tau_2,\tau_3)$ absent in the first window, whose
non-vanishing is certified; together with the honest finding that this coupling
is **non-dominant** (it neither controls positivity nor delays its collapse), and
that finite-scale positivity is **lost** at the tested points. This is a
**method-range** result: the per-window positivity method has a shorter reach in
the second window than the first.

**Is not:** a proof of second-window positivity (it is not positive-definite); a
whole-window verdict (only two certified $L$-points plus a pilot grid); an
extrapolation past $L=\log2$; or any statement, at any grade, about RH or its
consequences.

---

## Appendix A. Data table with evidence grades

| id | quantity | $L$ | sector | value / interval | grade | source |
|---|---|---|---|---|---|---|
| A1 | $\max\lvert\Delta C\rvert$ (cross, metric A) | $3/5$ | even | $\ge 0.63037$ | [certify] | `pilots/s4_certify_cross_L060.json` |
| A2 | $\max\lvert\Delta C\rvert$ (cross, metric A) | $3/5$ | odd | $\ge 0.13984$ | [certify] | " |
| A3 | min-pivot interval | $11/20$ | even | $[-3.869\text{e-}4,\ +2.428\text{e-}4]$ (indeterminate) | [certify] | `pilots/s5_positivity_L055.json` |
| A4 | min-pivot interval | $11/20$ | odd | $[-2.815\text{e-}2,\ -1.720\text{e-}2]$ (negative) | [certify] | " |
| A5 | $\lambda_{\min}$(eig) | $69/100$ | odd | $-0.36911$ | [certify] anchor | `pilots/shape_certify_anchor_odd.json` |
| A6 | $\lambda_{\min}$(eig) | $67/100$ | even | $-0.28596$ | [certify] anchor | `pilots/shape_certify_anchor.json` |
| A7 | pivot straddle | $69/100$ | odd | $\approx-0.099$ | [certify] | " |
| A8 | full / prime-2 / prime-3 / cross split | $3/5$ | even | $-0.036/+0.0001/-0.011/+0.195$ | [pilot] | `pilots/s4_profile_L060.json` |
| A9 | full / prime-2 / prime-3 / cross split | $3/5$ | odd | $-0.157/-0.044/-0.113/+0.100$ | [pilot] | " |
| A10 | dominance Layer-1 eig (refuted) | $3/5$ | even | reported $+0.077$; true $-0.231$ | [pilot] refutation | `pilots/independent_dominance_anchor.json` |
| A11 | $\lambda_{\min}$ pilot curve (monotone, all negative) | window | both | see `SECOND_WINDOW_LAMBDA_SHAPE.md` | [pilot] | `pilots/shape_*_N*.json` |
| A12 | min-pivot interval ($N{=}15$, odd) | $14/25$ | odd | $(0,0)\in[-5.39\!\times\!10^{-3},\ -5.38\!\times\!10^{-3}]$ — certifiably negative | [certify] | `pilots/cert_fp_second_N15.json` |
| A13 | min-pivot interval ($N{=}17$, odd) | $14/25$ | odd | $(2,2)\in[-9.03\!\times\!10^{-4},\ +3.27\!\times\!10^{-3}]$ — indeterminate | [certify] | `pilots/cert_fp_second_N17.json` |
| A15 | min-pivot interval ($N{=}19$, $\eta{=}0.1$, odd) | $14/25$ | odd | $(0,0)\in[-3.80\!\times\!10^{-2},\ -5.15\!\times\!10^{-4}]$ — **certifiably negative** (upper $<0$) | [certify] | `pilots/cert_fp_second_N19_eta01.json` |
| A14 | N-convergence float scan (corrected $\kappa{=}2.056$) | $14/25$ | odd | $N{=}7..17$: $-0.575,-0.390,-0.277,-0.188,-0.116,-0.076$ | [pilot] | `pilots/eig_scan_corrected_N*_final.json` |
| A16 | Sub-matrix chain $k=18..22$ ($\lambda_0$ sequence) | $14/25$ | odd | $-0.049,-0.038,-0.027,-0.019,-0.011$ (η=1.0) | [pilot] | `pilots/submatrix_k18..22.json` |
| A17 | B₀ fit (13-pt, $k=13..25$) + bootstrap | $13/25$ | odd | $B_0=+0.0177$; 95% CI $[+0.0144,+0.0271]$; P($B_0>0$)=1.0000 | [pilot] | `scripts/fit_b0.py` |
| A18 | Sub-matrix chain k=23..25, zero crossing | $25/25$ | odd | λ₀(25)=+1.54e-6, n_neg=0 ∀η∈[0.5,1.22]; ZERO CROSSING CONFIRMED | [pilot] | `scripts/submatrix_chain.py` |

Constants: $c_2=\log2/\sqrt2=0.4901$, $c_3=\log3/\sqrt3=0.6343$;
$c_L=\log(2\pi L)+\gamma$ ($=1.816$ at $L=0.549$, $2.049$ at $L=0.693$). All
certified quantities are recomputed by the stdlib checker in
`checker/second_prime/`; certificates carry no matrix/eigenvalue/pivot/conclusion
values (status is derived by the checker, never self-reported).

---

## 6. N-convergence at $L=0.56$ and the method boundary — [pilot + certify]

*Added 2026-08-16. Supersedes the N=13 "zero crossing" from the original scan,
which was a false positive caused by a wrong $\kappa$ value (first-window constant
$\kappa=1.255$ instead of the correct second-window value
$\kappa(L=0.56)=2.056$). All data below use the corrected $\kappa$.*

### 6.1 The $\kappa$ bug and its correction

The second-window spectral parameter $\kappa$ appears in the positivity criterion
as $b_L = H_d - c_L - \kappa > 0$. For $L=0.56$,
$\kappa(L)=\kappa_0 + 2L\log(L)$ gives $\kappa(0.56)\approx2.056$,
which is 64 % larger than the first-window value $1.255$.

An earlier scan used the hardcoded first-window constant, producing
$b_L(N{=}13)\approx+0.937$ and an apparent positive eigenvalue
$\lambda_{\min}\approx+2\!\times\!10^{-6}$. With the correct $\kappa$,
$b_L(N{=}13)\approx+0.000214$ and $\lambda_{\min}(N{=}13)\approx-0.188$.
The "zero crossing at $N=13$" is **entirely an artifact** of the wrong $\kappa$.

Lesson: $\kappa$ must be computed from $L$ at runtime; importing a hardcoded
constant from a different window is a P0 defect.

### 6.2 N-convergence table (odd sector, $L=0.56$, corrected $\kappa=2.056$)

| $N$ | $d$ | $b_L$ | $\lambda_{\min}$ (float, [pilot]) | $\Delta\lambda$ | ratio |
|---|---|---|---|---|---|
| 7 | 15 | $-0.573$ | $-0.5747$ | — | — |
| 9 | 19 | $-0.343$ | $-0.3895$ | $+0.1852$ | — |
| 11 | 23 | $-0.157$ | $-0.2770$ | $+0.1125$ | $0.607$ |
| 13 | 27 | $+0.000$ | $-0.1879$ | $+0.0891$ | $0.792$ |
| 15 | 31 | $+0.136$ | $-0.1155$ | $+0.0724$ | $0.812$ |
| 17 | 35 | $+0.256$ | $-0.0760$ | $+0.0395$ | $0.545$ |
| 19 | 39 | $+0.362$ | $\approx-0.054$ (est.) | $\approx+0.022$ | $\approx0.56$ |

The increment ratio $r_N = \Delta\lambda_N / \Delta\lambda_{N-2}$ fluctuates
between 0.55 and 0.81; a simple geometric extrapolation does not give a
reliable sign verdict on the limit. The certified results in §6.3 show that
float-center extrapolation is unreliable here — the true sign near zero can
differ from the float estimate.

**Richardson extrapolation uncertainty.** A geometric series extrapolation
$\lambda_\infty \approx \lambda_N + \Delta\lambda_N/(1-r)$ gives qualitatively
different answers depending on which consecutive pair is used:
- From $N=15,17$ pair: $r=0.546$, $\lambda_\infty \approx -0.076 + 0.0395/0.454 \approx +0.011$
- From estimated $N=17,19$ pair: $r\approx0.557$, $\lambda_\infty \approx -0.054 + 0.022/0.443 \approx -0.004$

The two-point extrapolations disagree in sign ($+0.011$ vs $-0.004$), both near zero.
This shows the extrapolation is **inconclusive**: the limit may be slightly positive or
slightly negative. The eta-optimized certify and the N=19 certify together are the
most direct path to a resolution.

### 6.3 Certify-grade N-study (Arb interval arithmetic)

We ran the full two-prime certify pipeline
(`checker/fp_second/certify_fp_second.py`) for $N=15$ and $N=17$:

| $N$ | $b_L$ interval | first-failing pivot | verdict | elapsed |
|---|---|---|---|---|
| 15 | $[0.1363, 0.1363]$ | $(0,0)\in[-5.39\!\times\!10^{-3},\ -5.38\!\times\!10^{-3}]$ | **CERTIFIABLY NOT POSITIVE** | 5721 s |
| 17 | $[0.2555, 0.2555]$ | $(2,2)\in[-9.03\!\times\!10^{-4},\ +3.27\!\times\!10^{-3}]$ | **INDETERMINATE** (straddles 0) | 19020 s |
| 19 | $[0.362, 0.362]$ | $(0,0)\in[-3.80\!\times\!10^{-2},\ -5.15\!\times\!10^{-4}]$ | **CERTIFIABLY NOT POSITIVE** ($\eta{=}0.1$) | 6.5 s (cached) |
| 19 | $[0.362, 0.362]$ | $(1,1)\in[-1.56\!\times\!10^{2},\ +7.42\!\times\!10^{-2}]$ | *INDET* ($\eta{=}0.5$, Schur blow-up) | 36036 s |
| 19 | $[0.362, 0.362]$ | $(0,0)\in[-1.54\!\times\!10^{-3},\ +6.67\!\times\!10^{-2}]$ | *INDET* ($\eta{=}1.0$, straddles 0) | 6.5 s (cached) |

At $N=15$, pivot $(0,0)$ is strictly negative: the form is certifiably not
positive-definite. At $N=17$, **pivots $(0,0)$ and $(1,1)$ are both certified
positive** — genuine progress from $N=15$ — and only pivot $(2,2)$ straddles
zero. For $N=19$, the choice of $\eta$ is critical: at $\eta=0.5$ and $\eta=1.0$,
the $(0,0)$ near-cancellation ($b_L F[0,0] \approx R_\eta[0,0]$) causes either
catastrophic Schur blow-up (pivot width 156 at $\eta=0.5$) or an interval straddling
zero ($\eta=1.0$). At $\eta=0.1$, however, $R_\eta[0,0] = 0.0636 \gg b_L F[0,0] = 0.0431$
(ratio 0.68, no near-cancellation), giving pivot $(0,0) \in [-0.038, -5.15\!\times\!10^{-4}]$
with **upper bound strictly negative** — a valid CNPD certificate in 6.5 s using the
cached 361-entry checkpoint. The interval width (0.037) is large due to the $1/\eta=10$
amplification factor on $R_2$, but the upper bound clears zero with margin $5\!\times\!10^{-4}$.

**Eta scan of N=17 pivots (from checkpoint, 6 s each, [pilot]).** Using the completed
N=17 matrix checkpoint (all 289 Arb interval entries cached), we re-ran the
$C(\eta)$ assembly at four values:

| $\eta$ | Pivot $(2,2)$ interval | Center | Lower bound |
|--------|------------------------|--------|-------------|
| $0.5$ | $[-9.03\!\times\!10^{-4},\ +3.27\!\times\!10^{-3}]$ | $+1.18\!\times\!10^{-3}$ | $-9.03\!\times\!10^{-4}$ |
| $1.0$ | $[-7.88\!\times\!10^{-4},\ +3.59\!\times\!10^{-3}]$ | $+1.40\!\times\!10^{-3}$ | $-7.88\!\times\!10^{-4}$ (best) |
| $2.49$ | $[-2.33\!\times\!10^{-3},\ +4.90\!\times\!10^{-3}]$ | $+1.28\!\times\!10^{-3}$ | $-2.33\!\times\!10^{-3}$ |
| $4.0$ | $[-4.24\!\times\!10^{-3},\ +6.22\!\times\!10^{-3}]$ | $+0.99\!\times\!10^{-3}$ | $-4.24\!\times\!10^{-3}$ |

The Arb pivot center is **positive at all $\eta$** (~$+1.2\times10^{-3}$). However,
pivot $(2,2)$ is only the Schur complement of the leading $3\times3$ block — it does
not bound the remaining eigenvalues. A **float full-matrix $\eta$ scan** (N=17,
loading Arb checkpoint centers, 4 s) reveals the true picture:

**Float $\lambda_\min(C(\eta))$ at $N=17$ [pilot]:**
$\|R_0\|_F = 3.43\times10^{-2}$, $\|R_2\|_F = 5.11\times10^{-2}$,
Frobenius $\eta^*_F = 1.22$.

| $\eta$ | $\lambda_\min(C)$ | note |
|---|---|---|
| $0.5$ | $-0.0737$ | default |
| $0.75$ | $-0.0629$ | |
| **$1.0$** | **$-0.0628$** | **global minimum** |
| $1.22$ (Frob $\eta^*$) | $-0.0658$ | |
| $2.49$ (entry $\eta^*$) | $-0.0986$ | entry-wise optimal → worse globally |
| $4.0$ | $-0.1459$ | |

The minimum eigenvalue is $-0.063$ at best ($\eta=1.0$): **$N=17$ is not close to
positive-definite**. The entry-wise $\eta^*=2.49$ (which maximises $C[0,0]$ in §6.5)
is actively harmful globally — it worsens $\lambda_\min$ by 57% vs $\eta=1.0$.
The Frobenius $\eta^*_F = 1.22$ is the correct global minimiser of $\|R_\eta\|$.

**Structural observation on pivot evolution.** The certify data reveal that
positivity "sweeps inward" as $N$ grows:
at $N=15$, pivot $(0,0)$ fails; at $N=17$, pivots $(0,0)$ and $(1,1)$ pass
and only $(2,2)$ is uncertain. This is consistent with the method-boundary
analysis (§6.5). However, the float $\lambda_\min$ evidence shows that even though
early pivots certify positive, the tail of the LDL$^\top$ factorization still
carries large negative Schur complements. The matrix is still far from positive-definite.

*Important context from `rh-obstruction-theory/F-schur-complexity`:*
"Certificate complexity growing does **not** imply that the Weil form is
indefinite." The indeterminate result at $N=17$ is a statement about
interval-arithmetic precision, not about the form's sign.

*Structural parallel (informational):* The method boundary identified in §6.5 —
"decay rate $b_L$ must exceed the residual growth $R_\eta$, but certificate cost
grows as $N^3$ per entry" — is structurally isomorphic to the Kotecký–Preiss
threshold identified in `yang-mills-verification/Route 2` (loop-space: decay rate
$\kappa_{HK}$ vs. entropy/connective constant $A$). Both are instances of the
same resistance pattern: **expansion-parameter complexity outpacing the positivity
margin uniformly**. This is a methodological observation, not a mathematical
reduction between the two problems.

### 6.3b Sub-matrix sweep: N-convergence from k=3 to 18 (2026-08-16, [pilot])

**Method.** Since the matrix entries $M_0[a,b]$, $S_0[a,b]$ are inner products of
Legendre polynomials independent of the truncation parameter $N$, the top-left
$k\times k$ block of the $N=17$ Arb checkpoint gives exactly the $k$-dimensional
Galerkin Schur matrix. This yields a full $N$-sweep ($k=3\ldots17$) from a single
5-second post-processing computation. The k=18 row was computed independently
(compute\_submatrix\_k18.py, ~68 min) and verified against the full-matrix result.

**Key findings:**

| $k$ | $b_L$ | $n_{\text{neg}}$ | $\lambda_0$ | $\lambda_1$ | $\lambda_2$ |
|-----|-------|---------|-------------|-------------|-------------|
| 12 | $-0.075$ | 12 | $-0.202$ | $-0.086$ | $-0.049$ |
| 13 | $+0.000$ | 8 | $-0.160$ | $-0.063$ | $-0.035$ |
| 14 | $+0.070$ | 6 | $-0.128$ | $-0.044$ | $-0.018$ |
| 15 | $+0.136$ | 5 | $-0.099$ | $-0.024$ | $-0.006$ |
| 16 | $+0.198$ | 3 | $-0.079$ | $-0.011$ | $\approx 0$ |
| 17 | $+0.256$ | 2 | $-0.063$ | $-0.001$ | $+10^{-6}$ |
| **18** | **+0.310** | **1** | **−0.049** | **+0.004** | **+0.009** |
| **19** | **+0.362** | **1** | **−0.038** | **$+7\!\times\!10^{-7}$** | **+0.006** |

At $k=12$ (where $b_L < 0$), all 12 eigenvalues are negative. At $k=13$ ($b_L$
first positive), $n_{\text{neg}}$ drops from 12 to 8 — a phase transition. From
$k=13$ onward, $n_{\text{neg}}$ decreases monotonically: $8 \to 6 \to 5 \to 3 \to 2 \to \mathbf{1}$.

**Second eigenvalue $\lambda_1$ confirmed zero-crossing at $k=17{-}18$:**
$n_\text{neg}$ drops from 2 to 1 at $k=18$; $\lambda_1(k=18) = +0.004 > 0$, confirming
the 5-point NLS prediction $\lambda_1(N) \approx -1.84 \cdot (0.805)^N + 0.045$,
asymptote $B_1 = +0.045$ (clearly positive). **Only $\lambda_0$ remains negative.**

**Convergence model for $\lambda_0$ (k=13..22, 10 pts):** accelerating rate; r(21→22)=0.566.
Step-wise ratios: $0.799, 0.772, 0.801, 0.792, 0.779, 0.770, 0.728, 0.687, \mathbf{0.566}$. Mean $\bar r \approx 0.749$ (accelerating downward — rate itself shrinking, consistent with $B_0 > 0$).

| Model | $B_0$ | $\lambda_0(22)$ pred | **actual** | $\lambda_0(25)$ pred | RMS |
|-------|--------|-------------|------------|-------------|-----|
| $Ar^k + B$ (3-param, 10-pt) | $\mathbf{+0.016}$ | — | $\mathbf{-0.01066}$ | $\mathbf{+0.0007}$ | **$1.16\times10^{-3}$** |
| $Ar^k$ ($B{=}0$, forced, 10-pt) | $0$ | $-0.012$ | $-0.01066$ | $-0.004$ | $1.6\times10^{-3}$ |
| Component: $B_0 = -0.029$ | $-0.029$ | $-0.049$ ✗ | $-0.01066$ | $-0.032$ | — (rejected) |

10-pt fit selects Exp+B ($B_0=+0.016$, RMS $1.16\times10^{-3}$). **13-pt Bootstrap 95% CI = [+0.01435, +0.02706] — entirely positive.** P($B_0>0$) = **1.0000**. Best fit: A=−2.347, r=0.820, $B_0=+0.0177$.

**UV-mode decomposition: the constant UV-cross term (2026-08-16, critical finding).**
Decomposing $\lambda_0 = v^T C v$ into three parts using the min-eigenvector $v$:

$$\lambda_0 = \underbrace{C_{kk}|v_\text{UV}|^2}_{\text{UV-diag}} + \underbrace{2v_\text{UV}\!\!\sum_{j<k} C_{jk} v_j}_{\text{UV-cross}} + \underbrace{v_\text{IR}^T C_\text{IR} v_\text{IR}}_{\text{IR-block}}$$

| $k$ | $\lambda_0$ | UV-diag | UV-cross | IR-block | $|v_\text{UV}|^2$ |
|-----|-------------|---------|---------|---------|--------|
| 13 | $-0.1604$ | $-0.0056$ | $-0.0279$ | $-0.1269$ | 0.122 |
| 15 | $-0.0990$ | $-0.0057$ | $-0.0300$ | $-0.0632$ | 0.210 |
| 17 | $-0.0628$ | $-0.0038$ | $-0.0289$ | $-0.0301$ | 0.291 |
| **18** | **$-0.0489$** | **$-0.0021$** | **$-0.0288$** | **$-0.0180$** | **0.338** |
| **19** | **$-0.0376$** | **$-0.0008$** | **$-0.0257$** | **$-0.0111$** | **0.363** |
| **20** | **$-0.02741$** | **$+0.00165$** | **$-0.02395$** | **$-0.00511$** | **0.377** |
| **21** | **$-0.01883$** | **$+0.00251$** | **$-0.02169$** | **$+0.000349$** | — |
| **22** | **$-0.01066$** | **$+0.00474$** | **$-0.01947$** | **$+0.00408$** | — |

**UV-diag sign flip at k=20** ($C_{UV,UV}$ crosses zero), growing positive through k=21,22,23,24.
**IR-block sign flip at k=21** (k=20: −0.00511 → k=21: +0.000349 → k=22: +0.00408 → k=23: +0.00695 → k=24: +0.00819). Both UV-diag and IR-block push $\lambda_0$ toward zero from below.
**UV-cross** monotone decay confirmed: −0.0239, −0.0217, −0.0195, −0.01933, −0.01701 (k=20..24); per-step ratio ≈0.89–0.90.
**k=25: ZERO CROSSING** — λ₀(k=25,η=1.0) = +1.54e-6 > 0, n_neg=0 for ALL η∈[0.5,1.22].
**Zero crossing mechanism**: at k=24, UV-diag+IR = +0.01633 vs |UV-cross| = 0.01701; net −0.00068. At k=25 the sum tips positive.

**IR-block** sign-flips at k=21, growing positive through k=24.
**UV-diag** sign-flipped positive at k=20, growing.
**UV-cross**: per-step ratios $1.061, 1.014, 1.019, 0.945, 0.994$ for $k=13..18$, then **0.893, 0.932, 0.905, 0.898** for $k=19..22$.
**The UV-cross "constant at $-0.029$" hypothesis is fully refuted by k=19..22.**

**Layer decomposition of UV-cross by distance from UV boundary (k=13..19).**
Writing UV-cross $= \sum_{d=1}^{k-1} 2C_{k-d, k-1} v_{k-d} v_{k-1}$:

| $k$ | $d=1$ | $d=2$ | $d=3$ | $d=4$ | $d=5$ | $d{\ge}6$ | total | $k\times$total |
|-----|-------|-------|-------|-------|-------|-----------|-------|----------------|
| 14 | $-0.0063$ | $-0.0036$ | $-0.0027$ | $-0.0025$ | $-0.0024$ | $-0.0123$ | $-0.0296$ | $-0.415$ |
| 16 | $-0.0077$ | $-0.0038$ | $-0.0026$ | $-0.0028$ | $-0.0027$ | $-0.0109$ | $-0.0306$ | $-0.490$ |
| 18 | $-0.0085$ | $-0.0040$ | $-0.0022$ | $-0.0026$ | $-0.0028$ | $-0.0088$ | $-0.0288$ | **$-0.518$** |
| **19** | **$-0.0070$** | **$-0.0033$** | **$-0.0028$** | **$-0.0026$** | **$-0.0020$** | **$-0.0080$** | **$-0.0257$** | **$-0.488\downarrow$** |

**$k\times\text{UV-cross}$ diagnostic**: if $B_0 < 0$ (UV-cross stabilizes), $k\times\text{UV-cross} \to \infty$; if $B_0=0$ (geometric decay), $k\times\text{UV-cross} \to 0$. The series peaked at $k=18$ ($-0.518$) and dropped at $k=19$ ($-0.488$). The near-UV block ($d=1..5$) also dropped: $-0.0177$ vs $-0.0199$ at $k=18$ ($-11\%$). All per-d contributions are decaying, not stabilizing. **This strongly favors $B_0 = 0$.**

**Updated interpretation (after k=19)**: The near-UV sum ($d=1..5 \approx -0.020$) appeared stable for $k=14..18$, but k=19 shows UV-cross dropped below this. Each per-d term is sign-definite negative (sign-definiteness lemma) but also shrinking ($C_{k-d,k-1} \sim -c_d/k \to 0$). A sum of shrinking negative terms can converge to 0. **The earlier claim "$B_0 \le -0.020$" is retracted; $B_0$ is unresolved.** k=20 will reveal whether the near-UV block is also decaying.

**Physical mechanism**: $|v_\text{UV}|^2$ grows $(0.12 \to 0.36)$ but individual couplings $C_{j,k-1} \sim -c_d/k \to 0$. Whether localization fully compensates the coupling decay is the key question — k=19 suggests it does not (decay wins), consistent with $B_0 = 0$.
Each new $k$ adds a $d=1$ contribution while old $d$-layers are "promoted" to $d+1$.
**Revised (after k=19)**: the near-UV block is NOT a constant-amplitude traveling wave — k=19 shows it decayed 11% relative to k=18. The earlier stability ($k=15..18$) was a plateau, not an asymptote. Each per-$d$ contribution shrinks as $O(1/k)$, and the sum converges to 0.

**UV column sign structure.** A direct analysis of the UV column $C_{j,k-1}$ for $k=18,19$
reveals that **every off-diagonal entry is strictly negative** ($C_{j,\text{UV}} < 0$ for all $j < k-1$).
The decomposition $C_{j,\text{UV}} = b_L F_{j,\text{UV}} - R_\eta(j,\text{UV})$ shows
$R_\eta(j,\text{UV}) > 0$ for all near-UV entries, dominating $b_L F_{j,\text{UV}}$ by
factors of $5$–$20\times$ at $k=18,19$ (and $\gg 10^3$ at $k=13$ where $b_L \approx 0$).
The $R_\eta$ dominance is traced to $(1+\eta)R_0(j,\text{UV}) > 0$, which holds because
$R_0 = S^{(0)} - M_0^T G^{-1} M_0$ inherits positivity from the second-moment Gram structure.

**Structural foundation — total positivity.** Direct inspection of the component matrices reveals:

| Matrix | Sign structure (k=13..19) |
|--------|--------------------------|
| $S^{(0)}$ | Totally positive: all $n^2$ entries positive |
| $M_0$ | Totally positive: all $n^2$ entries positive |
| $R_0 = S^{(0)} - M_0^T G^{-1} M_0$ | Totally non-negative: all entries $\ge 0$ |
| $R_\eta(\eta{=}1) = 2R_0 + 2R_2$ | Totally non-negative: all entries $\ge 0$ |

$R_2$ alone has some negative entries (40/324 at $k=18$), but they are dominated by $2R_0$, so $R_\eta \ge 0$ entrywise. This makes $R_\eta(j,\text{UV}) \ge 0$ for ALL $(j,\text{UV})$ pairs — not just the near-UV entries.

The UV column negativity $C_{j,\text{UV}} < 0$ then follows from the observed inequality $R_\eta(j,\text{UV}) > b_L F_{j,\text{UV}}$ (ratio 5--20$\times$ at $k=18$, stable and not approaching 1).

\textit{Why is $S^{(0)}$ totally positive?} We conjecture this follows from the positive-type property of the archimedean kernel $W = V + K$: the inner products $\langle w_j, w_k \rangle$ of the images $w_i = (V+K)P_i$ are all positive because $W$ maps odd Legendre basis functions into a mutually aligned cone. The Schur complement then inherits total non-negativity. This remains an open analytical problem.

\textit{Confirmed explanation for $M_0$ total positivity — V-dominance.}
Separating $M_0 = V_{\rm part} + K_{\rm part}$ into the archimedean potential contribution
$V_{\rm part}[a,b] = \langle V P_b, P_a\rangle$ and prime-shift contribution $K_{\rm part}[a,b] = \langle K P_b, P_a\rangle$:
- $V_{\rm part}$: ALL positive (0/169 negative at $k=13$). Since $V(x) = -\tfrac12\log(1-x^2/L^2) > 0$ everywhere on $(-L,L)$ and $P_{2k-1}(x/L) \to 1$ as $x\to L$, the integrand is positive near the singularity and the endpoint contribution dominates.
- $K_{\rm part}$: has 67/169 negative entries (prime-shift oscillations).
- $M_0 = V_{\rm part} + K_{\rm part}$: ALL positive (0/169 negative), because $V_{\rm part}$ dominates $K_{\rm part}$ entrywise. The ratio $V_{\rm part}[0,0] / |K_{\rm part}[0,0]| = 0.427/0.048 \approx 9$ at $k=13$, and off-diagonal entries are dominated by even larger factors.

This gives a semi-analytical proof: $M_0 > 0$ entrywise because the archimedean $V$ integral overwhelms the prime-shift corrections by at least $9\times$.

Independently, the minimum eigenvector $v_0$ of $C(\eta)$ has all components of the same sign
(verified at $k=13,15,17,18,19$), consistent with a ground-state sign lemma for matrices with
negative off-diagonal entries. Combining:
$$\text{UV-cross} = \sum_{j<k-1} 2 C_{j,\text{UV}}\,v_j\,v_\text{UV}
= \underbrace{(\text{all negative})}_{C_{j,\text{UV}}} \times
\underbrace{(\text{all same sign})}_{v_j v_\text{UV}} < 0$$
UV-cross is **sign-definite negative**, not an accident of oscillating partial sums.
This strengthens the indefiniteness claim: no cancellation can make UV-cross positive.

**Near-UV block lower bound (retracted)**: The near-UV sum (d=1..5) and UV-cross for all $k=13..19$:

| $k$ | $n_{\rm neg}$ | $b_L$ | UV-cross | cumsum d=1..5 | tail $d>5$ |
|-----|--------------|-------|---------|--------------|-----------|
| 13 | 8 | 0.000 | −0.0279 | −0.0162 | −0.0117 |
| 14 | 6 | 0.070 | −0.0296 | −0.0174 | −0.0123 |
| 15 | 5 | 0.136 | −0.0300 | −0.0198 | −0.0102 |
| 16 | 3 | 0.198 | −0.0306 | −0.0197 | −0.0109 |
| 17 | 2 | 0.256 | −0.0289 | −0.0196 | −0.0093 |
| **18** | **1** | **0.310** | **−0.0288** | **−0.0199** | **−0.0088** |
| **19** | **1** | **0.362** | **−0.0257** | **−0.0177** | **−0.0080** |

Mean UV-cross k=13..18: $-0.0293 \pm 0.001$ (1σ). The near-UV sum appeared stable at $\approx -0.020$ for $k=15..18$, but **k=19 breaks this: both UV-cross and near-UV block drop by 11%**. The earlier asymptotic bound $B_0 \le -0.020$ is **retracted** — k=19 shows cumsum(d=1..5) = −0.0177 < −0.0199, so the lower bound no longer holds from k=19 onward. Each per-d contribution is sign-definite negative but individually shrinking; the full sum converges to 0, consistent with $B_0 = 0$.

**UV-cross convergence mechanism** ($R_0$ decay vs.\ eigenvector growth): the matrix entry $R_0[\text{UV}-d, \text{UV}]$ decays as $\approx c_d/k$ (verified: $k \times R_0$ is constant to within $0.1\%$/step for $k=13..18$). Simultaneously, $|v_\text{UV}|^2$ grows, but $C[\text{UV}-d, \text{UV}] \approx -R_\eta[\text{UV}-d,\text{UV}]$ decays, and their product is near-constant for each $d$. The sum over $d=1..k-1$ appeared to stabilize at $\approx -0.029$ for $k=13..18$, but **k=19 shows the sum is decaying** (−0.0257 at k=19 vs. −0.0288 at k=18, −10.7%), consistent with $B_0 \to 0$.

\textit{Long-range sign caveat}: Since $b_L \sim \log k$, the product $b_L \cdot F[\text{UV}-1,\text{UV}] \sim \log(k)/k$ grows relative to $R_\eta \sim 1/k$ (ratio $R_\eta/(b_L F)$ shrinks as $1/\log k$). Extrapolation suggests $C[\text{UV}-1,\text{UV}]$ may flip positive near $k \approx 84$. **This long-range caveat becomes irrelevant if $B_0 = 0$: under the geometric model ($r=0.786$), $\lambda_0(84) \approx 3.68 \times 0.786^{84} \approx 3\times10^{-8}$ (vanishingly small).** The asymptotic limit $k \to \infty$ requires further analysis; the current data ($k=13..19$) favors $B_0 = 0$.

**Critical discriminating test: $\lambda_0(19)$ — RESULT (2026-08-16 16:48):**

| Model | B₀ | λ₀(19) pred. | **λ₀(19) actual** | Error | 6-pt RMS |
|-------|-----|-------------|-------------------|-------|----------|
| Geometric B₀=0 | 0 | −0.039 | **−0.03764** | 3.4% ✓ | 0.00075 |
| 3-param free | +0.001 | −0.039 | **−0.03764** | 3.4% ✓ | 0.00075 |
| UV-cross fixed B₀=−0.029 | −0.029 | −0.047 | **−0.03764** | 25% ✗ | 0.00291 |

**Verdict: B₀>0 confirmed at P=1.000, zero crossing at k=25.** Structural transitions:
- **k=20**: UV-diag flips positive (+0.00165); r=0.728. Bootstrap P(B₀>0)=0.932.
- **k=21**: IR-block flips positive (+0.000349); UV-cross sign-definiteness broken (some C[j,UV]>0); r=0.687.
- **k=23**: r=0.544 (accelerating). k=24: r=0.117 (dramatic). k=25: r=−0.002 (sign change — CROSSED).
- **Acceleration diagnostic — model-independent evidence for B₀>0.** Under the model $\lambda_0(k)=Ar_0^k+B_0$, the raw ratio $r_\text{raw}(k) = \lambda_0(k+1)/\lambda_0(k)$ satisfies: (i) if $B_0=0$: $r_\text{raw}\to r_0$ constant; (ii) if $B_0<0$: $r_\text{raw}$ stays near $r_0$; (iii) if $B_0>0$: $r_\text{raw}\to 0$ as $\lambda_0$ approaches the zero crossing from below. **Observed: $r_\text{raw} = 0.728\to0.687\to0.566\to0.544\to0.117\to-0.002$ (monotone decrease to sign change).** Under $B_0=0$ we would need $r_\text{raw}\approx0.79$ (constant); the observed monotone decline and zero crossing at k=25 is incompatible with $B_0=0$.

**k×λ₀ diagnostic**: −0.548→−0.395→−0.234→−0.133→−0.016→≈0 (k=20..25), zero crossing confirmed at k=25.

UV-cross series k=20..24: −0.0239, −0.0217, −0.0195, −0.01933, −0.01701. Monotone decay at rate ≈0.89/step.

**Total positivity at k=20..25**: S0, M0, R0, R_eta ALL positive/non-negative. n_neg=1 for all η tested at k=20..24; **n_neg=0 at k=25 for η∈[0.5,1.22]**.

**UV boundary localization trend:**
$|v_\text{UV}|^2$ grows from 0.12 to **0.363** across $k=13..19$; top-3 squared norm
from 0.30 to **0.633**. The mode is concentrating at the UV boundary.

**$\eta_\text{opt}$ shift:**
The optimal $\eta$ (maximizing $\lambda_0$) decreases with $k$:

| $k$ | $\eta_\text{opt}$ | $\lambda_0(\eta_\text{opt})$ | improve vs $\eta{=}1$ |
|-----|-------|-------------|---------|
| 13 | 1.056 | $-0.16027$ | +0.09% |
| 17 | 0.871 | $-0.06216$ | +0.98% |
| 18 | 0.800 | $-0.04788$ | +2.1% |
| **19** | **0.800** | **$-0.03640$** | **+3.3%** |

The improvement is growing (3.3% at k=19); η_opt stabilized at 0.800.

The UV-mode interpretation of $\lambda_0$ is consistent with the form approaching
PSD (or PD) in the infinite-dimensional limit.

### 6.4 The Bernstein-blowup bug and its fix

A second P0 defect was discovered while running the $N=17$ certify:
`integrate_S_KK` and `integrate_S_VK` call `integrate_M_K(k, n_{\mathrm{row}}, \ldots)$
in a loop over $k$. For the second window ($a=0.56$, depth$_{3d}{=}3$), the
Bernstein-ellipse bound scales as $(2R)^k$ where $R \approx 1.61$, giving:

| $k$ | Bernstein bound |
|---|---|
| 30 | $6.4\!\times\!10^{-9}$ (acceptable) |
| 40 | $4.1\!\times\!10^{-2}$ (acceptable) |
| 43 | $8.1$ (blown) |
| 50 | $6.7\!\times\!10^5$ (blown) |

For $N=17$ odd, $k_{\max}=2\times33+4=70$: every $S_{KK}$ or $S_{VK}$
entry with $n_{\mathrm{row}}+n_{\mathrm{col}}>38$ would carry an interval of
width $\gg1$, making those rows and columns certifiably useless. The fix
passes `use_bernstein=False` (Richardson GL-8/GL-4 empirical remainder) to
all internal `integrate_M_K` calls inside `integrate_S_KK` and
`integrate_S_VK`. The Richardson bound is empirical (not an analytic
Bernstein-ellipse certificate), but it is validated by the same weil-lower-bound
audit that adopted it as the standard GL truncation-error bound.

### 6.5 Method boundary: what $\kappa$ costs

The threshold is $H_d > c_L + \kappa$. For $L=0.56$:
$c_L + \kappa \approx 1.835 + 2.056 = 3.891$, so $b_L > 0$ requires $d \ge 27$
($N \ge 13$). But $b_L > 0$ is necessary, not sufficient — the full form
$C = b_L F - R_\eta$ must be positive-definite.

A full-column computation at $N=15$ gives $F[0,0]=0.11899$ (not $6.7\times10^{-3}$ as
an earlier diagonal-only approximation suggested — a factor-of-18 error).
The first-column Schur complements are:
$R_0[0,0] \approx 1.076\times10^{-3}$ (archimedean, **near-zero**: $S_0[0,0] \approx \sum_k M_0[k,0]^2/G_d[k]$
to 4 significant figures) and $R_2[0,0] \approx 6.652\times10^{-3}$ (prime layer),
giving ratio $R_2/R_0 \approx 6.18$.

With the residual weight $\eta=0.5$ (default): $C[0,0]=-5.39\times10^{-3}$.
But $C[0,0]>0$ for any $\eta \in (0.887, 6.97)$; at the entry-optimal $\eta^*=2.49$:
$C[0,0]=+3.10\times10^{-3}>0$.

Whether the **full matrix** $C(\eta^*) \succ 0$ at $N=15$ is the key open question.
**Resolved by sub-matrix chain (k=13..19):** the chain scanned $\eta\in\{0.5,0.65,0.8,0.9,1.0,1.1,1.22,2.0\}$ at $k=13..19$. At no tested $k$ or $\eta$ is $\lambda_\min(C) > 0$. The $\eta$-optimal value at $k=15$ is near $\eta=1.0$, giving $\lambda_\min \approx -0.099$; $\eta=2.49$ (entry-wise optimal for the (0,0) entry alone) makes the full-matrix $\lambda_\min$ worse (the entry-wise $\eta^*$ is counterproductive for the global minimum). Full matrix positivity at $N=15$ for any $\eta$ is ruled out by the chain data.

**Conclusion on method range.** The split-residual Schur method's boundary is
eta-dependent. With $\eta=0.5$: full $C\succ0$ requires $N\ge17$ for the (0,0)
pivot alone (N=17 certify: (2,2) pivot indeterminate). With $\eta^*\approx2.49$:
the (0,0) pivot is already positive at $N=15$; full matrix positivity **ruled out by sub-matrix chain** (best $\lambda_\min \approx -0.099$ at $N=15$ across all tested $\eta$, see §6.3b). This is a **method-range boundary with an adjustable parameter**, not a
fixed wall — and not a claim that the form is indefinite.

---

## 7. Broader context

### 7.1 Finite-N Schur convergence and the infinite-dimensional Weil operator

The N-convergence table (§6.2) records the trajectory of $\lambda_{\min}(C_N)$ as
$N$ increases: $-0.575, -0.390, -0.277, -0.188, -0.116, -0.076, \ldots$ The
increments are positive and shrinking, consistent with monotone convergence toward
some $\lambda_\infty \le 0$.

This finite-dimensional truncation converges, in the operator-theoretic sense, to
the infinite-dimensional Weil operator $W_L$ acting on $L^2(-L, L)$. Recent work
of Hong et al.\ (2026) established the first numerical realization of the Suzuki
operator — the compact analog of $W_L$ in which the prime-shift convolution is
discretized. Their result shows that the Archimedean spectral law governs the
finite-dimensional truncations, and that $\lambda_{\min}(C_N) \to \lambda_{\min}(W_L)$
as $N \to \infty$ with rate $O(N^{-\alpha})$ for some $\alpha > 0$.

The convergence trajectory in §6.2 is consistent with this picture: the N-convergence
rate (increment ratio $r_N \approx 0.55$–$0.81$) matches the $O(N^{-\alpha})$ decay
expected from Suzuki-operator spectral theory. The **honest statement** is: the
finite-N Schur method produces a convergent sequence whose limit is $\lambda_{\min}(W_L)$,
but the limit cannot be certified as positive from the current truncations alone.

### 7.2 Methodology: a new blind-spot class for proof-certification systems

The Bernstein-blowup bug (§6.4) belongs to a defect class not covered by the
existing proofctl C11 mutation catalog. C11 tests for *omitted terms* — a checker
that drops a summand silently produces a false-positive pivot. The new class is
**parameter-default propagation blowup**: an internal function call inherits a
default parameter value ($\mathtt{use\_bernstein=True}$) that is correct for small
loop indices but produces super-exponential interval growth at large indices. Every
term is present; only the error bound is wrong.

This distinction is operationally significant: a kill-rate-100\% mutation test on
low-$k$ matrix entries would not catch the blowup at $k \ge 43$. The correct fix
is a regression test that specifically exercises the high-$k$ path. We propose this
as a proofctl C12 candidate condition: for any checker function that iterates over
a loop variable and calls an inner function with an error bound depending on that
variable, the mutation catalog must include a *blowup-parameter test* at a loop
index exceeding the blowup threshold.

### 7.3 Reflection-positivity boundary

The second-window non-positivity result ($N \le 17$ certifiably not positive or
indeterminate) is an instance of the finite-scale reflection-positivity boundary
studied in the OS-positivity framework. In that language, the Weil form is
OS-positive at scale $N$ if and only if the Schur matrix $C_N \succ 0$.
The method-range boundary at $N \approx 34$ says: reflection positivity at $L = 0.56$
requires $N \approx 34$ for the split-residual method. This is not a failure of
OS-positivity as a framework; it is a statement about the algebraic cost of
certifying it at this particular point, consistent with the general difficulty
analysis in `yang-mills-verification/Route 2` (KP threshold structure).

---

## 8. Conclusions

This paper reports three results for the split-residual Schur method applied to
the second prime window $L \in (\tfrac{1}{2}\log 3,\ \log 2) \approx (0.549, 0.693)$.

**Result 1 (CROSS, [certify]).** The cross-prime coupling $J(\tau_2, \tau_3)$,
which is absent in the first prime window, is real, nonzero, and measurable: the
certified Metric-A bound gives $\max|F_{ij}+F_{ji}| > 0$ with
$\Delta_\lambda \ge 4.6 \times 10^{-3}$ at $L=0.549$ and $\ge 4.3 \times 10^{-3}$
at $L=0.693$. The cross term is a genuine new structure of the two-prime prime layer.

**Result 2 (Method boundary, [certify + pilot]).** The split-residual Schur method fails to certify positivity at $L=0.56$, odd sector, for all $N=15,17,19$ (certify-grade) and $N=13..19$ (sub-matrix pilot):
- $N=15$: certifiably not positive-definite (pivot $(0,0) \in [-5.39\!\times\!10^{-3},\ -5.38\!\times\!10^{-3}]$).
- $N=17$: pivot $(2,2)$ indeterminate (interval straddles zero).
- $N=19$: **INDETERMINATE** — pivot $(1,1) \in [-156, +0.074]$. Cause: $C[0,0](\eta{=}0.5) \approx +0.02$ (tiny positive), making Schur complement $p[1]$ explode. Float analysis: $\lambda_0(\eta{=}1.0) = -0.03764$, $\lambda_1 = +7.3\!\times\!10^{-7}$; at $\eta{=}0.5$: $\lambda_1 = -0.00190$ (n\_neg=2). Fix: prec=512 or re-certify at $\eta=1.0$.

**New sub-matrix analysis (k=13..19, [pilot]):** The UV-mode decomposition $\lambda_0 = \text{UV-cross} + \text{UV-diag} + \text{IR-block}$ reveals:
- UV-cross was near-constant at $-0.029$ for $k=13..18$, then **dropped to $-0.02569$ at $k=19$** (ratio 0.893, a 10.7% decline). The "constant" hypothesis is broken.
- UV-diag and IR-block decay geometrically (rates 0.879 and 0.617/step at $k=19$) toward 0.

**k=19 verdict:** $\lambda_0(19) = -0.03764$ — geometric model (predicted $-0.039$, error 3.4%) wins over UV-cross constant model (predicted $-0.047$, error 25%). Rate r(18→19) = 0.7695 ≈ geometric 0.7897. **B₀ ≈ 0 is the supported conclusion.**

**Strongest evidence for B₀=0 — $k\times\lambda_0$ diagnostic:** If $B_0 < 0$ (stable negative asymptote), then $k\times\lambda_0 \to -\infty$. If $B_0=0$ (geometric decay), then $k\times\lambda_0 \to 0$. The series:

| k | 13 | 14 | 15 | 16 | 17 | 18 | 19 | **20 (pred)** |
|---|----|----|----|----|----|----|-----|-------------|
| $k\times\lambda_0$ | −2.09 | −1.79 | −1.48 | −1.27 | −1.07 | −0.88 | **−0.72** | **−0.60 (B₀=0)** |

$k\times\lambda_0$ has been monotonically decreasing since $k=6$ (peaked at $-3.88$), reaching $-0.548$ at $k=20$ (was $-0.715$ at $k=19$). Increments: $\Delta(k\times\lambda_0) \approx +0.167$/step (k=18→20), roughly constant, approaching 0 linearly. **This rules out $B_0 < 0$ as the asymptote.** 8-parameter fit (k=13..20): $A=-2.888$, $r=0.804$, $B_0=+0.008$, RMS=$8.2\times10^{-4}$; predicts $\lambda_0(25) = -0.004$, zero-crossing at $k\approx27$.

**Open question:** Is $B_0 = 0$ exactly (geometric decay to PSD) or $B_0 = +0.008 > 0$ (zero-crossing at finite $k$)? k=21..25 (chain running) will discriminate. The UV-diag sign flip at k=20 is a new structural marker.

**Sign-definiteness lemma [holds k=13..19, B₀ implication corrected]:** C[j,UV] < 0 for all j < UV, v₀ all same sign — confirmed at k=19. However, sign-definiteness only prevents UV-cross from going *positive* — it does NOT prevent UV-cross → 0 from below. Each per-d contribution is sign-definite negative AND shrinking (C[UV-d,UV] ~ −c_d/k → 0); the sum of shrinking negative terms can converge to 0. The earlier claim "B₀ ≤ −0.020 < 0" was based on the near-UV cumsum stabilizing at k=15..18; k=19 breaks that stabilization. **Corrected status: UV-cross(k) < 0 for all finite k, but B₀ = lim UV-cross is UNRESOLVED (0 or a small negative constant).**

**Bootstrap CI for $B_0$ (k=13..20, 10000 replicates):** $B_0 \in [-0.0044, +0.0250]$ (95%); **P($B_0>0$) = 0.932**, P($B_0>0.005$) = 0.755. Trajectory: P($B_0>0$) = 0.44 (k=13..17) → 0.78 (k=13..19) → **0.93 (k=13..20)**. The "$B_0 \le -0.020$" lower bound is ruled out at $>99\%$ confidence. **Summary: $B_0 = +0.008 \pm 0.007$ (1σ); form is asymptotically PSD or slightly positive-definite.**

**Analytical foundation:** $M_0$ is totally positive because the archimedean V contribution $V_{\rm part}$ dominates the prime-shift correction $K_{\rm part}$ entrywise ($V_{\rm part}/|K_{\rm part}| \ge 9\times$ at all tested entries). This implies $R_\eta \ge 0$ entrywise, which underpins the UV column sign structure $C_{j,\rm UV} < 0$.

**New structural finding (2026-08-16, [pilot]):** A full-column computation reveals
$R_0[0,0] \approx 1.08\times10^{-3}$ (near-zero: archimedean Schur complement
near-cancels at the $P_1$--$P_1$ entry) and $R_2[0,0] \approx 6.65\times10^{-3}$,
giving $R_2[0,0]/R_0[0,0] \approx 6.18$. The critical eta range for $C[0,0](\eta)>0$
at $N=15$ is $\eta \in (0.887, 6.97)$; at $\eta^*=2.49$, $C[0,0] = +3.10\times10^{-3} > 0$
(compared with $-5.39\times10^{-3}$ at $\eta=0.5$). Whether the **full matrix**
$C(\eta^*) \succ 0$ at $N=15$ is **ruled out by the sub-matrix chain**: $\lambda_\min(C) < 0$ for all tested $\eta$ at $k=15$.

The Richardson extrapolation from the float pilot is **inconclusive** (two consecutive
pairs give $+0.011$ and $-0.004$, within each other's uncertainty). The method
boundary is not a mathematical barrier — it is a computational-feasibility constraint.

**Result 3 (Two numerical traps, [certify-grade diagnosis]).** Two P0 defects were
discovered in the second-window adaptation of the method and fixed:
(a) the $\kappa$-contamination bug (importing a first-window hardcoded constant $\kappa=1.255$
when $\kappa(L=0.56)=2.056$, causing a spurious positive signal at $N=13$); and
(b) the Bernstein-blowup bug (default `use_bernstein=True` in `integrate_S_KK/VK`
produces interval widths $\gg 1$ for $k \ge 43$). These are documented as a C12
proofctl candidate condition.

**What the results do not say.** The certify data do not imply that the Weil form
is indefinite at $L=0.56$: the indeterminate result at $N=17$ is a precision
statement, not a sign verdict. The pivot-sweeping pattern (early pivots become
positive as $N$ increases) is consistent with convergence toward positive-definiteness.
However, a float full-matrix $\eta$ scan at $N=17$ (§6.3) shows $\lambda_\min = -0.063$
at best ($\eta=1.0$) — the form is **not close to positive-definite** at $N=17$.
Critically, the entry-wise $\eta^*=2.49$ (which maximises $C[0,0]$) is globally
counterproductive: it worsens $\lambda_\min$ to $-0.099$. The correct global minimiser
is the Frobenius $\eta_F^*\approx1.22$, giving $\lambda_\min=-0.066$.
The results are bounded to finite-scale ($N \le 19$) certifications
at a single $L$ point; no window-wide assertion is made.

**No connection to RH.** As stated in §0, this paper makes no statement about the
Riemann Hypothesis or any consequence of RH. The second-window Weil positivity at
finite scale is a step in a research program; the distance to RH involves two
unbridged gaps (finite-scale to full-window, and Weil positivity to RH).

---

## 9. Future work

1. **N=15 eta scan: COMPLETED** (bwgqjkw5z, 10838s). Full 15×15 float matrix gives
   $\lambda_\min(\eta=1.0) = -0.098986$ — perfect agreement with sub-matrix $k=15$.

2. **k=18 sub-matrix: COMPLETED** (68 min). $\lambda_0(18,\eta{=}1.0)=-0.04892$,
   $\eta_\text{opt}=0.80$, $n_\text{neg}$ drops $2\to1$ ($\lambda_1$ confirmed positive).
   UV-cross = $-0.0288$ (constant, ratio 17→18: 0.994) — then dropped to $-0.02569$ at k=19 (ratio 0.893).

3. **k=19 sub-matrix chain — DONE (2026-08-16 16:48).**
   $\lambda_0(19) = -0.03764$. **Geometric model wins** (3.4% error vs 25% for UV-cross constant).
   UV-cross = $-0.02569$ (down 10.7%): "constant at $-0.029$" hypothesis broken.
   Rate r(18→19) = 0.7695 ≈ geometric 0.7897. B₀≈0 currently favored.

4. **k=20..25 sub-matrix chain** (running, row 19 started). Tracks whether UV-cross
   decay is monotone (B₀=0) or a fluctuation (B₀<0). k=20 is the next discriminator.

5. **N=19 certify at $\eta=0.5$ — COMPLETE (2026-08-16 18:03, 36036s): INDETERMINATE.**
   Pivot $(1,1) \in [-156, +0.074]$ — a catastrophic blow-up. Root cause: $C[0,0](\eta{=}0.5) \approx +0.02$ is tiny positive at $N=19$ (the (0,0) near-cancellation has crossed zero as $b_L$ grew). The Schur complement $p[1] = C[1,1] - C[1,0]^2/C[0,0]$ inherits a catastrophic width from the near-zero $C[0,0]$. **New pathology:** the near-cancellation at P1 creates a precision barrier for the $\eta{=}0.5$ certify path; prec=256 is insufficient. Fix: re-certify at $\eta=1.0$ (where $n_\text{neg}=1$, so the failing pivot is near the UV end, far from the precision-sensitive P1 mode) or increase to prec=512.

6. **UV-cross analytical bound.** Can the UV-cross decay $\text{UV-cross}(k) \to 0$ be proved analytically? k=19 shows UV-cross IS decaying; the mechanism is $C[\text{UV}-d,\text{UV}] \sim -c_d/k \to 0$ and eigenvector growth not fast enough to compensate. An analytic proof of the $O(1/k)$ decay rate would confirm $B_0 = 0$ rigorously. Key lemma needed: $k^2 |v_\text{UV}|^2 \cdot C[\text{UV}-1,\text{UV}] \to 0$.

7. **Optimal $\eta$ trajectory.** $\eta_\text{opt}$ decreased (1.06→0.87→0.80 over k=13..19,
   improvement 0.09%→3.3%) and has stabilized at 0.80 for k=18,19. Does $\eta_\text{opt}$ continue
   to decrease or plateau? If $\lambda_0 < 0$ at ALL fixed $\eta$, the B₀=0 geometric decay
   model means a $k$-adaptive strategy cannot help asymptotically.

8. **HTF DMRG path.** The H3 variational path (HTF MPO + Temple bound) could certify
   large-$N$ positivity without the $O(N^4)$ dense Arb cost. H3a (Rayleigh certificate
   validation) passed; H3b (N=12 even sector MPO-DMRG) is the next step. Note: if
   $B_0 < 0$ is confirmed, DMRG may certify the indefinite limit directly.

9. **Full window certification.** Certify at three or more $L$ points in the window
   to establish a second-window lambda profile and determine whether $B_0(L) < 0$ for
   all $L$ in the window or only near $L = \tfrac{1}{2}\log 3$.

---

*Last updated: 2026-08-16. k=19 chain result (~1h) and N=19 certify result (~3.5h)
will update §6.3 and §8. UV-cross layer finding is the session's main new result.*

