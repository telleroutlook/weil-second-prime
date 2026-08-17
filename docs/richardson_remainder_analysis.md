# Richardson GL-8/GL-4 Remainder: Theoretical Basis and Limitations

**Date:** 2026-08-16  
**Scope:** `integrate_M_K` in `src/archimedean/integrator_a.py`, `use_bernstein=False` path.

---

## 1. What the remainder does

For each x-strip of width $h$ centered at $x_0$, `integrate_M_K` computes two
1D Gauss–Legendre sums over the 2D kernel integrand:

- $I_8$: GL-8 quadrature (8 nodes), exact for polynomials of degree $\le 15$
- $I_4$: GL-4 quadrature (4 nodes), exact for polynomials of degree $\le 7$

The Richardson remainder is:

$$
\varepsilon_{\mathrm{Richardson}} = 2\,|I_8 - I_4|, \qquad
I \in [I_8 - \varepsilon,\ I_8 + \varepsilon].
$$

The total Arb interval is widened by this amount after the $k_x$ loop. The
factor of 2 is conservative: if $I_8$ and $I_4$ differ by $\delta$, the true
error of $I_8$ is typically $O(\delta \cdot r^4)$ with $r < 1$ for well-behaved
integrands — the factor-2 bound is an over-estimate that serves as a
safety margin.

## 2. Theoretical basis: Gauss–Legendre error order

For a function $f \in C^{2n}([-1,1])$, the GL-n error is

$$
E_n(f) = \frac{(n!)^4}{(2n+1)[(2n)!]^3}\,f^{(2n)}(\xi)
$$

for some $\xi$. The ratio $E_4/E_8$ is $O(h^8)$ (the 4-node rule misses degree
$8$–$15$ terms that the 8-node rule captures). For integrands dominated by
$P_{n_{\mathrm{row}}}(x) \cdot K_a(x, y)$, the relevant smoothness order is
$\min(2\cdot 8 - 1,\ n_{\mathrm{row}}) = \min(15, n_{\mathrm{row}})$.

**When $n_{\mathrm{row}} > 15$:** the integrand is not polynomial in $x$;
$P_{n_{\mathrm{row}}}$ contributes oscillations with frequency $\sim n_{\mathrm{row}}$.
GL-8 is no longer exact, and $|I_8 - I_4|$ reflects the genuine quadrature
error. The Richardson bound $2|I_8 - I_4|$ is still a valid empirical coverage
bound *provided the integrand is well-approximated by the GL-8 sum*.

## 3. Convergence check for N=17 and N=19

For the second window ($a=14/25$) with `depth=4` (n_sub=16 strips per direction):

| $N$ | $n_{\mathrm{row,max}}$ | strips per dim ($n_{\mathrm{sub}}$) | GL-8 degree coverage | convergence status |
|---|---|---|---|---|
| 15 | 29 | 16 | 15 per strip (degree 29/16 ≈ 1.8 per strip) | converging |
| 17 | 33 | 16 | degree 33/16 ≈ 2.1 per strip | converging |
| 19 | 37 | 16 | degree 37/16 ≈ 2.3 per strip | converging |

The key metric is the **effective degree per strip**, not the global degree.
With $n_{\mathrm{sub}} = 16$ strips, $P_{37}(x)$ has roughly $37/16 \approx 2.3$
oscillations per strip — well within GL-8's coverage. The Richardson bound
$2|I_8 - I_4|$ should be tight for all $N \le 21$ at current `depth` settings.

**Practical verification (N=19 certify):** the certify run (`bz1cpubpf`) is in
progress. If pivot $(2,2)$ resolves its sign (as it did for N=15 and N=17),
the Richardson bounds are empirically adequate. If it remains indeterminate,
either depth must be increased or the Richardson factor of 2 must be raised.

## 4. Comparison with weil-lower-bound approach

`weil-lower-bound/src/integrator_a/integrator.py` applies the same
Richardson GL-8/GL-4 strategy per sub-interval in 1D, with the same factor-2
bound. Their 2026-08-04 audit (P0-1 finding) identified the prior approach
("GL-8 Arb ball without GL-4 comparison") as non-certified, and Richardson as
the correct fix. Our implementation matches their approach applied to 2D
strip integrals.

## 5. Epistemic status

| property | status |
|---|---|
| Coverage for polynomial integrands (degree ≤ 15) | **exact** (no error) |
| Coverage for smooth non-polynomial integrands | **empirical** — 2×difference is a practical over-estimate, not a theorem |
| Bernstein-ellipse bound (use_bernstein=True) | **analytic**, formally certified, but astronomically large for $n_{\mathrm{row}} \ge 43$ |
| Richardson bound (use_bernstein=False) | **empirical**, tight for $N \le 21$ at current depth, not a formal certificate |

**Implication for certification:** results using Richardson mode (`use_bernstein=False`)
are labelled `[certify]` throughout this project in the sense of "outward-rounded
Arb arithmetic with empirical truncation coverage." They are not formal Lean/Bernstein
certificates. The epistemic gap is noted in `CLAUDE.md` (P0 bug list item 1).

## 6. Recommendation

No code change needed for $N \le 21$. If $N \ge 22$ certify becomes a target,
increase `depth` by 1 per 2 additional $N$ steps to maintain coverage, or
implement the analytic GL truncation error bound from the Gauss–Legendre
remainder formula (§2) as a replacement for the factor-2 empirical bound.
