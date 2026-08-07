# weil-second-prime

Certificate-first proof repository for the **second prime window** of the Weil
quadratic form: does the local Weil quadratic form have strictly positive
infimum for $L$ in the interval

$$\tfrac{1}{2}\log 3 \;<\; L \;<\; \log 2 \qquad (\approx 0.5493 < L < 0.6931)?$$

This window is the sibling of the *first* prime window $(\tfrac12\log2,\tfrac12\log3)$
studied in [`weil-first-prime`](https://github.com/telleroutlook/weil-first-prime)
(FP-0.35, published at doi:10.5281/zenodo.21807498). The second window is
genuinely new mathematics: **two primes** $p=2$ and $p=3$ both contribute to the
explicit formula (both in the single-hop regime), so the prime layer carries a
**cross-prime coupling** $J_{ij}(\tau_2,\tau_3)$ absent from the first window.

## Scope boundary (non-negotiable)

- The only conclusions are **finite-scale Weil positivity** at specific $L$ in
  the second window. This does **NOT** imply the Riemann Hypothesis, nor
  "near-RH", nor global positivity.
- The window is horizontal method extension, not vertical progress toward RH.
  There remain two gaps with no known path (finite-scale → full-interval
  positivity, then Weil equivalence); this repo does not touch them.
- **Hard boundary at $L=\log 2$**: prime $n=2$ exits the single-hop regime and
  $n=4$ enters. The Theorem-3 three-interval framework needs genuine extension
  there. Do not extrapolate beyond $L < \log 2$.

## Why this repo exists (dual purpose)

1. **New mathematics.** Extend the split-residual Schur certification method to a
   window where the prime layer is a genuine two-prime object with cross coupling.
   New algebra ($\mathbb{Q}[\tau_2,\tau_3]$), new numerical traps, new $c_L\approx1.82$.
2. **A live host for `proofctl`.** The first-window pilot forced two verification
   kernel conditions (C10 no-copy generator, C11 checker mutation coverage) into
   `proofctl` — because a *real* proof produced real bugs a tool cannot invent by
   self-testing. The second window is a fresh live host: its new bug classes
   (cross-prime omitted terms, incommensurable shifts) are expected to drive the
   next round of `proofctl` conditions. The tool co-evolves with the research;
   it is not a standalone line.

## Status

**Scaffold** (2026-08-08). No certificate yet. See `PLAN.md` for the phased plan
and `HANDOFF.md` for the operational entry point. The first concrete action is a
per-sector prime-influence profile (cheap mutation-style probe) to decide where
the cross-prime coupling actually moves the pivot — do **not** budget compute
symmetrically between sectors (first-window even-sector prime term was nearly
inert; the tension lives in the odd sector / larger $L$ / cross terms).

## Relationship to weil-first-prime

Shared machinery (archimedean integrators, interval arithmetic, LDL$^\mathsf{T}$,
Legendre algebra) is ported from `weil-first-prime`; the prime layer and domain
are new. The prototype `src/prime_layer/legendre_shift_2prime.py` in
`weil-first-prime` is the seed for this repo's two-shift algebra.
