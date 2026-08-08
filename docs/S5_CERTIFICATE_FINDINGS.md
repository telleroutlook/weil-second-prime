# S5 — Schema, domain, and the first pilot certificate

**Date: 2026-08-08.** S5 delivers the second-window certification infrastructure
and the first pilot certificate. The certificate asserts ONLY what is
certify-true; it does NOT claim second-window positivity (which is false at the
left end — see below).

## What is certified: the cross-term new-structure finding

`claim_id: thm-second-cross-structure`. Statement: in the window
$(\tfrac12\log3,\log2)$, the cross-prime term
$F_{ij}(\tau_2,\tau_3)=\langle C_{\tau_3,1}P_j, C_{\tau_2,1}P_i\rangle$ is a real,
nonzero contribution to $S^{(2)}$ — the genuine new structure the first window
lacks (S4 metric A: $\Delta C=-3c_2c_3(F_{ij}+F_{ji})$ nonzero at Arb grade).

The checker (`checker/second_prime/check_cross_structure.py`) recomputes the
cross term in exact `Fraction` (fast — the cross term is purely the prime layer,
no archimedean integrals) and verifies six obligations:
1. `second.cross-term-present-and-nonzero` — recompute $\max|F_{ij}+F_{ji}|>0$.
2. `second.both-shifts-present` — recomputed $J(\tau_2)$ AND $J(\tau_3)$ nonzero.
3. `second.window-bounds-hold` — $\tfrac12\log3<L<\log2$ by certified rational
   bounds on $\log2,\log3$ (arctanh series, computed in-checker, no float).
4. `second.four-term-S0-declared` — $S_{VV}+S_{VK}+S_{KV}+S_{KK}$ + cross term.
5. `second.positivity-not-claimed` — no positivity/pivot/conclusion fields.
6. `second.conclusion-bounded-and-no-rh` — $L<\log2$, two-prime method, primes=[2,3].

## What is NOT claimed: positivity (certifiably false at the left end)

S5 also ran the certify-grade positivity scan from the window LEFT END
(`scripts/scan_positivity_leftend.py`, `pilots/s5_positivity_L055.json`), per the
supervisor steer (c_L smallest there; first-window positivity lived only near the
left edge). At $L=11/20=0.55$ (edge+0.0007, $\tau_3=1.9975$, prime 3 barely
single-hop):

| sector | full min-pivot interval | verdict |
|---|---|---|
| even N=8 | [−3.869e-4, +2.428e-4] | **indeterminate** (straddles 0; interval width > distance to 0) |
| odd N=7  | [−2.815e-2, −1.720e-2] | **certifiably NEGATIVE** (upper < 0) |

The odd sector is certifiably negative even at the most-favorable left end, so
**the second window is NOT positive-definite there** — and the certificate does
not claim otherwise. Compared to the first window (even +0.0087, odd +0.053 at
L=7/20, both clearly positive), the second window sits at or past the boundary.
This is consistent with the honest prediction that the method's certifiable reach
is SHORTER in the second window (larger $c_L\approx1.82$ + prime-3 negative
contribution). This is an E2-type boundary result, not a failure.

## Mutation coverage (C11 + two-prime extensions)

`tests/mutation/mutation_catalog_second_cross.py`: baseline certifies,
**kill_rate 100% (8/8)**. Mutants (each killed by a sensible obligation):
zero τ₂ shift, zero τ₃ shift, drop cross term F, S_KK-only S0, L>log2 (window
high), L<½log3 (window low), wrong primes (single-prime masquerade), wrong method
(first-window). The window-low mutant is also caught by the cross-term recompute
(at L=1/2, τ₃>2 leaves the single-hop regime and F vanishes) — the recompute is
real, not a declaration check.

Honest gap: the CLAUDE.md "swap c₂↔c₃" mutant is not literally in this catalog —
in a *structure* checker (which certifies nonzero-ness, not a positivity value)
c₂/c₃ enter only symmetric recomputed quantities, so a pure swap need not flip
the verdict. The c-swap belongs to a future *positivity* checker (when the window
has a positive point to certify) and is recorded here rather than faked.

## Three propositions, strictly separate (PROOF_CONSTITUTION D3)

1. Cross term is significant new structure — **CERTIFIED** (this certificate).
2. Second window is positive-definite — **NOT established** (odd negative, even
   indeterminate at the left end). Independent question, not this certificate's.
3. Approaching RH — **no**. Finite-scale, $L<\log2$, no implication.

## Files

- `schemas/certificate-second-prime-v1.schema.json`
- `domains/fp_second/contracts/thm-second-cross-structure.json`
- `checker/second_prime/check_cross_structure.py`
- `pilots/cert_second_cross_structure.json` (CERTIFIED)
- `pilots/s5_positivity_L055.json` (left-end boundary)
- `tests/mutation/{mutation_catalog_second_cross.py,test_mutation_second_cross.py}`
- `tests/test_second_cert_endtoend.py`
