# CLAUDE.md — weil-second-prime development rules

/ This repo inherits the discipline of weil-first-prime. Where a rule is
identical, it is restated briefly; where the second window differs, the
difference is called out. The authoritative epistemic discipline lives in
`docs/PROOF_CONSTITUTION.md` (ported); CLAUDE.md points there, does not duplicate.

## Project identity

weil-second-prime is a **certificate-first proof repository** for second-prime-
window Weil positivity: the Weil quadratic form on $L^2(-L,L)$ has strictly
positive infimum for $L \in (\tfrac12\log3,\ \log2)$. Both primes $p=2,3$ are in
the single-hop regime here, so the prime layer is a genuine **two-shift** object
with cross-prime coupling $J_{ij}(\tau_2,\tau_3)$, $\tau_p=\log p / L$.

## Hard invariants (inherited, non-negotiable)

- **No PASS/RELEASED in certificate JSON.** Only `proofverify` derives status
  from obligations. A certificate that self-reports its conclusion is a P0 defect.
- **No floating-point conclusions.** All bounds entering a certificate are
  outward-rounded interval arithmetic (Arb balls). Float centers are pilot only.
- **Full second moment, no omitted terms.** $S^{(0)} = S_{VV}+S_{VK}+S_{KV}+S_{KK}$
  (four terms). The second window adds a two-prime $S^{(2)}$ with cross terms
  $S^{(2)} = \langle(V{+}K)P_j, C_{\tau_2}P_i\rangle + \langle(V{+}K)P_j, C_{\tau_3}P_i\rangle$
  — **every prime shift and every cross term must be present**. Omitting a
  second-moment or cross term shrinks the residual and yields a false-positive
  pivot: the exact bug class C11 exists to catch. This is the single most
  dangerous pattern in this repo.
- **min-pivot judge.** Positivity is the min LDL$^\mathsf{T}$ pivot, not the
  symmetrized min eigenvalue.
- **Window check mandatory.** Any certificate must carry
  $\log 3 \le 2L$ ... wait — the second-window condition is
  $\tfrac12\log3 < L < \log2$, i.e. certified rational bounds proving
  $\log 3 < 2L$ and $L < \log 2$. Verify by certified rationals, not enumeration.
- **Conclusion boundary.** Published conclusions are bounded to "finite-scale
  Weil positivity at $L$ in the second window, $L < \log 2$." Never write RH,
  "near RH", or extrapolate past $L=\log 2$.

## proofctl integration (inherited C01–C11 + expect new conditions)

- Every certify assertion must pass proofctl C01–C11. C10 (no copy-only
  generator) and C11 (checker mutation coverage, kill_rate 100% + catalog digest)
  are enforced; do not disable them or fake mutation numbers.
- **Co-evolution expectation.** The two-prime layer is a new bug host. When a new
  blind-spot class appears (e.g. a cross-prime omitted-term that C11's single-term
  catalog misses, or an incommensurable-shift numerical trap), the correct
  response is to (a) fix the checker honestly, (b) extend the mutation catalog to
  cover the new term class, and (c) if it is a *kernel* gap, propose a new
  proofctl condition upstream (`~/github/proofctl`) — exactly as the first-window
  pilot produced C10/C11. Record every such finding for the methodology paper.
- `runtime.class` is `scripted` for native Python checkers (reaches
  GLOBALLY_VERIFIED). Never use `native`, `native-dev`, or `wasi`.

## Long-running computation requirements (inherited)

Any computation > 30s MUST be: (1) **observable** — `print(..., flush=True)` per
unit of work, format `[sector] step N/total (elapsed Xs)`; (2) **pausable** —
catch KeyboardInterrupt, save a JSON checkpoint; (3) **resumable** — `--resume`
loads the latest checkpoint and skips completed work. The second-window
integrals are heavier ($c_L\approx1.82$, two shifts); certify-grade runs are
long. Use `~/.local/bin/run_and_wait.sh -t <sec> -- <cmd>`; never bare `&`.

## Python conventions (inherited)

- After any Python change: `python -m pytest tests/ -x` — zero failures is the bar.
- Numeric results entering a certificate must be `python-flint` Arb balls with
  outward rounding. Never pass `float()` into a certificate.
- `Fraction` arithmetic for $\mathbb{Q}[\tau_2,\tau_3]$ polynomial algebra — no
  mpmath/sympy shortcuts in the exact algebra.
- stdlib only in `checker/` and `schemas/` (no numpy/scipy in checker code).
- Type annotations required in `src/` and `checker/`.
- Checker exit codes: 0=certified, 1=uncertified, 2=malformed/resource, 3=blocked.

## Ported P0 bugs (never re-introduce)

1. `integrate_M_K` must call the 1-D integrator with a GL-8/GL-4 remainder;
   returning a raw GL-8 Arb ball without truncation-error coverage is a P0 defect.
2. Near-zero Taylor cubic coefficient for the RPP kernel is $7s^3/11520$, not
   $s^3/2880$. A remainder without an analytic domain + Bernstein-ellipse bound +
   theorem constant is not a certified remainder.

## Schema conventions (inherited)

- `additionalProperties: false` on every schema — unknown fields rejected.
- `format_version: "second-prime-1.0"` and `method: "exact_two_prime_split_v1"`
  are const fields.
- No matrix/eigenvalue/pivot/conclusion values in certificate JSON — recomputed
  by the checker or refused as unknown fields.

## Mutation / negative test requirements (extended for two primes)

Every push must pass the mutation suite in `tests/mutation/`. Beyond the
first-window mutants, the two-prime layer requires:
- Zeroing the $\tau_2$ shift alone → checker rejects.
- Zeroing the $\tau_3$ shift alone → checker rejects.
- Dropping a cross term of $S^{(2)}$ → checker rejects.
- Swapping $c_2=\log2/\sqrt2$ and $c_3=\log3/\sqrt3$ → checker rejects.
- Changing the window bounds ($\tfrac12\log3$, $\log2$) → checker rejects.
Kill-criterion nuance (PROOF_CONSTITUTION E1): a term with tiny true influence is
probed by a large-factor scaling mutant, not by demanding a sign flip it cannot
cause. Profile per-sector prime influence FIRST.

## Commit conventions (inherited)

- English commit messages.
- Do not commit floating-point discovery pilot values as proof artifacts.
- Do not commit certificates not independently replayed by the checker.
- `git status` before any commit; never stage `.proofctl/attestations/` or
  `.proofctl/keys/*.priv` (private keys are gitignored).

## What this project does NOT do

- Does not claim to prove RH or any consequence of RH.
- Does not extrapolate beyond $L = \log 2$.
- Does not supersede proofctl (proofctl is the orchestration layer, co-evolving).
- Scope is the second prime window only.
