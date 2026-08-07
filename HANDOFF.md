# Handoff — weil-second-prime (scaffold, 2026-08-08)

## 1. One-line status
S2 DONE (2026-08-08). Shared archimedean machinery + single-prime layer ported
verbatim from weil-first (82 upstream tests + 16 new self-check tests = 98 pass).
The S2 acceptance gate — single-prime-limit self-check — PASSES: with prime 3
switched off (c3=0), the two-prime layer reproduces weil-first's assembled Schur
matrix element-wise (max|dC| = 0.00e+00, exact). Next: S3 (cross-prime term F).
**Scope: finite-scale positivity, does NOT imply RH.**

## 2. What exists now
- `README.md`, `CLAUDE.md`, `PLAN.md`, this file — identity, rules, plan.
- `src/archimedean/{integrator_a,integrator_b,interval,ldlt,log_moments,kernel,bernstein}.py`
  — ported verbatim (trusted, four-term S0, P0 bugs already fixed upstream).
- `checker/archimedean/{integrate,check_archimedean,replay}.py` — shared checker machinery.
- `checker/fp035/recompute_schur.py` — trusted four-term S0 + min-pivot, single-prime
  ground truth for the S2 self-check.
- `src/prime_layer/legendre_shift.py` — single-prime J/E (ported).
- `src/prime_layer/legendre_shift_2prime.py` — two-prime layer: M2 both shifts complete;
  S2 cross term F NOT yet implemented, RAISES when c3!=0 (no silent F=0 — C11 discipline).
- `scripts/single_prime_limit_check.py` + `tests/prime_layer/test_single_prime_limit.py`
  — the S2 acceptance gate.
- `docs/PROOF_CONSTITUTION.md` (ported from weil-first — the epistemic discipline).

## 3. The next concrete steps (in order)
1. **S2 — port shared machinery** — ✅ DONE. Ported verbatim, tests pass,
   single-prime-limit self-check PASSES (max|dC|=0).
2. **S3 — two-shift prime layer**: implement the cross-prime term
   $F_{ij}(\tau_2,\tau_3) = \langle C_{\tau_3,1}P_j, C_{\tau_2,1}P_i\rangle$ in
   `legendre_shift_2prime.py` (currently RAISES for c3!=0). It is the exchange
   Gram between the two shifts — derive the exact $\mathbb{Q}[\tau_2,\tau_3]$
   integral, do NOT stub it to 0 (that is the C11 bug the prototype carried).
   Ground-truth check: c3=0 already reproduces weil-first (S2 gate); add a
   parity/symmetry check on F and a `tau_2=tau_3` consistency check (F must equal
   E there).
3. **S4 — profile prime influence FIRST** (the key steer): a cheap mutation-style
   probe of per-sector, per-prime, per-cross-term influence on the Schur pivot.
   Do NOT budget compute symmetrically. Spend hard compute where the pivot moves.
4. **S5 — schema + domain + first pilot cert** once S3/S4 say the method is sound.

## 4. Trust discipline (inherited, do not relax)
- Full four-term $S^{(0)}$; two-prime $S^{(2)}$ with ALL cross terms. Omitting any
  is the C11 bug class.
- Verdicts require certify grade (Arb interval), never pilot/float.
- A number supporting an exciting narrative demands MORE scrutiny (PROOF_CONSTITUTION D1).
- Diff artifacts element-wise before narrating a disagreement (D2).
- Process defect ≠ wrong conclusion (D3).

## 5. Environment
- proofctl at `~/github/proofctl` (modifiable + publishable — if a new kernel
  blind spot appears, fix upstream first, then continue here). `~/bin/proofctl` deployed.
- github push may need a proxy via `${HTTPS_PROXY:-}` (local environment).
- python-flint (Arb), numpy; tectonic for LaTeX.
- Long tasks: `~/.local/bin/run_and_wait.sh -t <sec> -- <cmd>`; no bare `&`.
- Second-window integrals are heavier than first window ($c_L\approx1.82$, two
  shifts) — expect long certify runs; checkpoint/resume is mandatory (CLAUDE.md).

## 6. One sentence to the next maintainer
This repo's value is that it is a *live* new-math host: it will produce real bugs
that neither the math nor proofctl's designers can invent by self-testing —
capture each one honestly (it advances both the mathematics and the verifier),
and never let the second window's novelty tempt a claim past $L=\log2$ or toward RH.
