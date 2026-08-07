# Handoff — weil-second-prime (scaffold, 2026-08-08)

## 1. One-line status
Scaffold only. No math computed, no certificate, no proofctl graph yet. This repo
is the second-prime-window ($\tfrac12\log3 < L < \log2$) sibling of
weil-first-prime, created to be (a) a new mathematics host and (b) a fresh live
host that co-evolves proofctl. **Scope: finite-scale positivity, does NOT imply RH.**

## 2. What exists now
- `README.md`, `CLAUDE.md`, `PLAN.md`, this file — identity, rules, plan.
- Empty skeleton dirs: `src/{archimedean,prime_layer,assemble}`, `checker/`,
  `schemas/`, `domains/fp_second/contracts/`, `docs/`, `tests/`, `scripts/`,
  `pilots/`, `lean4/`.
- `docs/PROOF_CONSTITUTION.md` (ported from weil-first — the epistemic discipline).

## 3. The next concrete steps (in order)
1. **S2 — port shared machinery** from `../weil-first-prime`:
   `src/archimedean/{integrator_a,integrator_b,interval,ldlt,log_moments,kernel,bernstein}.py`.
   These are the trusted, tested archimedean primitives (P0 bugs already fixed).
   Port verbatim, then run their tests. Do NOT rewrite.
2. **S3 — two-shift prime layer**: move `src/prime_layer/legendre_shift_2prime.py`
   (currently a prototype in weil-first) here and complete the cross-prime $S^{(2)}$
   terms. Ground-truth check: single-prime limit ($c_3=0$) must reproduce
   weil-first's `legendre_shift.py`.
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
