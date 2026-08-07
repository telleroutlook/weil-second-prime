# Proof Constitution — Non-Negotiable Discipline for Numerical Proof Work

Single authoritative source. CLAUDE.md points here; do not duplicate these rules
elsewhere. Every rule is bound to a concrete TRIGGER (when it must fire) and,
where possible, an automatable CHECK. Rules exist because each was VIOLATED at
real cost (see "Incident" per rule). A rule written but not enforced is theatre —
prefer an executable check over a markdown virtue.

Established 2026-08-07 after an all-night campaign that (a) discovered the FP-0.35
certificate was assembled with two compounding bugs and (b) killed seven
high-dimensional "isomorphism" mapping attempts. Both halves produced the rules
below.

---

## PART A — Numerical / Computation Discipline

### A1. Normalization-first
TRIGGER: before reusing any quantity (M_K, compute_J, V_matrix_entry, ...) to
assemble a new quantity.
RULE: independently verify its normalization/definition (factor, inner-product
convention) by cross-checking against direct numeric evaluation FIRST.
INCIDENT: compute_J carried a factor 2 vs the L2 shifted-overlap; M_K was plain
L2 (factor 1). Assembling S_cross without checking would have injected a
systematic factor error.
CHECK: a ratio test (assembled-primitive / independent-numeric) must be a clean
constant across >=5 sample indices before use.

### A2. Ground-truth reproduction before trust
TRIGGER: any new script or parallel implementation of an existing pipeline.
RULE: it must reproduce a KNOWN certified value before it may compute anything
new. If it cannot, STOP and debug the script — do not report its new numbers.
INCIDENT: a parallel numpy script printed "CERTIFIED lambda(0.42)>0" but could
not reproduce FP-0.35's +0.01494 ground truth (gave +0.001). The "certification"
was a self-issued lie.
CHECK: tests/test_certificate_integrity.py must reproduce every published cert
value; CI fails if any drifts.

### A3. Exploration precision never issues a verdict
TRIGGER: any number that will enter a conclusion, a doc, or a decision.
RULE: pilot/explorer precision (depth<=2, float) is for DIRECTION ONLY. Verdicts
require certify grade (depth 4/3 + Arb interval).
INCIDENT: depth=2 rendered L=0.42 min_eig as -0.0007; certify grade gave -0.022
(30x error). A "lean negative" verdict was nearly issued on the depth=2 value.

### A4. Boundary values require error bars
TRIGGER: a decision quantity within ~10x of zero (e.g. min_eig in +-0.01).
RULE: a point estimate is insufficient; require interval-arithmetic error bounds
before declaring sign.
INCIDENT: split min_eig = +0.0009 vs joint = -0.0006 — the 0.0015 gap is the
same order as the kmax truncation bias (9e-4). Neither can be called a verdict
without error control.

### A5. Sufficient != necessary
TRIGGER: any Schur / LDL / positivity criterion returning a negative value.
RULE: a failing sufficient criterion (min_pivot<0) proves ONLY that THIS
certificate failed, NOT that the underlying claim is false. State it that way.
INCIDENT: min_pivot=-0.043 at N=8 means the N=8 certificate fails; since RH is
verified to great height, lambda(7/20)>0 almost surely — the certificate, not
the truth, is what failed.

### A6. Non-monotone / self-contradictory series = numerical alarm
TRIGGER: a sequence over N (or any refinement parameter) that is non-monotone or
internally contradictory.
RULE: treat it as instability, NOT physics; escalate precision before
interpreting.
INCIDENT: min_eig(N): N=16 -0.0007, N=20 -0.022 (fast path) — contradiction that
exposed depth=2 unreliability.

### A7. Batch-cache expensive primitives before assembly
TRIGGER: assembling a quantity that loops over many calls to an expensive
integral (M_K, S_KK, ...).
RULE: pre-compute and cache ALL needed primitives once, THEN assemble in memory.
Never embed uncached expensive integrals inside an assembly loop.
INCIDENT: S_cross embedded live integrate_M_K(k,i) for k up to 60 (cache misses),
turning an estimated "few minutes" into 30+ minutes.

### A8. Theory-bound sanity gate on every new estimator
TRIGGER: computing a quantity that has a KNOWN theoretical ordering vs an existing
one.
RULE: check the ordering; a violation means the new estimator is wrong, not that
a discovery was made.
INCIDENT: joint Schur MUST be >= best split (no Young loss). It came out -0.0006
< split's +0.0009 — instantly proving the S_cross assembly was buggy, preventing
a false "joint also fails" conclusion.

---

## PART B — Proof / Epistemology Discipline

### B1. Conservation of difficulty
Any reframing equivalent to the original problem reduces difficulty by zero — it
only changes clothes. TEST: to which field with independent proven theorems does
it transfer the difficulty? No answer => decoration.

### B2. Isomorphism five-questions (for any "map to field X" proposal)
(1) non-equivalent to the target inequality; (2) target field has an
RH-independent PROVEN theorem to absorb it; (3) can express uniformity/full
spectrum; (4) non-circular (no RH assumed); (5) difficulty visibly transferred.
Any "no/unknown/circular" => rejected.

### B3. Counting quantities cannot capture the spectral sign
Rank, dimension, prime count, single-direction overlap all failed to predict the
sign of lambda(L) (five falsifications). Uniform-in-L needs the FULL spectrum; no
structural shortcut.

### B4. Finer data overrules coarser
TRIGGER: any "pattern" seen on a coarse grid.
RULE: re-verify on a finer grid before believing it. "No narrative survives
contact with a finer grid."
INCIDENT: node-transition == sign-flip synchrony (coarse 3 points) was falsified
by a fine grid (node onset L~0.40 leads sign flip L~0.43-0.45).

### B5. The executioner audits their own blade
TRIGGER: issuing a "kill"/verdict on any direction.
RULE: re-check your own kill argument (scale choice, basis-space vs real-space,
data precision) before it lands.
INCIDENT: the dir1/5 "nothing grows" kill was too fast (prime norm DID grow +84%,
reserve +16%); dir7's IPR was in coefficient space, its trend needed real-space
reinterpretation.

### B6. Honest negative result > false breakthrough
Finding a bug, or proving a path impossible, is real progress. Manufacturing a
"CERTIFIED" that collapses on inspection is damage. When in doubt, downgrade and
verify.

---

## Meta-rule: rules must have teeth
CLAUDE.md already contained "no duplicate logic" and fail-closed intent; both
were violated tonight under the excitement of a apparent breakthrough. Therefore:
prefer to encode each rule as an automated check (CI test) over trusting
self-discipline. A markdown rule is a reminder; a failing test is a wall.

---

## PART C — Diagnosis Discipline (added 2026-08-07, from a human correction)

### C1. Read the code before you compute
TRIGGER: two implementations / scripts disagree on a number, OR you are about to
run an expensive recomputation to "find out which is right".
RULE: FIRST diff the two code paths and inspect the key constants. Most
disagreements are located by reading source + checking a handful of scalars
(seconds), NOT by re-running the whole pipeline (minutes/hours). Only compute
when code review genuinely cannot resolve it.
INCIDENT: recompute_schur (+0.0088) vs o1b_gate (-0.043) for the same L=7/20
even sector. A 40-min full recompute was launched to "decide". Code review
instead pinned it in seconds: R_eta identical, F identical, kappa identical to
4e-9 — so the difference had to be in the LAST step (min-pivot implementation),
not the matrices. The 40-min job was redundant and was killed.

### C2. Judge implementation credibility by provenance, not by which answer you like
TRIGGER: conflicting numbers from two implementations.
RULE: rank the implementations by provenance BEFORE trusting either number:
production code with test coverage + verified numeric methods (e.g. mpmath LDL,
realtime Arb-certified constants) outranks a freshly-written, untested script
with hardcoded constants and a hand-rolled numeric kernel (e.g. float64 LDL with
no pivoting). When they conflict, the higher-provenance result is the prior.
INCIDENT: o1b_gate uses mpmath LDL + realtime kappa() (124 tests, proofctl-
integrated); recompute_schur used a hand-written float64 LDL with no symmetric
pivoting + hardcoded kappa. The hand-rolled kernel is the natural suspect for a
sign flip on a near-singular matrix — identifiable by reading it, without a
verdict-grade run.

### C3. Isolate the differing stage; re-run only that stage
TRIGGER: code review has narrowed a disagreement to one stage of a pipeline.
RULE: test ONLY that stage (e.g. feed both min-pivot implementations a known
ill-conditioned matrix), not the whole pipeline. Do not recompute upstream
stages that code review already proved identical.

---

## PART D — Narrative Resistance (added 2026-08-07, hardest-won)

### D1. A number supporting an exciting narrative demands MORE scrutiny, not less
TRIGGER: a computed value would confirm a dramatic story — "found a critical
bug", "the proof fails", "a breakthrough", "the tool caught a fraud".
RULE: raise the verification bar precisely when the number is narratively
attractive. Excitement is a bias multiplier, not evidence.
INCIDENT: -0.043 ("FP-0.35 fails!") was accepted three times because it fit the
story "we found a fatal bug". It was an ad-hoc-script error; the true value was
+0.0087 all along. Downstream conclusions had to be retracted.

### D2. Compare artifacts element-wise before narrating a disagreement
TRIGGER: two implementations give different scalar results.
RULE: dump and diff the full intermediate ARTIFACTS (matrices) element-wise to
locate the exact divergence BEFORE explaining "why". A scalar mismatch plus a
plausible story is not a diagnosis.
INCIDENT: +0.0087 vs -0.043 was "explained" twice with confident wrong stories.
The element-wise diff proved the matrices AGREE (max|C_A-C_B|=4e-3), both give
+0.0087; -0.043 never existed in any correct computation.

### D3. Distinguish "process defect" from "wrong conclusion"
TRIGGER: discovering a bug in how a result was produced.
RULE: a defective process (copy-generation, omitted term, hard-coded verdict)
does NOT automatically mean the conclusion is wrong. State precisely which is
broken.
INCIDENT: the FP-0.35 certificate had a copy generator AND a 16x-inflated
eigenvalue, yet the SIGN (lambda(7/20)>0) is correct. "Bug in the certificate"
was true; "the proof fails" was false.

---

## PART E — Mutation Testing Nuance (added 2026-08-07)

### E1. Kill criterion must match the mutated term's real magnitude of influence
TRIGGER: designing a mutation catalog to prove a checker is sensitive to each
asserted term (proofctl C11).
RULE: a mutant that zeros a term whose true influence is tiny will NOT flip the
verdict — this does not mean the checker is blind; that term genuinely barely
matters. Classify mutants: JUDGE-sensitivity (sign flip, zero a dominant term,
swap min-pivot->min-eig) MUST be killed; PRECISION mutants (zero a ~1e-4-
influence term) judged by a RELATIVE tolerance, or replaced by a large-factor
scaling mutant that DOES move the pivot.
INCIDENT: L=7/20 EVEN sector, zeroing S2 moved min_pivot +0.0088 -> +0.0087
(~1e-4). A naive "100% kill" rule would wrongly flag the checker insensitive to
S2. S2 has negligible even-sector influence; a large-factor S2 scaling mutant is
the correct sensitivity probe.

### E2. Prefer self-contained recomputation over certificate-reading checkers
TRIGGER: writing or auditing a checker.
RULE: a checker that RECOMPUTES the claim from primitives is structurally immune
to copy-cert and self-attestation attacks (unlike one that reads+validates a
supplied certificate). Recompute-style checkers still need mutation coverage for
OMITTED-TERM bugs, but cannot be fooled by a doctored certificate.
INCIDENT: recompute_schur.py recomputes C from scratch; the retired
check_fp035.py hard-coded all obligations True while its generator merely copied
the certificate. Recompute is the fix for both.
