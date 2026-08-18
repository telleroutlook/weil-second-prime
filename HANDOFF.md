# Handoff — weil-second-prime (2026-08-18)

## 1. One-line status

The k=18..28 odd-sector submatrix chain is complete. The apparent positive
minimum at k=25..27 was an IR-branch window, not a global positivity trend:
k=28 introduces a distinct negative frontier mode with
`lambda0(eta=1) = -1.8181565e-1`. The N=25 Arb attempt is running from its
resumable checkpoint; N=27 is the next queued finite-scale check. No RH claim,
no global positivity claim, and no extrapolation beyond `L < log 2`.

## 2. Newly completed work

- k=28 chain output:
  - `pilots/submatrix_k28.json`
  - `pilots/submatrix_k28_analysis.json`
  - `pilots/submatrix_row27.npz`
- Spectrum at `eta=1`:
  - frontier mode: `lambda0 = -1.8181565e-1`
  - IR mode: `lambda1 = +1.8655017e-6`
  - continuing UV mode: `lambda2 = +1.2972974e-2`
  - every scanned `eta in {0.5,...,2.0}` has one negative eigenvalue.
- Float Rayleigh split for the k=28 frontier eigenvector:
  - low block (`P1..P53`): `+2.0733604`
  - `P53`-`P55` cross: `-4.4755881`
  - `P55` diagonal: `+2.2204120`
  - total: `-0.1818157`
  - largest weights: `P53` (`0.8896`) and `P55` (`0.0956`).
- Corrected branch tracking:
  - continuing UV fit, k=13..28: `B0_UV=+0.02361`, RMS `1.33e-3`,
    bootstrap CI `[+0.02028,+0.02917]`, seed `20260818`;
  - IR 4-point fit: `B0_IR=+8.4655e-6`, but the bootstrap upper endpoint equals
    the imposed `1e-5` bound, so this is upper-censored discovery data;
  - the negative frontier mode is tracked separately and globally refutes the
    pilot-level `B0>0` narrative at the current finite scale.
- `scripts/fit_b0.py` fixes:
  - deterministic bootstrap seed and explicit discovery-grade label;
  - scale-normalized exponential fitting for IR-size data;
  - separate UV, IR, and frontier branches (k=28 `lambda0` is no longer miscoded
    as the continuing UV branch).
- Documentation synchronized in `docs/method_boundary.md`,
  `docs/SECOND_WINDOW_PAPER_DRAFT.md`, `PLAN.md`, and this handoff.

## 3. Running work

N=25 Arb attempt:

```text
PID 8606: run_and_wait wrapper
PID 8609: checker.fp_second.certify_fp_second
log: /tmp/cert_fp_second_N25_eta1_p512.log
checkpoint: pilots/cert_fp_second_N25_eta1_p512.ckpt.json
```

Command being resumed:

```bash
env PYTHONPATH=. python3 -m checker.fp_second.certify_fp_second \
  --L 56 100 --sector odd --N 25 --prec 512 \
  --eta 1/1 --no-bernstein \
  --out pilots/cert_fp_second_N25_eta1_p512.json \
  --resume
```

The checker writes its checkpoint after every completed matrix pair and supports
KeyboardInterrupt plus `--resume`. N=27 should be launched only after N=25 ends;
it cannot override the k=28 negative pilot by itself.

## 4. Epistemic status

- k=28 chain values are discovery-tier float pilot data. In particular,
  `submatrix_chain.py` calls `integrate_M_K(..., skip_remainder=True)` for speed.
  Do not present these values as Arb certificates.
- A positive N=25/N=27 result would be finite-scale evidence only. The k=28
  negative mode remains the controlling pilot obstacle.
- The direct min-pivot judge remains mandatory for any certificate. Float
  eigenvalues are for mode diagnosis only.
- Richardson mode (`--no-bernstein`) is Arb outward-rounded arithmetic with the
  documented empirical-remainder gap, not a formal Bernstein analytic certificate.
- The 4-point IR asymptote is upper-censored; do not quote `8.4655e-6` as a
  sharp bound.

## 5. Next concrete steps

1. Let N=25 finish or interrupt it cleanly; inspect both the JSON and checker
   exit status. If interrupted, resume from the pair checkpoint.
2. Launch N=27 with the same eta/precision/Richardson settings and resume
   discipline.
3. Design a certify-grade check of the k=28 frontier mode (preferably a focused
   `P53`-`P55` coupling certificate if the mathematics permits one), rather than
   relying on the float Rayleigh decomposition.
4. If the N=25/N=27 interval widths are decisive, record the result honestly as
   finite-scale indeterminate/positive/negative; do not broaden the conclusion.
5. Keep `fit_b0.py` branch tracking and scale normalization under regression
   tests before changing the model again.

## 6. Verification snapshot

- Targeted regression: `python3 -m pytest tests/test_fit_b0.py -q` — 2 passed.
- Full suite: `python3 -m pytest tests/ -x` — **142 passed in 50.52s**.

## 7. One sentence for the next maintainer

Resist the attractive k=25..27 window: k=28 exposed a new negative high-degree
coupling mode, so the honest task is now to certify that frontier, not to narrate
positivity from the IR branch.
