# Handoff — weil-second-prime (2026-08-18, raw-GL correction)

## 1. One-line status

The apparent k=28 negative frontier mode was a raw-GL8 truncation artifact, not
a spectral result. A focused Richardson-Arb audit of the integer witness
`v = 3*P53 + P55` gives `v^T C v ∈ [+0.2292033466, +0.3953780542]`, entirely
positive. This refutes that particular negative witness but does not prove the
full N=28 matrix positive definite. The legacy k=18..28 chain and all UV/IR/B0
fits built from it are quarantined. N=25 remains in progress; N=27 is queued.

## 2. Root cause

The legacy `submatrix_chain.py` called:

```python
integrate_M_K(..., use_bernstein=False, skip_remainder=True)
```

That returns a raw GL8 Arb ball with no Richardson or Bernstein truncation
remainder. Its center is not sign-safe at high degree.

Focused Arb comparison:

| M0 entry | Legacy raw-GL center | Richardson-Arb center ± radius |
|---|---:|---:|
| `P53,P53` | `0.050151` | `0.012955 ± 0.000144` |
| `P53,P55` | `-0.138160` | `0.009174 ± 0.000127` |
| `P55,P55` | `0.422735` | `0.012488 ± 0.000073` |

The raw `P55,P55` value is more than 33 times the focused center. The apparent
k=28 mode and its Rayleigh decomposition were artifacts of this error.

## 3. Focused k=28 audit

Artifacts:

- `scripts/certify_k28_frontier.py`
- `pilots/cert_fp_second_N28_frontier_eta1_p512.json`
- `pilots/cert_fp_second_N28_frontier_eta1_p512.frontier.ckpt.json`

The script computes only the 55 unique M0 entries needed for the `P53` and `P55`
columns plus three S0 entries. It writes a durable checkpoint after every unit.
Result:

```text
rayleigh_witness_status = nonnegative
v^T C v ∈ [+0.2292033466, +0.3953780542]
certified_not_positive_definite = false
```

Interpretation:

- It does **not** certify N=28 positive definite.
- It does refute the proposed `3*P53 + P55` negative witness.
- It uses Richardson GL-8/GL-4 remainder mode, with the documented empirical
  remainder boundary; it is not a formal Bernstein analytic certificate.

## 4. Corrections now in place

- `scripts/submatrix_chain.py` no longer requests `skip_remainder=True`.
- Corrected chain rows/spectra use a new namespace:
  - `pilots/submatrix_rich_row*.npz`
  - `pilots/submatrix_rich_k*.json`
- The script will not silently resume legacy raw-GL row files.
- `scripts/fit_b0.py`:
  - prefers `submatrix_rich_k*.json` when available;
  - refuses to mix corrected and legacy sequences;
  - prints an explicit quarantine warning while only legacy `submatrix_k*.json`
    files exist.
- Legacy UV/IR/frontier mode labels, zero crossings, B0 fits, and bootstrap
  probabilities are withdrawn from the active narrative.
- Documentation and PLAN now record a proofctl C13 candidate:
  **pilot-sign firewall** — raw-center outputs cannot support sign narratives
  unless their remainder mode/enclosure is explicit.

## 5. Running and queued work

N=25 Arb attempt:

```text
checker PID: 8609
log: /tmp/cert_fp_second_N25_eta1_p512.log
checkpoint: pilots/cert_fp_second_N25_eta1_p512.ckpt.json
progress at 2026-08-18 12:07 CST: 121/625 pairs
```

Command:

```bash
env PYTHONPATH=. python3 -m checker.fp_second.certify_fp_second \
  --L 56 100 --sector odd --N 25 --prec 512 \
  --eta 1/1 --no-bernstein \
  --out pilots/cert_fp_second_N25_eta1_p512.json \
  --resume
```

N=27 is queued by `/tmp/post_n27_after_n25.sh`. It starts only after N=25
produces its final JSON. If N=25 dies without a final JSON, the queue errors
rather than launching N=27 from unrelated state.

## 6. Next steps

1. Let N=25 finish; inspect the JSON and checker exit status, not a self-reported
   narrative.
2. Run N=27 only through the existing serial queue.
3. Do not restore the withdrawn k=18..28 mode/B0 claims.
4. If mode tracking is still needed, rerun the corrected chain into
   `submatrix_rich_*`; do not reuse legacy rows.
5. Propose C13 upstream in proofctl: every numeric artifact must carry a
   remainder mode, and raw-center values must be blocked from sign claims.

## 7. Verification snapshot

- Focused helper tests:
  `python3 -m pytest tests/test_certify_k28_frontier.py -q` — 3 passed
  (including replay from the completed interval checkpoint).
- Chain remainder regression:
  `python3 -m pytest tests/test_submatrix_chain_remainder.py -q` — 1 passed.
- Combined targeted run:
  `python3 -m pytest tests/test_submatrix_chain_remainder.py tests/test_certify_k28_frontier.py -q`
  — 3 passed.
- Full suite:
  `python3 -m pytest tests/ -x` — **146 passed in 56.03s**.
- The focused result was regenerated from its completed checkpoint after adding
  `rayleigh_witness_status=nonnegative`.

## 8. One sentence for the next maintainer

Do not narrate either positivity or a k=28 negative mode from the legacy chain:
the immediate task is the N=25/N=27 min-pivot results and, if needed, a fully
remainder-safe corrected chain.
