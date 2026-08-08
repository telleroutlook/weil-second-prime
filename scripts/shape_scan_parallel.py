"""
Parallel driver for the second-window lambda_min(L) SHAPE characterization.

Reuses eig_scan_second_window.scan_L (pilot grade, float center) but runs one
(sector, L) build per worker process across cores, since each build is heavy
(~9 min single-threaded) and independent. Incremental-durable: each point is
written to its own JSON shard as it finishes, then merged; a crash loses at most
one in-flight point (CLAUDE.md long-task rule). --resume skips shards already on
disk.

GRADE: pilot (float center). This locates the infimum; certify-grade anchoring
at the winning L is done separately (authoritative_eig_check.py). No window-wide
assertion is made from these finite samples (PROOF_CONSTITUTION A3).
"""
from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from scripts.eig_scan_second_window import scan_L
from src.prime_layer.legendre_shift_2prime import window_check


def _one(num: int, den: int, sector: str, N: int, shard: str) -> dict:
    t0 = time.time()
    r = scan_L(num, den, sector, N)
    r["elapsed_s"] = round(time.time() - t0, 1)
    Path(shard).write_text(json.dumps(r, indent=2))
    return r


def main() -> int:
    ap = argparse.ArgumentParser(description="Parallel second-window shape scan")
    ap.add_argument("--sector", required=True, choices=["even", "odd"])
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--Ls", required=True, help="comma-separated num/den, e.g. 56000/100000,62000/100000")
    ap.add_argument("--out", required=True, help="merged output json")
    ap.add_argument("--shard-dir", default="pilots/shape_shards")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()

    shard_dir = Path(args.shard_dir)
    shard_dir.mkdir(parents=True, exist_ok=True)

    jobs = []
    for tok in args.Ls.split(","):
        num, den = (int(x) for x in tok.split("/"))
        L = num / den
        if not window_check(L):
            print(f"[skip] L={tok}={L:.6f} outside window", flush=True)
            continue
        shard = shard_dir / f"{args.sector}_N{args.N}_{num}_{den}.json"
        if args.resume and shard.exists():
            print(f"[skip] {tok} {args.sector} N={args.N} shard exists", flush=True)
            continue
        jobs.append((num, den, args.sector, args.N, str(shard)))

    print(f"[shape] launching {len(jobs)} builds on {args.workers} workers "
          f"(sector={args.sector} N={args.N})", flush=True)
    results = []
    if jobs:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(_one, *j): j for j in jobs}
            for fut in as_completed(futs):
                r = fut.result()
                results.append(r)
                print(f"[done] L={r['L_float']:.5f} {r['sector']} N={r['N']} "
                      f"eig_full={r['eig_full_min']:+.6f}(gap {r['eig_full_gap']:.3f}) "
                      f"eig_off={r['eig_off_min']:+.6f} eig_p1={r['eig_p1_min']:+.6f} "
                      f"({r['elapsed_s']}s)", flush=True)

    # Merge all shards for this sector/N into the requested out file
    merged = []
    for shard in sorted(shard_dir.glob(f"{args.sector}_N{args.N}_*.json")):
        merged.append(json.loads(shard.read_text()))
    merged.sort(key=lambda r: r["L_float"])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(merged, indent=2))
    print(f"\n[shape] merged {len(merged)} points -> {args.out}", flush=True)
    for r in merged:
        print(f"  L={r['L_float']:.5f} eig_full={r['eig_full_min']:+.6f} "
              f"gap={r['eig_full_gap']:.4f}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
