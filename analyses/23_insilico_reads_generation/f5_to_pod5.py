#!/usr/bin/env python3
"""
fast5 -> pod5 converter that does NOT silently drop reads.

Why this exists
---------------
`pod5 convert fast5` calls get_read_from_fast5(), which wraps the h5py lookup
in `try: ... except KeyError:` and, on failure, issues a warning and returns
None.  The caller then does `if f5_read is None: continue`.  The read is
dropped and conversion carries on.  That path raises no exception, so
`--strict` does NOT catch it and the process still exits 0.

This script instead:
  * enumerates the reads in each fast5 up front (that is the expected count),
  * retries a failed read on a freshly-opened file handle, several times,
  * refuses to keep a pod5 that has fewer reads than its input fast5,
  * exits non-zero listing every file that could not be fully converted.

Usage
-----
    python3 f5_to_pod5.py --input  DIR_OF_MULTI_FAST5 \\
                          --output DIR_FOR_POD5 \\
                          --jobs 20 [--attempts 5] [--scratch $TMPDIR]

  --scratch  copy each fast5 to local disk before reading, convert there, then
             move the pod5 to --output.  Use if the shared filesystem is
             implicated.

Re-running skips inputs whose pod5 already exists and has the right count, so
it is resumable.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import h5py
import pod5 as p5
from pod5.tools.pod5_convert_from_fast5 import convert_fast5_read
from pod5.signal_tools import DEFAULT_SIGNAL_CHUNK_SIZE


def read_group_names(path: Path) -> list[str]:
    with h5py.File(str(path), "r") as f5:
        return [k for k in f5.keys() if k.startswith("read_")]


def convert_one(
    src: Path, dest: Path, attempts: int, chunk: int, scratch: str | None,
    allow_missing: int = 0
) -> dict:
    """Convert one multi-read fast5 to one pod5. Returns a result record."""
    res = {"src": str(src), "expected": 0, "written": 0, "kept": 0,
           "retried": 0, "tolerated": 0, "missing": [], "error": None}
    staged = None
    try:
        expected_names = read_group_names(src)
        res["expected"] = len(expected_names)

        # resume: a complete output already present is left alone
        if dest.exists():
            try:
                with p5.Reader(dest) as r:
                    if r.num_reads >= res["expected"] - allow_missing:
                        res["written"] = res["kept"] = r.num_reads
                        res["tolerated"] = res["expected"] - r.num_reads
                        res["skipped"] = True
                        return res
            except Exception:
                pass
            dest.unlink()

        read_path = src
        if scratch:
            staged = Path(scratch) / f"{os.getpid()}_{src.name}"
            shutil.copy2(src, staged)
            read_path = staged

        pending = list(expected_names)
        tmp_out = dest.with_suffix(".pod5.partial")
        if tmp_out.exists():
            tmp_out.unlink()

        with p5.Writer(tmp_out) as writer:
            cache: dict = {}
            for attempt in range(1, attempts + 1):
                if not pending:
                    break
                if attempt > 1:
                    res["retried"] += len(pending)
                failed: list[str] = []
                # a FRESH handle each attempt - stale HDF5 cache state is the
                # most likely reason a retry succeeds where the first try failed
                with h5py.File(str(read_path), "r") as f5:
                    for name in pending:
                        try:
                            grp = f5[name]
                            writer.add_read(
                                convert_fast5_read(grp, cache, chunk)
                            )
                            res["written"] += 1
                        except Exception:
                            failed.append(name)
                pending = failed

        res["missing"] = pending
        if len(pending) > allow_missing:
            # too damaged to keep: discard so a re-run retries the whole file
            tmp_out.unlink(missing_ok=True)
            res["error"] = (f"{len(pending)} reads unrecoverable after "
                            f"{attempts} attempts (limit {allow_missing})")
        else:
            tmp_out.replace(dest)
            res["kept"] = res["written"]
            res["tolerated"] = len(pending)

    except Exception:
        res["error"] = traceback.format_exc(limit=3)
    finally:
        if staged is not None:
            staged.unlink(missing_ok=True)
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--attempts", type=int, default=5)
    ap.add_argument("--scratch", default=None,
                    help="local dir to stage each fast5 through (e.g. $TMPDIR)")
    ap.add_argument("--allow-missing", type=int, default=0, metavar="N",
                    help="keep a file's pod5 even if up to N of its reads are "
                         "permanently unreadable (default 0 = discard the file)")
    ap.add_argument("--signal-chunk-size", type=int,
                    default=DEFAULT_SIGNAL_CHUNK_SIZE)
    args = ap.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    if args.scratch:
        Path(args.scratch).mkdir(parents=True, exist_ok=True)

    srcs = sorted(args.input.glob("*.fast5"))
    if not srcs:
        print(f"no *.fast5 under {args.input}", file=sys.stderr)
        return 2
    print(f"{len(srcs)} fast5 files, {args.jobs} workers, "
          f"{args.attempts} attempts/read", flush=True)

    tot_exp = tot_kept = tot_retry = tot_tol = 0
    bad: list[dict] = []
    done = 0

    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futs = {
            pool.submit(convert_one, s, args.output / (s.stem + ".pod5"),
                        args.attempts, args.signal_chunk_size, args.scratch,
                        args.allow_missing): s
            for s in srcs
        }
        for fut in as_completed(futs):
            r = fut.result()
            done += 1
            tot_exp += r["expected"]
            tot_kept += r["kept"]
            tot_retry += r["retried"]
            tot_tol += r["tolerated"]
            if r["error"]:
                bad.append(r)
                print(f"  FAIL {Path(r['src']).name}: {r['error']}", flush=True)
            elif r["tolerated"]:
                print(f"  ok*  {Path(r['src']).name}: {r['kept']} reads kept, "
                      f"{r['tolerated']} unreadable (tolerated)", flush=True)
            if done % 500 == 0:
                print(f"  ... {done}/{len(srcs)} files, {tot_kept:,} reads",
                      flush=True)

    lost_discarded = sum(r["expected"] for r in bad)
    print()
    print(f"files                  : {len(srcs)}  ({len(bad)} discarded)")
    print(f"reads expected         : {tot_exp:,}")
    print(f"reads ON DISK          : {tot_kept:,}")
    print(f"reads recovered by retry: {tot_retry:,}")
    print(f"reads unreadable, tolerated : {tot_tol:,}")
    if lost_discarded:
        print(f"reads lost to discarded files: {lost_discarded:,}")
    if tot_exp:
        print(f"completeness           : {100.0 * tot_kept / tot_exp:.4f}%")

    if bad:
        print(f"\nINCOMPLETE - {len(bad)} file(s) discarded. Re-run to retry them "
              f"(add --scratch $TMPDIR, or --allow-missing 1 to accept a few "
              f"unreadable reads).", file=sys.stderr)
        return 1
    if tot_tol:
        print(f"\nDone with {tot_tol} unreadable read(s) tolerated "
              f"({100.0 * tot_tol / tot_exp:.6f}% of the dataset).")
        return 0
    print("\nComplete: every read in every input fast5 is present in the pod5 output.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
