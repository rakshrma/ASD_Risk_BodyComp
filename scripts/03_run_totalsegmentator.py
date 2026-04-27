#!/usr/bin/env python3
"""
03_run_totalsegmentator.py
==========================
Run TotalSegmentator (`-ta total` and `-ta tissue_4_types`) on every .nii.gz
file under a NIfTI source root, preserving the <accession>/ subdirectory layout
in the output. Existing outputs are skipped.

Default folders:
  source : data/input/nifti/
  dest   : data/output/totalsegmentator/

Outputs per file:
  <dest>/<accession>/<stem>_ts_total.nii.gz
  <dest>/<accession>/<stem>_ts_tissue.nii.gz

Usage:
  python scripts/03_run_totalsegmentator.py
  python scripts/03_run_totalsegmentator.py --source <dir> --dest <dir> [--skip-total|--skip-tissue]

Requirements:
  pip install TotalSegmentator        (provides the TotalSegmentator CLI)
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import List

REPO_ROOT      = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = REPO_ROOT / "data" / "input" / "nifti"
DEFAULT_DEST   = REPO_ROOT / "data" / "output" / "totalsegmentator"


def check_totalsegmentator() -> None:
    if shutil.which("TotalSegmentator") is None:
        sys.exit("[ERROR] TotalSegmentator not found on PATH. "
                 "Install: pip install TotalSegmentator")


def find_nii(source: Path) -> List[Path]:
    return sorted(source.rglob("*.nii.gz"))


def run_task(label: str, cmd: List[str], output_file: Path,
             dry_run: bool, verbose: bool) -> str:
    if output_file.exists():
        print(f"    [{label}] exists; skipping.")
        return "skipped"
    print(f"    [{label}] {' '.join(cmd)}")
    if dry_run:
        return "ok"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        r = subprocess.run(cmd, capture_output=not verbose, text=True)
        if verbose and r.stdout:
            print(r.stdout)
        if r.stderr and verbose:
            print(f"    [STDERR] {r.stderr.strip()}")
        if r.returncode == 0:
            print(f"    [{label}] ok")
            return "ok"
        print(f"    [{label}] FAILED (exit {r.returncode})")
        if not verbose and r.stdout:
            print(f"    [STDOUT] {r.stdout.strip()}")
        return "failed"
    except Exception as e:
        print(f"    [{label}] ERROR: {e}")
        return "failed"


def run(source: Path, dest: Path, dry_run: bool, verbose: bool,
        skip_total: bool, skip_tissue: bool) -> None:
    if not source.exists():
        sys.exit(f"[ERROR] Source directory does not exist: {source}")

    files = find_nii(source)
    if not files:
        print(f"[WARNING] No .nii.gz files under {source}")
        return

    print(f"\n{'[DRY RUN] ' if dry_run else ''}TotalSegmentator")
    print(f"  Source : {source}")
    print(f"  Dest   : {dest}")
    print(f"  Files  : {len(files)}")
    print(f"  Tasks  : {'total ' if not skip_total else ''}"
          f"{'tissue_4_types' if not skip_tissue else ''}")
    print(f"  Time   : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 60)

    stats = {"ok": 0, "skipped": 0, "failed": 0}
    failed_files = []

    for i, nii in enumerate(files, 1):
        rel       = nii.relative_to(source)
        base_name = nii.name.replace(".nii.gz", "")
        out_dir   = dest / rel.parent

        print(f"\n[{i}/{len(files)}] {rel}")
        file_failed = False

        if not skip_total:
            out_total = out_dir / f"{base_name}_ts_total.nii.gz"
            cmd = ["TotalSegmentator", "-ml",
                   "-i", str(nii), "-o", str(out_total),
                   "-ta", "total"]
            r = run_task("total", cmd, out_total, dry_run, verbose)
            stats[r] += 1
            file_failed |= (r == "failed")

        if not skip_tissue:
            out_tissue = out_dir / f"{base_name}_ts_tissue.nii.gz"
            cmd = ["TotalSegmentator", "-ml",
                   "-i", str(nii), "-o", str(out_tissue),
                   "-ta", "tissue_4_types"]
            r = run_task("tissue_4_types", cmd, out_tissue, dry_run, verbose)
            stats[r] += 1
            file_failed |= (r == "failed")

        if file_failed:
            failed_files.append(str(rel))

    total_runs = stats["ok"] + stats["skipped"] + stats["failed"]
    print("\n" + "=" * 60)
    print(f"Summary  task-runs={total_runs}  ok={stats['ok']}  "
          f"skipped={stats['skipped']}  failed={stats['failed']}")
    for f in failed_files:
        print(f"  - {f}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", type=Path, default=DEFAULT_SOURCE,
                   help="NIfTI source root (default: data/input/nifti/)")
    p.add_argument("--dest",   type=Path, default=DEFAULT_DEST,
                   help="Segmentation output root "
                        "(default: data/output/totalsegmentator/)")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--skip-total",  action="store_true")
    p.add_argument("--skip-tissue", action="store_true")
    args = p.parse_args()

    if not args.dry_run:
        check_totalsegmentator()

    run(args.source.resolve(), args.dest.resolve(),
        args.dry_run, args.verbose, args.skip_total, args.skip_tissue)


if __name__ == "__main__":
    main()
