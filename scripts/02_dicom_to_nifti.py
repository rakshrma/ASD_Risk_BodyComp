#!/usr/bin/env python3
"""
02_dicom_to_nifti.py
====================
Run dcm2niix on every <accession>/ subdirectory under a DICOM source root,
producing matching <accession>/ subdirectories of NIfTI files.

Default folders (relative to the repo root):
  source : data/input/dicom_decompressed/   (output of script 01)
           Falls back to data/input/dicom/  if the decompressed folder is empty.
  dest   : data/input/nifti/

Usage:
  python scripts/02_dicom_to_nifti.py
  python scripts/02_dicom_to_nifti.py --source <dir> --dest <dir>

Requirements:
  dcm2niix on PATH (https://github.com/rordenlab/dcm2niix/releases)
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

REPO_ROOT          = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_DEC = REPO_ROOT / "data" / "input" / "dicom_decompressed"
DEFAULT_SOURCE_RAW = REPO_ROOT / "data" / "input" / "dicom"
DEFAULT_DEST       = REPO_ROOT / "data" / "input" / "nifti"


def _resolve_default_source() -> Path:
    """Prefer dicom_decompressed if it has accessions; else fall back to dicom/."""
    if DEFAULT_SOURCE_DEC.is_dir() and any(DEFAULT_SOURCE_DEC.iterdir()):
        return DEFAULT_SOURCE_DEC
    return DEFAULT_SOURCE_RAW


def check_dcm2niix() -> None:
    if shutil.which("dcm2niix") is None:
        sys.exit("[ERROR] dcm2niix not found on PATH. "
                 "Install: https://github.com/rordenlab/dcm2niix/releases")


def run(source: Path, dest: Path, dry_run: bool, verbose: bool) -> None:
    if not source.exists():
        sys.exit(f"[ERROR] Source directory does not exist: {source}")

    subdirs = sorted(d for d in source.iterdir() if d.is_dir())
    if not subdirs:
        print(f"[WARNING] No subdirectories found in: {source}")
        return

    print(f"\n{'[DRY RUN] ' if dry_run else ''}DCM2NIIX")
    print(f"  Source : {source}")
    print(f"  Dest   : {dest}")
    print(f"  Jobs   : {len(subdirs)}")
    print(f"  Time   : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 60)

    success, failed = [], []
    for i, subdir in enumerate(subdirs, 1):
        dest_subdir = dest / subdir.name
        cmd = ["dcm2niix",
               "-i", "y", "-z", "y", "-a", "y",
               "-f", "%i_%f_%s",
               "-o", str(dest_subdir),
               str(subdir)]

        print(f"\n[{i}/{len(subdirs)}] {subdir.name}")
        print(f"  CMD: {' '.join(cmd)}")
        if dry_run:
            success.append(subdir.name)
            continue

        dest_subdir.mkdir(parents=True, exist_ok=True)
        try:
            r = subprocess.run(cmd, capture_output=not verbose, text=True)
            if verbose and r.stdout:
                print(r.stdout)
            if r.stderr:
                print(f"  [STDERR] {r.stderr.strip()}")
            if r.returncode == 0:
                print("  [OK]")
                success.append(subdir.name)
            else:
                print(f"  [FAILED] exit {r.returncode}")
                failed.append(subdir.name)
        except Exception as e:
            print(f"  [ERROR] {e}")
            failed.append(subdir.name)

    print("\n" + "=" * 60)
    print(f"Summary  total={len(subdirs)}  ok={len(success)}  failed={len(failed)}")
    for n in failed:
        print(f"  - {n}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", type=Path, default=None,
                   help="DICOM root (default: data/input/dicom_decompressed/ "
                        "or data/input/dicom/ if the former is empty)")
    p.add_argument("--dest",   type=Path, default=DEFAULT_DEST,
                   help="NIfTI output root (default: data/input/nifti/)")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()

    if args.source is None:
        args.source = _resolve_default_source()
    if not args.dry_run:
        check_dcm2niix()

    run(args.source.resolve(), args.dest.resolve(), args.dry_run, args.verbose)


if __name__ == "__main__":
    main()
