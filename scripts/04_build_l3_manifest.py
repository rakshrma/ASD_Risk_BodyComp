#!/usr/bin/env python3
"""
04_build_l3_manifest.py
=======================
Walk the TotalSegmentator output directory, find every CT case, locate the
L3 vertebral centroid, and write a manifest CSV consumed by scripts 5/6.

Output CSV columns:
  l3_slice_index   axial slice index of the L3 centroid (label 29 in *_ts_total)
  nifti_path       path to original CT NIfTI
  ts_total_path    path to *_ts_total.nii.gz
  ts_tissue_path   path to *_ts_tissue.nii.gz
  accession        accession ID (folder name under both roots)

Default folders:
  ts_root    : data/output/totalsegmentator/
  nifti_root : data/input/nifti/
  out_csv    : data/output/l3_manifest.csv

Usage:
  python scripts/04_build_l3_manifest.py
  python scripts/04_build_l3_manifest.py --ts-root <dir> --nifti-root <dir> --out <csv>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Make src/ importable when run as a script
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from asd_bodycomp.manifest import build_manifest

DEFAULT_TS_ROOT    = REPO_ROOT / "data" / "output" / "totalsegmentator"
DEFAULT_NIFTI_ROOT = REPO_ROOT / "data" / "input"  / "nifti"
DEFAULT_OUT_CSV    = REPO_ROOT / "data" / "output" / "l3_manifest.csv"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ts-root",    type=Path, default=DEFAULT_TS_ROOT)
    p.add_argument("--nifti-root", type=Path, default=DEFAULT_NIFTI_ROOT)
    p.add_argument("--out",        type=Path, default=DEFAULT_OUT_CSV)
    args = p.parse_args()

    print(f"TS root    : {args.ts_root}")
    print(f"NIfTI root : {args.nifti_root}")
    print(f"Out CSV    : {args.out}")

    df = build_manifest(str(args.ts_root), str(args.nifti_root), str(args.out))
    print(f"\nWrote {len(df)} rows -> {args.out}")
    if len(df):
        print(df.head().to_string(index=False))


if __name__ == "__main__":
    main()
