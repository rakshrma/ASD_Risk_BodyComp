#!/usr/bin/env python3
"""
05_process_l3.py
================
Apply the L3 / L1-L5 measurement pipeline to every row in the manifest CSV.

For each row this script writes (under <output>/<accession>/):
  <stem>_l3_image.dcm        L3 CT slice as DICOM
  <stem>_l3_label.dcm        L3 label map as DICOM
  <stem>_l3_image.png        PNG of CT slice
  <stem>_l3_label.png        PNG of label overlay
  <stem>_l1_l5_image.nii.gz  L1-L5 CT sub-volume
  <stem>_l1_l5_label.nii.gz  L1-L5 label sub-volume (labels 1/3/4)

It also appends per-row measurement columns and writes the combined results to
  <output>/results_summary.csv

Default folders:
  input  : data/output/l3_manifest.csv   (output of script 04)
  output : data/output/l3_results/

Usage:
  python scripts/05_process_l3.py
  python scripts/05_process_l3.py --input <csv> --output <dir> \
                                  --bbox-ant-mm 55 --bbox-lr-mm 85
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from asd_bodycomp.pipeline import run_manifest

DEFAULT_INPUT  = REPO_ROOT / "data" / "output" / "l3_manifest.csv"
DEFAULT_OUTPUT = REPO_ROOT / "data" / "output" / "l3_results"

REQUIRED_COLS = {"l3_slice_index", "nifti_path", "ts_total_path",
                 "ts_tissue_path", "accession"}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input",        type=Path,  default=DEFAULT_INPUT)
    p.add_argument("--output",       type=Path,  default=DEFAULT_OUTPUT)
    p.add_argument("--bbox-ant-mm",  type=float, default=55.0,
                   help="Anterior bounding-box extent in mm (default: 55)")
    p.add_argument("--bbox-lr-mm",   type=float, default=85.0,
                   help="Left/right bounding-box extent in mm (default: 85)")
    p.add_argument("--limit",        type=int,   default=-1,
                   help="Max rows to process; -1 = all (default: -1)")
    args = p.parse_args()

    if not args.input.is_file():
        sys.exit(f"ERROR: manifest CSV not found: {args.input}")

    df = pd.read_csv(args.input)
    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        sys.exit(f"ERROR: manifest missing columns: {missing}")

    args.output.mkdir(parents=True, exist_ok=True)
    print(f"Rows         : {len(df)}")
    print(f"bbox ant/lr  : {args.bbox_ant_mm} / {args.bbox_lr_mm} mm")
    print(f"Output root  : {args.output}\n")

    result_df, errors = run_manifest(
        df, args.bbox_ant_mm, args.bbox_lr_mm,
        str(args.output), alt_sc_map=None, limit=args.limit,
    )

    out_csv = args.output / "results_summary.csv"
    result_df.to_csv(out_csv, index=False)
    print(f"\nResults  : {out_csv}")
    print(f"Succeeded: {len(df) - len(errors)} / {len(df)}")
    for i, msg in errors:
        print(f"  Row {i}: {msg}")


if __name__ == "__main__":
    main()
