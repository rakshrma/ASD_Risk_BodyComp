#!/usr/bin/env python3
"""
06_process_l3_with_alt_sc.py
============================
Same as scripts/05_process_l3.py but accepts a second CSV that supplies
an alternative spinal-cord segmentation per case. Useful when the
TotalSegmentator spinal-cord label (79) is missing or unreliable on a case.

Filter CSV columns:
  nifti_path   path matching a row in the manifest
  seg_path     path to a NIfTI containing the alternative SC mask (label 1)

Only manifest rows whose `nifti_path` appears in the filter CSV are processed,
and each one uses the matching `seg_path` for SC detection.

Default folders:
  input        : data/output/l3_manifest.csv          (output of script 04)
  filter_csv   : data/input/alt_sc_paths.csv          (user-provided)
  output       : data/output/l3_results_with_alt_sc/

Usage:
  python scripts/06_process_l3_with_alt_sc.py \\
      --input data/output/l3_manifest.csv \\
      --filter-csv data/input/alt_sc_paths.csv \\
      --output data/output/l3_results_with_alt_sc
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from asd_bodycomp.pipeline import run_manifest

DEFAULT_INPUT      = REPO_ROOT / "data" / "output" / "l3_manifest.csv"
DEFAULT_FILTER_CSV = REPO_ROOT / "data" / "input"  / "alt_sc_paths.csv"
DEFAULT_OUTPUT     = REPO_ROOT / "data" / "output" / "l3_results_with_alt_sc"

REQUIRED_COLS        = {"l3_slice_index", "nifti_path", "ts_total_path",
                        "ts_tissue_path", "accession"}
REQUIRED_FILTER_COLS = {"nifti_path", "seg_path"}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input",       type=Path,  default=DEFAULT_INPUT)
    p.add_argument("--filter-csv",  type=Path,  default=DEFAULT_FILTER_CSV)
    p.add_argument("--output",      type=Path,  default=DEFAULT_OUTPUT)
    p.add_argument("--bbox-ant-mm", type=float, default=55.0)
    p.add_argument("--bbox-lr-mm",  type=float, default=85.0)
    p.add_argument("--limit",       type=int,   default=-1)
    args = p.parse_args()

    if not args.input.is_file():
        sys.exit(f"ERROR: manifest CSV not found: {args.input}")
    if not args.filter_csv.is_file():
        sys.exit(f"ERROR: filter CSV not found: {args.filter_csv}")

    df      = pd.read_csv(args.input)
    fdf     = pd.read_csv(args.filter_csv)
    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        sys.exit(f"ERROR: manifest missing columns: {missing}")
    fmissing = REQUIRED_FILTER_COLS - set(fdf.columns)
    if fmissing:
        sys.exit(f"ERROR: filter CSV missing columns: {fmissing}")

    alt_map = {str(r["nifti_path"]).strip(): str(r["seg_path"]).strip()
               for _, r in fdf.iterrows()}
    df = df[df["nifti_path"].apply(lambda x: str(x).strip() in alt_map)].copy()
    if df.empty:
        sys.exit("ERROR: no manifest rows match the filter CSV.")

    args.output.mkdir(parents=True, exist_ok=True)
    print(f"Filtered rows : {len(df)}")
    print(f"Output root   : {args.output}\n")

    result_df, errors = run_manifest(
        df, args.bbox_ant_mm, args.bbox_lr_mm,
        str(args.output), alt_sc_map=alt_map, limit=args.limit,
    )

    out_csv = args.output / "results_summary.csv"
    result_df.to_csv(out_csv, index=False)
    print(f"\nResults  : {out_csv}")
    print(f"Succeeded: {len(df) - len(errors)} / {len(df)}")
    for i, msg in errors:
        print(f"  Row {i}: {msg}")


if __name__ == "__main__":
    main()
