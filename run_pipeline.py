#!/usr/bin/env python3
"""
run_pipeline.py
===============
Single-command end-to-end driver. Drop your data into one of these folders
and run this script:

  data/input/dicom/<accession>/...    (raw or compressed DICOMs)
              -> decompress -> NIfTI -> TotalSegmentator -> manifest -> L3
  data/input/nifti/<accession>/<x>.nii.gz
              -> TotalSegmentator -> manifest -> L3
  data/output/totalsegmentator/<accession>/<x>_ts_total.nii.gz + _ts_tissue.nii.gz
              -> manifest -> L3

The script auto-detects which entry point applies. Outputs land in
data/output/l3_results/ (DICOM, PNG, NIfTI sub-volumes, results_summary.csv).

Usage:
  python run_pipeline.py
  python run_pipeline.py --skip-decompress       # already decompressed DICOMs
  python run_pipeline.py --bbox-ant-mm 60 --bbox-lr-mm 90
  python run_pipeline.py --dry-run               # print steps without executing
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
SCRIPTS   = REPO_ROOT / "scripts"

DICOM_IN          = REPO_ROOT / "data" / "input"  / "dicom"
DICOM_DECOMP      = REPO_ROOT / "data" / "input"  / "dicom_decompressed"
NIFTI_IN          = REPO_ROOT / "data" / "input"  / "nifti"
TS_OUT            = REPO_ROOT / "data" / "output" / "totalsegmentator"
MANIFEST_CSV      = REPO_ROOT / "data" / "output" / "l3_manifest.csv"
RESULTS_OUT       = REPO_ROOT / "data" / "output" / "l3_results"


def _has_files(p: Path, glob: str = "*") -> bool:
    return p.is_dir() and any(p.glob(glob))


def _has_dicoms(p: Path) -> bool:
    return p.is_dir() and any(p.iterdir()) and \
           any(any(d.iterdir()) for d in p.iterdir() if d.is_dir())


def _has_niftis(p: Path) -> bool:
    return p.is_dir() and any(p.rglob("*.nii.gz"))


def _has_ts_outputs(p: Path) -> bool:
    return p.is_dir() and any(p.rglob("*_ts_tissue.nii.gz"))


def _step(label: str):
    print(f"\n{'='*72}\n  STEP: {label}\n{'='*72}")


def _run(cmd, dry_run: bool):
    print(f"  $ {' '.join(str(c) for c in cmd)}")
    if dry_run:
        return 0
    return subprocess.call([str(c) for c in cmd])


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--skip-decompress", action="store_true",
                   help="Skip step 1 even if data/input/dicom is populated.")
    p.add_argument("--bbox-ant-mm", type=float, default=55.0)
    p.add_argument("--bbox-lr-mm",  type=float, default=85.0)
    p.add_argument("--dry-run", action="store_true",
                   help="Print steps without executing.")
    args = p.parse_args()

    print("ASD Risk Body Composition pipeline")
    print(f"  Repo root: {REPO_ROOT}")

    # ── Decide entry point ──────────────────────────────────────────────────
    do_decompress = False
    do_dcm2niix   = False
    do_totalseg   = False

    # Precedence: prefer the latest existing artifact so we don't redo work.
    # Users with TS outputs already on disk skip directly to manifest + L3.
    if _has_ts_outputs(TS_OUT):
        pass  # straight to step 04 + 05
    elif _has_dicoms(DICOM_IN) and not args.skip_decompress:
        do_decompress = True
        do_dcm2niix   = True
        do_totalseg   = True
    elif _has_dicoms(DICOM_DECOMP):
        do_dcm2niix = True
        do_totalseg = True
    elif _has_niftis(NIFTI_IN):
        do_totalseg = True
    else:
        sys.exit(
            "ERROR: No input data found.\n"
            "  Place DICOMs under   data/input/dicom/<accession>/\n"
            "  OR     NIfTIs under  data/input/nifti/<accession>/<file>.nii.gz\n"
            "  OR     TS outputs under data/output/totalsegmentator/<accession>/"
        )

    # ── 01: decompress ──────────────────────────────────────────────────────
    if do_decompress:
        _step("01  decompress DICOM")
        rc = _run(["bash", SCRIPTS / "01_decompress_dicom.sh"], args.dry_run)
        if rc != 0:
            sys.exit(f"01_decompress_dicom.sh exited {rc}")

    # ── 02: dcm2niix ────────────────────────────────────────────────────────
    if do_dcm2niix:
        _step("02  DICOM -> NIfTI")
        rc = _run([sys.executable, SCRIPTS / "02_dicom_to_nifti.py"], args.dry_run)
        if rc != 0:
            sys.exit(f"02_dicom_to_nifti.py exited {rc}")

    # ── 03: TotalSegmentator ────────────────────────────────────────────────
    if do_totalseg:
        _step("03  TotalSegmentator (total + tissue_4_types)")
        rc = _run([sys.executable, SCRIPTS / "03_run_totalsegmentator.py"],
                  args.dry_run)
        if rc != 0:
            sys.exit(f"03_run_totalsegmentator.py exited {rc}")

    # ── 04: manifest ────────────────────────────────────────────────────────
    _step("04  build L3 manifest")
    rc = _run([sys.executable, SCRIPTS / "04_build_l3_manifest.py"], args.dry_run)
    if rc != 0:
        sys.exit(f"04_build_l3_manifest.py exited {rc}")

    if not args.dry_run and not MANIFEST_CSV.is_file():
        sys.exit(f"ERROR: manifest CSV not produced at {MANIFEST_CSV}")

    # ── 05: L3 measurements ────────────────────────────────────────────────
    _step("05  L3 + L1-L5 measurements")
    rc = _run([sys.executable, SCRIPTS / "05_process_l3.py",
               "--bbox-ant-mm", args.bbox_ant_mm,
               "--bbox-lr-mm",  args.bbox_lr_mm],
              args.dry_run)
    if rc != 0:
        sys.exit(f"05_process_l3.py exited {rc}")

    print(f"\n{'='*72}")
    print(f"  DONE. Outputs under: {RESULTS_OUT}")
    print(f"  Per-case files     : <accession>/<stem>_l3_image.{{dcm,png}}, "
          f"<stem>_l3_label.{{dcm,png}}, <stem>_l1_l5_{{image,label}}.nii.gz")
    print(f"  Combined CSV       : {RESULTS_OUT / 'results_summary.csv'}")
    print(f"{'='*72}")


if __name__ == "__main__":
    main()
