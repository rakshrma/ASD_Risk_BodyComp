"""Build the L3 manifest CSV consumed by step 5/6.

The manifest enumerates one row per CT case with these columns:

  l3_slice_index   axial slice index of the L3 vertebral centroid
  nifti_path       absolute path to original CT NIfTI
  ts_total_path    absolute path to TotalSegmentator -ta total output
  ts_tissue_path   absolute path to TotalSegmentator -ta tissue_4_types output
  accession        accession ID (taken from the TS subdir name)

The manifest's role is to decouple the search/discovery stage from the
measurement stage; downstream scripts only need to read this CSV.
"""

from __future__ import annotations

import glob
import os
from typing import List, Optional

import nibabel as nib
import numpy as np
import pandas as pd

from .labels import L3_VERTEBRA_TOTAL_LABEL
from .nifti_io import get_axial_axis


def _l3_slice_index(ts_total_path: str) -> Optional[int]:
    """Compute the axial slice index of the L3 vertebra centroid.

    Returns None if the L3 label is not present in the volume.
    """
    img    = nib.load(ts_total_path)
    data   = np.asarray(img.dataobj)
    affine = img.affine
    axial  = get_axial_axis(affine)

    mask = data == L3_VERTEBRA_TOTAL_LABEL
    if not mask.any():
        return None

    if axial == 0:
        occupied = np.where(mask.any(axis=(1, 2)))[0]
    elif axial == 1:
        occupied = np.where(mask.any(axis=(0, 2)))[0]
    else:
        occupied = np.where(mask.any(axis=(0, 1)))[0]
    return int(round(occupied.mean()))


def discover_cases(ts_root: str, nifti_root: str) -> List[dict]:
    """Walk `ts_root/<accession>/` for *_ts_tissue.nii.gz and pair with siblings.

    For each tissue file we expect:
      ts_root/<accession>/<stem>_ts_tissue.nii.gz
      ts_root/<accession>/<stem>_ts_total.nii.gz
      nifti_root/<accession>/<stem>.nii.gz   (or <stem>*.nii.gz)
    """
    cases: List[dict] = []
    if not os.path.isdir(ts_root):
        return cases

    for accession in sorted(os.listdir(ts_root)):
        acc_dir = os.path.join(ts_root, accession)
        if not os.path.isdir(acc_dir):
            continue

        for tissue_path in sorted(glob.glob(os.path.join(acc_dir, "*_ts_tissue.nii.gz"))):
            stem = os.path.basename(tissue_path).replace("_ts_tissue.nii.gz", "")

            total_matches = glob.glob(os.path.join(acc_dir, f"{stem}_ts_total.nii.gz"))
            if not total_matches:
                continue
            total_path = total_matches[0]

            nifti_dir = os.path.join(nifti_root, accession)
            nifti_matches = sorted(glob.glob(os.path.join(nifti_dir, f"{stem}*.nii.gz")))
            if not nifti_matches:
                continue
            nifti_path = nifti_matches[0]

            cases.append({
                "accession":      accession,
                "stem":           stem,
                "nifti_path":     nifti_path,
                "ts_total_path":  total_path,
                "ts_tissue_path": tissue_path,
            })
    return cases


def build_manifest(ts_root: str, nifti_root: str, out_csv: str) -> pd.DataFrame:
    """Write the manifest CSV; return it as a DataFrame."""
    cases = discover_cases(ts_root, nifti_root)
    rows = []
    for case in cases:
        slice_idx = _l3_slice_index(case["ts_total_path"])
        if slice_idx is None:
            print(f"  SKIP {case['accession']}/{case['stem']}: "
                  f"L3 label {L3_VERTEBRA_TOTAL_LABEL} absent in total seg.")
            continue
        rows.append({
            "l3_slice_index": slice_idx,
            "nifti_path":     case["nifti_path"],
            "ts_total_path":  case["ts_total_path"],
            "ts_tissue_path": case["ts_tissue_path"],
            "accession":      case["accession"],
        })

    cols = ["l3_slice_index", "nifti_path", "ts_total_path",
            "ts_tissue_path", "accession"]
    df = pd.DataFrame(rows, columns=cols)
    os.makedirs(os.path.dirname(os.path.abspath(out_csv)) or ".", exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df
