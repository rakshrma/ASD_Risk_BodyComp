"""Shared pytest fixtures.

Builds a small synthetic CT volume + matching TotalSegmentator outputs that
contain the structures the pipeline expects (spinal cord, psoas, L1-L5, fat,
muscle). All tests run against this fixture so the repo can ship without any
real patient data.
"""

from __future__ import annotations

import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))


# ── synthetic-volume layout (constants used by tests) ───────────────────────

VOL_SHAPE   = (32, 32, 16)
PIXDIM_MM   = (2.0, 2.0, 3.0)        # in-plane 2mm, 3mm slice thickness
ABDOMEN_Z   = range(4, 13)           # z indices that contain the abdomen
SC_LABEL    = 79                     # TotalSegmentator total label for SC
PSOAS_L     = 89                     # TS total label for left psoas
PSOAS_R     = 88                     # TS total label for right psoas
VERT_LABEL_PER_Z = {4: 27, 6: 28, 8: 29, 10: 30, 12: 31}  # L5..L1
L3_Z        = 8                      # axial slice that contains the L3 centroid
TISSUE_SUBQ_FAT = 1
TISSUE_MUSCLE   = 3


def _build_synthetic_volumes():
    """Return (ct, total_seg, tissue_seg, affine) numpy arrays + affine matrix."""
    nx, ny, nz = VOL_SHAPE
    ct     = np.full(VOL_SHAPE, -1000, dtype=np.int16)   # air outside body
    total  = np.zeros(VOL_SHAPE, dtype=np.int16)
    tissue = np.zeros(VOL_SHAPE, dtype=np.int16)

    for z in ABDOMEN_Z:
        # Body block (square cross-section, simulates abdomen)
        ct[4:28, 4:28, z] = 50                            # muscle-equivalent default

        # Subcutaneous fat ring (outermost pixels of body)
        for i in range(4, 28):
            for j in range(4, 28):
                if i in (4, 27) or j in (4, 27):
                    ct[i, j, z]     = -100
                    tissue[i, j, z] = TISSUE_SUBQ_FAT

        # Vertebral body (bone) — bright region (no tissue label, bone isn't in tissue_4_types)
        ct[6:12, 12:20, z] = 500

        # Iliopsoas pair, anterior to spinal cord
        # left psoas (label 89) and right psoas (label 88)
        ct[13:17, 10:14, z]     = 60
        total[13:17, 10:14, z]  = PSOAS_L
        tissue[13:17, 10:14, z] = TISSUE_MUSCLE
        ct[13:17, 18:22, z]     = 60
        total[13:17, 18:22, z]  = PSOAS_R
        tissue[13:17, 18:22, z] = TISSUE_MUSCLE

        # Add a strip of "lateral abdominal muscle" so muscle area is measurable
        ct[14:20, 5:8, z]     = 55
        tissue[14:20, 5:8, z] = TISSUE_MUSCLE
        ct[14:20, 24:27, z]     = 55
        tissue[14:20, 24:27, z] = TISSUE_MUSCLE

    # Vertebral labels (one z per L1..L5). Write these BEFORE the spinal cord
    # so the SC label is not overwritten where the two regions overlap.
    for z, vlabel in VERT_LABEL_PER_Z.items():
        total[6:12, 12:20, z] = vlabel

    # Spinal cord runs through every abdomen slice (in the vertebral foramen).
    # Write last so it always wins at the SC pixels.
    for z in ABDOMEN_Z:
        ct[8:11, 15:18, z]    = 30
        total[8:11, 15:18, z] = SC_LABEL

    affine = np.diag([PIXDIM_MM[0], PIXDIM_MM[1], PIXDIM_MM[2], 1.0])
    return ct, total, tissue, affine


@pytest.fixture(scope="session")
def synthetic_arrays():
    """Return raw numpy arrays + affine (no disk I/O)."""
    return _build_synthetic_volumes()


@pytest.fixture
def synthetic_case(tmp_path, synthetic_arrays):
    """Write the three NIfTI files into tmp_path and return a dict of paths."""
    ct, total, tissue, affine = synthetic_arrays
    accession = "ACC0001"
    nifti_dir = tmp_path / "nifti"   / accession
    ts_dir    = tmp_path / "tsout"   / accession
    nifti_dir.mkdir(parents=True)
    ts_dir.mkdir(parents=True)

    stem = "patient_001_ACC0001_001"
    nifti_path  = nifti_dir / f"{stem}.nii.gz"
    total_path  = ts_dir    / f"{stem}_ts_total.nii.gz"
    tissue_path = ts_dir    / f"{stem}_ts_tissue.nii.gz"

    nib.save(nib.Nifti1Image(ct,     affine), str(nifti_path))
    nib.save(nib.Nifti1Image(total,  affine), str(total_path))
    nib.save(nib.Nifti1Image(tissue, affine), str(tissue_path))

    return {
        "accession":      accession,
        "stem":           stem,
        "nifti_path":     str(nifti_path),
        "ts_total_path":  str(total_path),
        "ts_tissue_path": str(tissue_path),
        "nifti_root":     str(tmp_path / "nifti"),
        "ts_root":        str(tmp_path / "tsout"),
        "out_dir":        str(tmp_path / "out"),
        "affine":         affine,
        "l3_slice":       L3_Z,
    }
