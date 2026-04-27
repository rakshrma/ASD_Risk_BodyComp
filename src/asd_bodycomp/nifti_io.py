"""NIfTI loading helpers and 2-D slice extraction."""

from __future__ import annotations

import os

import nibabel as nib
import numpy as np


def nifti_stem(path: str) -> str:
    """Filename without the .nii or .nii.gz extension."""
    base = os.path.basename(path)
    for ext in (".nii.gz", ".nii"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def load_nifti(path: str):
    """Return (data, affine, header, pixdim_xyz_mm)."""
    img    = nib.load(path)
    data   = np.asarray(img.dataobj)
    pixdim = np.abs(img.header.get_zooms()[:3])
    return data, img.affine, img.header, pixdim


def get_axial_axis(affine: np.ndarray) -> int:
    """Index (0/1/2) of the volume axis closest to the patient inferior-superior direction."""
    return int(np.argmax(np.abs(affine[:3, :3])[2]))


def volume_n_slices(volume: np.ndarray, axis: int) -> int:
    return volume.shape[axis]


def extract_slice_2d(volume: np.ndarray, idx: int, axis: int = 2) -> np.ndarray:
    """Return a 2-D slice along the given axis."""
    if axis == 0:
        return volume[idx, :, :]
    if axis == 1:
        return volume[:, idx, :]
    return volume[:, :, idx]
