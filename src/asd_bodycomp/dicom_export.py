"""Write a single 2-D NumPy array as a CT-style DICOM file.

Geometry (ImagePositionPatient / ImageOrientationPatient / PixelSpacing /
SliceThickness) is recovered from the original NIfTI affine so the exported
slice can be loaded into any DICOM viewer with correct orientation.
"""

from __future__ import annotations

import datetime
import os
from typing import Optional

import numpy as np
import pydicom
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import ExplicitVRLittleEndian, generate_uid


COPY_TAGS = [
    "PatientID", "PatientName", "PatientBirthDate", "PatientSex",
    "StudyInstanceUID", "StudyDate", "StudyTime", "StudyDescription",
    "AccessionNumber", "Modality", "Manufacturer", "InstitutionName",
    "RescaleIntercept", "RescaleSlope", "WindowCenter", "WindowWidth",
    "KVP", "XRayTubeCurrent",
]


def affine_to_dicom_geometry(affine: np.ndarray, axial_axis: int):
    """Return (orientation_6, voxel_size_xyz, origin_xyz) from a NIfTI affine."""
    M         = affine[:3, :3]
    col_norms = np.linalg.norm(M, axis=0)
    plane_axes = [ax for ax in range(3) if ax != axial_axis]
    row_ax, col_ax = plane_axes[0], plane_axes[1]
    row_cos    = (M[:, row_ax] / col_norms[row_ax]).tolist()
    col_cos    = (M[:, col_ax] / col_norms[col_ax]).tolist()
    voxel_size = (col_norms[row_ax], col_norms[col_ax], col_norms[axial_axis])
    return row_cos + col_cos, voxel_size, affine[:3, 3].tolist()


def slice_image_position(affine: np.ndarray, slice_idx: int, axial_axis: int):
    vox = np.zeros(3)
    vox[axial_axis] = slice_idx
    return (affine[:3, :3] @ vox + affine[:3, 3]).tolist()


def build_dicom(pixel_array, series_uid, instance_number,
                series_description, series_number,
                image_position, image_orientation,
                pixel_spacing, slice_thickness,
                orig_meta=None):
    """Construct a CT image FileDataset from a 2-D int16 array and geometry."""
    arr     = pixel_array.astype(np.int16)
    sop_uid = generate_uid()

    file_meta = FileMetaDataset()
    file_meta.MediaStorageSOPClassUID    = "1.2.840.10008.5.1.4.1.1.2"
    file_meta.MediaStorageSOPInstanceUID = sop_uid
    file_meta.TransferSyntaxUID          = ExplicitVRLittleEndian

    ds = FileDataset("", {}, file_meta=file_meta, preamble=b"\0" * 128)
    ds.is_implicit_VR   = False
    ds.is_little_endian = True

    now = datetime.datetime.now()
    ds.ContentDate = now.strftime("%Y%m%d")
    ds.ContentTime = now.strftime("%H%M%S.%f")

    if orig_meta is not None:
        for tag in COPY_TAGS:
            if hasattr(orig_meta, tag):
                try:
                    setattr(ds, tag, getattr(orig_meta, tag))
                except Exception:
                    pass

    ds.SOPClassUID             = file_meta.MediaStorageSOPClassUID
    ds.SOPInstanceUID          = sop_uid
    ds.SeriesInstanceUID       = series_uid
    ds.SeriesDescription       = series_description
    ds.SeriesNumber            = str(series_number)
    ds.InstanceNumber          = str(instance_number)
    if not hasattr(ds, "Modality"):
        ds.Modality = "CT"

    ds.ImagePositionPatient    = [f"{v:.6f}" for v in image_position]
    ds.ImageOrientationPatient = [f"{v:.6f}" for v in image_orientation]
    ds.PixelSpacing            = [f"{pixel_spacing[0]:.6f}",
                                  f"{pixel_spacing[1]:.6f}"]
    ds.SliceThickness          = f"{slice_thickness:.4f}"
    if not hasattr(ds, "RescaleIntercept"):
        ds.RescaleIntercept = "0"
    if not hasattr(ds, "RescaleSlope"):
        ds.RescaleSlope = "1"

    ds.SamplesPerPixel           = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.Rows, ds.Columns          = arr.shape
    ds.BitsAllocated             = 16
    ds.BitsStored                = 16
    ds.HighBit                   = 15
    ds.PixelRepresentation       = 1
    ds.PixelData                 = arr.tobytes()
    return ds


def save_dicom(ds, path: str):
    pydicom.dcmwrite(path, ds)


def try_load_dicom_meta(nifti_path: str):
    """If the NIfTI's parent dir contains DICOM(s), return one dataset for tag-copying.

    Returns None when no DICOM is present (the typical NIfTI-only case).
    """
    d = os.path.dirname(nifti_path)
    if not os.path.isdir(d):
        return None
    for f in sorted(os.listdir(d)):
        if f.lower().endswith(".dcm"):
            try:
                return pydicom.dcmread(os.path.join(d, f), stop_before_pixels=True)
            except Exception:
                pass
    return None
