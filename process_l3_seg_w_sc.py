"""
L3 / L1-L5 Slice Segmentation Processor (Modified)
===================================================
Processes NIfTI segmentation files at the L3 vertebral level (and L1-L5 volume)
and exports results as DICOM, preserving original orientation and metadata.

MODIFICATIONS:
- Accepts optional --filter_csv to process only specific nifti_paths
- Uses alternative spinal cord segmentation (seg_path) when SC not found in _total.nii.gz
- Alternative SC segmentation has label 1 for spinal cord
- Reports L3 vertebral body cross-sectional area (cm²)
- Reports mean/min/max/median HU for each label (soft-tissue + vertebra) at L3

Usage:
    python process_l3_segmentation.py --input <csv_file> [--filter_csv <filter_csv>] [--bbox_ant_mm 70] [--bbox_lr_mm 100] [--output ./l3_output]

CSV columns required (--input):
    - l3_slice_index  : int, axial slice index for L3
    - nifti_path      : path to original *.nii.gz
    - ts_total_path   : path to *_total.nii.gz segmentation
    - ts_tissue_path  : path to *_tissue.nii.gz segmentation
    - accession       : folder name for this case's outputs

CSV columns required (--filter_csv):
    - nifti_path      : path to original *.nii.gz (used for filtering)
    - seg_path        : path to alternative spinal cord segmentation (label 1)

Label reference (_total.nii.gz):
    79  = spinal cord
    88  = iliopsoas_left
    89  = iliopsoas_right

Label reference (alternative seg_path):
    1   = spinal cord

Label reference (_tissue.nii.gz):
    1   = bone  -> retained as label 1
    3,4 = merged -> label 2, then split by HU into 3 (fat) and 4 (muscle)
    27  = L5 vertebra
    28  = L4 vertebra
    29  = L3 vertebra
    30  = L2 vertebra
    31  = L1 vertebra

Output per row  (inside <output>/<accession>/):
    <stem>_l3_image.dcm          L3 CT slice as DICOM
    <stem>_l3_label.dcm          L3 label map as DICOM
    <stem>_l3_image.png          PNG of CT slice
    <stem>_l3_label.png          PNG of label overlay
    <stem>_l1_l5_image.nii.gz    L1-L5 CT sub-volume (NIfTI)
    <stem>_l1_l5_label.nii.gz    L1-L5 label sub-volume (NIfTI, labels 1/3/4)

Also appended to input CSV (saved as <output>/results_summary.csv):
    l3_<name>_area_cm2           Cross-sectional area of each soft-tissue label on L3 slice
    l3_<name>_mean_hu            Mean HU of each soft-tissue label on L3 slice
    l3_<name>_min_hu             Min HU of each soft-tissue label on L3 slice
    l3_<name>_max_hu             Max HU of each soft-tissue label on L3 slice
    l3_<name>_median_hu          Median HU of each soft-tissue label on L3 slice
    l3_vertebra_area_cm2         Cross-sectional area of L3 vertebral body on L3 slice
    l3_vertebra_mean_hu          Mean HU of L3 vertebral body on L3 slice
    l3_vertebra_min_hu           Min HU of L3 vertebral body on L3 slice
    l3_vertebra_max_hu           Max HU of L3 vertebral body on L3 slice
    l3_vertebra_median_hu        Median HU of L3 vertebral body on L3 slice
    vol_<name>_volume_cm3        Volume of each label across L1-L5
    vol_<name>_mean_hu           Mean HU of each label across L1-L5

where <stem> is the original NIfTI filename without extension,
and <name> is one of: subq_fat, imat, muscle.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import nibabel as nib
import pydicom
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import generate_uid, ExplicitVRLittleEndian
import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ─────────────────────────────────────────────────────────────────────────────
# Label constants
# ─────────────────────────────────────────────────────────────────────────────
SPINAL_CORD_LABEL     = 79
SPINAL_CORD_ALT_LABEL = 1              # Label in alternative segmentation
ILIOPSOAS_LEFT_LABEL  = 88
ILIOPSOAS_RIGHT_LABEL = 89
TOTAL_MERGE_LABELS    = [88, 89]
TISSUE_MERGE_LABELS   = [3, 4]
TISSUE_RETAIN_LABEL   = 1
NEW_MERGED_LABEL      = 2
FAT_SUB_LABEL         = 3             # HU < -30 and > -190 within label 2
MUSCLE_SUB_LABEL      = 4             # HU >= -30 within label 2
HU_FAT_THRESHOLD_HIGH = -30
HU_FAT_THRESHOLD_LOW  = -190
BODY_HU_THRESHOLD     = -720
SC_SEARCH_RADIUS      = 5             # ±slices to search if SC / psoas missing

VERTEBRA_LABELS = {"L1": 31, "L2": 30, "L3": 29, "L4": 28, "L5": 27}
L3_VERTEBRA_TISSUE_LABEL = 29         # tissue label for L3 vertebral body

# Labels present in final output maps
OUTPUT_LABELS = [1, 3, 4]
LABEL_NAMES   = {1: "subq_fat", 3: "imat", 4: "muscle"}

COPY_TAGS = [
    "PatientID", "PatientName", "PatientBirthDate", "PatientSex",
    "StudyInstanceUID", "StudyDate", "StudyTime", "StudyDescription",
    "AccessionNumber", "Modality", "Manufacturer", "InstitutionName",
    "RescaleIntercept", "RescaleSlope", "WindowCenter", "WindowWidth",
    "KVP", "XRayTubeCurrent",
]

# ─────────────────────────────────────────────────────────────────────────────
# NIfTI / utility helpers
# ─────────────────────────────────────────────────────────────────────────────

def nifti_stem(path):
    """Return filename without .nii.gz or .nii extension."""
    base = os.path.basename(path)
    for ext in (".nii.gz", ".nii"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def load_nifti(path):
    img    = nib.load(path)
    data   = np.asarray(img.dataobj)
    pixdim = np.abs(img.header.get_zooms()[:3])
    return data, img.affine, img.header, pixdim


def get_axial_axis(affine):
    return int(np.argmax(np.abs(affine[:3, :3])[2]))


def volume_n_slices(volume, axis):
    return volume.shape[axis]


def extract_slice_2d(volume, idx, axis=2):
    if axis == 0: return volume[idx, :, :]
    if axis == 1: return volume[:, idx, :]
    return volume[:, :, idx]

# ─────────────────────────────────────────────────────────────────────────────
# Spinal cord / psoas with ±N-slice fallback
# ─────────────────────────────────────────────────────────────────────────────

def _sc_center_at(total_data, idx, axis, sc_label=SPINAL_CORD_LABEL):
    """Find spinal cord center at a specific slice."""
    slc  = extract_slice_2d(total_data, idx, axis).astype(np.int32)
    unique_in_slice = np.unique(slc[slc > 0])
    print(f"      DEBUG: slice {idx}, axis {axis}, looking for label {sc_label}, slice dtype={slc.dtype}")
    print(f"      DEBUG: unique labels in slice: {list(unique_in_slice)}")
    mask = slc == sc_label
    print(f"      DEBUG: mask sum = {mask.sum()} pixels")
    if not mask.any():
        raise ValueError(f"SC label {sc_label} absent at slice {idx}")
    return slc, mask, np.argwhere(mask).mean(axis=0)


def find_spinal_cord_center_with_fallback(total_data, slice_idx, axial_axis,
                                          alt_sc_data=None,
                                          search_radius=SC_SEARCH_RADIUS):
    """
    Find spinal cord center with fallback strategy:
    1. If alt_sc_data provided, use alternative segmentation (label 1) ONLY
    2. Otherwise, try primary segmentation (label 79) at slice_idx ± search_radius
    3. Raise ValueError if not found
    """
    total_n = volume_n_slices(total_data, axial_axis)
    offsets = [0]
    for d in range(1, search_radius + 1):
        offsets += [-d, +d]

    # If alternative segmentation provided, use it exclusively
    if alt_sc_data is not None:
        print(f"    Using alternative segmentation for SC detection (label {SPINAL_CORD_ALT_LABEL})...")
        for offset in offsets:
            idx = slice_idx + offset
            if not (0 <= idx < total_n):
                continue
            try:
                slc_alt, mask_alt, center = _sc_center_at(alt_sc_data, idx, axial_axis, SPINAL_CORD_ALT_LABEL)
                n_pixels = int(mask_alt.sum())
                if offset != 0:
                    print(f"    SC found in alternative segmentation at slice {idx} (offset {offset:+d})")
                else:
                    print(f"    SC found in alternative segmentation at slice {idx}")
                print(f"    Alternative seg: {n_pixels} pixels with label {SPINAL_CORD_ALT_LABEL}")
                return center, idx
            except ValueError:
                continue
        raise ValueError(
            f"Spinal cord (label {SPINAL_CORD_ALT_LABEL}) not found in alternative segmentation "
            f"within ±{search_radius} slices of slice {slice_idx}."
        )

    # Otherwise use primary segmentation
    for offset in offsets:
        idx = slice_idx + offset
        if not (0 <= idx < total_n):
            continue
        try:
            _, _, center = _sc_center_at(total_data, idx, axial_axis, SPINAL_CORD_LABEL)
            if offset != 0:
                print(f"    SC not found at slice {slice_idx}; "
                      f"using center from slice {idx} (offset {offset:+d})")
            return center, idx
        except ValueError:
            continue

    raise ValueError(
        f"Spinal cord (label {SPINAL_CORD_LABEL}) not found in primary segmentation "
        f"within ±{search_radius} slices of slice {slice_idx}."
    )


def _psoas_at(total_data, idx, axis, sc_center):
    slc   = extract_slice_2d(total_data, idx, axis).astype(np.int32)
    sc_r, sc_c = sc_center
    found = {}
    for name, lbl in [("left", ILIOPSOAS_LEFT_LABEL), ("right", ILIOPSOAS_RIGHT_LABEL)]:
        mask = slc == lbl
        if mask.any():
            found[name] = np.argwhere(mask).mean(axis=0)
    if not found:
        raise ValueError(f"Psoas absent at slice {idx}")
    centers   = np.array(list(found.values()))
    psoas_ctr = centers.mean(axis=0)
    vec       = psoas_ctr - np.array([sc_r, sc_c])
    if abs(vec[0]) >= abs(vec[1]):
        ant_axis = 0
        ant_sign = int(np.sign(vec[0])) or 1
    else:
        ant_axis = 1
        ant_sign = int(np.sign(vec[1])) or 1
    left_col  = found.get("left",  np.array([sc_r, sc_c]))[1]
    right_col = found.get("right", np.array([sc_r, sc_c]))[1]
    return ant_axis, ant_sign, left_col, right_col


def find_psoas_direction_with_fallback(total_data, slice_idx, axial_axis,
                                       sc_center, search_radius=SC_SEARCH_RADIUS):
    total_n = volume_n_slices(total_data, axial_axis)
    offsets = [0]
    for d in range(1, search_radius + 1):
        offsets += [-d, +d]
    for offset in offsets:
        idx = slice_idx + offset
        if not (0 <= idx < total_n):
            continue
        try:
            ant_axis, ant_sign, left_col, right_col = _psoas_at(
                total_data, idx, axial_axis, sc_center)
            if offset != 0:
                print(f"    Psoas not found at slice {slice_idx}; "
                      f"using direction from slice {idx} (offset {offset:+d})")
            return ant_axis, ant_sign, left_col, right_col, idx
        except ValueError:
            continue
    raise ValueError(
        f"Iliopsoas (labels 88/89) not found within "
        f"±{search_radius} slices of slice {slice_idx}."
    )

# ─────────────────────────────────────────────────────────────────────────────
# Bounding box & label map
# ─────────────────────────────────────────────────────────────────────────────

def create_bounding_box_mask(ct_slice, sc_center, ant_axis, ant_sign,
                              bbox_ant_mm, pixdim, bbox_lr_mm=100.0):
    sc_row, sc_col = sc_center
    nrows, ncols   = ct_slice.shape
    row_mm, col_mm = pixdim[0], pixdim[1]
    rows_idx, cols_idx = np.mgrid[0:nrows, 0:ncols]
    if ant_axis == 0:
        ant_px = bbox_ant_mm / row_mm
        lr_px  = bbox_lr_mm  / col_mm
        d_ant  = (rows_idx - sc_row) * ant_sign
        d_lr   = cols_idx - sc_col
    else:
        ant_px = bbox_ant_mm / col_mm
        lr_px  = bbox_lr_mm  / row_mm
        d_ant  = (cols_idx - sc_col) * ant_sign
        d_lr   = rows_idx - sc_row
    ant_mask  = (d_ant >= 0) & (d_ant <= ant_px)
    lr_mask   = np.abs(d_lr) <= lr_px
    body_mask = ct_slice >= BODY_HU_THRESHOLD
    post_mask = (d_ant < 0) & body_mask
    return lr_mask & (ant_mask | post_mask)


def build_label_map(ct_slice, total_slice, tissue_slice, bbox_mask):
    label_map = np.zeros_like(total_slice, dtype=np.int16)
    merged = (np.isin(total_slice,  TOTAL_MERGE_LABELS) |
              np.isin(tissue_slice, TISSUE_MERGE_LABELS)) & bbox_mask
    retain = (tissue_slice == TISSUE_RETAIN_LABEL) & bbox_mask
    label_map[merged] = NEW_MERGED_LABEL
    label_map[retain] = TISSUE_RETAIN_LABEL
    return label_map


def split_merged_label_by_hu(ct_slice, label_map):
    out       = label_map.copy()
    merged_px = label_map == NEW_MERGED_LABEL
    out[merged_px & (ct_slice <  HU_FAT_THRESHOLD_HIGH) &
                    (ct_slice >  HU_FAT_THRESHOLD_LOW)]  = FAT_SUB_LABEL
    out[merged_px & (ct_slice >= HU_FAT_THRESHOLD_HIGH)] = MUSCLE_SUB_LABEL
    return out

# ─────────────────────────────────────────────────────────────────────────────
# Measurement functions
# ─────────────────────────────────────────────────────────────────────────────

def _hu_stats(hu_values):
    """Return dict of mean/min/max/median for a 1-D array of HU values."""
    return {
        "mean_hu":   round(float(hu_values.mean()),              2),
        "min_hu":    round(float(hu_values.min()),               2),
        "max_hu":    round(float(hu_values.max()),               2),
        "median_hu": round(float(np.median(hu_values)),          2),
    }


def measure_l3_slice(ct_slice, label_map, pixdim, axial_axis):
    """
    Compute cross-sectional area (cm²) and mean/min/max/median HU for each
    output soft-tissue label on a single 2-D slice.

    Returns dict with keys:
        l3_<name>_area_cm2
        l3_<name>_mean_hu
        l3_<name>_min_hu
        l3_<name>_max_hu
        l3_<name>_median_hu
    for each name in LABEL_NAMES.  Missing labels get NaN.
    """
    # In-plane pixel area in cm²
    dims = [pixdim[i] for i in range(3) if i != axial_axis]
    pixel_area_cm2 = (dims[0] / 10.0) * (dims[1] / 10.0)

    stats = {}
    for lbl in OUTPUT_LABELS:
        name = LABEL_NAMES[lbl]
        mask = label_map == lbl
        n_px = int(mask.sum())
        if n_px > 0:
            hu_vals = ct_slice[mask]
            stats[f"l3_{name}_area_cm2"] = round(n_px * pixel_area_cm2, 4)
            for stat_key, stat_val in _hu_stats(hu_vals).items():
                stats[f"l3_{name}_{stat_key}"] = stat_val
        else:
            stats[f"l3_{name}_area_cm2"]   = np.nan
            stats[f"l3_{name}_mean_hu"]    = np.nan
            stats[f"l3_{name}_min_hu"]     = np.nan
            stats[f"l3_{name}_max_hu"]     = np.nan
            stats[f"l3_{name}_median_hu"]  = np.nan
    return stats


def measure_vertebra_l3_slice(ct_slice, total_slice, pixdim, axial_axis):
    """
    Compute cross-sectional area (cm²) and mean/min/max/median HU for the
    L3 vertebral body (tissue label 29) on the L3 axial slice.

    Returns dict with keys:
        l3_vertebra_area_cm2
        l3_vertebra_mean_hu
        l3_vertebra_min_hu
        l3_vertebra_max_hu
        l3_vertebra_median_hu

    All values are NaN when the L3 vertebra label is absent from the slice
    (e.g. the l3_slice_index points to a gap between labels).
    """
    dims = [pixdim[i] for i in range(3) if i != axial_axis]
    pixel_area_cm2 = (dims[0] / 10.0) * (dims[1] / 10.0)

    mask = total_slice == L3_VERTEBRA_TISSUE_LABEL
    n_px = int(mask.sum())

    if n_px > 0:
        hu_vals = ct_slice[mask]
        result = {"l3_vertebra_area_cm2": round(n_px * pixel_area_cm2, 4)}
        for stat_key, stat_val in _hu_stats(hu_vals).items():
            result[f"l3_vertebra_{stat_key}"] = stat_val
        print(f"  L3 vertebra: {n_px} px  area={result['l3_vertebra_area_cm2']} cm²  "
              f"mean={result['l3_vertebra_mean_hu']} HU  "
              f"min={result['l3_vertebra_min_hu']} HU  "
              f"max={result['l3_vertebra_max_hu']} HU  "
              f"median={result['l3_vertebra_median_hu']} HU")
    else:
        print(f"  WARNING: L3 vertebra label ({L3_VERTEBRA_TISSUE_LABEL}) absent "
              f"at this slice. All vertebra stats will be NaN.")
        result = {
            "l3_vertebra_area_cm2":  np.nan,
            "l3_vertebra_mean_hu":   np.nan,
            "l3_vertebra_min_hu":    np.nan,
            "l3_vertebra_max_hu":    np.nan,
            "l3_vertebra_median_hu": np.nan,
        }
    return result


def measure_volume(ct_volume, label_volume, pixdim):
    """
    Compute volume (cm³) and mean HU for each output label across a 3-D
    label volume (already cropped to L1-L5).

    pixdim : (x_mm, y_mm, z_mm) of the cropped volume.

    Returns dict with keys:
        vol_<name>_volume_cm3
        vol_<name>_mean_hu
    for each label in OUTPUT_LABELS.  Missing labels get NaN.
    """
    voxel_vol_cm3 = (pixdim[0] / 10.0) * (pixdim[1] / 10.0) * (pixdim[2] / 10.0)

    stats = {}
    for lbl in OUTPUT_LABELS:
        name = LABEL_NAMES[lbl]
        mask = label_volume == lbl
        n_vx = int(mask.sum())
        if n_vx > 0:
            stats[f"vol_{name}_volume_cm3"] = round(n_vx * voxel_vol_cm3, 4)
            stats[f"vol_{name}_mean_hu"]    = round(float(ct_volume[mask].mean()), 2)
        else:
            stats[f"vol_{name}_volume_cm3"] = np.nan
            stats[f"vol_{name}_mean_hu"]    = np.nan
    return stats

# ─────────────────────────────────────────────────────────────────────────────
# DICOM helpers
# ─────────────────────────────────────────────────────────────────────────────

def affine_to_dicom_geometry(affine, axial_axis):
    M         = affine[:3, :3]
    col_norms = np.linalg.norm(M, axis=0)
    plane_axes     = [ax for ax in range(3) if ax != axial_axis]
    row_ax, col_ax = plane_axes[0], plane_axes[1]
    row_cos    = (M[:, row_ax] / col_norms[row_ax]).tolist()
    col_cos    = (M[:, col_ax] / col_norms[col_ax]).tolist()
    voxel_size = (col_norms[row_ax], col_norms[col_ax], col_norms[axial_axis])
    return row_cos + col_cos, voxel_size, affine[:3, 3].tolist()


def slice_image_position(affine, slice_idx, axial_axis):
    vox = np.zeros(3)
    vox[axial_axis] = slice_idx
    return (affine[:3, :3] @ vox + affine[:3, 3]).tolist()


def build_dicom(pixel_array, series_uid, instance_number,
                series_description, series_number,
                image_position, image_orientation,
                pixel_spacing, slice_thickness,
                orig_meta=None):
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
    if not hasattr(ds, "Modality"): ds.Modality = "CT"
    ds.ImagePositionPatient    = [f"{v:.6f}" for v in image_position]
    ds.ImageOrientationPatient = [f"{v:.6f}" for v in image_orientation]
    ds.PixelSpacing            = [f"{pixel_spacing[0]:.6f}", f"{pixel_spacing[1]:.6f}"]
    ds.SliceThickness          = f"{slice_thickness:.4f}"
    if not hasattr(ds, "RescaleIntercept"): ds.RescaleIntercept = "0"
    if not hasattr(ds, "RescaleSlope"):     ds.RescaleSlope     = "1"
    ds.SamplesPerPixel           = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.Rows, ds.Columns          = arr.shape
    ds.BitsAllocated             = 16
    ds.BitsStored                = 16
    ds.HighBit                   = 15
    ds.PixelRepresentation       = 1
    ds.PixelData                 = arr.tobytes()
    return ds


def save_dicom(ds, path):
    pydicom.dcmwrite(path, ds)
    print(f"  Saved DICOM : {path}")


def try_load_dicom_meta(nifti_path):
    d = os.path.dirname(nifti_path)
    if os.path.isdir(d):
        for f in sorted(os.listdir(d)):
            if f.lower().endswith(".dcm"):
                try:
                    return pydicom.dcmread(os.path.join(d, f), stop_before_pixels=True)
                except Exception:
                    pass
    return None

# ─────────────────────────────────────────────────────────────────────────────
# PNG helpers
# ─────────────────────────────────────────────────────────────────────────────

_CMAP = mcolors.ListedColormap(["black", "cyan", "orange", "yellow", "red"])
_NORM = mcolors.BoundaryNorm([0, 0.5, 1.5, 2.5, 3.5, 4.5], 5)


def _window(arr, wc=40, ww=400):
    lo = wc - ww / 2
    hi = wc + ww / 2
    return np.clip((arr.astype(float) - lo) / (hi - lo), 0, 1)


def save_png_ct(ct_slice, path):
    fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
    ax.imshow(_window(ct_slice), cmap="gray", origin="upper")
    ax.set_title("L3 CT Slice (soft-tissue window)", fontsize=10)
    ax.axis("off")
    fig.tight_layout(pad=0)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved PNG   : {path}")


def save_png_label(label_map, ct_slice, path):
    bg = _window(ct_slice, wc=-20, ww=400)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), dpi=150)
    axes[0].imshow(bg, cmap="gray", origin="upper")
    axes[0].set_title("CT slice", fontsize=9)
    axes[0].axis("off")
    axes[1].imshow(bg, cmap="gray", origin="upper")
    rgba = _CMAP(_NORM(label_map))
    rgba[label_map == 0, 3] = 0
    axes[1].imshow(rgba, origin="upper", alpha=0.65)
    axes[1].set_title("Label overlay", fontsize=9)
    axes[1].axis("off")
    from matplotlib.patches import Patch
    legend = [
        Patch(facecolor="cyan",   label="1 – subcutaneous fat"),
        Patch(facecolor="yellow", label="3 – intramuscular adipose tissue (< -30 HU)"),
        Patch(facecolor="red",    label="4 – muscle (>= -30 HU)"),
    ]
    axes[1].legend(handles=legend, loc="lower right", fontsize=7, framealpha=0.7)
    fig.tight_layout(pad=0.5)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved PNG   : {path}")

# ─────────────────────────────────────────────────────────────────────────────
# L1-L5 volumetric processing  →  NIfTI output
# ─────────────────────────────────────────────────────────────────────────────

def find_vertebral_slice_range(total_data, axial_axis):
    mask = np.isin(total_data, list(VERTEBRA_LABELS.values()))
    if not mask.any():
        return None
    if axial_axis == 0:
        occupied = np.where(mask.any(axis=(1, 2)))[0]
    elif axial_axis == 1:
        occupied = np.where(mask.any(axis=(0, 2)))[0]
    else:
        occupied = np.where(mask.any(axis=(0, 1)))[0]
    return int(occupied.min()), int(occupied.max())


def process_volume_l1_l5(ct_data, total_data, tissue_data,
                          affine, pixdim, axial_axis, bbox_ant_mm, bbox_lr_mm,
                          out_dir, stem, orig_meta,
                          orientation_6, voxel_size, alt_sc_data=None,
                          l3_ant_axis=None, l3_ant_sign=None):
    """
    Process every axial slice in the L1-L5 span.
    Uses psoas direction from L3 slice for consistency; SC re-detected per slice.
    Saves two NIfTI files and returns volumetric measurement dict.

    alt_sc_data: Optional alternative spinal cord segmentation volume
    l3_ant_axis: Anterior axis determined from L3 slice (0 or 1)
    l3_ant_sign: Anterior sign determined from L3 slice (+1 or -1)
    """
    slc_range = find_vertebral_slice_range(total_data, axial_axis)
    if slc_range is None:
        print("  WARNING: No L1-L5 vertebral labels found. Skipping volumetric export.")
        return {f"vol_{LABEL_NAMES[l]}_volume_cm3": np.nan for l in OUTPUT_LABELS} | \
               {f"vol_{LABEL_NAMES[l]}_mean_hu":    np.nan for l in OUTPUT_LABELS}

    slc_min, slc_max = slc_range
    print(f"  L1-L5 slice range: {slc_min}-{slc_max}  ({slc_max-slc_min+1} slices)")

    # ── Use L3's psoas direction if provided ──────────────────────────────────
    if l3_ant_axis is not None and l3_ant_sign is not None:
        cached_ant_axis = l3_ant_axis
        cached_ant_sign = l3_ant_sign
        axis_name = "row" if cached_ant_axis == 0 else "col"
        print(f"  Using psoas direction from L3 slice: "
              f"axis={axis_name} sign={cached_ant_sign:+d}  "
              f"(applied to all L1-L5 slices)")
    else:
        # ── Fallback: Lock psoas direction from first detectable slice ───────
        print("  WARNING: L3 psoas direction not provided. Detecting from L1-L5 range...")
        cached_ant_axis = None
        cached_ant_sign = None
        for probe_idx in range(slc_min, slc_max + 1):
            try:
                probe_sc, _ = find_spinal_cord_center_with_fallback(
                    total_data, probe_idx, axial_axis, alt_sc_data=alt_sc_data)
                cached_ant_axis, cached_ant_sign, _, _, src = \
                    find_psoas_direction_with_fallback(
                        total_data, probe_idx, axial_axis, probe_sc)
                axis_name = "row" if cached_ant_axis == 0 else "col"
                print(f"  Psoas direction locked from slice {src}: "
                      f"axis={axis_name} sign={cached_ant_sign:+d}  "
                      f"(reused for all L1-L5 slices)")
                break
            except ValueError:
                continue

    if cached_ant_axis is None:
        print("  WARNING: Could not determine psoas direction in any L1-L5 slice. "
              "All label maps will be empty.")

    # ── Per-slice processing ──────────────────────────────────────────────────
    ct_vol  = np.zeros_like(ct_data, dtype=np.int16)
    lbl_vol = np.zeros_like(ct_data, dtype=np.int16)

    for idx in range(slc_min, slc_max + 1):
        ct_sl     = extract_slice_2d(ct_data,     idx, axial_axis).astype(np.float32)
        total_sl  = extract_slice_2d(total_data,  idx, axial_axis).astype(np.int32)
        tissue_sl = extract_slice_2d(tissue_data, idx, axial_axis).astype(np.int32)

        if cached_ant_axis is None:
            label_map = np.zeros(ct_sl.shape, dtype=np.int16)
        else:
            try:
                sc_ctr, _ = find_spinal_cord_center_with_fallback(
                    total_data, idx, axial_axis, alt_sc_data=alt_sc_data)
                bbox_mask = create_bounding_box_mask(
                    ct_sl, sc_ctr, cached_ant_axis, cached_ant_sign,
                    bbox_ant_mm, pixdim, bbox_lr_mm=bbox_lr_mm)
                label_map = build_label_map(ct_sl, total_sl, tissue_sl, bbox_mask)
                label_map = split_merged_label_by_hu(ct_sl, label_map)
            except ValueError as e:
                print(f"    Slice {idx}: SC fallback exhausted – zero label map ({e})")
                label_map = np.zeros(ct_sl.shape, dtype=np.int16)

        if axial_axis == 0:
            ct_vol[idx, :, :]  = ct_sl.astype(np.int16)
            lbl_vol[idx, :, :] = label_map
        elif axial_axis == 1:
            ct_vol[:, idx, :]  = ct_sl.astype(np.int16)
            lbl_vol[:, idx, :] = label_map
        else:
            ct_vol[:, :, idx]  = ct_sl.astype(np.int16)
            lbl_vol[:, :, idx] = label_map

    # ── Crop to L1-L5 range ───────────────────────────────────────────────────
    if axial_axis == 0:
        ct_crop  = ct_vol[slc_min:slc_max + 1, :, :]
        lbl_crop = lbl_vol[slc_min:slc_max + 1, :, :]
    elif axial_axis == 1:
        ct_crop  = ct_vol[:, slc_min:slc_max + 1, :]
        lbl_crop = lbl_vol[:, slc_min:slc_max + 1, :]
    else:
        ct_crop  = ct_vol[:, :, slc_min:slc_max + 1]
        lbl_crop = lbl_vol[:, :, slc_min:slc_max + 1]

    # Shift affine origin to first included slice
    cropped_affine = affine.copy()
    shift_vox = np.zeros(3)
    shift_vox[axial_axis] = slc_min
    cropped_affine[:3, 3] = affine[:3, :3] @ shift_vox + affine[:3, 3]

    ct_nii_path  = os.path.join(out_dir, f"{stem}_l1_l5_image.nii.gz")
    lbl_nii_path = os.path.join(out_dir, f"{stem}_l1_l5_label.nii.gz")
    nib.save(nib.Nifti1Image(ct_crop,  cropped_affine), ct_nii_path)
    print(f"  Saved NIfTI CT   : {ct_nii_path}  shape={ct_crop.shape}")
    nib.save(nib.Nifti1Image(lbl_crop, cropped_affine), lbl_nii_path)
    print(f"  Saved NIfTI Label: {lbl_nii_path}  shape={lbl_crop.shape}")

    # ── Volumetric measurements ───────────────────────────────────────────────
    vol_stats = measure_volume(ct_crop.astype(np.float32), lbl_crop, pixdim)
    for k, v in vol_stats.items():
        print(f"  {k}: {v}")
    return vol_stats

# ─────────────────────────────────────────────────────────────────────────────
# Per-row entry point
# ─────────────────────────────────────────────────────────────────────────────

def process_row(row, bbox_ant_mm, bbox_lr_mm, base_output_dir, sc_seg_path=None):
    """
    Process a single row from the input CSV.

    sc_seg_path: Optional path to alternative spinal cord segmentation
    """
    slice_idx   = int(row["l3_slice_index"])
    nifti_path  = str(row["nifti_path"]).strip()
    total_path  = str(row["ts_total_path"]).strip()
    tissue_path = str(row["ts_tissue_path"]).strip()
    accession   = str(row["accession"]).strip()

    stem    = nifti_stem(nifti_path)
    out_dir = os.path.join(base_output_dir, accession)
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"NIfTI stem : {stem}")
    print(f"Accession  : {accession}")
    print(f"Output     : {out_dir}")
    if sc_seg_path:
        print(f"Alt SC seg : {sc_seg_path}")
    print(f"{'='*60}")

    ct_data,     affine, _, pixdim = load_nifti(nifti_path)
    total_data,  *_                = load_nifti(total_path)
    tissue_data, *_                = load_nifti(tissue_path)

    # Load alternative SC segmentation if provided
    alt_sc_data = None
    if sc_seg_path and os.path.exists(sc_seg_path):
        print(f"  Loading alternative SC segmentation from {sc_seg_path}")
        alt_sc_data_raw, alt_affine, alt_header, alt_pixdim = load_nifti(sc_seg_path)
        print(f"  Alternative SC raw dtype: {alt_sc_data_raw.dtype}, shape: {alt_sc_data_raw.shape}")

        # Convert to int32 to ensure proper label comparison
        alt_sc_data = np.round(alt_sc_data_raw).astype(np.int32)
        print(f"  Alternative SC converted dtype: {alt_sc_data.dtype}")

        unique_labels = np.unique(alt_sc_data[alt_sc_data > 0])
        print(f"  Alternative SC seg unique labels: {unique_labels} (dtype: {unique_labels.dtype})")
        print(f"  Using label {SPINAL_CORD_ALT_LABEL} for spinal cord")
        print(f"  Total non-zero voxels in alt SC: {np.sum(alt_sc_data > 0)}")

    axial_axis = get_axial_axis(affine)
    print(f"  Axial axis: {axial_axis} | L3 slice: {slice_idx} | Voxel mm: {pixdim}")

    orig_meta = try_load_dicom_meta(nifti_path)
    orientation_6, voxel_size, _ = affine_to_dicom_geometry(affine, axial_axis)
    pixel_spacing = (voxel_size[0], voxel_size[1])
    slice_thick   = voxel_size[2]

    # ── L3 slice ──────────────────────────────────────────────────────────────
    ct_sl     = extract_slice_2d(ct_data,     slice_idx, axial_axis).astype(np.float32)
    total_sl  = extract_slice_2d(total_data,  slice_idx, axial_axis).astype(np.int32)
    tissue_sl = extract_slice_2d(tissue_data, slice_idx, axial_axis).astype(np.int32)

    sc_ctr, sc_src = find_spinal_cord_center_with_fallback(
        total_data, slice_idx, axial_axis, alt_sc_data=alt_sc_data)
    print(f"  Spinal cord: row={sc_ctr[0]:.1f} col={sc_ctr[1]:.1f}"
          + (f"  [fallback slice {sc_src}]" if sc_src != slice_idx else ""))

    ant_axis, ant_sign, left_col, right_col, psoas_src = \
        find_psoas_direction_with_fallback(
            total_data, slice_idx, axial_axis, sc_ctr)
    axis_name = "row" if ant_axis == 0 else "col"
    print(f"  Anterior axis: {axis_name} (sign={ant_sign:+d}) "
          f"| psoas L={left_col:.1f} R={right_col:.1f}"
          + (f"  [fallback slice {psoas_src}]" if psoas_src != slice_idx else ""))

    bbox_mask = create_bounding_box_mask(
        ct_sl, sc_ctr, ant_axis, ant_sign, bbox_ant_mm, pixdim,
        bbox_lr_mm=bbox_lr_mm)
    label_map = build_label_map(ct_sl, total_sl, tissue_sl, bbox_mask)
    label_map = split_merged_label_by_hu(ct_sl, label_map)
    print(f"  L3 labels present: {np.unique(label_map[label_map > 0])}")

    # ── L3 measurements (soft tissue) ─────────────────────────────────────────
    l3_stats = measure_l3_slice(ct_sl, label_map, pixdim, axial_axis)
    for k, v in l3_stats.items():
        print(f"  {k}: {v}")

    # ── L3 measurements (vertebra) ────────────────────────────────────────────
    vertebra_stats = measure_vertebra_l3_slice(ct_sl, total_sl, pixdim, axial_axis)

    # ── L3 DICOM + PNG ────────────────────────────────────────────────────────
    img_pos_l3 = slice_image_position(affine, slice_idx, axial_axis)
    geom_kw = dict(image_position=img_pos_l3, image_orientation=orientation_6,
                   pixel_spacing=pixel_spacing, slice_thickness=slice_thick,
                   orig_meta=orig_meta)

    ds_img = build_dicom(ct_sl.astype(np.int16), generate_uid(), 1,
                         f"L3 CT [{stem}]", 1, **geom_kw)
    save_dicom(ds_img, os.path.join(out_dir, f"{stem}_l3_image.dcm"))

    ds_lbl = build_dicom(label_map, generate_uid(), 1,
                         f"L3 Label [{stem}]", 2, **geom_kw)
    save_dicom(ds_lbl, os.path.join(out_dir, f"{stem}_l3_label.dcm"))

    save_png_ct(ct_sl,    os.path.join(out_dir, f"{stem}_l3_image.png"))
    save_png_label(label_map, ct_sl,
                   os.path.join(out_dir, f"{stem}_l3_label.png"))

    # ── L1-L5 volume ──────────────────────────────────────────────────────────
    print("\n  Processing L1-L5 volume ...")
    vol_stats = process_volume_l1_l5(
        ct_data, total_data, tissue_data,
        affine, pixdim, axial_axis, bbox_ant_mm, bbox_lr_mm,
        out_dir, stem, orig_meta, orientation_6, voxel_size,
        alt_sc_data=alt_sc_data,
        l3_ant_axis=ant_axis, l3_ant_sign=ant_sign,
    )

    # Return combined stats for this row
    return {**l3_stats, **vertebra_stats, **vol_stats}

# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="L3 / L1-L5 segmentation -> DICOM + NIfTI + PNG + measurements")
    parser.add_argument("--input",
                        default="/project/hipaa_shinjohnlab/Project56/l3_label_summary.csv",
                        help="Path to input CSV file")
    parser.add_argument("--filter_csv",
                        default="/project/hipaa_shinjohnlab/Project56/nifti_seg_path.csv",
                        help="Optional CSV file to filter which nifti_paths to process "
                             "(must have 'nifti_path' and 'seg_path' columns)")
    parser.add_argument("--bbox_ant_mm", type=float, default=55.0,
                        help="Anterior bounding box extent in mm (default: 55)")
    parser.add_argument("--bbox_lr_mm",  type=float, default=85.0,
                        help="Left/right bounding box extent in mm (default: 85)")
    parser.add_argument("--output",  default="./l3_output3",
                        help="Base output directory (default: ./l3_output)")
    parser.add_argument("--limit",   default=-1, type=int,
                        help="Max rows to process (-1 = all)")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f"ERROR: CSV not found: {args.input}")

    df = pd.read_csv(args.input)
    required = {"l3_slice_index", "nifti_path", "ts_total_path",
                "ts_tissue_path", "accession"}
    missing  = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: CSV missing columns: {missing}")

    # Load and process filter CSV if provided
    filter_map = {}  # nifti_path -> seg_path
    if args.filter_csv:
        if not os.path.isfile(args.filter_csv):
            sys.exit(f"ERROR: Filter CSV not found: {args.filter_csv}")

        filter_df = pd.read_csv(args.filter_csv)
        filter_required = {"nifti_path", "seg_path"}
        filter_missing = filter_required - set(filter_df.columns)
        if filter_missing:
            sys.exit(f"ERROR: Filter CSV missing columns: {filter_missing}")

        # Create mapping of nifti_path to seg_path
        for _, row in filter_df.iterrows():
            nifti = str(row["nifti_path"]).strip()
            seg = str(row["seg_path"]).strip()
            filter_map[nifti] = seg

        print(f"Filter CSV loaded: {len(filter_map)} entries")

        # Filter main dataframe to only include rows with matching nifti_paths
        df = df[df["nifti_path"].apply(lambda x: str(x).strip() in filter_map)].copy()
        print(f"Filtered to {len(df)} matching rows")

        if len(df) == 0:
            sys.exit("ERROR: No matching rows found after filtering")

    if args.limit > 0:
        df = df.iloc[:args.limit]

    os.makedirs(args.output, exist_ok=True)
    print(f"Rows        : {len(df)}")
    print(f"Bbox ant mm : {args.bbox_ant_mm}  |  Bbox L/R mm: {args.bbox_lr_mm}")
    print(f"Output root : {args.output}\n")

    # Collect per-row measurement dicts; None on failure
    all_stats = []
    errors    = []

    for i, row in df.iterrows():
        try:
            # Get seg_path if this nifti_path is in filter map
            nifti = str(row["nifti_path"]).strip()
            seg_path = filter_map.get(nifti, None)

            stats = process_row(row, args.bbox_ant_mm, args.bbox_lr_mm,
                                args.output, sc_seg_path=seg_path)
            all_stats.append(stats)
        except Exception as e:
            import traceback
            print(f"\n  !! ERROR on row {i}: {e}")
            traceback.print_exc()
            errors.append((i, str(e)))
            all_stats.append(None)   # keep alignment with df rows

    # ── Build results table ───────────────────────────────────────────────────
    # Align stats list with the (possibly filtered) df index
    stats_rows = []
    for i in range(len(df)):
        stats_rows.append(all_stats[i] if all_stats[i] is not None else {})

    stats_df = pd.DataFrame(stats_rows, index=df.index)
    result_df = pd.concat([df, stats_df], axis=1)

    out_csv = os.path.join(args.output, "results_summary2.csv")
    result_df.to_csv(out_csv, index=False)
    print(f"\n  Results table saved: {out_csv}")

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"Completed: {len(df) - len(errors)}/{len(df)} rows succeeded.")
    if errors:
        print("Failed rows:")
        for idx, msg in errors:
            print(f"  Row {idx}: {msg}")

    print("\nMeasurement columns added:")
    for col in stats_df.columns:
        print(f"  {col}")


if __name__ == "__main__":
    main()