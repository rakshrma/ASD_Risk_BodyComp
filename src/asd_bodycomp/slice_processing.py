"""Slice-level geometry: spinal cord, psoas, bbox, label-map building.

The bounding-box logic isolates the L3 abdominal-wall region by:
  1. Centering on the spinal cord.
  2. Determining "anterior" from the iliopsoas centroids.
  3. Keeping pixels within ±bbox_lr_mm laterally and 0..bbox_ant_mm anteriorly,
     plus all body-HU pixels posterior to the cord (back muscles).
"""

from __future__ import annotations

import numpy as np

from .labels import (
    BODY_HU_THRESHOLD,
    FAT_SUB_LABEL,
    HU_FAT_THRESHOLD_HIGH,
    HU_FAT_THRESHOLD_LOW,
    ILIOPSOAS_LEFT_LABEL,
    ILIOPSOAS_RIGHT_LABEL,
    MUSCLE_SUB_LABEL,
    NEW_MERGED_LABEL,
    SC_SEARCH_RADIUS,
    SPINAL_CORD_ALT_LABEL,
    SPINAL_CORD_LABEL,
    TISSUE_MERGE_LABELS,
    TISSUE_RETAIN_LABEL,
    TOTAL_MERGE_LABELS,
)
from .nifti_io import extract_slice_2d, volume_n_slices


# ── spinal cord ─────────────────────────────────────────────────────────────

def _sc_center_at(volume, idx, axis, sc_label):
    """Centroid (row, col) of `sc_label` at slice `idx`, or ValueError."""
    slc = extract_slice_2d(volume, idx, axis).astype(np.int32)
    mask = slc == sc_label
    if not mask.any():
        raise ValueError(f"SC label {sc_label} absent at slice {idx}")
    return slc, mask, np.argwhere(mask).mean(axis=0)


def find_spinal_cord_center_with_fallback(total_data, slice_idx, axial_axis,
                                          alt_sc_data=None,
                                          search_radius=SC_SEARCH_RADIUS):
    """Locate the spinal-cord centroid near `slice_idx`.

    Search order:
      • If `alt_sc_data` is given (script 06): use it exclusively, with label
        `SPINAL_CORD_ALT_LABEL`. Scan ±search_radius slices.
      • Otherwise: use `total_data` with label `SPINAL_CORD_LABEL`.
    Returns ((row, col), source_slice_idx).
    """
    target_data, target_label = (
        (alt_sc_data, SPINAL_CORD_ALT_LABEL)
        if alt_sc_data is not None
        else (total_data, SPINAL_CORD_LABEL)
    )
    n_slices = volume_n_slices(target_data, axial_axis)
    offsets  = [0] + [s for d in range(1, search_radius + 1) for s in (-d, +d)]

    for offset in offsets:
        idx = slice_idx + offset
        if not (0 <= idx < n_slices):
            continue
        try:
            _, _, center = _sc_center_at(target_data, idx, axial_axis, target_label)
            if offset != 0:
                src = "alternative" if alt_sc_data is not None else "primary"
                print(f"    SC: {src} seg used slice {idx} (offset {offset:+d})")
            return center, idx
        except ValueError:
            continue

    src = "alternative" if alt_sc_data is not None else "primary"
    raise ValueError(
        f"Spinal cord (label {target_label}) not found in {src} segmentation "
        f"within ±{search_radius} of slice {slice_idx}."
    )


# ── psoas direction ─────────────────────────────────────────────────────────

def _psoas_at(total_data, idx, axis, sc_center):
    slc = extract_slice_2d(total_data, idx, axis).astype(np.int32)
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
    """Determine the (axis, sign) describing the anterior direction relative to SC.

    The vector from the spinal cord to the psoas centroid points anteriorly;
    we collapse it to one of (±row, ±col) so downstream code can use simple
    inequalities.
    """
    n_slices = volume_n_slices(total_data, axial_axis)
    offsets  = [0] + [s for d in range(1, search_radius + 1) for s in (-d, +d)]

    for offset in offsets:
        idx = slice_idx + offset
        if not (0 <= idx < n_slices):
            continue
        try:
            ant_axis, ant_sign, left_col, right_col = _psoas_at(
                total_data, idx, axial_axis, sc_center)
            if offset != 0:
                print(f"    Psoas: used slice {idx} (offset {offset:+d})")
            return ant_axis, ant_sign, left_col, right_col, idx
        except ValueError:
            continue

    raise ValueError(
        f"Iliopsoas (labels 88/89) not found within "
        f"±{search_radius} slices of slice {slice_idx}."
    )


# ── bbox & label map ────────────────────────────────────────────────────────

def create_bounding_box_mask(ct_slice, sc_center, ant_axis, ant_sign,
                             bbox_ant_mm, pixdim, bbox_lr_mm=100.0):
    """Boolean mask isolating the L3 abdominal cross-section.

    Anterior half:  rectangle of bbox_lr_mm × bbox_ant_mm centered on SC.
    Posterior half: every pixel >= BODY_HU_THRESHOLD (back muscles).
    """
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


def build_label_map(total_slice, tissue_slice, bbox_mask):
    """Combine total + tissue segmentations within the bbox.

    Output labels (intermediate):
      1 = subcutaneous fat (tissue label 1)
      2 = merged muscle/imat region (tissue 3,4 + total 88,89 psoas)
    Label 2 is then split by HU in `split_merged_label_by_hu`.
    """
    label_map = np.zeros_like(total_slice, dtype=np.int16)
    merged = (np.isin(total_slice,  TOTAL_MERGE_LABELS) |
              np.isin(tissue_slice, TISSUE_MERGE_LABELS)) & bbox_mask
    retain = (tissue_slice == TISSUE_RETAIN_LABEL) & bbox_mask
    label_map[merged] = NEW_MERGED_LABEL
    label_map[retain] = TISSUE_RETAIN_LABEL
    return label_map


def split_merged_label_by_hu(ct_slice, label_map):
    """Split label 2 into IMAT (3) and muscle (4) by HU."""
    out = label_map.copy()
    merged_px = label_map == NEW_MERGED_LABEL
    out[merged_px & (ct_slice <  HU_FAT_THRESHOLD_HIGH) &
                    (ct_slice >  HU_FAT_THRESHOLD_LOW)]  = FAT_SUB_LABEL
    out[merged_px & (ct_slice >= HU_FAT_THRESHOLD_HIGH)] = MUSCLE_SUB_LABEL
    return out
