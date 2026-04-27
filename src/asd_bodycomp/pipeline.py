"""Per-row pipeline shared by scripts 5 and 6.

`process_row` reads one row of the L3 manifest (built by step 4), produces:
  - L3 DICOM image + L3 DICOM label
  - L3 PNG image + L3 PNG label-overlay
  - L1-L5 cropped CT NIfTI + L1-L5 cropped label NIfTI
  - dict of per-row measurements that the caller appends to a CSV

The optional `alt_sc_path` parameter switches script 6's behavior: when given,
spinal-cord detection uses the user-supplied SC-only segmentation (label 1)
instead of label 79 in the TotalSegmentator total volume.
"""

from __future__ import annotations

import os
from typing import Optional

import numpy as np
import nibabel as nib
from pydicom.uid import generate_uid

from .dicom_export import (
    affine_to_dicom_geometry,
    build_dicom,
    save_dicom,
    slice_image_position,
    try_load_dicom_meta,
)
from .labels import LABEL_NAMES, OUTPUT_LABELS, VERTEBRA_LABELS
from .measurements import (
    measure_l3_slice,
    measure_vertebra_l3_slice,
    measure_volume,
)
from .nifti_io import (
    extract_slice_2d,
    get_axial_axis,
    load_nifti,
    nifti_stem,
)
from .png_export import save_png_ct, save_png_label
from .slice_processing import (
    build_label_map,
    create_bounding_box_mask,
    find_psoas_direction_with_fallback,
    find_spinal_cord_center_with_fallback,
    split_merged_label_by_hu,
)


# ── L1-L5 vertebral span ────────────────────────────────────────────────────

def _find_vertebral_slice_range(total_data, axial_axis):
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


def _process_volume_l1_l5(ct_data, total_data, tissue_data, affine, pixdim,
                          axial_axis, bbox_ant_mm, bbox_lr_mm,
                          out_dir, stem,
                          alt_sc_data=None,
                          l3_ant_axis=None, l3_ant_sign=None):
    """Process every axial slice in the L1-L5 span and write cropped NIfTI outputs."""
    slc_range = _find_vertebral_slice_range(total_data, axial_axis)
    if slc_range is None:
        print("  WARNING: No L1-L5 vertebral labels found; skipping volumetric export.")
        return {f"vol_{LABEL_NAMES[l]}_volume_cm3": np.nan for l in OUTPUT_LABELS} | \
               {f"vol_{LABEL_NAMES[l]}_mean_hu":    np.nan for l in OUTPUT_LABELS}

    slc_min, slc_max = slc_range
    print(f"  L1-L5 slice range: {slc_min}-{slc_max}  ({slc_max-slc_min+1} slices)")

    if l3_ant_axis is not None and l3_ant_sign is not None:
        cached_ant_axis, cached_ant_sign = l3_ant_axis, l3_ant_sign
    else:
        cached_ant_axis = cached_ant_sign = None
        for probe_idx in range(slc_min, slc_max + 1):
            try:
                probe_sc, _ = find_spinal_cord_center_with_fallback(
                    total_data, probe_idx, axial_axis, alt_sc_data=alt_sc_data)
                cached_ant_axis, cached_ant_sign, *_ = \
                    find_psoas_direction_with_fallback(
                        total_data, probe_idx, axial_axis, probe_sc)
                break
            except ValueError:
                continue

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
                label_map = build_label_map(total_sl, tissue_sl, bbox_mask)
                label_map = split_merged_label_by_hu(ct_sl, label_map)
            except ValueError:
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

    if axial_axis == 0:
        ct_crop  = ct_vol[slc_min:slc_max + 1, :, :]
        lbl_crop = lbl_vol[slc_min:slc_max + 1, :, :]
    elif axial_axis == 1:
        ct_crop  = ct_vol[:, slc_min:slc_max + 1, :]
        lbl_crop = lbl_vol[:, slc_min:slc_max + 1, :]
    else:
        ct_crop  = ct_vol[:, :, slc_min:slc_max + 1]
        lbl_crop = lbl_vol[:, :, slc_min:slc_max + 1]

    cropped_affine = affine.copy()
    shift_vox = np.zeros(3)
    shift_vox[axial_axis] = slc_min
    cropped_affine[:3, 3] = affine[:3, :3] @ shift_vox + affine[:3, 3]

    nib.save(nib.Nifti1Image(ct_crop,  cropped_affine),
             os.path.join(out_dir, f"{stem}_l1_l5_image.nii.gz"))
    nib.save(nib.Nifti1Image(lbl_crop, cropped_affine),
             os.path.join(out_dir, f"{stem}_l1_l5_label.nii.gz"))

    return measure_volume(ct_crop.astype(np.float32), lbl_crop, pixdim)


# ── per-row entry point ─────────────────────────────────────────────────────

def process_row(row, bbox_ant_mm: float, bbox_lr_mm: float,
                base_output_dir: str,
                alt_sc_path: Optional[str] = None) -> dict:
    """Run the full L3 + L1-L5 pipeline for one manifest row.

    Parameters
    ----------
    row : Mapping with keys {l3_slice_index, nifti_path, ts_total_path,
                              ts_tissue_path, accession}.
    bbox_ant_mm, bbox_lr_mm : float
        Anterior and lateral extents of the L3 bounding box, in mm.
    base_output_dir : str
        Outputs are written to <base_output_dir>/<accession>/.
    alt_sc_path : str, optional
        When provided, spinal-cord detection uses this SC-only segmentation
        (label 1) instead of label 79 in the total volume.
    """
    slice_idx   = int(row["l3_slice_index"])
    nifti_path  = str(row["nifti_path"]).strip()
    total_path  = str(row["ts_total_path"]).strip()
    tissue_path = str(row["ts_tissue_path"]).strip()
    accession   = str(row["accession"]).strip()

    stem    = nifti_stem(nifti_path)
    out_dir = os.path.join(base_output_dir, accession)
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n{'='*60}\n{accession}/{stem}  ->  {out_dir}\n{'='*60}")

    ct_data,    affine, _, pixdim = load_nifti(nifti_path)
    total_data,  *_               = load_nifti(total_path)
    tissue_data, *_               = load_nifti(tissue_path)

    alt_sc_data = None
    if alt_sc_path and os.path.exists(alt_sc_path):
        alt_raw, *_ = load_nifti(alt_sc_path)
        alt_sc_data = np.round(alt_raw).astype(np.int32)

    axial_axis = get_axial_axis(affine)
    print(f"  axial_axis={axial_axis}  L3_slice={slice_idx}  voxel_mm={pixdim}")

    orig_meta = try_load_dicom_meta(nifti_path)
    orientation_6, voxel_size, _ = affine_to_dicom_geometry(affine, axial_axis)
    pixel_spacing = (voxel_size[0], voxel_size[1])
    slice_thick   = voxel_size[2]

    # ── L3 slice ────────────────────────────────────────────────────────────
    ct_sl     = extract_slice_2d(ct_data,     slice_idx, axial_axis).astype(np.float32)
    total_sl  = extract_slice_2d(total_data,  slice_idx, axial_axis).astype(np.int32)
    tissue_sl = extract_slice_2d(tissue_data, slice_idx, axial_axis).astype(np.int32)

    sc_ctr, _ = find_spinal_cord_center_with_fallback(
        total_data, slice_idx, axial_axis, alt_sc_data=alt_sc_data)
    ant_axis, ant_sign, *_ = find_psoas_direction_with_fallback(
        total_data, slice_idx, axial_axis, sc_ctr)

    bbox_mask = create_bounding_box_mask(
        ct_sl, sc_ctr, ant_axis, ant_sign, bbox_ant_mm, pixdim,
        bbox_lr_mm=bbox_lr_mm)
    label_map = build_label_map(total_sl, tissue_sl, bbox_mask)
    label_map = split_merged_label_by_hu(ct_sl, label_map)

    l3_stats       = measure_l3_slice(ct_sl, label_map, pixdim, axial_axis)
    vertebra_stats = measure_vertebra_l3_slice(ct_sl, total_sl, pixdim, axial_axis)

    # ── L3 DICOM + PNG ──────────────────────────────────────────────────────
    img_pos_l3 = slice_image_position(affine, slice_idx, axial_axis)
    geom = dict(image_position=img_pos_l3, image_orientation=orientation_6,
                pixel_spacing=pixel_spacing, slice_thickness=slice_thick,
                orig_meta=orig_meta)

    ds_img = build_dicom(ct_sl.astype(np.int16), generate_uid(), 1,
                         f"L3 CT [{stem}]", 1, **geom)
    save_dicom(ds_img, os.path.join(out_dir, f"{stem}_l3_image.dcm"))

    ds_lbl = build_dicom(label_map, generate_uid(), 1,
                         f"L3 Label [{stem}]", 2, **geom)
    save_dicom(ds_lbl, os.path.join(out_dir, f"{stem}_l3_label.dcm"))

    save_png_ct(ct_sl, os.path.join(out_dir, f"{stem}_l3_image.png"))
    save_png_label(label_map, ct_sl,
                   os.path.join(out_dir, f"{stem}_l3_label.png"))

    # ── L1-L5 volume ────────────────────────────────────────────────────────
    vol_stats = _process_volume_l1_l5(
        ct_data, total_data, tissue_data, affine, pixdim, axial_axis,
        bbox_ant_mm, bbox_lr_mm, out_dir, stem,
        alt_sc_data=alt_sc_data,
        l3_ant_axis=ant_axis, l3_ant_sign=ant_sign,
    )

    return {**l3_stats, **vertebra_stats, **vol_stats}


def run_manifest(df, bbox_ant_mm, bbox_lr_mm, output_dir,
                 alt_sc_map: Optional[dict] = None,
                 limit: int = -1):
    """Run process_row for every row in the manifest.

    Returns (results_df_with_stats, errors_list).
    """
    import pandas as pd

    if limit > 0:
        df = df.iloc[:limit]

    os.makedirs(output_dir, exist_ok=True)

    all_stats, errors = [], []
    for i, row in df.iterrows():
        try:
            alt_path = (alt_sc_map or {}).get(str(row["nifti_path"]).strip())
            stats = process_row(row, bbox_ant_mm, bbox_lr_mm,
                                output_dir, alt_sc_path=alt_path)
            all_stats.append(stats)
        except Exception as e:
            import traceback
            traceback.print_exc()
            errors.append((i, str(e)))
            all_stats.append({})

    stats_df  = pd.DataFrame(all_stats, index=df.index)
    return pd.concat([df, stats_df], axis=1), errors
