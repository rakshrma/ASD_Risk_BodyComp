"""Tests for slice-level geometry and label-map building."""

import numpy as np

from asd_bodycomp.nifti_io import extract_slice_2d, get_axial_axis
from asd_bodycomp.slice_processing import (
    build_label_map,
    create_bounding_box_mask,
    find_psoas_direction_with_fallback,
    find_spinal_cord_center_with_fallback,
    split_merged_label_by_hu,
)


def test_axial_axis_for_diagonal_affine(synthetic_arrays):
    *_, affine = synthetic_arrays
    assert get_axial_axis(affine) == 2


def test_find_spinal_cord_in_synthetic(synthetic_arrays):
    _, total, _, _ = synthetic_arrays
    center, src = find_spinal_cord_center_with_fallback(total, 8, 2)
    # SC is at rows 8..10, cols 15..17 → centroid roughly (9, 16)
    assert src == 8
    assert 8 <= center[0] <= 10
    assert 15 <= center[1] <= 17


def test_find_spinal_cord_falls_back_when_missing(synthetic_arrays):
    _, total, _, _ = synthetic_arrays
    # z=2 has no SC; z=4 does. With search_radius=5 we should find it at z=4 or so.
    center, src = find_spinal_cord_center_with_fallback(total, 2, 2)
    assert src in {3, 4, 5, 6, 7}


def test_find_psoas_direction(synthetic_arrays):
    _, total, _, _ = synthetic_arrays
    center, _ = find_spinal_cord_center_with_fallback(total, 8, 2)
    ant_axis, ant_sign, _, _, _ = \
        find_psoas_direction_with_fallback(total, 8, 2, center)
    # Psoas is at rows 13..16 (anterior of SC at rows 8..10) so ant_axis=0, sign=+1
    assert ant_axis == 0
    assert ant_sign == 1


def test_create_bbox_mask_includes_anterior_only(synthetic_arrays):
    ct, total, _, _ = synthetic_arrays
    sc_center = (9, 16)
    pixdim    = (2.0, 2.0, 3.0)
    ct_sl     = extract_slice_2d(ct, 8, 2)
    mask      = create_bounding_box_mask(
        ct_sl, sc_center, ant_axis=0, ant_sign=1,
        bbox_ant_mm=20.0, pixdim=pixdim, bbox_lr_mm=20.0)
    # Anterior slab: rows ≥ sc_row=9 and ≤ 9 + 20/2 = 19. Inside lateral cols.
    assert mask[15, 16] == True   # anterior pixel near center
    # Far posterior pixel outside body should be False (CT < BODY_HU_THRESHOLD)
    assert mask[0, 16] == False


def test_build_label_map_and_split(synthetic_arrays):
    ct, total, tissue, _ = synthetic_arrays
    ct_sl     = extract_slice_2d(ct,     8, 2).astype(np.float32)
    total_sl  = extract_slice_2d(total,  8, 2).astype(np.int32)
    tissue_sl = extract_slice_2d(tissue, 8, 2).astype(np.int32)

    bbox = np.ones_like(total_sl, dtype=bool)
    lm   = build_label_map(total_sl, tissue_sl, bbox)
    # Subq fat pixels become label 1
    assert (lm == 1).any()
    # Psoas / muscle pixels become label 2 (merged)
    assert (lm == 2).any()

    lm2 = split_merged_label_by_hu(ct_sl, lm)
    # After split, label 2 should be gone, replaced by 3 or 4
    assert not (lm2 == 2).any()
    # Muscle pixels (HU ~50, 60) should now be label 4
    assert (lm2 == 4).any()
