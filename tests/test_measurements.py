"""Tests for L3-slice and L1-L5 volume measurements."""

import numpy as np

from asd_bodycomp.measurements import (
    measure_l3_slice,
    measure_vertebra_l3_slice,
    measure_volume,
)


def test_measure_l3_slice_area_matches_pixel_count():
    label_map = np.zeros((10, 10), dtype=np.int16)
    label_map[0:2, 0:5] = 1            # 10 pixels of subq_fat
    label_map[3:5, 0:5] = 4            # 10 pixels of muscle
    ct = np.full((10, 10), 50, dtype=np.float32)
    ct[label_map == 1] = -100
    pixdim = (2.0, 2.0, 3.0)
    stats = measure_l3_slice(ct, label_map, pixdim, axial_axis=2)
    # 10 pixels × 0.2 cm × 0.2 cm = 0.4 cm²
    assert stats["l3_subq_fat_area_cm2"] == 0.4
    assert stats["l3_muscle_area_cm2"]   == 0.4
    assert np.isnan(stats["l3_imat_area_cm2"])
    assert stats["l3_subq_fat_mean_hu"]  == -100.0
    assert stats["l3_muscle_mean_hu"]    == 50.0


def test_measure_vertebra_l3_returns_nan_when_label_absent():
    total = np.zeros((10, 10), dtype=np.int32)
    ct    = np.zeros((10, 10), dtype=np.float32)
    stats = measure_vertebra_l3_slice(ct, total, (2.0, 2.0, 3.0), 2)
    assert np.isnan(stats["l3_vertebra_area_cm2"])


def test_measure_vertebra_l3_present():
    total = np.zeros((10, 10), dtype=np.int32)
    total[3:7, 3:7] = 29               # 16 pixels of L3 vertebra
    ct    = np.full((10, 10), 200.0, dtype=np.float32)
    stats = measure_vertebra_l3_slice(ct, total, (2.0, 2.0, 3.0), 2)
    # 16 pixels × 0.2 × 0.2 = 0.64 cm²
    assert stats["l3_vertebra_area_cm2"] == 0.64
    assert stats["l3_vertebra_mean_hu"]  == 200.0


def test_measure_volume_uses_voxel_size():
    label = np.zeros((4, 4, 5), dtype=np.int16)
    label[..., 0] = 1                   # 16 voxels of subq_fat in slice 0
    ct = np.full((4, 4, 5), -100.0, dtype=np.float32)
    pixdim = (2.0, 2.0, 3.0)            # voxel = 0.2 × 0.2 × 0.3 cm³ = 0.012 cm³
    stats  = measure_volume(ct, label, pixdim)
    assert round(stats["vol_subq_fat_volume_cm3"], 4) == round(16 * 0.012, 4)
    assert stats["vol_subq_fat_mean_hu"] == -100.0
    assert np.isnan(stats["vol_muscle_volume_cm3"])
