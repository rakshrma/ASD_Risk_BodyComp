"""Per-slice and per-volume body-composition statistics."""

from __future__ import annotations

import numpy as np

from .labels import L3_VERTEBRA_TOTAL_LABEL, LABEL_NAMES, OUTPUT_LABELS


def _hu_stats(hu_values):
    return {
        "mean_hu":   round(float(hu_values.mean()),     2),
        "min_hu":    round(float(hu_values.min()),      2),
        "max_hu":    round(float(hu_values.max()),      2),
        "median_hu": round(float(np.median(hu_values)), 2),
    }


def measure_l3_slice(ct_slice, label_map, pixdim, axial_axis):
    """Cross-sectional area (cm²) + HU stats for each soft-tissue label on the L3 slice."""
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
            for k, v in _hu_stats(hu_vals).items():
                stats[f"l3_{name}_{k}"] = v
        else:
            stats[f"l3_{name}_area_cm2"]  = np.nan
            stats[f"l3_{name}_mean_hu"]   = np.nan
            stats[f"l3_{name}_min_hu"]    = np.nan
            stats[f"l3_{name}_max_hu"]    = np.nan
            stats[f"l3_{name}_median_hu"] = np.nan
    return stats


def measure_vertebra_l3_slice(ct_slice, total_slice, pixdim, axial_axis):
    """Cross-sectional area (cm²) + HU stats for the L3 vertebral body on the L3 slice."""
    dims = [pixdim[i] for i in range(3) if i != axial_axis]
    pixel_area_cm2 = (dims[0] / 10.0) * (dims[1] / 10.0)

    mask = total_slice == L3_VERTEBRA_TOTAL_LABEL
    n_px = int(mask.sum())

    if n_px > 0:
        hu_vals = ct_slice[mask]
        result  = {"l3_vertebra_area_cm2": round(n_px * pixel_area_cm2, 4)}
        for k, v in _hu_stats(hu_vals).items():
            result[f"l3_vertebra_{k}"] = v
    else:
        result = {
            "l3_vertebra_area_cm2":  np.nan,
            "l3_vertebra_mean_hu":   np.nan,
            "l3_vertebra_min_hu":    np.nan,
            "l3_vertebra_max_hu":    np.nan,
            "l3_vertebra_median_hu": np.nan,
        }
    return result


def measure_volume(ct_volume, label_volume, pixdim):
    """Volume (cm³) + mean HU for each output label across a 3-D label volume."""
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
