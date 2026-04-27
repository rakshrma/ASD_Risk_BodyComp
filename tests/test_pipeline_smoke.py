"""End-to-end smoke test: run process_row on synthetic data and check outputs."""

import os
from pathlib import Path

import pandas as pd

from asd_bodycomp.pipeline import process_row, run_manifest


def test_process_row_on_synthetic(synthetic_case):
    row = {
        "l3_slice_index": synthetic_case["l3_slice"],
        "nifti_path":     synthetic_case["nifti_path"],
        "ts_total_path":  synthetic_case["ts_total_path"],
        "ts_tissue_path": synthetic_case["ts_tissue_path"],
        "accession":      synthetic_case["accession"],
    }
    out_dir = synthetic_case["out_dir"]
    stats = process_row(row, bbox_ant_mm=20.0, bbox_lr_mm=20.0,
                        base_output_dir=out_dir)

    expected = Path(out_dir) / synthetic_case["accession"]
    stem = synthetic_case["stem"]
    for suffix in (
        f"{stem}_l3_image.dcm", f"{stem}_l3_label.dcm",
        f"{stem}_l3_image.png", f"{stem}_l3_label.png",
        f"{stem}_l1_l5_image.nii.gz", f"{stem}_l1_l5_label.nii.gz",
    ):
        f = expected / suffix
        assert f.is_file(), f"missing output: {f}"
        assert f.stat().st_size > 0

    # All expected measurement keys are present
    for k in ("l3_subq_fat_area_cm2", "l3_muscle_area_cm2",
              "l3_vertebra_area_cm2",
              "vol_subq_fat_volume_cm3", "vol_muscle_volume_cm3"):
        assert k in stats, f"stats missing key: {k}"

    # Vertebra label is present at L3 → area should be a real number
    assert not (stats["l3_vertebra_area_cm2"] != stats["l3_vertebra_area_cm2"])  # not NaN


def test_run_manifest_writes_summary(synthetic_case, tmp_path):
    df = pd.DataFrame([{
        "l3_slice_index": synthetic_case["l3_slice"],
        "nifti_path":     synthetic_case["nifti_path"],
        "ts_total_path":  synthetic_case["ts_total_path"],
        "ts_tissue_path": synthetic_case["ts_tissue_path"],
        "accession":      synthetic_case["accession"],
    }])
    out = tmp_path / "results"
    result_df, errors = run_manifest(df, 20.0, 20.0, str(out))
    assert errors == []
    assert "l3_muscle_area_cm2" in result_df.columns
    assert len(result_df) == 1
