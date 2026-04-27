"""Tests for the L3 manifest builder (script 04)."""

import pandas as pd

from asd_bodycomp.manifest import build_manifest


def test_build_manifest_single_case(synthetic_case, tmp_path):
    out_csv = tmp_path / "manifest.csv"
    df = build_manifest(synthetic_case["ts_root"],
                        synthetic_case["nifti_root"],
                        str(out_csv))
    assert out_csv.is_file()
    assert len(df) == 1
    row = df.iloc[0]

    expected_cols = {"l3_slice_index", "nifti_path", "ts_total_path",
                     "ts_tissue_path", "accession"}
    assert expected_cols.issubset(df.columns)

    assert row["accession"]      == synthetic_case["accession"]
    assert row["nifti_path"]     == synthetic_case["nifti_path"]
    assert row["ts_total_path"]  == synthetic_case["ts_total_path"]
    assert row["ts_tissue_path"] == synthetic_case["ts_tissue_path"]
    # The fixture places L3 at z=8
    assert int(row["l3_slice_index"]) == synthetic_case["l3_slice"]


def test_manifest_is_empty_when_no_cases(tmp_path):
    (tmp_path / "tsout").mkdir()
    (tmp_path / "nifti").mkdir()
    out_csv = tmp_path / "manifest.csv"
    df = build_manifest(str(tmp_path / "tsout"),
                        str(tmp_path / "nifti"),
                        str(out_csv))
    assert out_csv.is_file()
    assert len(df) == 0
    # Even when empty, the CSV exists so downstream scripts get a clean failure
    parsed = pd.read_csv(out_csv) if out_csv.stat().st_size > 0 else pd.DataFrame()
    assert len(parsed) == 0
