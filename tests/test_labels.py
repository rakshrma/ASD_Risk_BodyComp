"""Sanity checks for label constants."""

from asd_bodycomp import labels as L


def test_output_label_names_consistent():
    assert set(L.OUTPUT_LABELS) == set(L.LABEL_NAMES.keys())
    assert L.LABEL_NAMES[1] == "subq_fat"
    assert L.LABEL_NAMES[3] == "imat"
    assert L.LABEL_NAMES[4] == "muscle"


def test_l3_label_in_vertebra_dict():
    assert L.VERTEBRA_LABELS["L3"] == L.L3_VERTEBRA_TOTAL_LABEL


def test_hu_thresholds_are_ordered():
    assert L.HU_FAT_THRESHOLD_LOW < L.HU_FAT_THRESHOLD_HIGH
    assert L.BODY_HU_THRESHOLD < L.HU_FAT_THRESHOLD_LOW


def test_vertebra_labels_cover_l1_l5():
    assert set(L.VERTEBRA_LABELS.keys()) == {"L1", "L2", "L3", "L4", "L5"}
