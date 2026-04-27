# Label reference

All numeric IDs are TotalSegmentator's own label IDs. The pipeline uses
two output volumes per CT case.

## `*_ts_total.nii.gz`  (task `total`)

The full TotalSegmentator label set (≈110 organs). The pipeline uses
only these labels:

| ID  | Structure                    | Used for                                |
| --- | ---------------------------- | --------------------------------------- |
| 27  | vertebra L5                  | L1–L5 slice-range detection (step 5/6)  |
| 28  | vertebra L4                  | "                                       |
| 29  | vertebra L3                  | L3 slice index (step 4) + vertebra area |
| 30  | vertebra L2                  | L1–L5 slice-range detection             |
| 31  | vertebra L1                  | "                                       |
| 79  | spinal cord                  | bounding-box origin (step 5)            |
| 88  | iliopsoas — right            | anterior-direction detection            |
| 89  | iliopsoas — left             | anterior-direction detection            |

Iliopsoas labels 88/89 are also folded into the merged muscle region
(see below) before HU-splitting.

## `*_ts_tissue.nii.gz`  (task `tissue_4_types`)

| ID  | Structure                    | Treatment                                                     |
| --- | ---------------------------- | ------------------------------------------------------------- |
| 1   | subcutaneous fat             | retained as output label **1**                                |
| 2   | visceral fat                 | not used (out of scope for L3 abdominal-wall analysis)        |
| 3   | skeletal muscle              | merged with label 4 and the psoas, then split by HU           |
| 4   | IMAT                         | "                                                             |

## Output label map

Final labels written into `*_l3_label.dcm` and `*_l1_l5_label.nii.gz`:

| ID  | Name in CSV / PNG legend  | Meaning                                          |
| --- | ------------------------- | ------------------------------------------------ |
| 1   | `subq_fat`                | subcutaneous fat (= tissue label 1)              |
| 3   | `imat`                    | intramuscular adipose tissue (HU < −30 & > −190) |
| 4   | `muscle`                  | skeletal muscle (HU ≥ −30)                       |

Label `2` exists transiently during processing (merged muscle/IMAT
region) and is split into 3 and 4 by HU thresholding before any output
is written.

## HU thresholds

| Constant                | Value | Used for                                      |
| ----------------------- | ----- | --------------------------------------------- |
| `HU_FAT_THRESHOLD_HIGH` | −30   | upper bound of IMAT HU window                 |
| `HU_FAT_THRESHOLD_LOW`  | −190  | lower bound of IMAT HU window (excludes air)  |
| `BODY_HU_THRESHOLD`     | −720  | "this pixel is body tissue" (posterior bbox)  |

All defined in [`src/asd_bodycomp/labels.py`](../src/asd_bodycomp/labels.py).

## CSV files in `labels/`

These three legacy CSVs are kept for downstream tools that key off
label IDs by name. They are **not** read by any script in this
repository — they are reference material.

- `ts_4_tissue_labels.csv` — TotalSegmentator `tissue_4_types` IDs.
- `jabba_label_list_fat.csv`, `jabba_label_list_muscle.csv` —
  label IDs from a separate JABBA segmentation system, retained for
  cross-reference.
