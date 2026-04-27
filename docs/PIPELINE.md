# Pipeline reference

The pipeline has six steps. They can be run individually (one CLI per
step under `scripts/`) or end-to-end (`run_pipeline.py`).

```
                              ┌──────────────────────────┐
data/input/dicom/  ──── 01 ─►│ data/input/dicom_decomp/ │
                              └────────────┬─────────────┘
                                           │ 02
                                           ▼
                              ┌──────────────────────────┐
data/input/nifti/  ──────────►│       NIfTI volumes      │
                              └────────────┬─────────────┘
                                           │ 03
                                           ▼
                          ┌──────────────────────────────────┐
                          │ data/output/totalsegmentator/    │
                          │   *_ts_total.nii.gz              │
                          │   *_ts_tissue.nii.gz             │
                          └──────────────┬───────────────────┘
                                         │ 04
                                         ▼
                          ┌──────────────────────────────────┐
                          │ data/output/l3_manifest.csv      │
                          └──────────────┬───────────────────┘
                                         │ 05  (or 06 with alt SC)
                                         ▼
                          ┌──────────────────────────────────┐
                          │ data/output/l3_results/          │
                          │   <accession>/<stem>_l3_*.dcm    │
                          │   <accession>/<stem>_l3_*.png    │
                          │   <accession>/<stem>_l1_l5_*.nii │
                          │   results_summary.csv            │
                          └──────────────────────────────────┘
```

---

## 01 — Decompress DICOM

`scripts/01_decompress_dicom.sh`

Some scanners write DICOM files using compressed transfer syntaxes
(JPEG 2000, JPEG-LS) that `dcm2niix` cannot read directly. This step
re-encodes them with `gdcmconv --raw` so step 02 will work. Files that
fail to decompress are copied as-is so the output set stays complete.

| In  | `data/input/dicom/<accession>/...`              |
| --- | ----------------------------------------------- |
| Out | `data/input/dicom_decompressed/<accession>/...` |
| Log | `data/input/dicom_decompressed/decompression.log` |

Skip this step if your DICOMs are already uncompressed.

---

## 02 — DICOM → NIfTI

`scripts/02_dicom_to_nifti.py`

Runs `dcm2niix -i y -z y -a y -f %i_%f_%s` on each accession. The flags
mean: ignore derived/2D images, gz-compress output, treat the directory
as adjacent DICOMs, and name the file `PatientID_FolderName_SeriesNumber`.

| In  | `data/input/dicom_decompressed/` (or `dicom/`)    |
| --- | ------------------------------------------------- |
| Out | `data/input/nifti/<accession>/<stem>.nii.gz`      |

---

## 03 — TotalSegmentator

`scripts/03_run_totalsegmentator.py`

Runs two TotalSegmentator tasks per CT volume:

| Task             | Output suffix          | Labels we need                       |
| ---------------- | ---------------------- | ------------------------------------ |
| `total`          | `*_ts_total.nii.gz`    | spinal cord (79), psoas L/R (88/89), vertebrae L1–L5 (27–31) |
| `tissue_4_types` | `*_ts_tissue.nii.gz`   | subq fat (1), visceral fat (2), muscle (3), IMAT (4) |

Files that already exist are skipped, so this step is safely re-runnable.

| In  | `data/input/nifti/`                                                         |
| --- | --------------------------------------------------------------------------- |
| Out | `data/output/totalsegmentator/<accession>/<stem>_ts_{total,tissue}.nii.gz`  |

---

## 04 — Build L3 manifest

`scripts/04_build_l3_manifest.py`

Walks the segmentation output, locates the **L3 vertebra centroid** (label
29 in the total volume), and writes one row per case to a manifest CSV.
This step replaces the old `generate_l3_results_ts_models.py`.

Output CSV (`data/output/l3_manifest.csv`):

| Column           | Description                                       |
| ---------------- | ------------------------------------------------- |
| `l3_slice_index` | axial slice index of L3 centroid                  |
| `nifti_path`     | absolute path to original CT NIfTI                |
| `ts_total_path`  | absolute path to `*_ts_total.nii.gz`              |
| `ts_tissue_path` | absolute path to `*_ts_tissue.nii.gz`             |
| `accession`      | accession ID (subdirectory name)                  |

Cases with no L3 label in the total volume are skipped with a console
warning.

---

## 05 — L3 + L1–L5 measurements (TS-only spinal cord)

`scripts/05_process_l3.py`

For each manifest row:

1. Load CT and the two segmentation volumes.
2. Find the spinal-cord centroid (label 79) at the L3 slice, with a
   `±5` slice fallback if missing.
3. Determine the anterior direction from the iliopsoas centroids
   (labels 88, 89).
4. Build a bounding box centered on the SC: 0..`bbox_ant_mm` mm
   anteriorly + `±bbox_lr_mm` mm laterally + every body-HU pixel
   posteriorly.
5. Build a label map within the bbox by combining `_ts_tissue.nii.gz`
   (label 1 → subq fat) with `_ts_total.nii.gz` (psoas labels merged
   into a transient label 2). Label 2 is then split by HU into IMAT
   (label 3, HU < −30 and > −190) and muscle (label 4, HU ≥ −30).
6. Compute area / mean / median / min / max HU per output label, plus
   the L3 vertebral body's area + HU stats.
7. Repeat the bbox/label-map logic for every slice in the L1–L5 span,
   write the cropped sub-volume as NIfTI, and compute volume + mean HU.
8. Write DICOMs (`*_l3_image.dcm`, `*_l3_label.dcm`) and PNGs.

CLI options:

| Flag             | Default | Meaning                            |
| ---------------- | ------- | ---------------------------------- |
| `--bbox-ant-mm`  | 55      | anterior bbox extent in mm         |
| `--bbox-lr-mm`   | 85      | left/right bbox extent in mm       |
| `--limit`        | -1      | max rows to process (-1 = all)     |

Combined results are written to
`data/output/l3_results/results_summary.csv`.

### Measurement columns

| Column                       | Unit | Source                                  |
| ---------------------------- | ---- | --------------------------------------- |
| `l3_<name>_area_cm2`         | cm²  | L3 slice, `name ∈ {subq_fat,imat,muscle}` |
| `l3_<name>_{mean,min,max,median}_hu` | HU | L3 slice                          |
| `l3_vertebra_area_cm2`       | cm²  | L3 slice (label 29)                     |
| `l3_vertebra_{mean,…}_hu`    | HU   | L3 slice (label 29)                     |
| `vol_<name>_volume_cm3`      | cm³  | L1–L5 sub-volume                        |
| `vol_<name>_mean_hu`         | HU   | L1–L5 sub-volume                        |

---

## 06 — L3 measurements with alternative spinal-cord seg

`scripts/06_process_l3_with_alt_sc.py`

Same as step 05, but accepts a second CSV
(`data/input/alt_sc_paths.csv`) listing per-case paths to a
user-supplied SC-only segmentation (label 1 = spinal cord). Only
manifest rows that appear in this CSV are processed; for those cases
the alternative seg is used **exclusively** for SC detection (it is not
combined with the TotalSegmentator label).

This is useful when the `total` task fails to label the SC reliably for
specific cases (rare, but happens with unusual scanner settings or
abnormal anatomy).

---

## End-to-end driver

`run_pipeline.py` chains 01 → 02 → 03 → 04 → 05 and chooses the right
entry point based on which input folder is populated. See
[DATA_LAYOUT.md](DATA_LAYOUT.md) for the four supported entry points.
