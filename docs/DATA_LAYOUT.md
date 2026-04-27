# Data layout

The pipeline uses a fixed folder layout under `data/`. You only need to
populate one of the input subfolders — `run_pipeline.py` figures out
which steps to run based on what's present.

## Input folders (you populate these)

### Option A: starting from raw DICOMs

```
data/input/dicom/
├── ACC0001/                           ← one folder per accession
│   ├── 0001.dcm
│   ├── 0002.dcm
│   └── ...
├── ACC0002/
│   └── ...
```

`run_pipeline.py` will then run **01 → 02 → 03 → 04 → 05** end-to-end.

### Option B: starting from already-decompressed DICOMs

```
data/input/dicom_decompressed/
├── ACC0001/
│   └── ...
```

Use `run_pipeline.py --skip-decompress` (or just don't put anything in
`data/input/dicom/`). The pipeline runs **02 → 03 → 04 → 05**.

### Option C: starting from NIfTIs

```
data/input/nifti/
├── ACC0001/
│   └── patient_001_ACC0001_001.nii.gz
├── ACC0002/
│   └── ...
```

Pipeline runs **03 → 04 → 05**.

### Option D: starting from existing TotalSegmentator outputs

```
data/output/totalsegmentator/
├── ACC0001/
│   ├── patient_001_ACC0001_001_ts_total.nii.gz
│   └── patient_001_ACC0001_001_ts_tissue.nii.gz
data/input/nifti/
├── ACC0001/
│   └── patient_001_ACC0001_001.nii.gz
```

Pipeline runs **04 → 05** only.

## Output folders (the pipeline writes these)

```
data/input/dicom_decompressed/         ← step 01 (auto-created)
data/input/nifti/                      ← step 02
data/output/totalsegmentator/          ← step 03
data/output/l3_manifest.csv            ← step 04
data/output/l3_results/
├── ACC0001/
│   ├── <stem>_l3_image.dcm
│   ├── <stem>_l3_label.dcm
│   ├── <stem>_l3_image.png
│   ├── <stem>_l3_label.png
│   ├── <stem>_l1_l5_image.nii.gz
│   └── <stem>_l1_l5_label.nii.gz
└── results_summary.csv                ← combined per-row CSV
```

## File-name convention

`<stem>` is the original NIfTI base name (without `.nii.gz`). The
default `dcm2niix` format string used in step 02 is `%i_%f_%s` →
`PatientID_FolderName_SeriesNumber`. So a typical stem looks like:

```
ANON12345_ACC0001_3
```

The pipeline does not depend on this format; it simply re-uses
whatever stem comes out of step 02 to name all downstream outputs.

## Optional: alternative spinal-cord segmentations (script 06)

When the TotalSegmentator spinal-cord label (79) is missing or
unreliable on certain cases, you can supply your own SC segmentation
(label 1) for those cases via a CSV:

```
data/input/alt_sc_paths.csv
```

Required columns:

| Column      | Description                                     |
| ----------- | ----------------------------------------------- |
| `nifti_path`| Path to original CT NIfTI (must match manifest) |
| `seg_path`  | Path to alternative SC NIfTI (label 1 = SC)     |

Then run:

```bash
python scripts/06_process_l3_with_alt_sc.py \
    --filter-csv data/input/alt_sc_paths.csv
```
