# ASD Risk — L3 Body Composition Pipeline

End-to-end pipeline that converts abdominal CT scans into L3-slice and
L1–L5 volumetric body-composition measurements, built on
[TotalSegmentator](https://github.com/wasserth/TotalSegmentator).

![Pipeline overview — original CT, TotalSegmentator output, and the L3 ROI with the final subcutaneous-fat / IMAT / muscle label map](docs/Spine%20CT%20Image%20Analysis%20Overview.png)

*From left to right: the original CT (axial + sagittal), the raw TotalSegmentator
labels, and the bounding-box ROI with the final three-class label map
(subcutaneous fat, IMAT, muscle) used for all area / volume / HU statistics.
Top row: spine CT. Bottom row: full abdominopelvic CT.*

For each CT case the pipeline produces:

- **L3 axial slice** as DICOM (`*_l3_image.dcm`) and PNG (`*_l3_image.png`)
- **L3 label map** as DICOM (`*_l3_label.dcm`) and PNG overlay (`*_l3_label.png`)
- **L1–L5 cropped sub-volumes** as NIfTI (`*_l1_l5_image.nii.gz`,
  `*_l1_l5_label.nii.gz`)
- **Combined CSV** with cross-sectional area (cm²), volume (cm³),
  and HU statistics for subcutaneous fat, IMAT, muscle, and the L3
  vertebral body

---

## Quick start

### 1. Install

```bash
git clone <this-repo-url> ASD_Risk_BodyComp
cd ASD_Risk_BodyComp
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

You also need two external binaries on your `PATH`:

| Binary         | Used by    | Install                                                   |
| -------------- | ---------- | --------------------------------------------------------- |
| `gdcmconv`     | step 01    | `brew install gdcm` / `apt-get install libgdcm-tools`     |
| `dcm2niix`     | step 02    | <https://github.com/rordenlab/dcm2niix/releases>          |

`TotalSegmentator` is installed from `requirements.txt` (it brings the
`TotalSegmentator` CLI and PyTorch).

See [docs/INSTALLATION.md](docs/INSTALLATION.md) for full details.

### 2. Drop your data into one of the input folders

```
data/input/dicom/<accession>/<dicom files...>          ← raw DICOMs
            OR
data/input/nifti/<accession>/<file>.nii.gz             ← already converted
```

See [docs/DATA_LAYOUT.md](docs/DATA_LAYOUT.md) for the full folder spec.

### 3. Run the whole pipeline with one command

```bash
python run_pipeline.py
```

Outputs land in `data/output/l3_results/`, with the master CSV at
`data/output/l3_results/results_summary.csv`.

That's it. The orchestrator auto-detects whether to start at decompression,
NIfTI conversion, segmentation, or just the L3 measurements depending on
which folders you populated.

---

## Project layout

```
ASD_Risk_BodyComp/
├── README.md                 ← this file
├── LICENSE                   ← MIT
├── requirements.txt          ← pip install -r
├── run_pipeline.py           ← one-command end-to-end driver
│
├── data/
│   ├── input/
│   │   ├── dicom/            ← user input: DICOMs go here
│   │   └── nifti/            ← OR pre-converted NIfTIs go here
│   ├── output/               ← generated outputs (gitignored)
│   └── sample/               ← synthetic example data
│
├── labels/                   ← CSVs describing label IDs
│   ├── jabba_label_list_fat.csv
│   ├── jabba_label_list_muscle.csv
│   └── ts_4_tissue_labels.csv
│
├── scripts/                  ← per-step CLI entry points
│   ├── 01_decompress_dicom.sh
│   ├── 02_dicom_to_nifti.py
│   ├── 03_run_totalsegmentator.py
│   ├── 04_build_l3_manifest.py
│   ├── 05_process_l3.py
│   └── 06_process_l3_with_alt_sc.py
│
├── src/asd_bodycomp/         ← shared library (importable as `asd_bodycomp`)
│   ├── labels.py
│   ├── nifti_io.py
│   ├── slice_processing.py
│   ├── measurements.py
│   ├── dicom_export.py
│   ├── png_export.py
│   ├── manifest.py
│   └── pipeline.py
│
├── tests/                    ← pytest suite (synthetic CT fixture)
└── docs/
    ├── PIPELINE.md
    ├── INSTALLATION.md
    ├── DATA_LAYOUT.md
    └── LABEL_REFERENCE.md
```

---

## Running individual steps

If you want to control the pipeline manually instead of using
`run_pipeline.py`:

```bash
bash    scripts/01_decompress_dicom.sh                # DICOMs → decompressed
python  scripts/02_dicom_to_nifti.py                  # → data/input/nifti/
python  scripts/03_run_totalsegmentator.py            # → data/output/totalsegmentator/
python  scripts/04_build_l3_manifest.py               # → data/output/l3_manifest.csv
python  scripts/05_process_l3.py                      # → data/output/l3_results/
```

Use `scripts/06_process_l3_with_alt_sc.py` instead of step 5 when the
TotalSegmentator spinal-cord label is unreliable on some cases and you
want to supply your own SC mask. See [docs/PIPELINE.md](docs/PIPELINE.md).

---

## Tests

```bash
pytest -v
```

The test suite generates a small synthetic CT volume (no patient data
needed) and exercises the full pipeline including DICOM/PNG/NIfTI export
and the manifest builder.

---

## Citation / license

Code is released under the MIT License — see [LICENSE](LICENSE).

If you use this in academic work, please also cite
[TotalSegmentator](https://github.com/wasserth/TotalSegmentator).
