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

## First-time setup (≈ 5 minutes)

This guide assumes you have **never used Python or a terminal before**. If
you're already set up, skip to [Quick start](#quick-start).

### Step 1 — Open a terminal

You'll be typing commands into a terminal (also called a "command line" or
"shell"). It's a built-in app on every computer.

- **macOS:** press `⌘ + Space`, type `Terminal`, hit Enter.
- **Windows:** press the `Windows` key, type `PowerShell`, hit Enter.
  (Throughout this README, "terminal" means PowerShell on Windows.)

Everything that looks like `this` in a gray box is a command you type into
the terminal and then press Enter.

### Step 2 — Install Python (3.9 or newer)

- **macOS:** download the latest installer from
  <https://www.python.org/downloads/macos/>, open the `.pkg` file, and
  click through. Once it's done, verify by running:
  ```bash
  python3 --version
  ```
  It should print something like `Python 3.12.x`.

- **Windows:** download the latest installer from
  <https://www.python.org/downloads/windows/>, open the `.exe` file, and
  **on the very first screen check the box "Add Python to PATH"** before
  clicking Install. Then verify with:
  ```powershell
  python --version
  ```

> Throughout this README, when you see `python` use whichever name works
> for you (`python` on Windows, `python3` on macOS).

### Step 3 — Download this project

You have two options. **Option A** is easier if you've never used git.

**Option A — download as a ZIP (no git needed):**

1. Open this link in your browser:
   <https://github.com/rakshrma/ASD_Risk_BodyComp/archive/refs/heads/main.zip>
2. Unzip the downloaded file.
3. Move the unzipped folder somewhere easy to find — for example, your
   `Documents` folder. Rename it to `ASD_Risk_BodyComp` if it isn't
   already.

**Option B — clone with git** (only if you have git installed):

```bash
git clone https://github.com/rakshrma/ASD_Risk_BodyComp.git
```

### Step 4 — Navigate into the project folder

In the terminal, type `cd ` (the letters c-d followed by a space), then
drag the project folder from your file explorer onto the terminal window
and press Enter. That fills in the full path for you. It should look like:

```bash
# macOS example
cd /Users/yourname/Documents/ASD_Risk_BodyComp

# Windows example
cd C:\Users\yourname\Documents\ASD_Risk_BodyComp
```

You can confirm you're in the right place by running `ls` (macOS) or
`dir` (Windows) — you should see `README.md`, `run_pipeline.py`, `data/`,
etc.

### Step 5 — Create a Python virtual environment and install dependencies

A *virtual environment* is an isolated copy of Python just for this
project. It keeps these packages from interfering with other Python work
on your computer.

- **macOS:**
  ```bash
  python3 -m venv .venv
  source .venv/bin/activate
  pip install --upgrade pip
  pip install -r requirements.txt
  ```

- **Windows (PowerShell):**
  ```powershell
  python -m venv .venv
  .venv\Scripts\Activate.ps1
  pip install --upgrade pip
  pip install -r requirements.txt
  ```
  If PowerShell blocks `Activate.ps1` with an execution-policy error, run
  this once and try again:
  `Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned`.

After this, your terminal prompt will show `(.venv)` at the start — that
means the virtual environment is active. Every time you open a new
terminal to work on this project, `cd` back into the folder and run the
`activate` line again.

> Installing Python + the dependencies takes about **5 minutes** on a
> typical broadband connection. The first time you run TotalSegmentator
> (step 3 of the pipeline), it will additionally download ~3 GB of model
> weights and PyTorch.

### Step 6 — Activate a free TotalSegmentator license

This pipeline uses TotalSegmentator's `tissue_4_types` model to label
subcutaneous fat, visceral fat, muscle, and IMAT. That model (along with
several other body-composition models) is **not bundled with the
open-source release** — it requires a free non-commercial / academic
license. Without it, step 03 of the pipeline will fail to produce the
`*_ts_tissue.nii.gz` file and everything downstream will be incomplete.

1. Request a license — it's free for academic / non-commercial use:
   <https://backend.totalsegmentator.com/license-academic/>

   Fill in the form with your name, institution, and email. You'll
   receive a license key that looks like `aca_XXXXXXXXXXX` by email
   (usually within a few minutes).

2. With your virtual environment still active, activate the key:
   ```bash
   totalseg_set_license -l aca_XXXXXXXXXXX
   ```
   (Replace `aca_XXXXXXXXXXX` with the key you received.) This writes
   the license to your user profile, so you only need to do it once per
   machine.

3. Confirm it works on a small test file:
   ```bash
   TotalSegmentator -i data/sample/<some_file>.nii.gz -o /tmp/ts_test -ta tissue_4_types
   ```
   If you see a license error, recheck the key you pasted. Any other
   warning about model download is normal — TotalSegmentator pulls the
   weights on first use.

> The `total` task (used for the spine + organ labels in this pipeline)
> is open-source and does **not** need a license. The license is only
> required for `tissue_4_types`.

### Step 7 — Install two external tools

These two command-line tools are not Python packages, so they need
separate installers.

| Tool        | What it does                              | macOS                                          | Windows                                                                                  |
| ----------- | ----------------------------------------- | ---------------------------------------------- | ---------------------------------------------------------------------------------------- |
| `gdcmconv`  | Decompresses compressed DICOMs (step 01)  | `brew install gdcm`                            | Download installer: <https://sourceforge.net/projects/gdcm/files/gdcm%20binaries/>       |
| `dcm2niix` | Converts DICOM → NIfTI (step 02)          | `brew install dcm2niix`                        | Download pre-built `.zip`: <https://github.com/rordenlab/dcm2niix/releases>              |

On macOS, `brew` is [Homebrew](https://brew.sh) — install it once with the
one-line command on that page.

On Windows, after downloading `dcm2niix`, unzip it and either (a) move
`dcm2niix.exe` into a folder already on your PATH, or (b) add the folder
containing it to your PATH via *System Properties → Environment
Variables*. Verify with:

```bash
gdcmconv --version
dcm2niix -h
```

If your DICOMs are already uncompressed, you don't need `gdcmconv` — you
can skip step 01 of the pipeline.

### Step 8 — Install a viewer for the results

The pipeline outputs NIfTI files and DICOMs. To **look at the
segmentations on top of the CT** you'll want one of these free viewers:

- **3D Slicer** — most full-featured, works on Windows and macOS:
  <https://download.slicer.org/>
- **ITK-SNAP** — lighter weight, very good for inspecting / editing
  segmentation masks: <http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.SNAP3>

Either one will load the original CT (`*.nii.gz`) and let you overlay the
label map (`*_l1_l5_label.nii.gz` or the L3 DICOMs) on top.

---

## Quick start

Once setup is done, the whole pipeline runs with one command. Make sure
your terminal is in the project folder and the virtual environment is
active (you should see `(.venv)` in your prompt).

### 1. Drop your data into one of the input folders

```
data/input/dicom/<accession>/<dicom files...>          ← raw DICOMs
            OR
data/input/nifti/<accession>/<file>.nii.gz             ← already converted
```

See [docs/DATA_LAYOUT.md](docs/DATA_LAYOUT.md) for the full folder spec
and [File-naming rules](#file-naming-rules) below for how to name files.

### 2. Run the pipeline

```bash
python run_pipeline.py
```

Outputs land in `data/output/l3_results/`, with the master CSV at
`data/output/l3_results/results_summary.csv`.

The orchestrator auto-detects whether to start at decompression, NIfTI
conversion, segmentation, or just the L3 measurements depending on which
folders you populated.

### 3. Inspect the results

Open `data/output/l3_results/results_summary.csv` in Excel or any
spreadsheet program. To check segmentation quality visually, open 3D
Slicer or ITK-SNAP, load the original NIfTI in `data/input/nifti/` and
overlay the matching `*_l1_l5_label.nii.gz` from the output folder.

---

## File-naming rules

**Every file for a different patient must have a unique name.** Two
patients with the same file name will overwrite each other's outputs and
silently corrupt your results CSV.

Recommended layout — one folder per accession (or per patient), with a
descriptive file name inside:

```
data/input/dicom/
├── ACC0001/
│   └── (DICOM slices for one CT series)
├── ACC0002/
│   └── ...
```

For pre-converted NIfTIs:

```
data/input/nifti/
├── ACC0001/
│   └── patient001_ACC0001.nii.gz
├── ACC0002/
│   └── patient002_ACC0002.nii.gz
```

Guidelines:

- Use the accession number (or another **unique** identifier) in both the
  folder name *and* the file name.
- Avoid spaces and special characters — stick to letters, digits,
  underscores, and dashes.
- If a single patient has multiple CTs, give each one a distinct suffix
  (e.g. `patient001_ACC0001_pre.nii.gz`, `patient001_ACC0002_post.nii.gz`).
- **Never** name two NIfTIs the same thing, even in different
  subfolders — the pipeline derives output names from the input file
  stem, so duplicates will collide.

---

## Supplying a manual spinal-cord segmentation (myelograms etc.)

TotalSegmentator's spinal-cord label (79) can be unreliable on
**myelograms** and other CTs with intrathecal contrast, unusual scanner
settings, or atypical anatomy. In those cases, segment the cord yourself
and tell the pipeline to use your mask instead.

### 1. Create the manual SC segmentation

1. Open the CT (`data/input/nifti/<accession>/<stem>.nii.gz`) in ITK-SNAP
   or 3D Slicer.
2. Segment the **spinal cord only**, anywhere from L1 through L5.
3. Save the mask as a NIfTI in the same geometry as the input CT (same
   shape, spacing, and affine). Label value **1 = spinal cord**, label 0
   = background.
4. Place the manual mask in a clearly named folder — for example:
   ```
   data/input/manual_sc/<accession>/<stem>_sc.nii.gz
   ```
   Keep one mask per accession, using the **same stem** as the source
   NIfTI so it's easy to match later.

### 2. Tell the pipeline where the manual masks live

Create `data/input/alt_sc_paths.csv` with these two columns:

| Column        | Description                                                          |
| ------------- | -------------------------------------------------------------------- |
| `nifti_path`  | Absolute path to the original CT NIfTI (must match the manifest row) |
| `seg_path`    | Absolute path to your manual SC NIfTI                                |

Example:

```csv
nifti_path,seg_path
/Users/you/.../data/input/nifti/ACC0001/patient001_ACC0001.nii.gz,/Users/you/.../data/input/manual_sc/ACC0001/patient001_ACC0001_sc.nii.gz
```

### 3. Run step 06 instead of step 05

First make sure you're back in the project folder with the virtual
environment active:

```bash
cd /path/to/ASD_Risk_BodyComp        # macOS / Linux
# or
cd C:\path\to\ASD_Risk_BodyComp      # Windows
```

Then run:

```bash
python scripts/06_process_l3_with_alt_sc.py --filter-csv data/input/alt_sc_paths.csv
```

Only manifest rows listed in your CSV are processed; for those cases the
manual segmentation is used **exclusively** for SC detection (it is not
combined with the TotalSegmentator label). All other cases should still
be processed with `scripts/05_process_l3.py` or `python run_pipeline.py`.

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
`run_pipeline.py`, first make sure your terminal is in the project
folder:

```bash
cd /path/to/ASD_Risk_BodyComp        # macOS / Linux
# or
cd C:\path\to\ASD_Risk_BodyComp      # Windows
```

Then run each step in order:

```bash
bash    scripts/01_decompress_dicom.sh                # DICOMs → decompressed
python  scripts/02_dicom_to_nifti.py                  # → data/input/nifti/
python  scripts/03_run_totalsegmentator.py            # → data/output/totalsegmentator/
python  scripts/04_build_l3_manifest.py               # → data/output/l3_manifest.csv
python  scripts/05_process_l3.py                      # → data/output/l3_results/
```

Use `scripts/06_process_l3_with_alt_sc.py` instead of step 5 when you're
supplying your own SC mask (see
[Supplying a manual spinal-cord segmentation](#supplying-a-manual-spinal-cord-segmentation-myelograms-etc)).
For more on each step, see [docs/PIPELINE.md](docs/PIPELINE.md).

> On Windows, `bash` is not available by default. To run step 01 you can
> either install [Git for Windows](https://git-scm.com/download/win)
> (which includes a `bash` shell) or use Windows Subsystem for Linux
> (WSL). Alternatively, decompress your DICOMs manually with `gdcmconv`
> and skip step 01.

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
