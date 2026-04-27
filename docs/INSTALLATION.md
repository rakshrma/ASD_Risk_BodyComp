# Installation

## 1. Python environment

Python 3.9 or newer is required (TotalSegmentator's minimum).

```bash
git clone <this-repo-url> ASD_Risk_BodyComp
cd ASD_Risk_BodyComp
python -m venv .venv
source .venv/bin/activate           # on Windows: .venv\Scripts\activate
pip install --upgrade pip
pip install -r requirements.txt
```

`pip install -r requirements.txt` brings in:

| Package          | Used by                                               |
| ---------------- | ----------------------------------------------------- |
| `numpy`          | numerical work everywhere                             |
| `pandas`         | manifest CSV / results CSV                            |
| `nibabel`        | reading and writing NIfTI                             |
| `pydicom`        | DICOM export (steps 5/6)                              |
| `matplotlib`     | PNG renderers (steps 5/6)                             |
| `TotalSegmentator` | the segmentation CLI **and** PyTorch (transitive) |
| `pytest`         | tests                                                 |

`TotalSegmentator` will install ~3 GB of PyTorch + segmentation weights on
first run. A GPU is **strongly** recommended for step 3.

## 2. External binaries

Two CLI tools must be on your `PATH`. They are not pip-installable.

### `gdcmconv` (used by step 01)

Decompresses DICOMs whose pixel data uses transfer syntaxes that `dcm2niix`
cannot read (notably JPEG 2000 and JPEG-LS).

```bash
# macOS
brew install gdcm

# Debian / Ubuntu
sudo apt-get install libgdcm-tools

# Verify
gdcmconv --version
```

If you already have decompressed DICOMs, you can skip step 01 entirely
and place them in `data/input/dicom_decompressed/` (or just let
`run_pipeline.py --skip-decompress` do it).

### `dcm2niix` (used by step 02)

Converts DICOM series to NIfTI.

```bash
# macOS
brew install dcm2niix

# Debian / Ubuntu — package may be outdated; prefer the GitHub release
sudo apt-get install dcm2niix
# or download a pre-built binary:
#   https://github.com/rordenlab/dcm2niix/releases

# Verify
dcm2niix -h | head -1
```

## 3. Quick verification

Run the test suite — it builds a synthetic CT and exercises the whole
pipeline (no external binaries needed for this, only the Python deps):

```bash
pytest -v
```

All tests should pass. If you see import errors, double-check that you
activated the virtualenv before running `pytest`.

## 4. Optional: GPU setup for TotalSegmentator

CPU-only inference works but is slow (~30 min/case). To use a GPU:

1. Install a CUDA-enabled PyTorch wheel matching your CUDA version
   *before* running `pip install -r requirements.txt`. See
   <https://pytorch.org/get-started/locally/>.
2. Verify with:
   ```bash
   python -c "import torch; print(torch.cuda.is_available())"
   ```
   If this prints `True`, TotalSegmentator will use the GPU automatically.
