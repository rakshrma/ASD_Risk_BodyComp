# Sample data

This folder is intentionally empty. The repository ships **no real CT
or DICOM data**.

For testing, the test suite generates a small synthetic CT volume
on-the-fly (see `tests/conftest.py`). To produce one outside the test
suite, run:

```python
from tests.conftest import _build_synthetic_volumes
import nibabel as nib

ct, total, tissue, affine = _build_synthetic_volumes()
nib.save(nib.Nifti1Image(ct,     affine), "sample_ct.nii.gz")
nib.save(nib.Nifti1Image(total,  affine), "sample_ts_total.nii.gz")
nib.save(nib.Nifti1Image(tissue, affine), "sample_ts_tissue.nii.gz")
```

The synthetic volume is 32 × 32 × 16 voxels (≈4 cm × 4 cm × 4.8 cm). It
contains a simulated abdomen with subcutaneous fat, a vertebral body,
spinal cord, iliopsoas pair, and L1–L5 vertebra labels — just enough
structure for the pipeline's geometry detection to succeed.
