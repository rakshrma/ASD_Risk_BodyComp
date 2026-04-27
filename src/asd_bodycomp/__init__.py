"""
asd_bodycomp
============
L3 / L1-L5 body composition analysis from CT scans, built on TotalSegmentator.

Public modules:
  - labels             label-ID constants (TotalSegmentator + project-specific)
  - nifti_io           NIfTI loading and slice extraction
  - slice_processing   spinal-cord / psoas detection, bbox, label-map building
  - measurements       per-slice and per-volume statistics
  - dicom_export       write a single 2-D slice as a DICOM
  - png_export         CT and label-overlay PNG renderers
  - manifest           build the L3 manifest CSV consumed by step 5/6
  - pipeline           per-row orchestration shared by scripts 5 and 6
"""

__version__ = "0.1.0"
