"""Label-ID constants used across the pipeline.

All numeric IDs come from the TotalSegmentator output volumes:
  *_ts_total.nii.gz   ← TotalSegmentator -ta total
  *_ts_tissue.nii.gz  ← TotalSegmentator -ta tissue_4_types

The single derived label map (written to *_l3_label.dcm and
*_l1_l5_label.nii.gz) uses three IDs only:
  1 = subcutaneous fat   3 = intramuscular adipose tissue (IMAT)   4 = muscle
"""

# ── total task labels ───────────────────────────────────────────────────────
SPINAL_CORD_LABEL     = 79
ILIOPSOAS_LEFT_LABEL  = 88
ILIOPSOAS_RIGHT_LABEL = 89

# Vertebrae in the total task
VERTEBRA_LABELS = {"L1": 31, "L2": 30, "L3": 29, "L4": 28, "L5": 27}
L3_VERTEBRA_TOTAL_LABEL = VERTEBRA_LABELS["L3"]

# ── tissue_4_types task labels ──────────────────────────────────────────────
# 1=subcutaneous fat  2=visceral fat  3=skeletal muscle  4=imat
TISSUE_RETAIN_LABEL = 1            # kept as label 1 (subcutaneous fat)
TISSUE_MERGE_LABELS = [3, 4]       # muscle/imat merged, then split by HU

# ── alternative spinal-cord seg (script 06) ─────────────────────────────────
SPINAL_CORD_ALT_LABEL = 1          # label 1 in user-supplied SC-only segmentation

# ── derived output labels ───────────────────────────────────────────────────
TOTAL_MERGE_LABELS = [ILIOPSOAS_LEFT_LABEL, ILIOPSOAS_RIGHT_LABEL]
NEW_MERGED_LABEL   = 2             # transient ID before HU split
FAT_SUB_LABEL      = 3             # final IMAT
MUSCLE_SUB_LABEL   = 4             # final muscle
OUTPUT_LABELS      = [1, 3, 4]
LABEL_NAMES        = {1: "subq_fat", 3: "imat", 4: "muscle"}

# ── HU thresholds ───────────────────────────────────────────────────────────
HU_FAT_THRESHOLD_HIGH = -30        # IMAT must be < this
HU_FAT_THRESHOLD_LOW  = -190       # and > this
BODY_HU_THRESHOLD     = -720       # everything >= this counts as body tissue

# ── search-range fallback ───────────────────────────────────────────────────
SC_SEARCH_RADIUS = 5               # ±slices to scan when SC/psoas missing
