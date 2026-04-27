#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────
# 01_decompress_dicom.sh
#
# Decompress every DICOM file under data/input/dicom/<accession>/ using
# gdcmconv --raw, writing the result to data/input/dicom_decompressed/.
# Files that fail to decompress are copied as-is so the output set is complete.
#
# Usage:
#   scripts/01_decompress_dicom.sh                       # default roots
#   scripts/01_decompress_dicom.sh <input_dicom_root>    # custom root
#
# Requirements:
#   gdcmconv  (brew install gdcm  /  apt-get install libgdcm-tools)
# ──────────────────────────────────────────────────────────────────────────

set -uo pipefail   # NOT -e: continue past per-file failures

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INDIR="${1:-$REPO_ROOT/data/input/dicom}"

if [[ ! -d "$INDIR" ]]; then
  echo "ERROR: Input directory does not exist: $INDIR" >&2
  exit 1
fi

if ! command -v gdcmconv >/dev/null 2>&1; then
  echo "ERROR: gdcmconv not found in PATH." >&2
  echo "       Install: brew install gdcm   (macOS)" >&2
  echo "                apt-get install libgdcm-tools  (Debian/Ubuntu)" >&2
  exit 1
fi

PARENT_DIR="$(dirname "$INDIR")"
OUTROOT="$PARENT_DIR/dicom_decompressed"
mkdir -p "$OUTROOT"

LOG="$OUTROOT/decompression.log"
: > "$LOG"

log() {
  local level="$1"; shift
  local ts; ts="$(date '+%Y-%m-%d %H:%M:%S')"
  echo "[$ts] [$level] $*" | tee -a "$LOG" >/dev/null
}

echo "Input  (dicom root) : $INDIR"
echo "Output root         : $OUTROOT"
echo "Log                 : $LOG"
echo

tmpcounts="$(mktemp)"
echo "0 0 0 0" > "$tmpcounts"   # total ok fail accessions

accession_count=$(find "$INDIR" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')
log "INFO" "Found $accession_count accession directories under $INDIR"

find "$INDIR" -mindepth 1 -maxdepth 1 -type d -print0 | while IFS= read -r -d '' ACCDIR; do
  read -r total ok fail acc < "$tmpcounts"
  acc=$((acc+1))
  echo "$total $ok $fail $acc" > "$tmpcounts"

  accession="$(basename "$ACCDIR")"
  OUTDIR="$OUTROOT/$accession"
  mkdir -p "$OUTDIR"
  log "INFO" "Processing accession=$accession"

  find "$ACCDIR" -type f -print0 | while IFS= read -r -d '' f; do
    read -r total ok fail acc < "$tmpcounts"
    total=$((total+1))

    rel="${f#$ACCDIR/}"
    out="$OUTDIR/$rel"
    mkdir -p "$(dirname "$out")"

    if gdcmconv --raw "$f" "$out" >>"$LOG" 2>&1; then
      ok=$((ok+1))
    else
      fail=$((fail+1))
      log "FAIL" "accession=$accession file=$rel (copying original as fallback)"
      cp -p "$f" "$out" >>"$LOG" 2>&1 || true
    fi

    echo "$total $ok $fail $acc" > "$tmpcounts"
  done
done

read -r total ok fail acc < "$tmpcounts"
rm -f "$tmpcounts"

echo
echo "Decompression complete"
echo "  Accessions processed : $acc / $accession_count"
echo "  Total files          : $total"
echo "  Success              : $ok"
echo "  Failed (copied raw)  : $fail"
echo "  Log                  : $LOG"
[[ "$fail" -gt 0 ]] && echo "WARNING: Some files failed; see log."
exit 0
