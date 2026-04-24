#!/usr/bin/env bash
#
# Rebuild the Mathematica audit summary from saved output files.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/output"
SUMMARY_FILE="$OUTPUT_DIR/_summary.txt"

mkdir -p "$OUTPUT_DIR"
: > "$SUMMARY_FILE"

total=0
pass=0
fail=0

for out_file in "$OUTPUT_DIR"/stage_v2_*.txt; do
  [ -f "$out_file" ] || continue
  basename="$(basename "$out_file" .txt)"
  total=$((total + 1))
  elapsed="$(sed -n 's/^# Elapsed: //p' "$out_file" | head -n 1)"
  if grep -q '^EXIT_CODE: 0$' "$out_file"; then
    echo "PASS  $basename  (${elapsed:-unknown})" >> "$SUMMARY_FILE"
    pass=$((pass + 1))
  else
    exit_code="$(sed -n 's/^# Exit code: //p' "$out_file" | head -n 1)"
    echo "FAIL  $basename  (exit ${exit_code:-unknown}, ${elapsed:-unknown})" >> "$SUMMARY_FILE"
    fail=$((fail + 1))
  fi
done

echo "" >> "$SUMMARY_FILE"
echo "TOTAL: $total  PASS: $pass  FAIL: $fail  SKIPPED: 0" >> "$SUMMARY_FILE"

cat "$SUMMARY_FILE"
