#!/usr/bin/env bash
#
# Rebuild the Mathematica audit summary from saved output files.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/output"
SCRIPTS_DIR="$SCRIPT_DIR"
SUMMARY_FILE="$OUTPUT_DIR/_summary.txt"

: > "$SUMMARY_FILE"

total=0
pass=0
fail=0

for out_file in "$OUTPUT_DIR"/moving_throat_pde_stage*_mathematica_audit.txt; do
  [ -f "$out_file" ] || continue

  total=$((total + 1))
  basename="$(basename "$out_file" .txt)"
  if [ ! -f "$SCRIPTS_DIR/${basename}.wl" ]; then
    total=$((total - 1))
    continue
  fi
  elapsed="$(sed -n 's/^# Elapsed: //p' "$out_file" | head -n 1)"

  if grep -q '^EXIT_CODE: 0$' "$out_file"; then
    if [ -n "$elapsed" ]; then
      echo "PASS  $basename  ($elapsed)" >> "$SUMMARY_FILE"
    else
      echo "PASS  $basename" >> "$SUMMARY_FILE"
    fi
    pass=$((pass + 1))
  else
    exit_code="$(sed -n 's/^# Exit code: //p' "$out_file" | head -n 1)"
    if [ -n "$elapsed" ] && [ -n "$exit_code" ]; then
      echo "FAIL  $basename  (exit $exit_code, $elapsed)" >> "$SUMMARY_FILE"
    elif [ -n "$exit_code" ]; then
      echo "FAIL  $basename  (exit $exit_code)" >> "$SUMMARY_FILE"
    else
      echo "FAIL  $basename" >> "$SUMMARY_FILE"
    fi
    fail=$((fail + 1))
  fi
done

echo "" >> "$SUMMARY_FILE"
echo "TOTAL: $total  PASS: $pass  FAIL: $fail  SKIPPED: 0" >> "$SUMMARY_FILE"

cat "$SUMMARY_FILE"
