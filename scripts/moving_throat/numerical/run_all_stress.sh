#!/usr/bin/env bash
#
# Run or verify the moving-throat Python numerical stress harnesses.
#
# Usage:
#   bash scripts/moving_throat/numerical/run_all_stress.sh
#   bash scripts/moving_throat/numerical/run_all_stress.sh --update
#   bash scripts/moving_throat/numerical/run_all_stress.sh stage137
#
# Default mode is --check: rerun each harness and compare stdout against the
# saved artifact under scripts/moving_throat/numerical/output/.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
SCRIPTS_DIR="$REPO_ROOT/scripts/moving_throat/numerical"
OUTPUT_DIR="$SCRIPTS_DIR/output"
SUMMARY_FILE="$OUTPUT_DIR/_summary.txt"

mkdir -p "$OUTPUT_DIR"

MODE="check"
FILTER=""
for arg in "$@"; do
  case "$arg" in
    --check) MODE="check" ;;
    --update) MODE="update" ;;
    *) FILTER="$arg" ;;
  esac
done

pass=0
fail=0
skip=0
total=0

: > "$SUMMARY_FILE"

echo "Moving-Throat Python Numerical Stress Runner"
echo "==========================================="
echo "Mode: $MODE"
echo "Output dir: $OUTPUT_DIR"
echo ""

for script in "$SCRIPTS_DIR"/stage*_stress.py; do
  [ -f "$script" ] || continue

  basename="$(basename "$script" .py)"
  out_file="$OUTPUT_DIR/${basename}.txt"

  # stage138_139_stress.py is a superseded exploratory harness; keep the
  # maintained fixed-point regression instead.
  if [ "$basename" = "stage138_139_stress" ]; then
    skip=$((skip + 1))
    continue
  fi

  if [ -n "$FILTER" ] && [[ "$basename" != *"$FILTER"* ]]; then
    continue
  fi

  total=$((total + 1))
  tmp_file="$(mktemp "$OUTPUT_DIR/${basename}.XXXXXX.tmp")"

  printf "  RUN   %s ... " "$basename"
  start_time=$(date +%s)
  set +e
  python "$script" > "$tmp_file" 2>&1
  exit_code=$?
  set -e
  elapsed=$(( $(date +%s) - start_time ))

  if [ "$exit_code" -ne 0 ]; then
    printf "FAIL (exit %d, %ds)\n" "$exit_code" "$elapsed"
    echo "FAIL  $basename  (exit $exit_code, ${elapsed}s)" >> "$SUMMARY_FILE"
    fail=$((fail + 1))
    rm -f "$tmp_file"
    continue
  fi

  if [ "$MODE" = "update" ]; then
    mv "$tmp_file" "$out_file"
    printf "UPDATED (%ds)\n" "$elapsed"
    echo "PASS  $basename  (updated, ${elapsed}s)" >> "$SUMMARY_FILE"
    pass=$((pass + 1))
    continue
  fi

  if [ ! -f "$out_file" ]; then
    printf "FAIL (missing saved output, %ds)\n" "$elapsed"
    echo "FAIL  $basename  (missing saved output)" >> "$SUMMARY_FILE"
    fail=$((fail + 1))
    rm -f "$tmp_file"
    continue
  fi

  if cmp -s "$tmp_file" "$out_file"; then
    printf "PASS (%ds)\n" "$elapsed"
    echo "PASS  $basename  (${elapsed}s)" >> "$SUMMARY_FILE"
    pass=$((pass + 1))
  else
    printf "DRIFT (%ds)\n" "$elapsed"
    echo "FAIL  $basename  (output drift, ${elapsed}s)" >> "$SUMMARY_FILE"
    echo "----- DIFF: $basename -----" >> "$SUMMARY_FILE"
    diff -u "$out_file" "$tmp_file" >> "$SUMMARY_FILE" || true
    echo "" >> "$SUMMARY_FILE"
    fail=$((fail + 1))
  fi

  rm -f "$tmp_file"
done

echo ""
echo "==========================================="
echo "Total: $total  |  Pass: $pass  |  Fail: $fail  |  Skipped: $skip"
echo ""
echo "Summary written to: $SUMMARY_FILE"

echo "" >> "$SUMMARY_FILE"
echo "TOTAL: $total  PASS: $pass  FAIL: $fail  SKIPPED: $skip  MODE: $MODE" >> "$SUMMARY_FILE"

if [ "$fail" -ne 0 ]; then
  exit 1
fi
