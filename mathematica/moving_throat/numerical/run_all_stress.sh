#!/usr/bin/env bash
#
# Run or verify the moving-throat Mathematica numerical stress harnesses.
#
# Usage:
#   bash mathematica/moving_throat/numerical/run_all_stress.sh
#   bash mathematica/moving_throat/numerical/run_all_stress.sh --update
#   bash mathematica/moving_throat/numerical/run_all_stress.sh stage137
#
# Default mode is --check: rerun each harness and compare stdout against the
# saved artifact under mathematica/moving_throat/numerical/output/.
#
# Licensing note:
#   Mathematica harnesses must run one at a time. This runner enforces a cooldown
#   between launches and never parallelizes kernels.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
SCRIPTS_DIR="$REPO_ROOT/mathematica/moving_throat/numerical"
OUTPUT_DIR="$SCRIPTS_DIR/output"
SUMMARY_FILE="$OUTPUT_DIR/_summary.txt"
COOLDOWN_SEC=10

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

echo "Moving-Throat Mathematica Numerical Stress Runner"
echo "================================================"
echo "Mode: $MODE"
echo "Output dir: $OUTPUT_DIR"
echo "Cooldown: ${COOLDOWN_SEC}s"
echo ""

for script in "$SCRIPTS_DIR"/stage*_stress.wl; do
  [ -f "$script" ] || continue

  basename="$(basename "$script" .wl)"
  out_file="$OUTPUT_DIR/${basename}.txt"

  # stage138_139_stress.wl is a superseded exploratory harness; keep the
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
  sleep "$COOLDOWN_SEC"
  start_time=$(date +%s)
  set +e
  math -script "$script" > "$tmp_file" 2>&1
  exit_code=$?
  set -e
  elapsed=$(( $(date +%s) - start_time ))

  # The sandboxed Mathematica kernel occasionally emits a non-mathematical
  # OpenMP shared-memory warning; strip it so drift checks track proof output.
  sed -i '/^OMP: Warning #179:/d' "$tmp_file"

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
echo "================================================"
echo "Total: $total  |  Pass: $pass  |  Fail: $fail  |  Skipped: $skip"
echo ""
echo "Summary written to: $SUMMARY_FILE"

echo "" >> "$SUMMARY_FILE"
echo "TOTAL: $total  PASS: $pass  FAIL: $fail  SKIPPED: $skip  MODE: $MODE" >> "$SUMMARY_FILE"

if [ "$fail" -ne 0 ]; then
  exit 1
fi
