#!/usr/bin/env bash
#
# Run all moving-throat SymPy audit scripts and save outputs.
# Safe to re-run — only runs scripts whose output is missing or older
# than the source script (use --force to re-run everything).
#
# Usage:
#   bash research/pde_ledger/scripts/run_all_audits.sh           # incremental
#   bash research/pde_ledger/scripts/run_all_audits.sh --force    # re-run all
#   bash research/pde_ledger/scripts/run_all_audits.sh stage059   # run one stage

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS_DIR="$SCRIPT_DIR"
OUTPUT_DIR="$SCRIPTS_DIR/output"
TIMEOUT_SEC=0  # 0 = no timeout

mkdir -p "$OUTPUT_DIR"

FORCE=0
FILTER=""
for arg in "$@"; do
  case "$arg" in
    --force) FORCE=1 ;;
    *) FILTER="$arg" ;;
  esac
done

pass=0
fail=0
skip=0
total=0

summary_file="$OUTPUT_DIR/_summary.txt"
: > "$summary_file"

echo "Moving-Throat SymPy Audit Runner"
echo "================================"
echo "Timeout per script: ${TIMEOUT_SEC}s"
echo "Output dir: $OUTPUT_DIR"
echo ""

for script in "$SCRIPTS_DIR"/moving_throat_pde_*_sympy_audit*.py; do
  [ -f "$script" ] || continue

  basename=$(basename "$script" .py)
  out_file="$OUTPUT_DIR/${basename}.txt"

  # Filter if a specific stage was requested
  if [ -n "$FILTER" ] && [[ "$basename" != *"$FILTER"* ]]; then
    continue
  fi

  total=$((total + 1))

  # Skip if output exists and is newer than source (unless --force)
  if [ "$FORCE" -eq 0 ] && [ -f "$out_file" ] && [ "$out_file" -nt "$script" ]; then
    # Check if previous run passed or failed
    if grep -q "^EXIT_CODE: 0$" "$out_file" 2>/dev/null; then
      printf "  SKIP (cached PASS) %s\n" "$basename"
      echo "PASS  $basename" >> "$summary_file"
      pass=$((pass + 1))
    else
      printf "  SKIP (cached FAIL) %s\n" "$basename"
      echo "FAIL  $basename  (cached)" >> "$summary_file"
      fail=$((fail + 1))
    fi
    skip=$((skip + 1))
    continue
  fi

  printf "  RUN   %s ... " "$basename"

  start_time=$(date +%s)

  # Run script, capture stdout+stderr
  set +e
  if [ "$TIMEOUT_SEC" -gt 0 ]; then
    timeout "$TIMEOUT_SEC" python3 "$script" > "$out_file.tmp" 2>&1
  else
    python3 "$script" > "$out_file.tmp" 2>&1
  fi
  exit_code=$?
  set -e

  end_time=$(date +%s)
  elapsed=$((end_time - start_time))

  # Build output file with metadata header
  {
    echo "# SymPy Audit Output"
    echo "# Script: $basename.py"
    echo "# Date: $(date -Iseconds)"
    echo "# Elapsed: ${elapsed}s"
    echo "# Exit code: $exit_code"
    if [ "$exit_code" -eq 124 ]; then
      echo "# Status: TIMEOUT (${TIMEOUT_SEC}s limit)"
    elif [ "$exit_code" -eq 0 ]; then
      echo "# Status: PASS"
    else
      echo "# Status: FAIL"
    fi
    echo "#"
    echo ""
    cat "$out_file.tmp"
    echo ""
    echo "EXIT_CODE: $exit_code"
  } > "$out_file"

  rm -f "$out_file.tmp"

  if [ "$exit_code" -eq 124 ]; then
    printf "TIMEOUT (%ds)\n" "$elapsed"
    echo "TIMEOUT  $basename  (${elapsed}s)" >> "$summary_file"
    fail=$((fail + 1))
  elif [ "$exit_code" -eq 0 ]; then
    printf "PASS (%ds)\n" "$elapsed"
    echo "PASS  $basename  (${elapsed}s)" >> "$summary_file"
    pass=$((pass + 1))
  else
    printf "FAIL (exit %d, %ds)\n" "$exit_code" "$elapsed"
    echo "FAIL  $basename  (exit $exit_code, ${elapsed}s)" >> "$summary_file"
    fail=$((fail + 1))
  fi

done

echo ""
echo "================================"
echo "Total: $total  |  Pass: $pass  |  Fail: $fail  |  Skipped: $skip"
echo ""
echo "Summary written to: $summary_file"
echo "Individual outputs in: $OUTPUT_DIR/"

# Append totals to summary
echo "" >> "$summary_file"
echo "TOTAL: $total  PASS: $pass  FAIL: $fail  SKIPPED: $skip" >> "$summary_file"
