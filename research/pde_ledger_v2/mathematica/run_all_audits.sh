#!/usr/bin/env bash
#
# Run ledger Mathematica stage audit scripts and save outputs.
#
# Usage:
#   bash research/pde_ledger/mathematica/run_all_audits.sh
#   bash research/pde_ledger/mathematica/run_all_audits.sh --force
#   bash research/pde_ledger/mathematica/run_all_audits.sh stage106

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS_DIR="$SCRIPT_DIR"
OUTPUT_DIR="$SCRIPTS_DIR/output"

mkdir -p "$OUTPUT_DIR"

FORCE=0
FILTER=""
COOLDOWN_SEC=10
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

echo "Ledger Mathematica Audit Runner"
echo "==============================="
echo "Output dir: $OUTPUT_DIR"
echo ""

for script in "$SCRIPTS_DIR"/ledger_stage*_mathematica_audit.wl; do
  [ -f "$script" ] || continue

  basename=$(basename "$script" .wl)
  out_file="$OUTPUT_DIR/${basename}.txt"

  if [ -n "$FILTER" ] && [[ "$basename" != *"$FILTER"* ]]; then
    continue
  fi

  total=$((total + 1))

  if [ "$FORCE" -eq 0 ] && [ -f "$out_file" ] && [ "$out_file" -nt "$script" ]; then
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
  sleep "$COOLDOWN_SEC"
  start_time=$(date +%s)
  tmp_file=$(mktemp "$OUTPUT_DIR/${basename}.XXXXXX.tmp")

  set +e
  math -script "$script" > "$tmp_file" 2>&1
  exit_code=$?
  set -e

  end_time=$(date +%s)
  elapsed=$((end_time - start_time))

  {
    echo "# Mathematica Audit Output"
    echo "# Script: $basename.wl"
    echo "# Date: $(date -Iseconds)"
    echo "# Elapsed: ${elapsed}s"
    echo "# Exit code: $exit_code"
    if [ "$exit_code" -eq 0 ]; then
      echo "# Status: PASS"
    else
      echo "# Status: FAIL"
    fi
    echo "#"
    echo ""
    cat "$tmp_file"
    echo ""
    echo "EXIT_CODE: $exit_code"
  } > "$out_file"

  rm -f "$tmp_file"

  if [ "$exit_code" -eq 0 ]; then
    printf "PASS (%ds)\n" "$elapsed"
    echo "PASS  $basename  (${elapsed}s)" >> "$summary_file"
    pass=$((pass + 1))
  else
    printf "FAIL (exit %d, %ds)\n" "$exit_code" "$elapsed"
    echo "FAIL  $basename  (exit $exit_code, ${elapsed}s)" >> "$summary_file"
    fail=$((fail + 1))
  fi
done

if [ "$total" -eq 0 ]; then
  echo "No ledger Mathematica stage audits found; nothing to do."
  echo "NO_STAGES: nothing to do" >> "$summary_file"
fi

echo ""
echo "==============================="
echo "Total: $total  |  Pass: $pass  |  Fail: $fail  |  Skipped: $skip"
echo ""
echo "Summary written to: $summary_file"
echo "Individual outputs in: $OUTPUT_DIR/"

echo "" >> "$summary_file"
echo "TOTAL: $total  PASS: $pass  FAIL: $fail  SKIPPED: $skip" >> "$summary_file"
