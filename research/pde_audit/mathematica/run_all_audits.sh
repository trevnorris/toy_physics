#!/usr/bin/env bash
#
# Run all PDE V2 Mathematica audit mirrors and save outputs.
#
# Usage:
#   bash research/pde_audit/mathematica/run_all_audits.sh
#   bash research/pde_audit/mathematica/run_all_audits.sh stage_v2_22c

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
AUDIT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/output"
FILTER="${1:-}"

mkdir -p "$OUTPUT_DIR"

summary_file="$OUTPUT_DIR/_summary.txt"
: > "$summary_file"

pass=0
fail=0
total=0
environment_status=0

echo "PDE V2 Mathematica Audit Runner"
echo "==============================="
echo "Output dir: $OUTPUT_DIR"
echo ""

env_file="$OUTPUT_DIR/_environment.txt"
env_tmp="$OUTPUT_DIR/_environment.tmp"
printf "  RUN   mathematica_environment_probe ... "
env_start_time="$(date +%s)"
math -script "$SCRIPT_DIR/mathematica_environment_probe.wl" > "$env_tmp" 2>&1
env_exit_code=$?
env_elapsed=$(( $(date +%s) - env_start_time ))
{
  echo "# Mathematica Environment Output"
  echo "# Script: mathematica_environment_probe.wl"
  echo "# Date: $(date -Iseconds)"
  echo "# Elapsed: ${env_elapsed}s"
  echo "# Exit code: $env_exit_code"
  if [ "$env_exit_code" -eq 0 ]; then
    echo "# Status: PASS"
  else
    echo "# Status: FAIL"
  fi
  echo "#"
  echo ""
  cat "$env_tmp"
  echo ""
  echo "EXIT_CODE: $env_exit_code"
} > "$env_file"
rm -f "$env_tmp"
if [ "$env_exit_code" -eq 0 ]; then
  printf "PASS (%ds)\n" "$env_elapsed"
else
  printf "FAIL (exit %d, %ds)\n" "$env_exit_code" "$env_elapsed"
  environment_status=1
fi

for script in "$SCRIPT_DIR"/stage_v2_*_mathematica_audit.wl; do
  [ -f "$script" ] || continue
  basename="$(basename "$script" .wl)"
  if [ -n "$FILTER" ] && [[ "$basename" != *"$FILTER"* ]]; then
    continue
  fi

  total=$((total + 1))
  out_file="$OUTPUT_DIR/${basename}.txt"
  tmp_file="$OUTPUT_DIR/${basename}.tmp"
  printf "  RUN   %s ... " "$basename"
  start_time="$(date +%s)"
  math -script "$script" > "$tmp_file" 2>&1
  exit_code=$?
  elapsed=$(( $(date +%s) - start_time ))

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

echo "" >> "$summary_file"
echo "TOTAL: $total  PASS: $pass  FAIL: $fail  SKIPPED: 0" >> "$summary_file"

python3 "$AUDIT_DIR/scripts/write_summary_json.py" \
  --suite mathematica \
  --output-dir "$OUTPUT_DIR" \
  --manifest "$AUDIT_DIR/scripts/fixtures/MANIFEST.json" \
  --filter "$FILTER"
summary_json_status=$?

echo ""
echo "==============================="
echo "Total: $total  |  Pass: $pass  |  Fail: $fail  |  Skipped: 0"
echo "Summary written to: $summary_file"

if [ "$fail" -ne 0 ] || [ "$summary_json_status" -ne 0 ] || [ "$environment_status" -ne 0 ]; then
  exit 1
fi
