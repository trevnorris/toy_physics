#!/usr/bin/env bash
#
# Referee-facing reproducibility harness for the PDE V2 audit stack.
#
# Usage:
#   bash research/pde_audit/run_all.sh

set -u

AUDIT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_DIR="$AUDIT_DIR/output"
mkdir -p "$OUTPUT_DIR"

run_top_check() {
  local name="$1"
  shift
  local out_file="$OUTPUT_DIR/${name}.txt"
  local tmp_file="$OUTPUT_DIR/${name}.tmp"
  local start_time elapsed exit_code

  printf "  RUN   %s ... " "$name"
  start_time="$(date +%s)"
  "$@" > "$tmp_file" 2>&1
  exit_code=$?
  elapsed=$(( $(date +%s) - start_time ))

  {
    echo "# PDE Audit Top-Level Output"
    echo "# Script: $name"
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
  else
    printf "FAIL (exit %d, %ds)\n" "$exit_code" "$elapsed"
  fi
  return "$exit_code"
}

check_root_json_clean() {
  local found
  found="$(find "$AUDIT_DIR" -maxdepth 1 -type f -name '*.json' -print | sort)"
  if [ -n "$found" ]; then
    echo "Root-level generated JSON files are not allowed:"
    echo "$found"
    return 1
  fi
  echo "PASS: no root-level JSON files in $AUDIT_DIR"
}

echo "PDE V2 Referee Reproducibility Harness"
echo "======================================"
echo "Audit dir: $AUDIT_DIR"
echo "Output dir: $OUTPUT_DIR"
echo ""

overall_status=0

run_top_check fixture_manifest python3 "$AUDIT_DIR/scripts/verify_fixtures.py" \
  --fixture-dir "$AUDIT_DIR/scripts/fixtures" || overall_status=1

run_top_check python_audits bash "$AUDIT_DIR/scripts/run_all_audits.sh" || overall_status=1

run_top_check root_json_clean check_root_json_clean || overall_status=1

run_top_check mathematica_audits bash "$AUDIT_DIR/mathematica/run_all_audits.sh" || overall_status=1

python3 "$AUDIT_DIR/scripts/write_referee_summary.py" \
  --audit-dir "$AUDIT_DIR" \
  --output-dir "$OUTPUT_DIR"
summary_status=$?
if [ "$summary_status" -ne 0 ]; then
  overall_status=1
fi

echo ""
echo "======================================"
cat "$OUTPUT_DIR/_summary.txt"

exit "$overall_status"
