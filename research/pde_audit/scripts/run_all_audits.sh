#!/usr/bin/env bash
#
# Run all PDE V2 Python audit scripts and save durable stdout captures.
#
# Usage:
#   bash research/pde_audit/scripts/run_all_audits.sh
#   bash research/pde_audit/scripts/run_all_audits.sh stage_v2_22c

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
AUDIT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
FIXTURE_DIR="$SCRIPT_DIR/fixtures"
OUTPUT_DIR="$SCRIPT_DIR/output"
ARTIFACT_DIR="$OUTPUT_DIR/artifacts"
FILTER="${1:-}"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$ARTIFACT_DIR"

summary_file="$OUTPUT_DIR/_summary.txt"
: > "$summary_file"

pass=0
fail=0
total=0

run_capture() {
  local label="$1"
  local script="$2"
  shift 2

  local basename
  basename="$(basename "$label" .py)"
  local out_file="$OUTPUT_DIR/${basename}.txt"
  local tmp_file="$OUTPUT_DIR/${basename}.tmp"
  local start_time end_time elapsed exit_code

  total=$((total + 1))
  printf "  RUN   %s ... " "$basename"
  start_time="$(date +%s)"
  python3 "$script" "$@" > "$tmp_file" 2>&1
  exit_code=$?
  end_time="$(date +%s)"
  elapsed=$((end_time - start_time))

  {
    echo "# Python Audit Output"
    echo "# Script: $(basename "$script")"
    if [ "$basename" != "$(basename "$script" .py)" ]; then
      echo "# Scenario: $basename"
    fi
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
}

args_for_script() {
  local basename="$1"
  SCRIPT_ARGS=()

  case "$basename" in
    stage_v2_21_branch_extraction_fixture)
      SCRIPT_ARGS=(
        --manifest "$FIXTURE_DIR/stage_v2_21_sample_branch_manifest.json"
        --out-json "$ARTIFACT_DIR/stage_v2_21_branch_extraction_packet.json"
      )
      ;;
    stage_v2_22a_profile_to_coefficient_adapter)
      SCRIPT_ARGS=(
        --profile-manifest "$FIXTURE_DIR/stage_v2_22a_profile_input_manifest.json"
        --out-v21-manifest "$ARTIFACT_DIR/stage_v2_22a_generated_v21_manifest.json"
        --out-json "$ARTIFACT_DIR/stage_v2_22a_observable_packet.json"
        --strict-profiles
      )
      ;;
    stage_v2_22b_solver_handoff_validator)
      SCRIPT_ARGS=(
        --solver-output "$FIXTURE_DIR/stage_v2_22b_sample_solver_output_valid.json"
        --out-profile-manifest "$ARTIFACT_DIR/stage_v2_22b_generated_profile_manifest.json"
        --out-report "$ARTIFACT_DIR/stage_v2_22b_validation_report.json"
      )
      ;;
    stage_v2_22c_end_to_end_smoke_pipeline)
      SCRIPT_ARGS=(
        --solver-output "$FIXTURE_DIR/stage_v2_22c_valid_solver_packet.json"
        --invalid-solver-output "$FIXTURE_DIR/stage_v2_22c_invalid_hardcap_solver_packet.json"
        --out-report "$ARTIFACT_DIR/stage_v2_22c_pipeline_report.json"
        --out-profile-manifest "$ARTIFACT_DIR/stage_v2_22c_generated_profile_manifest.json"
        --out-v21-manifest "$ARTIFACT_DIR/stage_v2_22c_generated_v21_manifest.json"
        --out-observable-packet "$ARTIFACT_DIR/stage_v2_22c_observable_packet.json"
        --out-tolerance-budget "$ARTIFACT_DIR/stage_v2_22c_tolerance_budget.json"
        --write-valid-solver "$ARTIFACT_DIR/stage_v2_22c_valid_solver_packet.json"
        --write-invalid-solver "$ARTIFACT_DIR/stage_v2_22c_invalid_hardcap_solver_packet.json"
      )
      ;;
    stage_v2_23_mesh_convergence_audit)
      SCRIPT_ARGS=(
        --out-json "$ARTIFACT_DIR/stage_v2_23_mesh_convergence_report.json"
      )
      ;;
    stage_v2_24_negative_controls)
      SCRIPT_ARGS=(
        --fixture-dir "$FIXTURE_DIR/negative_controls"
        --out-json "$ARTIFACT_DIR/stage_v2_24_negative_controls_report.json"
      )
      ;;
  esac
}

echo "PDE V2 Python Audit Runner"
echo "=========================="
echo "Output dir: $OUTPUT_DIR"
echo "Fixture dir: $FIXTURE_DIR"
echo ""

for script in "$SCRIPT_DIR"/stage_v2_*.py; do
  [ -f "$script" ] || continue
  basename="$(basename "$script" .py)"
  if [ -n "$FILTER" ] && [[ "$basename" != *"$FILTER"* ]]; then
    continue
  fi
  args_for_script "$basename"
  run_capture "$basename" "$script" "${SCRIPT_ARGS[@]}"
done

if [ -z "$FILTER" ] || [[ "stage_v2_22b_solver_handoff_validator_invalid_hardcap" == *"$FILTER"* ]]; then
  run_capture \
    "stage_v2_22b_solver_handoff_validator_invalid_hardcap" \
    "$SCRIPT_DIR/stage_v2_22b_solver_handoff_validator.py" \
    --solver-output "$FIXTURE_DIR/stage_v2_22b_sample_solver_output_invalid_hardcap.json" \
    --out-report "$ARTIFACT_DIR/stage_v2_22b_invalid_hardcap_validation_report.json"
fi

echo "" >> "$summary_file"
echo "TOTAL: $total  PASS: $pass  FAIL: $fail  SKIPPED: 0" >> "$summary_file"

python3 "$SCRIPT_DIR/write_summary_json.py" \
  --suite python \
  --output-dir "$OUTPUT_DIR" \
  --manifest "$FIXTURE_DIR/MANIFEST.json" \
  --artifact-dir "$ARTIFACT_DIR" \
  --filter "$FILTER"
summary_json_status=$?

echo ""
echo "=========================="
echo "Total: $total  |  Pass: $pass  |  Fail: $fail  |  Skipped: 0"
echo "Summary written to: $summary_file"

if [ "$fail" -ne 0 ] || [ "$summary_json_status" -ne 0 ]; then
  exit 1
fi
