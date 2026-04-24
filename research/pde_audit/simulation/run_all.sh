#!/usr/bin/env bash
#
# Run the target-blind reduced simulation bundle and evaluate frozen outputs.

set -u

SIM_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_DIR="$SIM_DIR/output"
mkdir -p "$OUTPUT_DIR"

summary_file="$OUTPUT_DIR/_summary.txt"
: > "$summary_file"

pass=0
fail=0
total=0

run_capture() {
  local name="$1"
  shift
  local out_file="$OUTPUT_DIR/${name}.txt"
  local tmp_file="$OUTPUT_DIR/${name}.tmp"
  local start_time elapsed exit_code

  total=$((total + 1))
  printf "  RUN   %s ... " "$name"
  start_time="$(date +%s)"
  "$@" > "$tmp_file" 2>&1
  exit_code=$?
  elapsed=$(( $(date +%s) - start_time ))

  {
    echo "# PDE Audit Simulation Output"
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
    echo "PASS  $name  (${elapsed}s)" >> "$summary_file"
    pass=$((pass + 1))
  else
    printf "FAIL (exit %d, %ds)\n" "$exit_code" "$elapsed"
    echo "FAIL  $name  (exit $exit_code, ${elapsed}s)" >> "$summary_file"
    fail=$((fail + 1))
  fi
  return "$exit_code"
}

echo "PDE Audit Simulation Runner"
echo "==========================="
echo "Simulation dir: $SIM_DIR"
echo "Output dir: $OUTPUT_DIR"
echo ""

overall_status=0

run_capture verify_reduced_fem python3 "$SIM_DIR/verify_reduced_fem.py" \
  --output-dir "$OUTPUT_DIR" || overall_status=1

run_capture verify_nonlinear_solver python3 "$SIM_DIR/verify_nonlinear_solver.py" \
  --output-dir "$OUTPUT_DIR" || overall_status=1

run_capture verify_physical_model python3 "$SIM_DIR/verify_physical_model.py" \
  --output-dir "$OUTPUT_DIR" || overall_status=1

run_capture generate_nonlinear_packets python3 "$SIM_DIR/generate_nonlinear_packets.py" \
  --output-dir "$OUTPUT_DIR" || overall_status=1

run_capture verify_nonlinear_target_blind python3 "$SIM_DIR/verify_target_blind.py" \
  --output-dir "$OUTPUT_DIR" \
  --mode nonlinear \
  --manifest-name nonlinear_manifest.json \
  --report-prefix nonlinear_target_blind_report || overall_status=1

run_capture evaluate_nonlinear_frozen_sweep python3 "$SIM_DIR/evaluate_frozen_sweep.py" \
  --output-dir "$OUTPUT_DIR" \
  --manifest-name nonlinear_manifest.json \
  --details-dir nonlinear_candidate_reports \
  --summary-prefix nonlinear_evaluation_summary || overall_status=1

run_capture diagnose_nonlinear_obstruction python3 "$SIM_DIR/diagnose_obstruction.py" \
  --output-dir "$OUTPUT_DIR" \
  --candidate-dir nonlinear_candidate_reports \
  --report-prefix nonlinear_obstruction_report || overall_status=1

run_capture diagnose_nonlinear_required_deformation python3 "$SIM_DIR/diagnose_required_deformation.py" \
  --output-dir "$OUTPUT_DIR" \
  --candidate-dir nonlinear_candidate_reports \
  --report-prefix nonlinear_required_deformation_report || overall_status=1

run_capture generate_reduced_sweep python3 "$SIM_DIR/generate_reduced_sweep.py" \
  --output-dir "$OUTPUT_DIR" || overall_status=1

run_capture verify_target_blind python3 "$SIM_DIR/verify_target_blind.py" \
  --output-dir "$OUTPUT_DIR" || overall_status=1

run_capture evaluate_frozen_sweep python3 "$SIM_DIR/evaluate_frozen_sweep.py" \
  --output-dir "$OUTPUT_DIR" || overall_status=1

run_capture diagnose_obstruction python3 "$SIM_DIR/diagnose_obstruction.py" \
  --output-dir "$OUTPUT_DIR" || overall_status=1

run_capture diagnose_required_deformation python3 "$SIM_DIR/diagnose_required_deformation.py" \
  --output-dir "$OUTPUT_DIR" || overall_status=1

run_capture diagnose_mechanism_gap python3 "$SIM_DIR/diagnose_mechanism_gap.py" \
  --output-dir "$OUTPUT_DIR" || overall_status=1

echo "" >> "$summary_file"
echo "TOTAL: $total  PASS: $pass  FAIL: $fail  SKIPPED: 0" >> "$summary_file"

python3 - "$OUTPUT_DIR" "$total" "$pass" "$fail" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

output_dir = Path(sys.argv[1])
total = int(sys.argv[2])
passed = int(sys.argv[3])
failed = int(sys.argv[4])

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

outputs = []
for path in sorted(output_dir.glob("*.txt")):
    if path.name.startswith("_"):
        continue
    outputs.append({"path": str(path), "sha256": sha256_file(path)})

artifacts = []
for path in sorted(output_dir.rglob("*.json")):
    artifacts.append({"path": str(path), "sha256": sha256_file(path)})

evaluation_path = output_dir / "evaluation_summary.json"
evaluation = {}
if evaluation_path.exists():
    evaluation = json.loads(evaluation_path.read_text(encoding="utf-8"))
nonlinear_evaluation_path = output_dir / "nonlinear_evaluation_summary.json"
nonlinear_evaluation = {}
if nonlinear_evaluation_path.exists():
    nonlinear_evaluation = json.loads(nonlinear_evaluation_path.read_text(encoding="utf-8"))
obstruction_path = output_dir / "obstruction_report.json"
obstruction = {}
if obstruction_path.exists():
    obstruction = json.loads(obstruction_path.read_text(encoding="utf-8"))
nonlinear_obstruction_path = output_dir / "nonlinear_obstruction_report.json"
nonlinear_obstruction = {}
if nonlinear_obstruction_path.exists():
    nonlinear_obstruction = json.loads(nonlinear_obstruction_path.read_text(encoding="utf-8"))
deformation_path = output_dir / "required_deformation_report.json"
deformation = {}
if deformation_path.exists():
    deformation = json.loads(deformation_path.read_text(encoding="utf-8"))
nonlinear_deformation_path = output_dir / "nonlinear_required_deformation_report.json"
nonlinear_deformation = {}
if nonlinear_deformation_path.exists():
    nonlinear_deformation = json.loads(nonlinear_deformation_path.read_text(encoding="utf-8"))
mechanism_gap_path = output_dir / "mechanism_gap_report.json"
mechanism_gap = {}
if mechanism_gap_path.exists():
    mechanism_gap = json.loads(mechanism_gap_path.read_text(encoding="utf-8"))
target_blind_path = output_dir / "target_blind_report.json"
target_blind = {}
if target_blind_path.exists():
    target_blind = json.loads(target_blind_path.read_text(encoding="utf-8"))
nonlinear_target_blind_path = output_dir / "nonlinear_target_blind_report.json"
nonlinear_target_blind = {}
if nonlinear_target_blind_path.exists():
    nonlinear_target_blind = json.loads(nonlinear_target_blind_path.read_text(encoding="utf-8"))
reduced_fem_path = output_dir / "reduced_fem_verification.json"
reduced_fem = {}
if reduced_fem_path.exists():
    reduced_fem = json.loads(reduced_fem_path.read_text(encoding="utf-8"))
nonlinear_path = output_dir / "nonlinear_readiness.json"
nonlinear = {}
if nonlinear_path.exists():
    nonlinear = json.loads(nonlinear_path.read_text(encoding="utf-8"))
physical_model_path = output_dir / "physical_model_status.json"
physical_model = {}
if physical_model_path.exists():
    physical_model = json.loads(physical_model_path.read_text(encoding="utf-8"))

summary = {
    "schema": "pde_audit_simulation_runner_summary/v1",
    "output_dir": str(output_dir),
    "total": total,
    "passed": passed,
    "failed": failed,
    "outputs": outputs,
    "artifacts": artifacts,
    "evaluation": {
        "candidate_count": evaluation.get("candidate_count"),
        "validation_pass_count": evaluation.get("validation_pass_count"),
        "open_stable_count": evaluation.get("open_stable_count"),
        "target_pass_count": evaluation.get("target_pass_count"),
        "score_distribution": evaluation.get("score_distribution"),
        "dominant_score_component_counts": evaluation.get("dominant_score_component_counts"),
        "pass_flag_counts": evaluation.get("pass_flag_counts"),
        "summary_hash": evaluation.get("summary_hash"),
        "best_post_hoc_candidate": evaluation.get("best_post_hoc_candidate", {}).get("branch_id") if evaluation.get("best_post_hoc_candidate") else None,
        "best_post_hoc_score": evaluation.get("best_post_hoc_candidate", {}).get("score") if evaluation.get("best_post_hoc_candidate") else None,
    },
    "reduced_fem_verification": {
        "pass": reduced_fem.get("pass"),
        "failed_checks": reduced_fem.get("failed_checks"),
        "report_hash": reduced_fem.get("report_hash"),
        "checks": [
            {"name": item.get("name"), "pass": item.get("pass")}
            for item in reduced_fem.get("checks", [])
        ],
    },
    "nonlinear_readiness": {
        "pass": nonlinear.get("pass"),
        "failed_checks": nonlinear.get("failed_checks"),
        "report_hash": nonlinear.get("report_hash"),
        "protocol_hash": nonlinear.get("protocol_hash"),
        "packets_emitted": nonlinear.get("packets_emitted"),
        "target_residuals_computed": nonlinear.get("target_residuals_computed"),
        "checks": [
            {"name": item.get("name"), "pass": item.get("pass")}
            for item in nonlinear.get("checks", [])
        ],
    },
    "physical_model_status": {
        "pass": physical_model.get("pass"),
        "failed_checks": physical_model.get("failed_checks"),
        "report_hash": physical_model.get("report_hash"),
        "status_hash": (physical_model.get("physical_model_status") or {}).get("status_hash"),
        "strict_parent_dynamic_wall_pass": physical_model.get("strict_parent_dynamic_wall_pass"),
        "effective_wall_closure_available": physical_model.get("effective_wall_closure_available"),
        "physical_export_permitted": physical_model.get("physical_export_permitted"),
        "packets_emitted": physical_model.get("packets_emitted"),
        "target_residuals_computed": physical_model.get("target_residuals_computed"),
        "blocker_count": physical_model.get("blocker_count"),
        "required_before_physical_export_count": physical_model.get("required_before_physical_export_count"),
        "checks": [
            {"name": item.get("name"), "pass": item.get("pass")}
            for item in physical_model.get("checks", [])
        ],
    },
    "nonlinear_export": {
        "candidate_count": nonlinear_evaluation.get("candidate_count"),
        "validation_pass_count": nonlinear_evaluation.get("validation_pass_count"),
        "open_stable_count": nonlinear_evaluation.get("open_stable_count"),
        "target_pass_count": nonlinear_evaluation.get("target_pass_count"),
        "score_distribution": nonlinear_evaluation.get("score_distribution"),
        "dominant_score_component_counts": nonlinear_evaluation.get("dominant_score_component_counts"),
        "pass_flag_counts": nonlinear_evaluation.get("pass_flag_counts"),
        "summary_hash": nonlinear_evaluation.get("summary_hash"),
        "best_post_hoc_candidate": nonlinear_evaluation.get("best_post_hoc_candidate", {}).get("branch_id") if nonlinear_evaluation.get("best_post_hoc_candidate") else None,
        "best_post_hoc_score": nonlinear_evaluation.get("best_post_hoc_candidate", {}).get("score") if nonlinear_evaluation.get("best_post_hoc_candidate") else None,
        "target_blind_guard": {
            "pass": nonlinear_target_blind.get("pass"),
            "issue_count": nonlinear_target_blind.get("issue_count"),
            "report_hash": nonlinear_target_blind.get("report_hash"),
            "source_imports_pass": (nonlinear_target_blind.get("source_imports") or {}).get("pass"),
            "manifest_and_packets_pass": (nonlinear_target_blind.get("manifest_and_packets") or {}).get("pass"),
            "packet_count": ((nonlinear_target_blind.get("manifest_and_packets") or {}).get("details") or {}).get("packet_count"),
        },
        "obstruction": {
            "report_hash": nonlinear_obstruction.get("report_hash"),
            "open_stable_one_pole_ratios_below_one": (nonlinear_obstruction.get("open_stable_candidates") or {}).get("all_one_pole_ratios_below_one"),
            "open_stable_C_shortfalls_positive": (nonlinear_obstruction.get("open_stable_candidates") or {}).get("all_C_shortfalls_positive"),
            "open_stable_near_one_pole_count_abs_ratio_minus_1_lt_0p1": (nonlinear_obstruction.get("open_stable_candidates") or {}).get("near_one_pole_count_abs_ratio_minus_1_lt_0p1"),
            "open_stable_one_pole_ratio_distribution": (nonlinear_obstruction.get("open_stable_candidates") or {}).get("one_pole_ratio_distribution"),
        },
        "required_deformation": {
            "report_hash": nonlinear_deformation.get("report_hash"),
            "candidate_count": nonlinear_deformation.get("candidate_count"),
            "target_pass_count": nonlinear_deformation.get("target_pass_count"),
            "open_stable_count": nonlinear_deformation.get("open_stable_count"),
            "dominant_score_component_counts": nonlinear_deformation.get("dominant_score_component_counts"),
            "local_continuation_looks_promising": (nonlinear_deformation.get("mechanism_assessment") or {}).get("local_continuation_looks_promising"),
            "suggested_missing_mechanism": (nonlinear_deformation.get("mechanism_assessment") or {}).get("suggested_missing_mechanism"),
            "open_stable_required_C_multiplier_distribution": (nonlinear_deformation.get("open_stable_candidates") or {}).get("required_C_multiplier_distribution"),
            "open_stable_required_P0_multiplier_distribution": (nonlinear_deformation.get("open_stable_candidates") or {}).get("required_P0_multiplier_distribution"),
            "open_stable_one_pole_ratio_distribution": (nonlinear_deformation.get("open_stable_candidates") or {}).get("one_pole_ratio_distribution"),
        },
    },
    "target_blind_guard": {
        "pass": target_blind.get("pass"),
        "issue_count": target_blind.get("issue_count"),
        "report_hash": target_blind.get("report_hash"),
        "source_imports_pass": (target_blind.get("source_imports") or {}).get("pass"),
        "manifest_and_packets_pass": (target_blind.get("manifest_and_packets") or {}).get("pass"),
        "packet_count": ((target_blind.get("manifest_and_packets") or {}).get("details") or {}).get("packet_count"),
    },
    "obstruction": {
        "report_hash": obstruction.get("report_hash"),
        "all_candidates": obstruction.get("all_candidates"),
        "open_stable_candidates": obstruction.get("open_stable_candidates"),
        "all_one_pole_ratios_below_one": obstruction.get("all_one_pole_ratios_below_one"),
        "all_C_shortfalls_positive": obstruction.get("all_C_shortfalls_positive"),
        "near_one_pole_count_abs_ratio_minus_1_lt_0p1": obstruction.get("near_one_pole_count_abs_ratio_minus_1_lt_0p1"),
        "one_pole_ratio_distribution": obstruction.get("one_pole_ratio_distribution"),
        "C_shortfall_distribution": obstruction.get("C_shortfall_distribution"),
        "P0_over_target_distribution": obstruction.get("P0_over_target_distribution"),
        "open_stable_one_pole_ratios_below_one": (obstruction.get("open_stable_candidates") or {}).get("all_one_pole_ratios_below_one"),
        "open_stable_C_shortfalls_positive": (obstruction.get("open_stable_candidates") or {}).get("all_C_shortfalls_positive"),
        "open_stable_near_one_pole_count_abs_ratio_minus_1_lt_0p1": (obstruction.get("open_stable_candidates") or {}).get("near_one_pole_count_abs_ratio_minus_1_lt_0p1"),
        "open_stable_one_pole_ratio_distribution": (obstruction.get("open_stable_candidates") or {}).get("one_pole_ratio_distribution"),
        "open_stable_C_shortfall_distribution": (obstruction.get("open_stable_candidates") or {}).get("C_shortfall_distribution"),
        "open_stable_P0_over_target_distribution": (obstruction.get("open_stable_candidates") or {}).get("P0_over_target_distribution"),
    },
    "required_deformation": {
        "report_hash": deformation.get("report_hash"),
        "candidate_count": deformation.get("candidate_count"),
        "target_pass_count": deformation.get("target_pass_count"),
        "open_stable_count": deformation.get("open_stable_count"),
        "dominant_score_component_counts": deformation.get("dominant_score_component_counts"),
        "local_continuation_looks_promising": (deformation.get("mechanism_assessment") or {}).get("local_continuation_looks_promising"),
        "suggested_missing_mechanism": (deformation.get("mechanism_assessment") or {}).get("suggested_missing_mechanism"),
        "open_stable_required_C_multiplier_distribution": (deformation.get("open_stable_candidates") or {}).get("required_C_multiplier_distribution"),
        "open_stable_required_P0_multiplier_distribution": (deformation.get("open_stable_candidates") or {}).get("required_P0_multiplier_distribution"),
        "open_stable_one_pole_ratio_distribution": (deformation.get("open_stable_candidates") or {}).get("one_pole_ratio_distribution"),
    },
    "mechanism_gap": {
        "pass": mechanism_gap.get("pass"),
        "report_hash": mechanism_gap.get("report_hash"),
        "post_hoc_only": mechanism_gap.get("post_hoc_only"),
        "target_residuals_used_to_generate_candidates": mechanism_gap.get("target_residuals_used_to_generate_candidates"),
        "candidate_generation_mutated": mechanism_gap.get("candidate_generation_mutated"),
        "physical_export_permitted": (mechanism_gap.get("overall") or {}).get("physical_export_permitted"),
        "local_continuation_recommended": (mechanism_gap.get("overall") or {}).get("local_continuation_recommended"),
        "primary_missing_channel": (mechanism_gap.get("mechanism_conclusion") or {}).get("primary_missing_channel"),
        "reduced_gap_class": (mechanism_gap.get("reduced_gap") or {}).get("gap_class"),
        "reduced_C_or_D0_multiplier_min": (mechanism_gap.get("reduced_gap") or {}).get("open_stable_required_C_or_D0_multiplier_min"),
        "reduced_C_or_D0_multiplier_median": (mechanism_gap.get("reduced_gap") or {}).get("open_stable_required_C_or_D0_multiplier_median"),
        "nonlinear_gap_class": (mechanism_gap.get("nonlinear_manufactured_gap") or {}).get("gap_class"),
        "nonlinear_C_or_D0_multiplier_min": (mechanism_gap.get("nonlinear_manufactured_gap") or {}).get("open_stable_required_C_or_D0_multiplier_min"),
        "next_physical_requirements": (mechanism_gap.get("mechanism_conclusion") or {}).get("next_physical_requirements"),
    },
}
(output_dir / "_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(f"Wrote {output_dir / '_summary.json'}")
PY
summary_status=$?
if [ "$summary_status" -ne 0 ]; then
  overall_status=1
fi

echo ""
echo "==========================="
echo "Total: $total  |  Pass: $pass  |  Fail: $fail  |  Skipped: 0"
echo "Summary written to: $summary_file"

exit "$overall_status"
