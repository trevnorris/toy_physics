#!/usr/bin/env bash
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
artifacts="$here/reports/u1_body_dynamics_artifacts"
stage1="$artifacts/stage1"
report="$here/reports/u1_body_dynamics.md"
results="$here/reports/u1_body_dynamics_results.yaml"
stage1_code=42

run_stage1() {
  rm -f "$artifacts/b1_mutation_results.yaml" "$artifacts/b1_engine_agreement.yaml"
  rm -rf "$stage1"
  mkdir -p "$stage1"

  timeout 600 python3 "$here/u1_body_dynamics_sympy.py" --output "$stage1/sympy_phase_a.json" \
    >"$stage1/sympy_phase_a.log" 2>&1
  timeout 600 wolframscript -file "$here/u1_body_dynamics_dual.wl" --output "$stage1/mathematica_phase_a.json" \
    >"$stage1/mathematica_phase_a.log" 2>&1

  timeout 600 python3 "$here/u1_body_mechanics_sympy.py" \
    --phase-artifact "$stage1/sympy_phase_a.json" --output "$stage1/sympy_phase_b1.yaml" \
    >"$stage1/sympy_phase_b1.log" 2>&1
  timeout 600 wolframscript -file "$here/u1_body_mechanics_dual.wl" \
    --phase-artifact "$stage1/mathematica_phase_a.json" --output "$stage1/mathematica_phase_b1.yaml" \
    >"$stage1/mathematica_phase_b1.log" 2>&1

  timeout 600 python3 "$here/u1_body_mechanics_compare.py" --input "$here/u1_body_mechanics_inputs.yaml" \
    --sympy-artifact "$stage1/sympy_phase_b1.yaml" --math-artifact "$stage1/mathematica_phase_b1.yaml" \
    --phase-a-artifact "$stage1/sympy_phase_a.json" --verify-only \
    --stage1-complete "$stage1/stage1_complete.yaml" >"$stage1/comparator.log" 2>&1

  test -s "$stage1/stage1_complete.yaml"
  timeout 60 python3 - "$stage1/stage1_complete.yaml" "$stage1/phase_a_amendment_agreement.yaml" <<'PY'
import sys, yaml
from pathlib import Path
source, target = map(Path, sys.argv[1:])
row = yaml.safe_load(source.read_text())
agreement = {
    "schema_version": "U1_PHASE_A_AMENDMENT_AGREEMENT_V3",
    "digest_gate": row["digest_gate"],
    "semantic_diff_gate": row["semantic_diff_gate"],
    "phase_a_acceptance_recheck": row["phase_a_acceptance_recheck"],
    "correction_finding": row["correction_finding"],
}
target.write_text(yaml.safe_dump(agreement, sort_keys=False, allow_unicode=True, width=220))
PY
  exit "$stage1_code"
}

run_stage2() {
  rm -f "$artifacts/b1_mutation_results.yaml" "$artifacts/b1_engine_agreement.yaml"
  test -s "$stage1/stage1_complete.yaml"
  timeout 60 python3 - "$stage1/stage1_complete.yaml" <<'PY'
import hashlib, json, sys, yaml
from pathlib import Path
path = Path(sys.argv[1]); row = yaml.safe_load(path.read_text())
expected = row.pop("resume_validation_sha256")
actual = hashlib.sha256(json.dumps(row, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()
if actual != expected or not row["digest_gate"]["agreement"] or row["final_b1_outputs_emitted"]:
    raise SystemExit(1)
PY

  core_batches=5
  core_batch_paths=()
  for batch in 0 1 2 3 4; do
    path="$stage1/b1_mutation_core_${batch}.yaml"
    core_batch_paths+=("$path")
    if [[ ! -s "$path" ]]; then
      timeout 600 python3 "$here/u1_body_mechanics_mutations.py" --artifacts "$stage1" --core-only \
        --core-batch "$batch/$core_batches" --output "$path" >"$stage1/mutations_core_${batch}.log" 2>&1
    fi
  done
  if [[ ! -s "$stage1/b1_mutation_core.yaml" ]]; then
    timeout 600 python3 "$here/u1_body_mechanics_mutations.py" --artifacts "$stage1" --core-only \
      --core-batch-results "${core_batch_paths[@]}" \
      --output "$stage1/b1_mutation_core.yaml" >"$stage1/mutations_core.log" 2>&1
  fi
  for batch in 0 1 2; do
    if [[ ! -s "$stage1/b1_liveness_${batch}.yaml" ]]; then
      timeout 600 python3 "$here/u1_body_mechanics_mutations.py" --artifacts "$stage1" \
        --core-results "$stage1/b1_mutation_core.yaml" --liveness-batch "$batch/3" \
        --output "$stage1/b1_liveness_${batch}.yaml" >"$stage1/mutations_liveness_${batch}.log" 2>&1
    fi
  done
  timeout 600 python3 "$here/u1_body_mechanics_mutations.py" --artifacts "$stage1" \
    --core-results "$stage1/b1_mutation_core.yaml" \
    --liveness-results "$stage1/b1_liveness_0.yaml" "$stage1/b1_liveness_1.yaml" "$stage1/b1_liveness_2.yaml" \
    --output "$stage1/b1_mutation_results.yaml" >"$stage1/mutations.log" 2>&1
  timeout 600 python3 "$here/u1_body_mechanics_compare.py" --artifacts "$stage1" \
    --sympy-artifact "$stage1/sympy_phase_b1.yaml" --math-artifact "$stage1/mathematica_phase_b1.yaml" \
    --phase-a-artifact "$stage1/sympy_phase_a.json" --mutation-results "$stage1/b1_mutation_results.yaml" \
    --report "$report" --results "$results" >"$stage1/stage2_comparator.log" 2>&1
  timeout 60 python3 "$here/u1_body_mechanics_compare.py" --results "$results" --print-summary
}

case "${1:-}" in
  "") run_stage1 ;;
  --resume-stage2) run_stage2 ;;
  *) printf 'usage: %s [--resume-stage2]\n' "$0" >&2; exit 64 ;;
esac
