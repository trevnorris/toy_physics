#!/usr/bin/env bash
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
artifacts="$here/reports/u1_body_dynamics_artifacts"
sym_out="$artifacts/sympy_phase_a.json"
math_out="$artifacts/mathematica_phase_a.json"
report="$here/reports/u1_body_dynamics.md"
results="$here/reports/u1_body_dynamics_results.yaml"
teeth=(
  INPUT_LEDGER SOURCE_COMPLETENESS DIMENSIONAL COMOVING_CONTINUITY
  COMOVING_MOMENTUM BASE_BALANCE TAIL_ODE ZERO_MODE PROJECTOR
  ENDPOINT_RESPONSE MOMENT_INTEGRALS RECONSTRUCTION CANONICAL_VARIATION
  CHANNEL_UNIQUENESS TYPED_DATAFLOW PROVENANCE_FORBIDDEN ANCESTRY PARITY
  NATIVE_PADDING OUTCOME_REACHABILITY
)

asserted_tooth() {
  local log="$1"
  if awk -F: '$1=="ASSERT_FAIL" && $2=="MUTATION_NOOP" {found=1} END {exit !found}' "$log"; then
    return 1
  fi
  awk -F: '$1=="ASSERT_FAIL" && $2!="MUTATION_NOOP" {print $2; exit}' "$log"
}

mkdir -p "$artifacts/sympy_ablation_logs" "$artifacts/math_ablation_logs"
rm -f "$sym_out" "$math_out" "$artifacts/engine_agreement.log" \
  "$artifacts/sympy_run.log" "$artifacts/mathematica_run.log" \
  "$artifacts/comparator_run.log" "$artifacts"/comparator_ablation_*.log \
  "$artifacts/sympy_ablation_logs"/*.log "$artifacts/math_ablation_logs"/*.log \
  "$report" "$results"

# R6 harness self-probe: a noop sentinel must never be accepted as a tooth.
r6_probe="$artifacts/r6_mutation_noop_rejection.log"
printf 'ASSERT_FAIL:MUTATION_NOOP:TAIL_ODE:synthetic harness probe\n' >"$r6_probe"
if [[ "$(asserted_tooth "$r6_probe" || true)" == "TAIL_ODE" ]]; then
  printf 'R6_HARNESS_FAILURE: mutation noop accepted\n' >&2
  exit 70
fi
printf 'R6_MUTATION_NOOP_REJECTION: PASS\n'

timeout 600 python3 "$here/u1_body_dynamics_sympy.py" --output "$sym_out" \
  >"$artifacts/sympy_run.log" 2>&1

run_sympy_ablation() {
  local tooth="$1"
  log="$artifacts/sympy_ablation_logs/$tooth.log"
  set +e
  timeout 600 python3 "$here/u1_body_dynamics_sympy.py" --mutation "$tooth" >"$log" 2>&1
  code=$?
  set -e
  asserted="$(asserted_tooth "$log" || true)"
  if [[ "$code" -ne 1 || "$asserted" != "$tooth" ]]; then
    printf 'SYMPY_ABLATION_FAILURE:%s expected_exit=1 actual=%s\n' "$tooth" "$code" >&2
    exit 71
  fi
  printf 'SYMPY_ABLATION:%s PASS\n' "$tooth"
}

# Python controls use no Wolfram seats; four-way batching keeps the complete
# runner inside the fixed outer 600-second process limit.
pids=()
for tooth in "${teeth[@]}"; do
  run_sympy_ablation "$tooth" &
  pids+=("$!")
  if [[ "${#pids[@]}" -eq 4 ]]; then
    for pid in "${pids[@]}"; do wait "$pid"; done
    pids=()
  fi
done
for pid in "${pids[@]}"; do wait "$pid"; done

# The Wolfram engine and all of its ablations are sequential: one seat maximum.
timeout 600 wolframscript -file "$here/u1_body_dynamics_dual.wl" --output "$math_out" \
  >"$artifacts/mathematica_run.log" 2>&1

for tooth in "${teeth[@]}"; do
  log="$artifacts/math_ablation_logs/$tooth.log"
  set +e
  timeout 600 wolframscript -file "$here/u1_body_dynamics_dual.wl" --mutation "$tooth" >"$log" 2>&1
  code=$?
  set -e
  asserted="$(asserted_tooth "$log" || true)"
  if [[ "$code" -ne 1 || "$asserted" != "$tooth" ]]; then
    printf 'MATHEMATICA_ABLATION_FAILURE:%s expected_exit=1 actual=%s\n' "$tooth" "$code" >&2
    exit 72
  fi
  printf 'MATHEMATICA_ABLATION:%s PASS\n' "$tooth"
done

for tooth in ENGINE_CANONICAL ENGINE_DEPENDENCIES; do
  log="$artifacts/comparator_ablation_${tooth}.log"
  set +e
  timeout 600 python3 "$here/u1_body_dynamics_compare.py" --artifacts "$artifacts" --mutation "$tooth" >"$log" 2>&1
  code=$?
  set -e
  asserted="$(asserted_tooth "$log" || true)"
  if [[ "$code" -ne 1 || "$asserted" != "$tooth" ]]; then
    printf 'COMPARATOR_ABLATION_FAILURE:%s expected_exit=1 actual=%s\n' "$tooth" "$code" >&2
    exit 73
  fi
  printf 'COMPARATOR_ABLATION:%s PASS\n' "$tooth"
done

timeout 600 python3 "$here/u1_body_dynamics_compare.py" --artifacts "$artifacts" \
  --report "$report" --results "$results" >"$artifacts/comparator_run.log" 2>&1

# The runner owns no summary literals: all ten lines are read from results.yaml.
timeout 600 python3 "$here/u1_body_dynamics_compare.py" --results "$results" --print-summary
