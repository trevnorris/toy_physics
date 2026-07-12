#!/usr/bin/env bash
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="$here/reports/native_p_gate_artifacts"
report="$here/reports/native_p_constraint_gate.md"
mkdir -p "$out"
rm -f "$out/sympy_results.json" "$out/mathematica_results.json" \
  "$out/engine_agreement.log" "$out/comparator_run.log"

timeout 600 python3 "$here/native_p_gate_sympy.py" --out-dir "$out" \
  >"$out/sympy_run.log" 2>&1

declare -A expected=(
  [maxwell]=41
  [gauged_hard_unit]=42
  [bare_sigma]=43
  [nonconserved_current]=44
  [gauge_fixed_maxwell]=45
  [global_only]=46
)

for tooth in maxwell gauged_hard_unit bare_sigma nonconserved_current gauge_fixed_maxwell global_only; do
  set +e
  timeout 600 python3 "$here/native_p_gate_sympy.py" --ablate-tooth "$tooth" \
    >"$out/ablation_${tooth}.log" 2>&1
  code=$?
  set -e
  if [[ "$code" -ne "${expected[$tooth]}" ]]; then
    printf '%s\n' "ABLATION_RUNNER_FAILURE: $tooth expected=${expected[$tooth]} actual=$code" >&2
    exit 70
  fi
done

# Exactly one Wolfram kernel is launched by a runner invocation.
timeout 600 wolframscript -file "$here/native_p_gate_dual.wl" \
  >"$out/mathematica_run.log" 2>&1

timeout 600 python3 "$here/native_p_gate_compare.py" --out-dir "$out" --report "$report" \
  >"$out/comparator_run.log" 2>&1

printf '%s\n' "RUNNER_STATUS: PASS"
printf '%s\n' "SymPy genuine quadratic Hamiltonian/Dirac engine: exit 0"
printf '%s\n' "Six per-tooth ablations: fired at expected own-point exit codes"
printf '%s\n' "Mathematica independent quadratic Hamiltonian/Dirac engine: exit 0"
printf '%s\n' "Computed guards: GUARD-COUPLINGS-ENTER and GUARD-SEARCH-CAPABLE PASS"
printf '%s\n' "Computed hardening guard: HARDENING-TUNED-DESCENDANT-REJECTION PASS"
printf '%s\n' "Comparator: ENGINE_AGREE theories and six controls"
printf '%s\n' "Report: $report"
