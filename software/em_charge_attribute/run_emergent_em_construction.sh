#!/usr/bin/env bash
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="$here/reports/artifacts"
mkdir -p "$out"

timeout 600 python3 "$here/emergent_em_sympy.py" --out-dir "$out" \
  >"$out/sympy_run.log" 2>&1
timeout 600 wolframscript -file "$here/emergent_em_dual.wl" \
  >"$out/mathematica_run.log" 2>&1
timeout 600 python3 "$here/compare_engines.py" --out-dir "$out" \
  >"$out/engine_compare_run.log" 2>&1

printf '%s\n' "RUNNER_STATUS: PASS"
printf '%s\n' "SymPy: exit 0"
printf '%s\n' "Mathematica: exit 0"
printf '%s\n' "Comparator: exit 0"
