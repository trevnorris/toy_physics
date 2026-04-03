#!/usr/bin/env bash
#
# Run one Mathematica moving-throat audit script and save its output.
#
# Usage:
#   bash mathematica/moving_throat/run_one_audit.sh stage118
#   bash mathematica/moving_throat/run_one_audit.sh moving_throat_pde_stage118_outlet_consistent_mouth_closure_mathematica_audit
#   bash mathematica/moving_throat/run_one_audit.sh /abs/path/to/script.wl

set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: bash mathematica/moving_throat/run_one_audit.sh <stage-token|basename|script-path>" >&2
  exit 2
fi

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SCRIPTS_DIR="$REPO_ROOT/mathematica/moving_throat"
OUTPUT_DIR="$SCRIPTS_DIR/output"
ARG="$1"

mkdir -p "$OUTPUT_DIR"

if [[ "$ARG" = /*.wl ]]; then
  script="$ARG"
else
  matches=()
  while IFS= read -r line; do
    matches+=("$line")
  done < <(find "$SCRIPTS_DIR" -maxdepth 1 -type f -name "moving_throat_pde_stage*_mathematica_audit.wl" | sort | grep -F "$ARG" || true)

  if [ "${#matches[@]}" -eq 0 ]; then
    echo "no Mathematica audit matched: $ARG" >&2
    exit 2
  fi
  if [ "${#matches[@]}" -gt 1 ]; then
    echo "multiple Mathematica audits matched: $ARG" >&2
    printf '  %s\n' "${matches[@]}" >&2
    exit 2
  fi
  script="${matches[0]}"
fi

basename="$(basename "$script" .wl)"
out_file="$OUTPUT_DIR/${basename}.txt"
tmp_file="$(mktemp "$OUTPUT_DIR/${basename}.XXXXXX.tmp")"

printf "RUN   %s\n" "$basename"
sleep 2
start_time=$(date +%s)

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

cat "$out_file"
rm -f "$tmp_file"

exit "$exit_code"
