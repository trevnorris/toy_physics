#!/usr/bin/env bash
set -euo pipefail
exec 9>"${B2_WOLFRAM_LOCK:-/tmp/u1_body_b2_wolframscript.lock}"
flock -x 9

# WolframKernel may open the -file payload with a direct syscall that bypasses
# the runner's libc preload hook.  Read every regular-file argument once in
# this runner process before dispatch so the same hook records its content
# digest at first use.  This also covers phase/input artifacts uniformly.
for argument in "$@"; do
  if [[ -f "$argument" ]]; then
    digest="$(/usr/bin/sha256sum -- "$argument")"
    metadata="$(/usr/bin/stat -Lc '%d %i' -- "$argument")"
    read -r device inode <<<"$metadata"
    printf '%s\tOPEN\t%s\t%s\t%s\t%s\n' \
      "$$" "${digest%% *}" "$device" "$inode" "$argument" \
      >>"${B2_FIRST_USE_LOG:?B2_FIRST_USE_LOG is required}"
  fi
done

# The protected SymPy phase-A baseline is a fixed cross-engine input computed
# inside u1_body_mechanics_dual.wl, rather than a native command-line input.
# Make that one implicit dependency runner-visible for both the protected
# script and the mutation harness's temporary source copies.
for argument in "$@"; do
  case "$(basename -- "$argument")" in
    u1_body_mechanics_dual.wl|.mutation_*_u1_body_mechanics_dual.wl)
      baseline="$(dirname -- "$argument")/reports/u1_body_dynamics_artifacts/sympy_phase_a.json"
      if [[ -f "$baseline" ]]; then
        digest="$(/usr/bin/sha256sum -- "$baseline")"
        metadata="$(/usr/bin/stat -Lc '%d %i' -- "$baseline")"
        read -r device inode <<<"$metadata"
        printf '%s\tOPEN\t%s\t%s\t%s\t%s\n' \
          "$$" "${digest%% *}" "$device" "$inode" "$baseline" \
          >>"${B2_FIRST_USE_LOG:?B2_FIRST_USE_LOG is required}"
      fi
      ;;
  esac
done

if [[ "${B2_CAPTURE_WOLFRAM_ASSERT_FROM_ARTIFACT:-0}" == 1 ]]; then
  capture="/tmp/u1_b2_wolfram_stdout.$$"
  : >"$capture"
  set +e
  "${WOLFRAMSCRIPT_REAL:?WOLFRAMSCRIPT_REAL is required}" -local "${WOLFRAM_KERNEL_REAL:?WOLFRAM_KERNEL_REAL is required}" "$@" >"$capture" 2>&1
  code=$?
  set -e
  while IFS= read -r line || [[ -n "$line" ]]; do
    printf '%s\n' "$line"
  done <"$capture"
  if [[ "$code" == 1 && ! -s "$capture" ]]; then
    output=""
    previous=""
    for argument in "$@"; do
      if [[ "$previous" == "--output" ]]; then output="$argument"; break; fi
      previous="$argument"
    done
    if [[ -f "$output" ]]; then
      /usr/bin/python3 - "$output" <<'PY'
import json
import sys

with open(sys.argv[1], "r", encoding="utf-8") as stream:
    artifact = json.load(stream)
for tooth, record in artifact.get("checks", {}).items():
    if record.get("status") != "PASS":
        print(f"ASSERT_FAIL:{tooth}:{record.get('evidence_digest', 'artifact_recovered_failure')}")
        break
PY
    fi
  fi
  /bin/rm -f "$capture"
  exit "$code"
fi

exec "${WOLFRAMSCRIPT_REAL:?WOLFRAMSCRIPT_REAL is required}" -local "${WOLFRAM_KERNEL_REAL:?WOLFRAM_KERNEL_REAL is required}" "$@"
