#!/bin/bash
# Installs the measurements gate (CLAUDE.md E1, orchestrator half) in this repository.
#
#   bash .claude/hooks/install_measurements_gate.sh
#
# 1. Records the baseline, once: the commit every ref points at now. The gate treats
#    those commits, and everything they reach, as already judged, so history from
#    before the gate is not re-judged. An existing baseline is never overwritten,
#    because refreshing it would grandfather whatever arrived since.
# 2. Links git's reference-transaction hook to require_measurements_gate.sh. Linked
#    worktrees share the main checkout's hooks, so run this once per clone.

set -euo pipefail

common=$(git rev-parse --path-format=absolute --git-common-dir)
main=${common%/.git}                 # the main checkout; linked worktrees use its gate
hooks=$(git rev-parse --path-format=absolute --git-path hooks)
baseline="$common/measurements-gate-baseline"
gate="$main/.claude/hooks/require_measurements_gate.sh"

[ -x "$gate" ] || { echo "not executable or missing: $gate" >&2; exit 1; }

if [ -e "$baseline" ]; then
  echo "baseline kept: $baseline ($(wc -l < "$baseline") commits)"
else
  tmp=$(mktemp "$common/measurements-gate-baseline.XXXXXX")
  git for-each-ref --format='%(objectname)' \
    | while read -r o; do git rev-parse -q --verify "$o^{commit}" 2>/dev/null || true; done \
    | sort -u > "$tmp"
  mv "$tmp" "$baseline"
  echo "baseline written: $baseline ($(wc -l < "$baseline") commits)"
fi

link="$hooks/reference-transaction"
if [ -L "$link" ] && [ "$(readlink -f "$link")" = "$(readlink -f "$gate")" ]; then
  echo "hook already linked: $link"
elif [ -e "$link" ] || [ -L "$link" ]; then
  echo "refusing to replace an existing hook: $link -> $(readlink -f "$link")" >&2
  exit 1
else
  ln -s "$gate" "$link"
  echo "hook linked: $link -> $gate"
fi
