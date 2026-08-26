#!/bin/bash
# Regenerates the §3c-clarification fidelity block appended to _measurements/S11c_a_SHARED_PHYSICS.md.
# rule 2: the §3c note asserts structural PREMISES; this twin carries the spec sources they are faithful to.
set -uo pipefail
cd "${CLAUDE_PROJECT_DIR:-/var/projects/toy_physics}"
SP=research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
claim () { printf '## %s\n\n```\n$ %s\n' "$1" "$2"; eval "$2"; printf '```\n\n'; }
cat <<'HDR'

# --- §3c shifted-trace clarification (appended 2026-08-26) ---

The §3c note added this round states only structural premises; it asserts no computed result ("Which trace
terms then survive is computed ... not stated here"). Below: the note, and the spec sources each premise
rests on. Verified independently by two residual-0 CAS derivations (see S11c_a_py_shifted_trace_fix_directive.md).
HDR
claim "the §3c clarification note" "sed -n '385,392p' $SP"
claim "premise: background face-normal bulk vel/flux vanish (§2d)" "sed -n '261p' $SP"
claim "premise: the drain is an inert rest-frame scope limit (§1)" "sed -n '48,49p' $SP"
claim "premise: perturbation pressure has zero background (§3b)" "sed -n '336p' $SP"
claim "premise: background face carries eta via W_bg (§2a); mixed-grade retained (§2a truncation)" "sed -n '176p' $SP; echo '---'; sed -n '195,198p' $SP"
