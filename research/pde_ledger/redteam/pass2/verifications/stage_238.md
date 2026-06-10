---
unit_id: 238
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
overall_verdict: verified
material_change: false
findings_reviewed: 0
edits_made: 0
diff_empty: true
relgate:
  sympy: { exit: 0, end: "All Stage 238 symbolic checks passed." }
  mathematica: { exit: 0, end: "All Stage 238 Mathematica checks passed." }
---

# Stage 238 verification

Zero-script-correction stage: audit (`reports/stage_238.md`) returned 0 findings,
Codex made no edits, captured diff `exec_logs/stage_238_diff.patch` is empty (0 bytes).
Reliability-gate re-run passed both engines, byte-identical to committed.

Re-checked the iter-2 reframe holds in BOTH fresh relgate outputs:

- SymPy (`exec_logs/relgate/sympy_238.txt`): negative control present and live
  ("live support channel derivative with respect to zeta/M_mix" + "live support
  channel depends on zeta and M_mix", lines 11-13), leak detector "support leak
  detector for R_tr" (14), structural exclusion "structural support exclusion for
  reduced observables and packet" (15). The vacuous support-blindness diffs follow
  and are load-bearing only as a set with the three guards. Ends
  "All Stage 238 symbolic checks passed."
- Mathematica (`exec_logs/relgate/mma_238.txt`): negative control witnesses are
  genuine nonzero rationals (zeta witness `-(((-1+eps)*Mmix)/(-1+eps*zeta)^2)`,
  Mmix witness `(1+zeta-2*eps*zeta)/(1-eps*zeta)`), leak-detector witness nonzero
  `-(((1+chi0+deltaU)*(-1+eps))/((1+chi0)*(1+deltaU)*(-1+eps*zeta)^2))`, six
  structural exclusions all True (lines 6-23). Ends "All Stage 238 Mathematica
  checks passed."

The negative control / leak detector would genuinely FIRE (raise / Exit[1]) if
support-blindness were faked (constant M_tr ⇒ derivative 0), confirming the
support-blindness evidence is not vacuous. No tautology survives without an
exercised counterpart; fully symbolic, no hardcoded numerics; paper↔script
alignment exact (only cosmetic naming differences).

note: clean; iter-2 reframe fix holds — neg-control + leak-detector + exclusion present/load-bearing both engines.
