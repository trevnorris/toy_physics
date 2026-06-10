---
unit_id: 239
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
overall_verdict: verified
material_change: false
checkpoint_bar: met
findings_reviewed: 0
relgate:
  sympy_exit: 0
  mathematica_exit: 0
  outputs_match_committed: true
diff_empty: true
---

# Verify unit 239 (CHECKPOINT) — verified

Zero-script-correction checkpoint batch. Audit report (`reports/stage_239.md`) carries 0 findings;
captured diff `exec_logs/stage_239_diff.patch` is empty (0 bytes); Codex made no edits.

Checkpoint bar MET: both engines present and substantive. The `.wl` is a genuinely independent
route — it DERIVES the compiler matrix via a `D[]`-Jacobian, the left inverse via native
`PseudoInverse`, and the orbit lock via `Reduce`/`Equivalent`, whereas the `.py` POSITS the compiler
matrix/left-inverse as literals and uses `sp.solve`. Orchestrator already performed the ground-truth
`.wl`-vs-`.py` read and confirmed the pass-1 re-author (away from a transliteration) is SUFFICIENT;
report's independent-derivation section corroborates the distinct information flows.

Reliability-gate re-run clean: `exec_logs/relgate/sympy_239.txt` (52 [ok], "All Stage 239 symbolic
checks passed.") and `mma_239.txt` (M1–M9 all PASS, banner "STAGE 239 — RIGID-MOUTH PHYSICAL NORMAL
FORM", "All Stage 239 symbolic checks passed.") both exit 0, byte-identical to committed.

note: independent route confirmed (Jacobian / PseudoInverse / Reduce-Equivalent vs. posited literals + sp.solve), re-author sufficient.
