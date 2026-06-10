---
unit_id: 248
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-10T00:00:00Z
overall_verdict: verified
material_change: false
checkpoint_bar: met
findings_reviewed: 0
relgate:
  sympy_exit: 0
  mma_exit: 0
  byte_identical_to_committed: true
diff_patch: empty
---

# Stage 248 verification (pass-2, checkpoint)

Confirmed: zero-script-correction batch. Diff patch is empty (0 bytes); Codex made no edits. Both relgate outputs end "All checks passed" and are byte-identical to the committed `output/*.txt`.

Checkpoint bar MET — both engines present and substantive; `.wl` is genuinely independent (satisfaction route): it POSITS `vcritNew = Sqrt[2(Vpeak-V0)/ms]` (L97) and verifies it SATISFIES the defining energy equality `EAtVcrit - Vpeak == 0` / `deltaNew /. v0 -> vcritNew == 0` (L121-122), gated by the non-vacuity guard `FreeQ[deltaNew, v0] → fail` (L109-112) — output L15 prints the v0-dependent pre-substitution gap, so the guard is real and the satisfaction non-vacuous. The `.py` instead `sp.solve`s for `v0` (L70-75) — opposite information flow, not a transliteration. ZERO surviving 168 (grep clean across scripts+outputs; ×168%→×100% fix holds); banner canonical (`STAGE 248`, `231→248` landed). Orchestrator ground-truth read corroborated.
