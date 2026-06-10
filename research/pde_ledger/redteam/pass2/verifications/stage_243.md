---
unit_id: 243
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-10T00:00:00Z
overall_verdict: verified
material_change: false
checkpoint_bar: met
findings_reviewed: 0
scripts_corrected: false
---

# Verification — unit 243 (checkpoint, zero-script-correction)

Report `stage_243.md` is clean (0 findings, checkpoint). Captured diff `stage_243_diff.patch` is empty (0 bytes) — Codex made no edits. Reliability-gate outputs `sympy_243.txt`/`mma_243.txt` both end all-checks-passed and are byte-identical to the committed outputs.

Checkpoint bar MET: both engines present and substantive. The `.wl` is a genuine independent re-derivation, not a transliteration — it adds the IBP-closure cross-check `Sleak + ibpInterior - boundary == 0` (wl:101), derives U/V via BOTH `Solve` and `LinearSolve[{{kU,-chiLam},{-chiLam,kV}},{fU,0}]` and asserts agreement (wl:131-132), and runs a `Series[…,{x,Infinity,0}]` asymptotic route alongside `Limit` for the short-range span (wl:213-217, 239-251) — vs the `.py`'s `sp.integrate`/`sp.solve`/`sp.limit`. Banner canonical: `STAGE 243 — …` (`.wl` source line 59, output line 3; the 226→243 fix landed). All 13 deliverable values reconcile to `.tex`/`.md`. Confirmed the orchestrator's ground-truth `.wl`-vs-`.py` read finding the re-author SUFFICIENT.

overall_verdict: verified.
