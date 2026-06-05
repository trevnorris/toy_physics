---
unit_id: 078
batch: III.4
verdict: verified
overall: verified
material_change: false
verified_at: 2026-06-05
---

# Verification — unit 078 (pass 2)

**Verdict: verified — clean.** Audit verdict was `clean` (0 findings). Both engines
reproduce the four boxed Pe windows (`96.5285…`, `11220.544…`, `22.0062…`, `2558.019…`)
and the ordering/overlap inequalities; 8 deliverable/input values reconcile exactly with
`stage_078.tex` and notes, 0 misaligned. Mathematica is genuinely independent on the
threshold side (symbolic cosh/sinh re-derivation of the Stage-075 coefficients, not a float
transliteration).

Orchestrator independent re-run: SymPy exit 0, Mathematica exit 0; committed outputs
byte-identical to HEAD; arbiter grep clean (no self-epoch 61 banner/closing). The SymPy
docstring's stale `Stage-60`/`Stage-58` CROSS-references are label-only and owned by the
deferred SCRIPT/OUTPUT-band numbering pass (not self-labels) — intentionally left untouched,
not a defect. `material_change: false`.
