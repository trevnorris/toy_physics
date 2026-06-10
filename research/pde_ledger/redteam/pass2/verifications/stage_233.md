---
unit_id: 233
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
overall_verdict: verified
material_change: false
audit_findings_count: 0
codex_edits: 0
diff_empty: true
relgate:
  sympy: pass   # exit 0, "All Stage 233 symbolic and numerical checks passed"
  mathematica: pass   # exit 0, "All Stage 233 Mathematica checks passed"
  byte_identical_to_committed: true
engines_agree: true
deferred_numbering_band:
  present: true
  comment_only: true   # .py L20/L32 "Stage 239" + L113/L129 "Stage 240/241" print+comment labels; canonical = 188 (compiler) / 224 (budgets, compat point 223)
  affects_computed_value: false
note: >-
  Clean math; zero-script-correction batch confirmed (audit 0 findings, empty
  diff, relgate re-run byte-identical, exit 0 both engines). The .py stale
  forward-renumber self-labels (188->239 compiler, 224->240/241 budgets) are
  print/comment text only — every computed value (Xi1, P1/P0 prefactor identity,
  transported ceiling, robust 0.367930328492646 / nonempty 0.737619063660757
  budget recovery) matches canonical and is corroborated by the independent .wl
  (canonical labels, exact-rational arithmetic). The .wl is genuinely
  independent (Coefficient-extraction route vs .py substitution), not a
  transliteration. Label drift deferred to redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md.
---

# Stage 233 verification (pass-2)

## Summary

`overall_verdict: verified`, `material_change: false`.

Zero-script-correction batch member. Independently re-confirmed:

1. **Report**: `redteam/pass2/reports/stage_233.md` — verdict `clean`, `findings_count: 0`,
   `paper_alignment: aligned`. All 9 reconciliation deliverable values MATCH; 0 misaligned.
2. **Empty diff**: `exec_logs/stage_233_diff.patch` is empty (Codex made no edits) — consistent
   with 0 findings.
3. **Reliability-gate re-run**: both fresh outputs end clean and exit 0 —
   - `relgate/sympy_233.txt` → "All Stage 233 symbolic and numerical checks passed."
   - `relgate/mma_233.txt` → "All Stage 233 Mathematica checks passed." (all M1–M8 PASS;
     M8 numeric residuals ≈ -4.3e-18, -3.2e-16, both below budgetTol 1e-15).
   Both byte-identical to committed.
4. **Engine agreement**: every symbolic residual 0 in both engines; numeric budgets agree
   (SymPy float vs Mathematica exact rational).

## Deferred numbering-band item (NOT a math finding)

The `.py` carries stale forward-renumber self-labels in print banners and comments:
"Stage 239" (the compiler, canonically **Stage 188**) and "Stage 240 / Stage 241" (the
carried budgets, canonically **Stage 224**; compat point **Stage 223**). Confirmed these
appear only in print/comment text — they do not enter any computed quantity, and the
canonical `.wl` (correct Stage-224 attribution) reproduces every number byte-for-byte.
This is the known deferred script/output-band numbering work
(redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md); left as-is per policy.
