---
unit_id: 247
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-10T00:00:00Z
overall_verdict: verified
material_change: false
findings_reviewed: 0
relgate:
  sympy_exit: 0
  mathematica_exit: 0
  diff_bytes: 0
---

# Verification unit 247 (pass-2)

Zero-script-correction stage: audit report carries 0 findings, captured diff
`stage_247_diff.patch` is empty (0 bytes), Codex made no edits. Trivially
verified. Sanity-confirmed the fresh relgate outputs both end clean and exit 0:
SymPy "All symbolic checks passed." (sympy_247.txt L71); Mathematica
"All Mathematica checks passed." (mma_247.txt L88).

Note: 142.17750000 holds in both engines (210.17750000 purged everywhere) and
the published card body is correctly abstract (no Delta numeric) — card clean.
