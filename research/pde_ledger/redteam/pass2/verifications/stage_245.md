---
unit_id: 245
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-10T00:00:00Z
overall_verdict: verified
material_change: false
findings_reviewed: 0
diff_reviewed: empty
relgate:
  sympy: exit 0, all symbolic checks passed
  mathematica: exit 0, M1-M8 all PASS
  engines_agree: true
---

# Stage 245 verification (pass-2)

Zero-script-correction unit: audit returned `clean` / 0 findings, captured diff `exec_logs/stage_245_diff.patch` is empty (0 bytes, Codex made no edits), so there is no material change to verify. Reliability-gate reruns are fresh and pass: `relgate/sympy_245.txt` ends "All symbolic checks passed." and `relgate/mma_245.txt` ends "Stage 245 Mathematica audit passed." with every M1-M8 check PASS, both exit 0. Spot-confirmed the load-bearing F1 self-test-trap guard is live in both engines — the M7/Section-6 support-blindness derivatives w.r.t. `Lam,varrho` are zero by genuine variable absence, backstopped by the positive control on `f_U+Lam` that asserts a nonzero derivative `= a_V/Delta_UV` (sympy out L74 / mma out L113-117). Nothing to apply; trivially verified.
