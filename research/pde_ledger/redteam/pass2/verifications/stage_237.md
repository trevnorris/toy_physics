---
unit_id: 237
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
overall_verdict: verified
material_change: false
codex_edits: none
diff_empty: true
relgate:
  sympy_exit: 0
  mathematica_exit: 0
  sympy_byte_identical_to_committed: true
  mathematica_byte_identical_to_committed: true
findings_count: 1
findings_disposition:
  - id: F1
    subtype: paper_missing_script_claim
    class: card-text-lag
    disposition: deferred
    deferred_to: P4-51
    user_gated: true
    script_change: none
self_test_trap_guard:
  present: true
  load_bearing: true
  both_engines: true
---

# Verification — unit 237 (pass-2, ZERO-SCRIPT-CORRECTION batch)

## Summary

Stage 237 (`actual-branch dressing compiler, finite static-blind curve, support-blind post-static orbit-lock theorem`) is **verified** with **no material change**. Codex made zero edits — the captured diff `redteam/pass2/exec_logs/stage_237_diff.patch` is empty (0 bytes). The reliability-gate re-run passed in both engines and is byte-identical to the committed outputs.

## Confirmations

1. **Empty diff** — `exec_logs/stage_237_diff.patch` is 0 bytes; no `.py`/`.wl`/output/doc edits. Consumes 0 corrective work.

2. **Relgate fresh outputs** —
   - `relgate/sympy_237.txt` ends `All Stage 237 symbolic checks passed.` (exit 0); byte-identical to committed `scripts/output/...stage237..._sympy_audit.txt`.
   - `relgate/mma_237.txt` ends `All Stage 237 Mathematica symbolic checks passed.` (exit 0); byte-identical to committed `mathematica/output/...stage237..._mathematica_audit.txt`.
   - All M1–M7 checks report PASS in both engines.

3. **Lone finding (F1) disposition — card-text-lag → DEFERRED to P4-51.**
   F1 is `paper_misalignment` / subtype `paper_missing_script_claim`: card line 11 of `paper/stages/stage_237.tex` says *"Mathematica audit: none yet"* while a complete, independent, all-passing `.wl` audit exists. This is a paper-side documentation lag, **not a script fix**. Directive `redteam/pass2/directives/stage_237.md` marks it `needs_user_resolution: true` (expected resolution (a): name the existing `.wl` on card line 11). Routed to standing user decision / P4-51; no Codex script edit warranted.

4. **Variable-independence self-test-trap fix present & load-bearing in BOTH engines** (re-checked this turn; holds):
   - **SymPy** (`scripts/...stage237..._sympy_audit.py`): explicit presence guard `set(support_args).issubset(q_eta_support_frame.free_symbols)` raises `AssertionError` if any support var is absent (py 252–253); `impose_support_independence` zeros ONLY support-function `Derivative` atoms (py 261–267), so the four `assert_zero` support-blindness checks (py 283–286) can pass only because the diffs first produced live `Derivative` terms that the physical-blindness rule kills; negative-control comment py 272–273.
   - **Mathematica** (`mathematica/...stage237..._mathematica_audit.wl`): TWO guards — presence guard `Not[FreeQ[qEtaSupport, #]]` over `{zeta, Mtr, lambdaPhi, KPhiEff}` (wl 231–234) AND the leak-detector `Not[FreeQ[#, Derivative]]` over `rawSupportDerivatives` (wl 243–246) which `fail`s if any raw derivative is trivially Derivative-free (dead channel); `supportBlindRules` (wl 248–252) zero only support-function derivatives; negative-control comment wl 254–255.
   - CONCLUSION: every `D[]`/`diff` differentiates w.r.t. a variable that genuinely appears in the abstract support functions; no vacuous derivative. The pass-1 fix (live-channel negative control + `Not[FreeQ[#, Derivative]]` leak detector) is intact, not removed or weakened.

## Verdict

`overall_verdict: verified`, `material_change: false`. Math sound and exactly aligned across both engines + Part VII appendix; the only defect is the card-text-lag F1 (deferred P4-51, user-gated, no script change). Self-test-trap guard confirmed present and load-bearing in both engines.
