---
unit_id: 244
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-10T00:00:00Z
overall_verdict: verified
material_change: false
codex_edits: none
diff_empty: true
reliability_gate:
  sympy_exit: 0
  mathematica_exit: 0
  sympy_byte_identical_to_committed: true
  mathematica_byte_identical_to_committed: true
findings_count: 1
findings_disposition:
  - id: F1
    subtype: paper_missing_script_claim
    severity: low
    class: card-text-lag
    disposition: deferred
    routed_to: PAPER_CLEANUP P4-51
    rationale: >
      Card \stagefield{Verification} line 4 says "Mathematica audit: none yet"
      despite a present, passing .wl (M1-M9 all PASS). Cosmetic/prose only; no
      math affected. Per standing user decision the card-text-lag class is a
      paper-side cleanup, not a script fix — the stage still verifies dual-engine.
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Stage 244 — pass-2 verification

## Summary

Zero-script-correction unit. Codex applied no edits; the captured diff
`redteam/pass2/exec_logs/stage_244_diff.patch` is empty (0 bytes). The lone audit
finding (F1) is the card-text-lag class and is deferred to PAPER_CLEANUP P4-51 per
the standing user decision — it is a paper-side prose lag, not a script defect.

## Reliability gate (re-run, fresh)

- `redteam/pass2/exec_logs/relgate/sympy_244.txt`: exit 0, ends "All symbolic
  checks passed.", no FAIL/error/traceback markers; byte-identical to committed
  `scripts/output/...stage244..._sympy_audit.txt`.
- `redteam/pass2/exec_logs/relgate/mma_244.txt`: exit 0, ends "All Stage 244
  Mathematica checks passed.", no FAIL/$Failed/error markers; byte-identical to
  committed `mathematica/output/...stage244..._mathematica_audit.txt`.

Both engines agree (S_leak, W_bulk, W_sess, Pi_tr, compiled (1-eps)/varrho forms,
support/orbit split, parity, recovery at eta=0).

## Findings disposition

### F1 — paper_missing_script_claim (card-text-lag) → DEFERRED P4-51

The stage_244.tex Verification line states "Mathematica audit: none yet" while an
independent, passing `.wl` (M1-M9) is present. Cosmetic prose under-reporting of
coverage; no math impact. Routed to PAPER_CLEANUP P4-51 (paper-side edit, deferred)
per standing user decision on the card-text-lag class. Stage still verifies.

## Re-checks confirmed (hold)

- The 128√2 notes correction holds (notes line 366 reads `128√2`, matching script
  line 106 / output line 55; consistent across the (1-eps) and varrho forms).
- The F1 variable-independence self-test-trap fix from pass-1 is present and
  load-bearing: the support/orbit split asserts `orbit_syms.isdisjoint(free)`
  (py 141) guarded by the positive anti-vacuity control
  `support_syms.issubset(free)` (py 143), with non-empty support coverage
  `[Lam, eta_leak, varrho]` in the output (no vacuous derivative trap). The .wl
  M7 mirrors with `Not[FreeQ[...]]`.

## Verdict

overall_verdict: verified — math holds dual-engine under attack, no Codex edits,
empty diff, fresh relgate outputs both exit 0 and byte-identical to committed; the
sole finding is the deferred card-text-lag (P4-51), so material_change: false.
