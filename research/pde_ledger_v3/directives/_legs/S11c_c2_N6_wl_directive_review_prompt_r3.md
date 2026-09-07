# Decision review (ROUND 3) — the blind Wolfram N6 build directive, after the round-2 fold

## Artifact
`research/pde_ledger_v3/directives/S11c_c2_N6_wl_build_directive.md` (orchestrator-written; physics-bearing build
directive for a NEW blind Mathematica engine). **Round 3.** Round 1 (12 issues) and round 2 (4 issues, incl. 2
regressions bred by round-1 folds) both returned FOLD-REQUIRED; all folded. Review the CURRENT directive **until clear**:
(a) a fresh pass against the sources with full rigor (⛔ a fold can breed a new defect — do not assume the fold is
correct); (b) confirm each of the four round-2 folds landed and bred no new regression. DECISION review of a directive —
no `.wl` exists, so ⛔ no fictional-script ablation. ⛔ You are reviewing the DIRECTIVE, ⛔ not re-adjudicating the settled
covariance verdict (Reading B).

## The four round-2 folds to confirm (recorded in `_measurements/S11c_c2_N6_wl_directive_review_adjudication.md`)
1. **Reviewer-side reclassification** — the two SymPy N6 build directives (`reconcile`, `covariance`) and `RESOLVED.md`
   are now REVIEWER-SIDE only, ⛔ not builder-facing authorities. Confirm the builder-facing authority list is exactly
   `S11c_c2_SHARED_PHYSICS.md` + route-2 spec + sibling SHARED_PHYSICS specs + this directive's translation, and that the
   directive text itself no longer NAMES a computed outcome anywhere (not even in the rationale for the exclusion).
2. **`R_cov` velocity channel** — the "commutes iff Φ θ-independent" sentence is gone; `R_cov` is stated to measure the
   complete combined-source square incl. independently-computed `V_E`,`V_M`; μ-isolation is downstream only. Confirm this
   is faithful to the covariance directive `:40-46` and route-2 `:157-159`, and answer-neutral.
3. **Slot-guard namespace** — the guards are now `WL_S11CC2_N6_SLOT_GUARD_{NATIVE,CARRIER,RESIDUAL}` +
   `CLOSURE_GUARD_*` (the diagnostic namespace matching `S11c_c2_N6_diagnostic_sympy.py:837-849`), `G_cross`/denominator
   data in the payload. Confirm these can T7-join and nothing else regressed.
4. **Pressure-symbol census** — `FROZEN_RELATIONS` now includes the four pressure-symbol identity/assumption census
   (route-2 `:89-97`; reconcile_sympy `:60-73`). Confirm it is complete.

## Source-of-truth (read to judge; quote both sides with `file:line`)
- `research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md` §5c, §§1–2, §6, §7.
- `research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py`, `.../S11c_c2_N6_covariance_sympy.py`, `.../S11c_c2_N6_diagnostic_sympy.py` (the cleared instruments — the objects/names the WL engine reproduces blind).
- `research/pde_ledger_v3/_measurements/S11c_c2_N6_route2_spec_astra.md`.
- `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`, `S11c_a_SHARED_PHYSICS.md`, `S11c_c1_SHARED_PHYSICS.md`, `S11b_SHARED_PHYSICS.md`, `S11c_c1_wl_build_directive.md`.
- The two SymPy N6 build directives (`S11c_c2_N6_reconcile_directive.md`, `S11c_c2_N6_covariance_directive.md`) are
  REVIEWER-SIDE sources for you to check faithfulness against — ⛔ they are correctly NOT builder authorities.

## Required method
For each finding, quote the directive AND the contradicting/omitted source with `file:line`; a physics claim without a
source `file:line` is discarded. Report a finding only if it catches a way the WL engine, built to this directive, could
compute the wrong thing or let a wrong claim be made. Cover: faithful engine-neutral translation; anything still missing
for the carrier bridge + `R_cov`; blindness by absence (no computed outcome named anywhere in the directive; "re-derived"
not "imported"; no designed-to-agree); rule-17 freezes; controls able-to-fail and one-sided; the two round-2 regressions
not recurring; no leaked value / residual-zero exit / VERDICT; the three caveats open. ⛔ Do not propose making the WL
engine agree with SymPy.

## Output
Findings (directive `file:line`, source `file:line`, why it changes what is computed/claimed, minimal fix); flag any
round-2 fold that did NOT land or bred a new defect. Sound sections briefly. End with: CLEAR-TO-BUILD, or FOLD-REQUIRED
with the blocking items.
