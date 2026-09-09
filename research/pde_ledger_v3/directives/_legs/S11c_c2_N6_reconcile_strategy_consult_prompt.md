# Strategy consult (advisory — NOT a blind review): S11c-c2 N6 cross-engine reconcile disposition

You are advising the orchestrator on a STRATEGIC decision. This is not a blind review leg — you have the full
picture on purpose. Read the state, then give a reasoned recommendation and, most importantly, tell us if we
missed something. Be adversarial about the approach, not just the details.

## The situation (read these to ground yourself)
- The reconcile QUESTION (vetted): `research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md`
- The comparator RUN result: `research/pde_ledger_v3/_measurements/S11c_c2_N6_comparator_run_data.md` (+ `…run_tally.txt`)
- The per-engine N6 result that ALREADY STANDS (committed): `research/pde_ledger_v3/_measurements/S11c_c2_N6_RESOLVED.md`
- The collapse-instrument directive saga (round 1→2→3, four rounds of map-4 findings):
  `research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_collapse_directive_gate.md`
- The v3 directive itself: `research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md`
- The engines: `research/pde_ledger_v3/scripts/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.py`,
  `research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl`; comparator
  `research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py`.

## The state in brief
The N6 cross-engine comparator ran clean. Two blind engines (SymPy imports the slab; WL re-derives blind) AGREE
on matched objects that are 0 within each engine (`R_cov`, the carrier bridge `C_E−C_M`) — but those are
`0 − 0 = 0`, TRIVIAL, and say nothing about whether the underlying OPERANDS match. The surfaced-nonzero
operands are: carrier (40 leaves), constitutive source ACTUAL/PREDICTED (76 each), field-map Φ (18).

The proposed "collapse instrument" applies a frozen bridge dictionary of justified representational identities to
both engines' operands and prints whether the residual vanishes (three-valued; decides nothing). The load-bearing
entry is **map 4**: the two engines use different thickness conventions (WL uses `WBg` directly; SymPy uses a
rescaled local thickness `E = W₀·e_W/WBg`), so their energy-basis coefficients differ by `R_W = W₀/WBg` factors.
Map 4 is a 19-row table translating WL coefficients → SymPy with those scales. It has now drawn THREE distinct
review findings across four rounds: (1) not frozen, (2) needs explicit scales, (3) the `R_W` factors are correct
UNGRADED but are applied to already-GRADED leaves (both engines expand `WBg→W₀(1+η·profile)` before grading), so
applying `R_W=W₀/WBg` post-grading re-inserts a spurious `WBg` — a stage artifact.

The per-engine result (`R_cov=0` ⇒ operator covariance holds, "Reading B", user-adopted) is COMMITTED and
independent of this; the collapse test only determines how STRONG the cross-engine corroboration of channel (b)
(the constitutive source/Φ) is.

## The three paths on the table
- **A. Continue full:** fix map-4 stage (Codex v4) → re-review → astra build → run → adjudicate + disposition
  legs → step record. Definitive channel-(b) answer, but ~8-12 more agent launches / several hours, and the
  collapse may still land partly UNDECIDED (e.g. an η¹ profile-jet term WL carries and SymPy cancelled).
- **B. Surface as UNDECIDED:** stop; declare the surfaced operands representational-difference-unadjudicated (the
  c1 precedent), write the step record. N6 cross-engine then rests on the trivial matched-zeros — weak.
- **C. R_W rows only:** test only the load-bearing R_W-normalization rows + contraction-structure match, surface
  the rest. Still needs the map-4 stage fix + a build.

## What we need from you
1. **Which path do you recommend, and why?** Weigh correctness value vs effort honestly.
2. **Did we miss something?** In particular:
   - Is the graded coefficient-table bridge the RIGHT instrument for channel (b), or does map-4's recurring
     `R_W`/staging difficulty signal a mis-decomposition — is there a cleaner construction (compare at a
     different stage; a numeric cross-engine PIT channel; absorb the thickness convention upstream; compare a
     convention-invariant object) that sidesteps the thickness-convention problem entirely?
   - Is the per-engine covariance result already SUFFICIENT for what c2 needs downstream (S11c-d binds the
     closed operators), making the cross-engine operand collapse a nice-to-have rather than load-bearing?
   - Is there a correctness TRAP in any path (e.g. does surfacing-as-UNDECIDED hide a real cross-engine
     disagreement that only the collapse test would catch — the two-errors-cancel risk)?
3. **If A or C:** is the map-4 stage fix (per-grade images, or consistent `WBg→W₀` expansion on both operands
   before the coefficientwise compare) sound, or is there a deeper reason it will keep breaking?

Ground claims in the files (cite where you can). Give a clear bottom-line recommendation. Being contrarian is
welcome — the point of this consult is to catch a wrong turn before we spend the hours.
