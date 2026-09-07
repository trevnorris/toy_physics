# Adjudication review — does the reconcile tally SUPPORT "R_N6 vanishes modulo the defining relations"?

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_adjudication.md` — the
orchestrator's adjudication of the per-engine N6 reconcile. Working dir `/var/projects/toy_physics`. This is a
**physics-bearing prose review**: does the disposition tally SUPPORT the conclusion, or is it an **over-clear**?
⚠ The orchestrator over-cleared N6 once before (the "L-CAS" incident: claimed invariance in a σ_W→0 projection, not
the retained order) — scrutinize specifically for that failure mode.

## What to check (form your OWN view from the sources FIRST, then read the adjudication)
Sources: the instrument `scripts/S11c_c2_N6_reconcile_sympy.py`; the build-clearance
`_measurements/S11c_c2_N6_reconcile_build_clearance.md` (both build legs CLEAR + the carrier-knife caveat); §5c of
`directives/S11c_c2_SHARED_PHYSICS.md`; the cleared route-2 spec `_measurements/S11c_c2_N6_route2_spec_astra.md`
(esp. §2 pullback, §4 what survives pressure projection, §5 combined source). The disposition **tally is in the
adjudication record** (the 4-case table); to re-verify it, run the committed parser
`_measurements/S11c_c2_N6_reconcile_disposition_tally.py` on the `.out` at
`/tmp/S11c_c2_N6_reconcile_sympy.<ANCH>.<DENS>.out` (regenerate a case with
`python3 scripts/S11c_c2_N6_reconcile_sympy.py --anchoring LAB_HELD --density RHOBR_CONSTANT` if `/tmp` is unavailable —
per-case PIT is ~2 min).

Scrutinize each adjudication claim:
1. **Carrier bridge = geometry, reconciles.** Is `CARRIER_BRIDGE_RESIDUAL = C_E − C_M` really the geometric/mechanical
   carrier, and is its `0/80` (all 4 cases) a GENUINE reconciliation `C_E = C_M` — given the build legs found the
   identity-map knife too weak but the material-NORMAL FORM knife DOES move it (so it is a live control, not a
   false-negative dead zero)? Or could `0/80` be a false negative from an inadequate PIT / a degenerate construction?
2. **R_N6 = SOURCE_CHANNEL = "sanctioned N4 channel" ⇒ vanishes modulo the defining relations.** This is the crux —
   is it too glib? `R_N6 = B(C_M, ΔS)` with `ΔS = es − ms`. Is `ΔS` ENTIRELY the sanctioned material↔Eulerian
   defining-relation content (the `material_pullback` μ + covector-mapped velocity difference, route-2 §2/§5), or
   could `ΔS` — or `B(C_M, ΔS)` — harbor content NOT derivable from the prescribed maps (which would be a real
   remainder, not "sanctioned")? Does "R_N6 lies in the ideal of the defining relations" actually follow, and is that
   the correct operationalization of N6 (per the question-vet), or a convenient tautology (is the conclusion forced
   for ANY residual that localizes to the source channel, making the test vacuous)?
3. **The a_ρ + h_α cross-check.** Is `R_N6 = 0` in exactly `MATERIAL_ADVECTED.RHO4_CONSTANT` (both pullback
   ingredients absent: `a_ρ=0` since `g_i=0`, `h_α=0` at MATERIAL_ADVECTED) genuinely evidence that the residual is
   DRIVEN by the sanctioned content and nothing else — or is it also consistent with a real defect that happens to
   vanish when the pullback is trivial? Verify `h_α=0` at MATERIAL_ADVECTED and `h_α≠0` at LAB_HELD (route-2 :45).
4. **Over-clear check.** Is the conclusion "per-engine N6 SATISFIED / representation invariance holds" claimed at the
   RETAINED order `(η^{≤1},σ_W^{≤1})` (where R_N6 is nonzero), ⛔ not in a projection/limit? Are the PIT caveats
   honest (all-zero = "no nonzero found" at conditional δ, ⛔ not absolute certified zero)? Is anything stated as
   settled that is actually only per-engine / still owed (the blind Wolfram N6)?

## Physics filter
Report a finding only if it catches a way the ADJUDICATION could be wrong, over-claimed, vacuous, or an over-clear —
⛔ not a re-review of the instrument (the build legs cleared that).

## Output
Findings each with the specific claim + why it is wrong/over-claimed + the corrected statement. End with
**ADJUDICATION SOUND** (the tally supports "R_N6 vanishes modulo the defining relations", per-engine, retained order)
or the exact correction. Brief, evidence-first.
