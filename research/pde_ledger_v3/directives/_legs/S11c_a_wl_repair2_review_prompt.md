# Clarity/correctness review — S11c-a WL engine repair-round-2 directive

## Artifact
`/var/projects/toy_projects/research/pde_ledger_v3/directives/S11c_a_wl_engine_repair2_directive.md`
(if that path is wrong, use `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_wl_engine_repair2_directive.md`)

This is an orchestrator-written repair directive for the blind Wolfram engine
`research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` (baseline `a15bc69c`).
Physics authority: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (`2926c71c`), §2a (first-jet
map `∂W_bg = σ_W ∂w1`), §2b/§4 T-0, §3c/§4 T-f (dynamic window + projection). It fixes exactly two
pre-existing defects and must preserve the prior round's verified §5a/T-h fixes.

## Your job — is this directive CLEAR, CORRECT, and SAFE to hand to the builder?
1. **T-f fix (Defect 1).** The directive says the shape derivative of `∫ window·field dw` must keep every
   window factor inside `Inactive[Integrate]` so no `Inactive[Integrate][1, {w,…}]` (divergent unit
   integrand) appears. Is that a correct, unambiguous instruction? Read the engine's `projectionTermsSource`
   (656–685) and `shapeDerivative`/`applyProfileJets` (933–972). Is the prescribed fix (differentiate into
   the integrand / take the η-derivative before expanding across the held integral) actionable and correct?
   Is the confirm criterion (no emitted PROJECTION tag contains `Inactive[Integrate][1,…]`) right?
2. **T-0 grading fix (Defect 2).** The directive says to grade `∇Σ_E^0` via the §2a first-jet map rather
   than the fragile literal substitution `etaBg W0/LW -> sigmaW` (which no-ops for RHO4_CONSTANT). Read
   engine lines 705–720. Is the instruction clear and correct (RHO4_CONSTANT should be `σ_W` first-jet
   `{0,0,1}`, RHOBR_CONSTANT `{0,0,0}`), and does it correctly forbid hard-coding the answer?
3. **Preservation.** Does the directive correctly fence off the verified §5a/T-h work (route-2, one-sided
   source corruption, T-c′, sigmaEulerianSource, primary derivations) so the builder won't disturb them?
4. **Blindness + rule 5.** The engine is blind (imports nothing, re-derives). Does the directive leak any
   SIBLING-engine (SymPy) construction or any expected/target value? It must quote only the WL engine's own
   defect and the spec. Flag any leak.
5. **Anything unclear or wrong** that would send the builder in the wrong direction, or any missing
   instruction needed to fix these two defects cleanly.

## Output
A short list: is the directive clear and correct on each of the two fixes? Is the preservation fence
adequate? Any blindness/rule-5 leak? Any ambiguity to tighten before launch? End with: safe to hand to the
builder as-is, or list the exact edits needed.
