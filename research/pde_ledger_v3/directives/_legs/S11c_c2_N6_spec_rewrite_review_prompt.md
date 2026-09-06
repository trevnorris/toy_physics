# Spec rewrite review — corrected c2 §5c (N6 rep-invariance). STATIC only (⛔ run nothing, load no `.out`).

You are ONE of two independent legs reviewing a **spec correction**. A prior 2-leg adjudication established that
c2 §5c was wrong: it compared the two physical **anchorings** and demanded their residual vanish, but N6 (per
parent `S11c_a_SHARED_PHYSICS.md` §5a + `S11c_decisions.md` N4/N6, and sibling `S11c_c1` §5a) is the **Eulerian vs
material-coordinate construction of the SAME object, with the anchoring HELD FIXED** — a different axis from the
anchoring choice. The orchestrator has rewritten §5c. Verify the rewrite is correct and introduces no new defect;
a real problem is a finding, ⛔ do not rubber-stamp.

## Artifact
The current `directives/S11c_c2_SHARED_PHYSICS.md` §5c (≈ lines 303–340) and its diff vs HEAD:
`git -C /var/projects/toy_physics --no-pager diff -- research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md`.
Read the full current §5c, and for consistency: c2 §3c (≈203–222, the increment definition), §4 (≈264–276), §5d
(≈314–319), §5e (≈321–331); parents `S11c_a_SHARED_PHYSICS.md` §2c (≈236–244) + §5a (≈460–475), `S11c_decisions.md`
N4/N6 (≈68–100); sibling `S11c_c1_SHARED_PHYSICS.md` §5a.

## Checks (derive from the sources; quote lines)
1. **Correct + complete fix.** Does §5c now define N6's two routes as **Eulerian** (`I_E`) vs **material-coordinate
   mapped back to Eulerian** (`I_{M→E}`) of the SAME self-energy increment, with **α AND ρ held fixed**; residual
   `= I_E − I_{M→E}`; emitting the two operands + residual keyed by (α,ρ); the one-sided corruption of **one
   representation route** (not one anchoring); and does it ⛔ forbid using `Δρ = δρ_E + u·∇ρ⁰` as an anchoring
   bridge? Anything missing or mis-stated?
2. **Cross-anchoring reclassified correctly.** Is the `LAB_HELD − MATERIAL_ADVECTED` difference now a **separate,
   non-N6 contract** (`S11CC2_ANCHORING_L_MINUS_M`), explicitly with **no prescribed zero target**, on the footing
   of `DENSITY_LIVE_MINUS_FROZEN`? Any residual place in the spec that still treats the anchoring pair as N6 /
   still says the cross-anchoring residual "must vanish"?
3. **No leaked acceptance criterion (M2).** The physics claim ("same operator in two representations") is fine, but
   the builder must emit operands + residual with the **diff adjudicated on our side** — ⛔ no supplied numeric
   target, ⛔ no "iterate until residual==0" exit. Is the wording clean on this? (Sibling S11c-a §5a is the model.)
3b. **No manufactured derivation path / recipe.** Does it name the OBJECTS/routes without dictating a specific
   algebraic path that pre-answers the result?
4. **Buildability of route 2.** Is the material-coordinate route named precisely enough that a build can implement
   it (an S11c-b material / face-flattening construction + the S11c-a flattening for the SAME anchoring + the Δρ
   map-back), or does it silently require an object the build does not have? Flag underspecification — but ⛔ do not
   demand build mechanics that belong in the build directive, not the spec.
5. **Consistency + no new defect.** Does the rewrite contradict §3c/§4/§5d/§5e or the parents/sibling anywhere?
   Does §5e's uniform-limit smoke test (which also cites N6) still read correctly alongside the corrected §5c? Is
   §4's object list still consistent (it references §5 control outputs generically)?

## Output
For each of 1–5: finding + quoted lines. Separate CONFIRMED defects from unsettled. End: is the corrected §5c
**correct to commit**, or the exact remaining wording that must change.
