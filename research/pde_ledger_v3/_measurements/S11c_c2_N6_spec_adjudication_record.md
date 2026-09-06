# Gate record — S11c-c2 N6 rep-invariance is MIS-SPECIFIED (c2 §5c defect) — 2026-09-06

**Trigger:** re-grounding STEP A's E/N6 the corrected way (orchestrator never authors the CAS instrument; CLAUDE.md
`6f8dbd34`). The biased `verify_EG.py` (retired) computed only `residual.subs(σ_W→0)==0` and I over-cleared E.

## Method (the corrected process, exercised end-to-end)
1. **I framed the question**, then **ran it by Codex** (question-vet, `directives/_legs/S11c_c2_N6_question_vet_prompt.md`
   → `scratchpad/codex_N6_question_vet.txt`). Codex-sol found my `(0,1)`-vs-`(1,1)` framing was itself a **proxy**
   (localization, not discrimination), settled the scope subtlety (`(1,1)=O(εησ_W)` is RETAINED, not deferrable —
   §2c forbids the `ησ_W~η²` rescaling), and surfaced a **deeper possible spec conflict**: the residual compares the
   two *anchorings*, but S11c-a says they are distinct physical setups, not Eulerian/material representations.
2. **Spec adjudication** (orchestrator-written spec → **Codex-sol + Grok**, neutral prompt, no leaked conclusion):
   `directives/_legs/S11c_c2_N6_spec_adjudication_prompt.md` → `scratchpad/{codex,grok}_N6_spec_adjudication.txt`.
   **Both legs independently converged: c2 §5c is wrong.**

## Finding (2-leg converged; citations verified by orchestrator reading — rule 13)
**c2 §5c ("the two anchorings ARE the representation-invariance pair … residual must vanish") is a substantive
definition error.** The correct N6 is already specified by the parents and the correct sibling:
- **S11c-a §2c** (`S11c_a_SHARED_PHYSICS.md:239`): LAB_HELD / MATERIAL_ADVECTED are "**distinct physical
  anchorings, not alternate Eulerian/material representations of one branch**."
- **S11c-a §5a** (`:461`, "Representation-invariance routes (N4/N6)"): "For **each physical anchoring**, compute by
  two independent routes — (1) direct Eulerian, (2) material-coordinate mapped back — … **The anchoring branch is
  held fixed across the two routes.**" Emits `S11CA_REP_INVARIANCE_{EULERIAN,MATERIAL}_OPERAND` + residual.
- **parent N4** (`S11c_decisions.md:72`): "Eulerian and material variables … **must agree after that field
  redefinition** (`Δρ = δρ_E + u·∇ρ⁰`)"; **N6** (`:94`): direct level-set vs material-coordinate, transform &
  compare. Δρ relates two *descriptions of one perturbation*, ⛔ not two backgrounds.
- **Sibling c1 §5a** implements N6 correctly (Eulerian-direct vs Hanzawa, anchoring fixed) and **explicitly forbids
  using Δρ as the anchoring map / treating MATERIAL_ADVECTED as the second N6 route** — so c2 §5c is an isolated
  defect, not the house pattern.

**Three axes, conflated by §5c:** Representation (Eulerian vs material — N6, agrees after Δρ) · Anchoring (L vs M —
distinct physics, N4) · Density (live vs frozen — §3d.1). §5c collapsed Representation onto Anchoring.

## Adjudication (CONFIRMED)
- The current `S11CC2_REP_INVARIANCE_RESIDUAL` = a **cross-anchoring** difference `I^L − I^M` — the **wrong object**.
  A nonzero value is **expected physics** (distinct anchorings differ), ⛔ NOT an N6 defect. Using Δρ to bridge L↔M
  is a category error (a same-state rewrite applied to two different states).
- **c2's real N6 (Eulerian vs material-coordinate construction of the increment, at a FIXED anchoring, mapped back
  by Δρ) was never computed — it is UNTESTED.**
- Both prior positions were mis-framed: my over-clear (wrong proxy AND wrong object) and Grok's original E flag
  ("σ-channel doesn't vanish" — but that residual isn't supposed to vanish). ⇒ the biased `verify_EG` E conclusion
  is **retired**; the "possible `representation_pullback` build defect" framing is superseded (the pullback
  *incompleteness* — trial-field-only, missing the test-covector `T*·O·T` transform — is a separate, real issue but
  is downstream of the correct-object question).

## Fix owed (routed through review; NOT folded here)
1. **Spec correction — c2 §5c** (review-until-clear, Codex-sol + Grok): N6's two routes for the self-energy
   increment are **Eulerian graph/level-set** vs **material-coordinate flattening mapped back by Δρ**, with **α and
   ρ held fixed**; emit `S11CC2_REP_INVARIANCE_{EULERIAN,MATERIAL}_OPERAND[α,ρ]` + residual (must vanish);
   one-sided corruption of **one representation route**. Reclassify the current cross-anchoring object as a
   **separate non-N6 anchoring-difference contract** (`S11CC2_ANCHORING_L_MINUS_M`, no zero target), like
   `DENSITY_LIVE_MINUS_FROZEN`.
2. **Build correction** (after the spec clears): add the material-coordinate fold route so the real N6 is computed;
   the current build lacks it. (Codex-written; G1 legs; then adjudicate.)

## Status of STEP A re-grounding
- **E/N6:** re-grounded to the correct object — result is a **spec defect + untested control** (above). The
  increment VALUES themselves (per-anchoring) are unaffected; only the N6 *control* was mis-defined.
- **F, G:** still owed — they rest on the biased `verify_F.py`/`verify_EG.py` and must be re-grounded the corrected
  way (Codex-written diagnostic, G1 legs) after E's path is set.
