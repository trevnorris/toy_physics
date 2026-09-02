# Independent physics review — S11c-b #89b WL slab-operator un-freeze (Codex build)

You are one of two independent build-review legs (fresh Claude agent + Grok) on a **Codex-written** engine
change. Whatever writes does not review; derive independently and demand computation. Report a finding only if
it catches a way the physics could be wrong.

## Artifact — the ENGINE SOURCE (there is no fresh full `.out` to read)
`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` (Codex build, +~750 lines vs
the pre-#89b engine). Build directives (what it was asked to do): `directives/S11c_b_89b_wl_operator_directive.md`
(the operator un-freeze + tractability) and `directives/S11c_b_89b_wl_defer_heavy_controls.md` (the deferral gate).

⚠⚠ **THE FULL IN-BAND `.out` IS DEFERRED TO A ≥64 GB BOX — do NOT expect to read the operator/kernel emits from
a committed `.out`.** Measured 2026-09-02: the correct un-frozen operator must hold its full jet tower un-reduced
(~16 GB/case) until the final reduction (reducing earlier is the re-freeze this build removes), so even with the
two heavy controls gated off, the operator loop alone peaks ~16 GB/case and the 4-case run does not fit 30 GB.
The whole regen (operator/kernel/μ_θ emits AND the two controls) is bigger-box work; see
`research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md`. ⇒ **Your review runs on this 30 GB box, so you verify by (A)
SOURCE FIDELITY of the un-freeze + tractability + gate, and (B) your own REDUCED-SCALE ablations that fit
`timeout 600` — build the operator at a TRUNCATED jet depth / low-order profile so a row is small — plus (C)
cross-check against the OUT-OF-BAND records.** Full-depth, all-cases in-band content is explicitly OUT OF SCOPE
here (deferred); flag it as such, do not attempt a full run.

Out-of-band records you may cross-check (⛔ never as your only evidence — your own reduced-scale computation is
the evidence): `directives/_measurements/S11c_b_89b_wl_operator_decision_review.md`,
`…/S11c_b_89b_wl_tractability_decision_review.md`, and (PY side) `_measurements/S11c_b_88_blast_radius_result.md`.

The two heavy controls (`CONTROL_TRACTABLE_ACTIVATION_EQUIVALENCE`, `CONTROL_SURFACE_BACKGROUND_JET_DIFFERENCE_ATOMS`)
are gated behind `S11CB_SKIP_HEAVY_CONTROLS`; verifying that gate is a clean, physics-neutral deferral is item 7.

## What the build did (the claim to test)
The §3b variable-coefficient slab operator was **un-frozen**: the background coefficients `widthBase`/
`modulusBase` are now kept LIVE (functions) through the whole Euler–Lagrange variation + in-plane divergence
handling + §3c weak-block kernel extraction + face EL, and the jet-retaining reduction is applied LAST — so
the emitted operator retains the full background jet tower, including the retained-grade curvature (2nd/3rd
spatial jets of the profile, at single-σ_W grade), which the previous FROZEN construction zeroed. The
computation was made tractable (the naive form did not terminate) by (i) activating divergences so the
generated derivatives actually evaluate — no `If`-on-`Association` re-firing, no held `D` — and (ii) exploiting
`Series` linearity per top-level summand. Both are claimed **algebraically equivalent** (no jet term dropped).

## What you must independently establish (each: derive it yourself, don't trust the engine's self-report)
1. **The operator genuinely retains the retained-grade curvature (no re-freeze) — at REDUCED scale.** The
   corrected operator (`SLAB_OPERATOR`, and the `MU_THETA_OPERATOR`, `COUPLING_KERNEL`) must depend on the 2nd
   (and, in the kernel, 3rd) background-jet atoms (`w1Jet*`/`m1Jet*`/`widthProfileJet[…]`/`modulusProfileJet[…]`
   at order ≥2) that the frozen operator lacks. ⚠ Full-depth is ~16 GB/case (deferred) — build a corrected row
   at a TRUNCATED jet depth / low-order profile that fits `timeout 600`, and the matching frozen row (the
   `applyProfile`-before-EL route, = `CONTROL_OPERATOR_REFREEZE_REGRESSION_OPERAND`), and show the difference set
   carries the order-≥2 jets at reduced depth. ⛔ A corrected row whose retained-grade atoms are ABSENT even at
   reduced depth is a silent re-freeze — the whole point of the build. Save script + literal stdout.
2. **The tractability speed-ups dropped NO term (the decisive check).** ⚠ `CONTROL_TRACTABLE_ACTIVATION_EQUIVALENCE`
   is DEFERRED in this `.out` (marker only), so you must ESTABLISH this yourself on a minimal harness — do not
   read it. Extract the two activation routes (`activateSpatialDivergences` and
   `activateOperatorPerRowTopDownReference`) and the per-summand `Series` path; on a representative row at a
   TRUNCATED jet depth that fits `timeout 600` (full depth ~16 GB, deferred) INCLUDING a NESTED-Div row (e.g.
   `U_INTERNAL`), compute (i) `tractable − top_down_reference` (must be 0),
   (ii) the derivative-normal counts (leftover `Inactive[Div]` AND held `D`/`Sum`/`If` all 0), (iii) the
   per-summand-vs-global `Series` residual (0). ⛔ Then ABLATE to prove the check is decisive: re-introduce a
   term-dropping activation (e.g. drop a generated jet, or a naive `ReplaceAll` that leaves a held `Sum[D]`) and
   confirm the residual goes NONZERO. If a residual stays 0 under that ablation, the check is decorative. Save
   script + literal stdout to a named absolute path. You MAY cross-check your numbers against the out-of-band
   record `directives/_measurements/S11c_b_89b_wl_tractability_decision_review.md`, but your own run is the
   evidence.
3. **FORM ablation (MANDATORY — the only thing that catches the worst defect).** On a `/tmp` copy: (a)
   Hessian-zero — set the order-2/3 background-jet atoms → 0 on the corrected operator; it must MOVE back toward
   the frozen construction. (b) Re-freeze — `applyProfile`-before-EL must reproduce the committed frozen
   operator. (c) Tower-depth — truncate the retaining depth one order below the generated depth (an emitted
   object MOVES) / extend one above (does NOT move). (d) Grade-support — show σ_W² is absent (including across
   an `Inactive[Div]`) while ησ_W survives. Report the literal before/after for each. A COEFFICIENT rescale is
   not a form ablation.
4. **Cross-engine CONTENT, not a rank — full check DEFERRED, reduced check here.** The corrected operator
   differs from the frozen one by non-absorbable Hessian structure (per `_measurements/S11c_b_88_blast_radius_result.md`,
   the PY side). ⚠ The full-depth all-cases cross-engine content check needs the bigger box (deferred — flag it).
   On this box, verify at reduced depth that the WL corrected operator carries that STRUCTURE — the kind of
   jet-atom content (2nd/3rd-order background jets entering non-absorbably) — not merely that some count changed.
   ⛔ Do not compare to a hard-coded rank; the corrected rank/termcounts are withheld — a matching integer is not
   evidence. If you cannot separate real structure from a reduced-depth artifact, say so and defer rather than
   guess.
5. **§3c kernel + §3d admissibility.** The `COUPLING_KERNEL` must be a weak block of the corrected §3b operator
   (its `Inactive[Div]` split preserved), not a parallel route; its `TERM_ORIGINS_SUM_RESIDUAL` should be 0.
   The `ADMISSIBILITY_OPERATOR_OPERAND` (§3d, a separate already-live route) should be UNCHANGED by this build
   — confirm.
6. **§5.E dimension fix.** The new dimension residual must use an INDEPENDENT primitive-atom route (dimensions
   of `uOne`, `D[uOne,x]`, `Derivative[widthBase]/WZero`, …), ⛔ not a second copy of the `FACTOR_NAMES` sum.
   Confirm it is no longer a structural zero (mutate a primitive atom's dimension on a copy → residual moves).
7. **Deferral integrity — the two skipped controls are a physics-neutral gate, not a silent drop.** On the
   engine source confirm: EXACTLY the two named controls are wrapped in `If[!skipHeavyControls, …]` around ONLY
   their heavy build (the `tractableOperator`/`referenceOperator` build inside the main `Do` ~L2238; the
   surface-Hessian-witness build inside the `Do` ~L2515), and NOTHING ELSE is gated — the core operator / kernel
   / μ_θ / origins / activation-postconditions assignments in those same loops, the
   `CONTROL_OPERATOR_REFREEZE_REGRESSION_OPERAND` streaming, and every other control all stay unguarded.
   Structurally confirm both branches of each gate: when set, the payload is the `DEFERRED_HEAVY_CONTROL` marker
   and the downstream `emitShared` still fires (tag present, ⛔ never silently absent); when unset, the payload
   is `<||>` and the loop's `If[!skipHeavyControls, …AssociateTo…]` fills it with the real per-case build.
   Decisive check WITHOUT a full run (a full case is ~16 GB/50 min — ⛔ do NOT run it): extract a MINIMAL harness
   that exercises the gate wrapper itself — feed the two `If[!skipHeavyControls, …]` bodies a tiny stand-in
   payload and toggle `skipHeavyControls` True/False — and confirm (marker vs real association) toggles correctly
   and that the surrounding must-keep assignments (`mainModels`/`mainKernels`/origins, the refreeze streaming)
   are OUTSIDE the gated body in the source. ⛔ A gate that also suppresses a non-deferred payload, a marker
   absent when the var is set, or the real-build path unreachable when it is unset, is a BLOCKER.

## Method + hygiene
- Derive independently; where you dispute equivalence or a control, **write your own script, run it, and save
  the script + its literal stdout to a named absolute path you report.** A prose re-derivation is discarded.
- Report any `assert` that precedes the value it guards; report any emitted residual that is zero for ANY input
  (tautology).
- ⛔ This is a Mathematica artifact: wrap EVERY kernel run in `timeout 600`; run ONE kernel at a time (2-seat
  licence); copy the engine to `/tmp` and ablate the COPY, never the working tree; after each run kill any
  leftover `WolframKernel` by its EXACT pid (never `pkill -f`). ⛔ Do not run the FULL engine — extract the
  functions you need into a minimal harness.
- Report a numbered list: `severity — file:line — finding — the physics it gets wrong — the computation that
  shows it`. Severity ∈ {BLOCKER, SHOULD-FIX, NIT}. Then `VERDICT: N findings (B blockers)` or `VERDICT: CLEAR`.
- A leg that finds nothing is weak evidence: if you clear it, name the two or three things you actively tried to
  break (with the script you ran) so the orchestrator can weigh the clearance.
