# S11c-b #89b WL operator un-freeze — build-review measurements (2 legs + orchestrator verification)

Codex-written engine change (operator un-freeze + tractability + deferral gate). Two independent build-review
legs on the then-working-tree engine `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
(+~750 lines vs #89a HEAD; NOW committed at `a1be8d8f` after the repair recorded below — line numbers in the
findings below are PRE-REPAIR, and the shift to HEAD VARIES by location: ~+1 near the Finding-1 emit (now
`:1405–1406`) but ~+25–28 for Findings 2–3, which sit after the repair's ~+28 inserted lines (dimension walker now
`:2986–2988`; independence marker now `:2489`, package `:2469–2496`). The authoritative post-repair line numbers
are in the "REPAIR + re-review outcome" section below). Both worked at REDUCED SCALE (full-depth operator ~16 GB/case,
deferred to a bigger box). Prompt: `directives/_legs/S11c_b_89b_wl_operator_build_review.md`.

## Leg verdicts

- **Leg 1 — Grok** (`~/.s11_build/S11c_b_89b_buildlegs/grok_leg.txt`): VERDICT 2 findings (1 blocker).
  Scripts+stdout under `/tmp/s11cb_89b_grok_leg/` (H1/H2/H3, gate).
- **Leg 2 — fresh Agent (opus)** (transcript task `a11e1025298eb00d4`): VERDICT 2 findings (0 blockers).
  Scripts+stdout under the session scratchpad `…/scratchpad/test_*.wl|out`.

The legs DISAGREED on the primary item (the emitted operator). That disagreement is the measurement; the
orchestrator resolved it by its own computation (below). Rule 13.

## Finding 1 — BLOCKER — emitted SLAB_OPERATOR / SLAB_OPERATOR_TERM_ORIGINS re-freeze (U-momentum rows)

`…audit.wl:1404` `operator = finalBackgroundReduction[operatorLive, source]` and `:1405`
`origins = finalBackgroundReduction[originsLive, source]` reduce the operator whose OUTER constrained
`Inactive[Div]` are STILL HELD. The correctly-activated `activatedOperator = activateSpatialDivergences[operatorLive]`
(`:1376`) and `activatedOrigins` (`:1378`) are computed and then DISCARDED. The sibling emits use the activated
route (`muTheta` `:1403`, `facesData` `:1406`, `divergenceSource` `:1407`). The L1372–1375 comment claims the
held-Div emit is needed "so the same live operator supplies the §3c weak split" — FALSE: the weak split consumes
`KERNEL_SOURCE_OPERATOR -> operatorLive` (`:1412`, un-reduced), not the emitted `operator`. The activation
postcondition (`:1388-1389`) is even computed on `activatedOperator`, mismatching the un-activated emit.

Effect: reducing before activating the outer Div freezes `widthBase`→constants INSIDE the held Div, so
activating afterward cannot regenerate the mixed/higher jets — the rule-17 re-freeze — on the emitted operator
that the cross-engine row comparison reads.

**Grok evidence** (`/tmp/s11cb_89b_grok_leg/H2.stdout`):
`H2_EMITTED_U_COEFF_widthProfileJet[1,1,0]: {0, 0, 0}` vs
`H2_FULLY2_U_COEFF_widthProfileJet[1,1,0]: {0, -((sigmaW*Derivative[1,0,0,0][thetaWave])/LWidth), 0}`;
`H2_FULLY2_MINUS_EMITTED_THEN_ACTIVATE_LISTZEROQ: False`.

**Orchestrator verification** (own script, rule 13):
`/tmp/claude-1000/-…/scratchpad/orch_refreeze_verify.wl` →
- `SANITY reduce(uLive) matches emitted OPERATOR U rows? True` (I tested the REAL emit).
- Coefficient of `widthProfileJet[1,1,0]` on the U row: CORRECT (activate→reduce) has 5 field-derivative terms;
  EMITTED-then-activated (reduce→activate) has only 3 → they differ.
- `correct − emittedActivated` NOT identically zero; missing higher jets up to order 3
  (`{widthProfileJet[2,0,0],[2,1,0],[3,0,0],[1,2,0],…}`). ⇒ activate/reduce DO NOT commute on the U row.

Leg 2 (Agent) CLEARED this because its test only checked jet-atom PRESENCE (`backgroundJetAtomsIn`) — the atoms
sit inside the held `Inactive[Div]` (`test_refreeze2.out` raw InputForm) — never the ACTIVATED coefficient. An
insufficient test; superseded by the coefficient/commutativity test above.

**Fix:** emit `finalBackgroundReduction[activatedOperator, source]` / `finalBackgroundReduction[activatedOrigins, source]`;
correct the L1372–1375 comment. (Kernel path already correct — both legs agree — so the fix is localized.)

## Finding 2 — SHOULD-FIX — §5.E "independent primitive-atom" dimension walker is non-functional (BOTH legs)

`…audit.wl:2958-2969`. `primitiveExpressionDimension[Times[factors__]]` / `[Plus[terms__]]` do not register as
usable DownValues (Times/Plus are Flat+Orderless+OneIdentity: the Times rule is dropped; the Plus rule is stored
below the generic `expression_` catch-all). Every composite invariant falls through to
`UnassignedPrimitiveDimension[HoldForm[…]]`; the residual `derived − Unassigned` is non-numeric and does NOT
depend on the primitive atom dimensions. Confirmed independently: Grok `H3.stdout`
(`H2_DIM_RESIDUAL_MOVED_AFTER_UONE_MUTATION: False`); Agent `test_dim_iso.out` + fix confirmation
`test_dim_fix.out` (`expr_Times`/`expr_Plus` head-patterns → `{0,0,0}` for `Times[1/2,theta]`, and residual MOVES
under atom-dim mutation). **Fix:** use `expression_Times` / `expression_Plus` head-patterns.

## Finding 3 — SHOULD-FIX — CONTROL_INDEPENDENCE_RESIDUAL tautology for MATERIAL_ADVECTED (Agent; diff-confirmed)

`…audit.wl:2464-2468`. #89b changed `corrupted = evaluatedModel[…,True]` (always corrupted) to
`corrupted = If[branch === First[branches], evaluatedModel[…,True], base]` — so for non-first branches
(the two MATERIAL_ADVECTED cases) the "corrupted" operand IS `base`, and the independence residual is
`base − base = 0` regardless of physics (a one-sided-corruption control turned vacuous). Diff-confirmed:
`git diff HEAD` line 700 (`-  corrupted = evaluatedModel…`) → 702 (`+  corrupted = If[branch === First[branches]…`).
Not authorized by either directive. Likely a memory shortcut (the corrupted model is another ~16 GB build).
**Fix:** don't emit a silent `base−base`; where the corruption is skipped for tractability, emit an explicit
`VALIDATED_ON_REPRESENTATIVE_BRANCH` marker (like the §3 streaming fallback), or restore the per-branch
corruption if it fits.

## Consolidated verdict

BLOCKER (finding 1) + 2 SHOULD-FIX (findings 2, 3). One Codex repair round handles all three (each fix is small
and localized; the correct object for finding 1 is already computed). Re-review the repair diff (2 legs) before
commit. Rule 9: no commit before the re-review legs report.

## REPAIR + re-review outcome (2026-09-02) — ALL THREE FIXED, cleared by 2 legs

Repair directive: `directives/S11c_b_89b_wl_operator_repair.md`. It got **2 decision legs** first (rule 7) —
Grok (`grok_repair_decision.txt`, 7 issues) + Codex (`codex_repair_decision.log`, 4 issues) = **11 directive-level
gaps** folded into v3 (the fix ripples to the tower-depth operands, the frozen Hessian-witness across ALL its
slots, the uniform-limit `SLAB_OPERATOR` **and** `TRANSVERSE_DISPERSION`, plus a rule-5 leak in the first
acceptance that could have passed the frozen object). Codex built the repair (+28 net lines; scope-verified =
only the 3 fixes + ripples, no streaming/scope creep).

- **Fix A** (`:1405–1406`): emit `finalBackgroundReduction[activatedOperator/activatedOrigins]` (activate-then-reduce);
  tower operands, frozen-witness comparison operands (all slots), and uniform-S11b comparison operand (SLAB +
  TRANSVERSE_DISPERSION) brought to derivative-normal; `KERNEL_SOURCE_*` and `LIVE_DIVERGENCE_FORM_OPERAND` left
  un-reduced-live.
- **Fix B** (`:2986–2988`): `primitiveExpressionDimension` now uses `expression_Times`/`expression_Plus`.
- **Fix C** (`:2469–2496`): `LAB_HELD` fresh corrupted build per slot; MATERIAL_ADVECTED whole package →
  `VALIDATED_ON_REPRESENTATIVE_BRANCH` marker.

**Re-review (2 legs, reduced-scale, rebuilt from `operatorLive` not the emit — Grok `grok_rereview.txt`, fresh
Agent opus):** BOTH **VERDICT CLEAR**, agreeing on every check:
- A: `EMIT − from-operatorLive reference = {0,0,0}`; one-sided swap to reduce-then-activate MOVES it; U-row
  `widthProfileJet[1,1,0]` coeff 5 terms (emit/ref) vs 3 (swap); higher-jet atoms 16 vs 6 (order-3 restored).
- Ripples: frozen witness both sides 0 leftover Div (like-with-like) on every slot; uniform SLAB + TRANSVERSE
  both derivative-normal; S11b reconstruction byte-identical to HEAD (comparison-only surface added); tower depth
  measured on the activated operator (truncate moves / extend does not).
- B: head-patterns register + fire; 30 invariants resolve; mutating `uOne`'s dimension moves EXACTLY the 22
  U-invariants, leaves the 8 non-U ones (selective, not blanket).
- C: whole-package marker on MATERIAL_ADVECTED; `LAB_HELD` corruption still live (a real defect drives the
  operator off zero on a θ-bearing representative). Agent noted RHOBR_CONSTANT corruption is a pre-existing no-op
  (σ_bg width-independent), NOT a repair defect.
- No regression: `KERNEL_SOURCE_*` un-reduced live (activatable to depth 3), gate toggles marker↔real, inner
  live-EL un-freeze + tractable activation unchanged.

Orchestrator (rule 13): independently reproduced the DEFECT before the repair
(`…/scratchpad/orch_refreeze_verify.wl`) and judged the repair self-check
(`/tmp/s11cb_89b_repair_harness/defect_{a,b,c}.stdout`); both agree with the legs. ⇒ committed.
