# S11c-b #89b WL directive — TRACTABILITY amendment decision-review (2 legs, folded once)

⚠ Context record for the orchestrator. No withheld production rank/termcount here.

Amendment reviewed: §3bis (Tractability) + §5 control 8, added to
`directives/S11c_b_89b_wl_operator_directive.md` after the naive live-operator regen failed to terminate (2h+).
Legs (rule 7: orchestrator-written ⇒ Codex + Grok), prompt
`_legs/S11c_b_89b_wl_tractability_decision_review.md`. Both COMPUTED (didn't argue). Evidence:
- Codex `~/.s11_build/S11c_b_89b_tractability_review_leg/` (div_order_probe, control8_decisiveness_probe,
  series_linearity_probe .wl+.stdout). VERDICT 1 finding (1 blocker), 133k tokens.
- Grok `~/.s11_build/S11c_b_89b_tractability_leg_grok/` (div_confluence, assoc_hang, if_rhs, series_linear
  .wl+.out). VERDICT 4 findings (1 blocker).

## What CLEARED (both legs computed it — do not reopen)
- **Div-activation confluence is sound.** top-down = innermost-first = by-level depth-walk, `PossibleZeroQ=True`
  on nested Div, live coefficient, variational wrapping, first/mixed derivative commutators, `List`/`Join`
  containers (Codex div_order_probe; Grok div_confluence: `zeroQ[topDown,inner]=zeroQ[inner,byLevel]=True`).
  `D[Inactive[Div][v,c],x] = Inactive[Div][D[v,x],c]` — activation commutes with `D`.
- **Fix 2 (per-summand `Series`) is bulletproof.** Zero residual vs global `Series` on polynomials, analytic
  rationals, shared denominators, Laurent-pole cancellation across summands, `ConditionalExpression`, `Hold`,
  `Log`, `O[]`, `Inactive[Div]` (both legs). `Series` linearity over `Plus` is exact.

## The converged BLOCKER — the fix mechanism + control 8 are BOTH under-nailed (verified by me, rule 13)
The precise root cause (Grok assoc_hang.out / if_rhs.out, verified): the 2h hang is **`ReplaceAll` of the
`If[SameQ[c,spatialCoordinates], Sum[D[…]], leftover Div]` RHS inside an `Association`** — the `If` condition
does not evaluate in the Association context (`hasUnevalIf=True`, `nDiv` stuck at 1), so the Div is never
removed and `ReplaceAll` re-fires unboundedly (assoc `nDiv 2→3→4→5→6→7`; the same expression as a `List`
terminates `2→0` in 2 passes). ⚠ Adding a Div-free `FreeQ` guard but **keeping the `If`** STILL hangs on an
Association — so a builder "adapting the mechanism, not pasting" re-creates the 2h hang.

And control 8 as written is NOT decisive (Codex control8_decisiveness_probe, verified): a naive innermost-first
`ReplaceRepeated` over an `Association` can drive leftover `Inactive[Div]` to 0 while leaving the generated
`Sum[D[…]]` **held/unevaluated** inside the container → 8(a)=0, 8(b) Div-count=0, 8(c)=0,
`ALL_CONTROL8_PASS=True`, YET the full profiled residual = `-sigmaW·wJetTwo·u - 2·sigmaW·wJetOne·∂u ≠ 0`: the
retained 1st/2nd jets are silently DROPPED when profiling substitutes the background before the held `D`
evaluates. A re-freeze that greens every control. 8(a) is also a nest-FREE row (EW_INTERNAL, `nNested=0`) so it
never exercises nested activation.

## Findings folded (all verified; all into the amendment)
1. **BLOCKER (Codex 1 / Grok 1) — name the ACTUAL trap + require derivative-normal.** §3bis must state: the
   hang is the `If`-RHS `ReplaceAll` on an `Association`; the fix must (i) use a `/; SameQ[c,spatialCoordinates]`
   pattern constraint with a fully-EVALUATING bare `Sum` RHS (no `If`), OR (ii) never run the rewrite on an
   `Association` — `Map` over list-valued rows / `Join` first (`elRowVector` ~L1205 already does this and
   terminates + equals innermost). Name the object = fully-activated Div `Σ_i ∂v_i/∂x_i` with derivatives
   **evaluated** (derivative-normal: no held `D`/`Sum`/`If`, no leftover `Inactive[Div]`).
2. **BLOCKER/SHOULD-FIX (Codex 1 / Grok 2) — control 8 decisive.** 8(a) → **per-row equality vs a top-down
   reference on the actual production container, including ≥1 nested-Div row (U_INTERNAL, `nNested=3`, top-down
   terminates 0.047 s)**: build the reference by top-down-activating each row (as List/Join) and reassembling;
   emit `bottom_up_full − reassembled_top_down` residual (=0). 8(b) → **derivative-normal postcondition**: emit
   the count of leftover `Inactive[Div]` AND of held/unevaluated `D`/`Sum`/`If` (both =0) — Div-count alone is
   not equality. Apply across every activation consumer (`elRowVector`).
3. **SHOULD-FIX (Grok 3) — do NOT make wall-clock a pass condition.** Acceptance = equality-to-naive (the
   per-row reference), not speed. The run is silence/RSS-bounded (§7). ⛔ Dropping retained-grade content to hit
   a clock is a FAILED build even if 8 residuals are 0. Partial-jet-drop teeth stay on §5.1 (Hessian-zero) +
   §5.5 (per-surface atom set) + new 8(a). Soften "a few minutes per case" from a target to "must complete
   (bounded, not the 2h+ non-termination)."
4. **NIT (Grok 4) — product-variant 8(c).** If the stronger grade-truncated-PRODUCT variant of Fix 2 is used,
   8(c) must include a product + re-truncation with a `thetaRule`-like rational coefficient (`Series[f]Series[g]
   ≠ Series[fg]` until re-truncated). Plain per-summand needs no change.

## Disposition
Two legs, one decision pass on the amendment (rule 7). The fix APPROACH is validated; the mechanism wording and
control 8 are sharpened. **Fold once and go** — ⛔ do not re-review the fold (that is iterating to green). NEXT:
re-leak-gate → Codex build (WL, `--sandbox danger-full-access`, xhigh) → 2 build legs (fresh Claude + Grok,
Mathematica SERIALIZED). The build legs verify the operator retains all jets + the controls bite.
