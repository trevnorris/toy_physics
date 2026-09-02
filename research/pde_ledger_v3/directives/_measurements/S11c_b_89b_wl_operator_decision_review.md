# S11c-b #89b WL operator directive — decision-review outcome (2 legs, folded once)

⚠ ANSWER-BEARING context (mechanism). For the orchestrator + directive record only — ⛔ NOT handed to the
builder. No withheld production rank/termcount appears here.

Directive reviewed: `directives/S11c_b_89b_wl_operator_directive.md` (v1). Legs (rule 7: orchestrator-written
⇒ Codex + Grok), prompt `directives/_legs/S11c_b_89b_wl_operator_decision_review.md` (byte-identical to both).
Both computed the jet-depth/freeze question rather than arguing it, and saved script + stdout.

- **Leg A — Codex (xhigh).** VERDICT 4 findings (2 blockers). 146,686 tokens. Evidence
  `~/.s11_build/S11c_b_89b_decision_codex/` (jet_depth_check.py+stdout, reduction_before_kernel_check.wl+stdout,
  inactive_grade_projection_check.wl+stdout).
- **Leg B — Grok (grok-4.6, high).** VERDICT 9 findings (2 blockers). Evidence
  `~/.s11_build/S11c_b_89b_wl_dirleg_grok/` (jet_depth_and_freeze.py+stdout). Clean launch (no path typo).

## The converged root cause (both legs, computation-backed — verified by me, rule 13)
The retained-grade Hessian curvature the spec requires (`∂²W_bg`, `∂²μ_R,bg` at grade η,σ_W ≤ 1) is the
**product-rule content of a LIVE divergence/variation** (`c∇·F ≡ −(∇c)·F`, §1d). Once the profile rule
replaces `W^{(n)}(x)` by a **constant** jet-symbol (`profileRulesRetainingGeneratedJets`, engine L339–357),
that product rule is gone. Grok computed: `RESIDUAL_UNEXPANDED_FROZEN_MINUS_LIVE_REDUCED_FLUX = 0` (a
live-then-reduce that does NOT differentiate through the coefficient reproduces the **frozen** operator);
`EXPANDED_LIVE_CARRIES_A2_HESSIAN = True` while `FROZEN_EXPANDED_CARRIES_A2_HESSIAN = False`. Codex computed
the same for the kernel: reduce-then-differentiate loses the 3rd jet (`wJet3` absent), differentiate-then-
reduce keeps it. ⇒ The fix is NOT "reduce the operator, then extract"; it is **keep `widthBase`/`modulusBase`
(and the density constraint) LIVE FUNCTIONS through EVERY order-raising differentiation — bulk EL, in-plane
divergence handling, the §3c weak-block extraction (`extractCoupling` IBP split + Helmholtz `curl`/`div`), and
face EL — with the jet-retaining reduction as the ABSOLUTE LAST step.** My v1 directive under-specified this
(said "reduce the resulting rows/operator, then `extractCoupling` consumes it") and offered `rawModel`+reduce
as a valid realisation — which Grok showed reproduces the frozen operator (BLOCKER).

## Findings and disposition (all verified against code/stdout; all FOLDED into v2 unless noted)
1. **BLOCKER (Codex 1 / Grok 1+2) — reduction ordering / partial-freeze.** Fold: state the construction
   INVARIANT (coefficients live through every order-raising differentiation incl. the kernel's `curl`/`div`;
   reduction last). ⛔ Withdraw `rawModel`+reduce as a sufficient recipe. Preserve the `Inactive[Div]` split
   INTO `extractCoupling` so the §3c weak-block split (`splitDivergenceRows`: IBP for divergence terms,
   Helmholtz for the remainder) is taken on the un-reduced operator; reduce the kernel AFTER the pairing —
   Grok showed fully activating the divergences first mis-classifies everything as remainder and changes the
   kernel split. Add a **tower-depth control** (truncate the retaining depth one order below generated depth
   ⇒ an emitted object must MOVE; extend one order above ⇒ must NOT), ⛔ no numeric depth in the pass condition
   (rule 5).
2. **BLOCKER (Codex 2) — grade leak across `Inactive[Div]`.** `truncateScalar` (engine L99–117) protects
   inactive operators as tokens and truncates exterior/interior SEPARATELY, so `σ_W·Inactive[Div][σ_W·f]`
   survives while plain `σ_W²·f`→0. Reachable via `muTheta·thetaRule` (L1031) and the face path. Fold:
   require a **global multigrade projection** that accounts σ_W/η factors jointly across `Inactive[Div]`
   boundaries, plus a grade-support control proving σ_W² absent (incl. across Div) while ησ_W survives.
3. **SHOULD-FIX (Grok 3) — constraint/`thetaRule` omitted from the live set.** `evaluatedModel` `applyProfile`s
   the constraint at L1141 before `constrainedRows`; the constraint is the density relation (ρ_4D,bg⁰/ρ_br,bg⁰,
   §3b L276–278). Fold: put the constraint / `THETA_VARIATION_RULE` in the same live-through-EL set.
4. **SHOULD-FIX (Codex 3 / Grok 5) — `DIVERGENCE_FORM_SOURCE` (+ `FACE_SHAPE_SUBSTRATE`) coverage.**
   `processModel` (L1120–1121) skips reducing `DIVERGENCE_FORM_SOURCE`, which `modelRecord` (L1290–1292)
   emits on the `SLAB_OPERATOR` payload. Fold: name these surfaces; require live-then-reduce consistency (or
   explicitly pre-reduction geometry); emit an origins-sum-minus-production residual for slab + kernel.
5. **SHOULD-FIX (Codex 4 / Grok 4) — controls must re-enter at the energy; Control 1 isn't the form ablation.**
   The real FORM ablation = **Hessian-zero of the corrected LIVE operator** (set the order-2/3 background-jet
   atoms → 0 on a copy, show it moves toward the frozen construction) — that is what PY used. Replaying the
   frozen path is the REGRESSION check, not the form control. The one-sided corruption (v1 Control 3) must
   inject at one anchoring branch's **energy/background input** and rerun the whole constructor, not "drop a
   jet in one branch's reduction" (a result-stage corruption that violates §7). Clarify "routes" = anchoring
   branches. Fold: rewrite the controls; per-surface Hessian-atom witness taken FROM the EL.
6. **SHOULD-FIX (Grok 6) — §5.E replacement independence.** The replacement is independent only if the second
   operand comes from **primitive** atom dimensions of the emitted invariant EXPRESSION (`uOne`, `D[uOne,x]`,
   `Derivative[widthBase]/WZero`, …), NOT from walking `FACTOR_NAMES`/`dofFactorSpecifications` (else
   tautological again; Grok: mutating the shared table leaves residual 0). Fold: forbid `FACTOR_NAMES`/
   `dofFactorSpecifications` as that operand; else delete + emit the derived dimension with a note.
7. **SHOULD-FIX (Grok 7) — tower-depth control absent.** Folded into #1's tower-depth control.
8. **NIT (Grok 8) — citation.** v1's "~L1979 formOnly" is wrong: `formOnly` early-return is L1165; parametric
   `formOnly->True` is L2064/L2203; L1979 is the `corrupted` (6th-arg) call. Fold: fix the citation.
9. **NIT (Grok 9) — CONFIRMATION, no change.** Keeping `operatorFreezeRankDiagnostic` is right (abstract
   per-record rank, a different object); ⛔ do NOT add a production-rank residual against it (would turn a
   withheld magnitude into a tuning target). The withholding is clean. No change needed.

## Confirmations (both legs, "not a finding" — recorded so #89b doesn't re-open them)
- Jet depth: bulk EL ≤ order 2; one further derivative → order 3; order 4 NOT generated. The retaining table
  covers 2 & 3. The residual freeze is the reduction-before-a-later-derivative, NOT the table depth.
- Faces are not a second Hessian factory (Eulerian LAB_HELD virtual work `= −traction·vζeta`, no W jets); the
  real face coupling is inherited frozen `muTheta` = the bulk-EL fix (#1) + the constraint fix (#3).
- §3d admissibility (`backgroundBalanceFromModel`, L1845–1859) already does expanded EL then live reduction —
  verify-unchanged is the right call (do not double-fix).
- §3b/§3a quotes faithful; withholding clean (public 40 = basis count; no production rank/termcount stated).

## Disposition
Two legs, one decision pass (rule 7). The findings SHARPEN a sound directive (the un-freeze goal + the
live-jet operator object are correct; jet-depth math and §3d/face scope confirmed) — not a wrong shape, so
**fold once and go**, ⛔ do not re-run the decision legs (never iterate to green). NEXT after the fold: re-leak-
gate → Codex build (WL, `--sandbox danger-full-access`, xhigh) → 2 BUILD legs (fresh Claude + Grok,
Mathematica SERIALIZED).
