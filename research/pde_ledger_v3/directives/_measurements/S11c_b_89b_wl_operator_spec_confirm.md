# S11c-b #89b WL side — spec-confirm (§3b binds WL; operator un-freeze = spec-compliance, no spec change)

⚠ ANSWER-BEARING (records the frozen 26 / live 40 diagnostic and the #88 PY-measured corrected-operator
numbers). For the orchestrator + directive-drafting + review record only — ⛔ NOT handed to the builder.

Prereq to drafting the #89b WL operator directive (resume step 1: "spec-confirm §3b WL-side → draft").
Records the claims the directive's §1–§2 rest on, each with the command/lines behind it (rule 2). No engine
changed; no spec changed. Every WL line number was re-located by NAME on the post-#89a file (HEAD `60e31504`).

## Claim 1 — §3b (272–282) is ENGINE-NEUTRAL, names the OBJECT, and binds WL
`research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md:272–282` (read verbatim): compute the
divergence-form first-order EOM "so their first jets appear explicitly. **Retain the full spatial dependence
of every background coefficient (`μ_R,bg`, `W_bg`, `ρ_4D,bg⁰`, `ρ_br,bg⁰`, and the `Σ_E⁰` map) and its first
jet; do not freeze a coefficient at its constant binding before differentiation.**" §3a:248–254 adds that the
first-background-*jet* order "bounds the number of independent background-**amplitude** factors — powers of
`η`/`σ_W` — **not the spatial-derivative order**", and "a second spatial derivative of `W_bg` is still first
order in background amplitude, `O(η)`/`O(σ_W)`, and is **not dropped**." §1d:163–168: the variable-coefficient
representative difference (`c∇·F ≡ −(∇c)·F`) is "physics in the operator/kernel, not a representational
identity."
⇒ The spec names an OBJECT — the variable-coefficient divergence-form slab operator whose construction keeps
the whole background jet tower LIVE through the EL variation, retaining every higher spatial jet at its
originating background-amplitude grade (η,σ_W ≤ 1). It prescribes NO mechanism (neither PY's differentiate-
then-reduce nor a WL rule-order). Both engines owe the SAME operator object by their own route.
**No spec change is required; the WL fix is a spec-compliance repair of a freeze-before-differentiate defect.**

## Claim 2 — PARITY: PY #89 already un-froze BOTH §3a-quotient AND §3b operator (verbatim)
`directives/S11c_b_89_sympy_3a_repair_directive.md:5–13` (the PY #89 deliverable): "Modify in place the SymPy
engine … so that every total in-plane divergence in its **§3a energy-basis quotient AND §3b operator
construction** is the variable-coefficient one — the background profile and its spatial-jet tower `{W_bg,
∂W_bg, ∂²W_bg, …}`, `{μ_R,bg, ∂μ_R,bg, ∂²μ_R,bg, …}` **are differentiated, not held constant** — to the jet
depth each consumer needs." L15: "a spec-compliance repair of an existing implementation defect, not new
physics and not a spec change."
⇒ On the PY side, both surfaces were un-frozen together in #89 (PY CLEARED `9f40c18e`). WL's #89a repaired
ONLY the §3a basis (emits 40, cleared `d4adbd99`); the WL §3b **operator is still frozen**. #89b brings the WL
operator to parity — a **DIFFERENT mechanism** (WL freezes via a jet-zeroing profile rule applied before EL;
PY froze via a constant lookup-map), the SAME object.

## Claim 3 — WL's PRODUCTION operator FREEZES the profile BEFORE differentiating (code-verified)
`mathematica/S11c_b_brane_operator_mathematica_audit.wl`, `evaluatedModel` (L1130–1195):
- L1136 builds `energyData` on the LIVE functions `widthBase`, `modulusBase`.
- **L1137–1139 `selectedRecords = applyProfile[…]` FREEZES the profile, then L1140 `energy = Total[… ENERGY_TERM]`** — the energy fed downstream is already frozen.
- L1143 `rows = constrainedRows[energy, constraint]` and L1181–1184 `variationalSource[energy, …]` take the
  EL variation on the ALREADY-FROZEN energy. `constraint` (L1141), `width` (L1144), `rhoBrValue` (L1146),
  `facesData` (L1148, via `faceSources`), `evolutionTerms` (L1152, via `evolutionSource`) are each
  `applyProfile`'d before entering the operator too.

`applyProfile` (L363) `= truncateBackground[expr /. profileRules[…]]`; `profileRules` (L312) →
`profileDerivativeRules` (L297), whose **`higherRules` (L304–307) map every jet with `2 ≤ i+j+k ≤ 3` to `0`.**
⇒ `applyProfile` **zeroes every 2nd and 3rd spatial jet** of `W_bg`/`μ_R,bg` — precisely the retained-grade
Hessian `∂²W_bg`/`∂²μ_R,bg` that §3a:251–254 says is "not dropped." Freezing before the divergence/variation
is the §3b violation.

## Claim 4 — the freeze CHANGES the operator (MEASURED, non-absorbable) — not spec-neutral
- #89a operator-freeze diagnostic (`operatorFreezeRankDiagnostic`, L1219–1244; emits
  `WL_S11CB_OPERATOR_BACKGROUND_FREEZE_DIAGNOSTIC`): `PRODUCTION_FROZEN_EL_RANK = 26` vs
  `LIVE_BACKGROUND_EL_RANK = 40`, difference 14, NOT asserted equal. `grep -a -o
  'WL_S11CB_OPERATOR_BACKGROUND_FREEZE_DIAGNOSTIC[^}]*}' mathematica/out/…audit.out`.
- #88 blast radius (`_measurements/S11c_b_88_blast_radius_result.md:4–22`) characterizes the Δ: the correction
  "retains the background curvature — the Hessian `∂²W_bg`, `∂²μ_R,bg` — … at retained grade (η,σ_W ≤ 1)", with
  `RANK_GAIN > 0` on every strong row ⇒ "no constant-coefficient reparametrisation of the frozen row can
  absorb" it, and the family verdicts (KINETIC, θ) "rested on the frozen operator" and must be re-adjudicated.
⇒ The freeze is a genuine loss of retained-grade operator content, not a representational choice. **These
numbers (26, 40, Δ14, the #88 termcounts/RANK_GAIN) are the answer — WITHHELD from the #89b builder (rule 5).**

## Claim 5 — the correct live path ALREADY EXISTS in the file; the fix is EL-before-freeze then live-reduce
`profileRulesRetainingGeneratedJets` (L339–357) retains every jet `i+j+k ≥ 1` as
`sigmaW · profileJetSymbol/LWidth^(n−1)` (a SINGLE σ_W factor — §3a retained grade), and
`applyBackgroundProfileWithGeneratedJets` (L360) reduces with it. The diagnostic's `liveRows` (L1231–1233)
demonstrates the exact fix pattern: **`applyBackgroundProfileWithGeneratedJets[elRowVector[ENERGY_TERM,
constraint]]`** — take EL (`elRowVector` → `constrainedRows`) on the LIVE energy first, THEN reduce with the
jet-retaining path. The production operator must follow this pattern instead of `applyProfile`-before-EL.
⛔ rule 1: the builder builds its OWN production integration (faces/flux/evolution/kernel), ⛔ NOT a
transliteration of the minimal abstract diagnostic.

## Claim 6 — blast radius of the WL fix (what the directive must cover)
The frozen `applyProfile` calls that feed the operator, all in `evaluatedModel`: energy records
(L1137–1140), `constraint` (L1141), `width`/`rhoBrValue` (L1144/1146), `facesData` via `faceSources`
(L1148–1150), `evolutionTerms` via `evolutionSource` (L1152). Downstream consumers of the operator:
`constrainedRows` (L1022; already live-capable — the defect is the frozen INPUT, not the routine),
`variationalSource` (L275), `extractCoupling` (L1376; the §3c kernel is a weak block of the §3b operator),
and the operator emits (`S11CB_SLAB_OPERATOR`, `S11CB_SLAB_OPERATOR_TERM_ORIGINS`, `S11CB_MU_THETA_OPERATOR`,
`S11CB_COUPLING_KERNEL`, `S11CB_COUPLING_KERNEL_TERM_ORIGINS`).
⚠ The §3d ADMISSIBILITY ε⁰ operand (`admissibility_operator` route) is per #88:42–44 a SEPARATE, already-
Hessian-aware route — the directive must CHECK it stays consistent, ⛔ not double-fix it into the wave operator.

## Claim 7 — §5.E dimension-residual tautology (folded control nit from #89a)
`…_audit.wl:2269–2290`: `RESIDUAL_DERIVED_MINUS_STORED` differences two copies of `dimensionGradient + Σ
(factor dims)` over the same `dofFactorSpecifications` ⇒ structurally 0, cannot catch a wrong invariant
dimension ([[feedback_a_check_cannot_audit_its_own_input]]). Benign (dims validated OUT-OF-BAND by both #89a
legs). Directive obligation: replace with an INDEPENDENT dimension route (dimensional analysis of the emitted
invariant EXPRESSION vs the factor-metadata sum), or DELETE the tautological residual.

## Rule-5 ledger for the #89b directive (what is public vs withheld)
- PUBLIC (hand to builder): the §3b/§3a/§1d spec obligation (Claim 1); the freeze diagnosis (Claim 3, code
  lines); the fix PATTERN (Claim 5); the basis count **40** (already emitted + cleared in #89a — public); the
  frozen public regression **26**; that the frozen and live EL ranks differ (the diagnostic already emits both).
- WITHHELD (rule 5, acceptance-value): the corrected PRODUCTION operator's expected rank / termcounts /
  RANK_GAIN per row (the #88 PY-measured values, e.g. termcounts 98/98/98/108/186). ⛔ Do not state a target
  rank. Acceptance = spec-compliant CONSTRUCTION (EL on live energy, retained-grade jets) + compute-and-emit;
  ⛔⛔ verify the unfrozen operator's cross-engine CONTENT/SPANS, ⛔ not a rank number (resume-prompt standard,
  [[feedback_basis_independence_must_not_freeze_spurion]]).

## Disposition
Spec-confirm PASSES: §3b binds WL (engine-neutral, object named), the WL operator unfreeze is a **spec-
compliance repair + parity with PY #89**, DIFFERENT mechanism, **no spec change**. The correct live path
already exists in the file (Claim 5). NEXT = draft `S11c_b_89b_wl_operator_directive.md` (orchestrator-written)
→ 2 DECISION legs (Codex + Grok, rule 7) BEFORE any builder → Codex build (WL, `--sandbox danger-full-access`)
→ 2 BUILD legs (fresh Claude + Grok, Mathematica SERIALIZED). Fold the §5.E fix (Claim 7). Withhold per the
rule-5 ledger above.
