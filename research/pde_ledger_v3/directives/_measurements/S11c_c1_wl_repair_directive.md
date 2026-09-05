# _measurements — S11c_c1_wl_repair_directive.md (the repair fold + its 2-leg re-review)

> **⛔⛔ CORRECTION (2026-09-04) — this repair DIRECTIVE skipped its rule-7 decision legs; run retroactively they
> found it NOT sound, and the repair-2 decision gate then relocated the parity finding. Read with these overrides:**
> 1. **The repair directive was built WITHOUT its 2 decision legs** (rule 7 gap; [[feedback_directive_design_review]]).
>    The remediation is `_measurements/S11c_c1_wl_remediation_plan.md`; the retroactive decision-leg record is
>    `_measurements/S11c_c1_wl_repair_directive_review.md`.
> 2. **The "byte-identical protected core" below (lines ~16–22) OVER-PROTECTED two objects that were actually
>    defective (emit-only):** `NONINVERTIBILITY_CONDITION`/`fredholmFunctionSpaceOperator` still carried the
>    single-leg freeze, and `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` carried a dead parity axis. Both are fixed by
>    repair-2 (baseline `13f0bd2c`).
> 3. **The "Residual NIT — R4b parity LABELLING" section below is REFUTED.** `PERMEABLE_PORT_HERMITIAN` is
>    **correct**: its diagonal blocks are the congruence `Aᵀ·diag(P₊,P₋)·A` (both even combinations, coupling odd),
>    correct given `V_s=s∂_tζ_s`. The claim "the `DELTA_W` diagonal block should be `(plus−minus)`" is a category
>    error (a coordinate's parity ≠ the parity of its diagonal block in a quadratic form). See the plan's
>    CORRECTION note and `_measurements/S11c_c1_wl_repair2_directive_review.md`.
> 4. **"c1's WOLFRAM per-engine side is DONE" (Verdict) OVERSTATES** — per-engine review is complete; c1 is not
>    done until the remediation (repair-2 + records + both `.out`) and T7. See [[project_s11c_c_state]].

The repair directive folds the three MUST-FIX of the WL build review (`_measurements/S11c_c1_wl_build_review.md`;
2 legs, Codex-written engine ⇒ fresh Claude Agent + Grok) plus four biting NITs. The repair was itself 2-leg
re-reviewed (prompt `directives/_legs/S11c_c1_wl_repair_review.md`, serialized: fresh-Claude Agent then Grok, both
ablating Mathematica). **Both verdicts: all four repaired controls now BITE; the 2-leg-sound core is unchanged;
NO MUST-FIX.** Raw Grok report `…/scratchpad/grok_wl_repair_review.log`; leg-1 report returned inline. Both legs'
saved ablation scripts + stdout are under `/tmp/{...}/wlrepair1_*` and `/tmp/s11cc1_wl_rereview/`.

## Commands (literal)
```bash
# Repair build (Codex, detached; danger-full-access)
codex exec -c model_reasoning_effort=xhigh --sandbox danger-full-access "$(<…/S11c_c1_wl_repair_directive.md)" \
  > …/scratchpad/codex_wl_repair.log 2>&1
# Re-review leg 1 — fresh Claude Agent (in-process), serialized FIRST
# Re-review leg 2 — Grok, serialized SECOND (both ablate Mathematica; 2-seat)
grok --prompt-file …/_legs/S11c_c1_wl_repair_review.md --cwd /var/projects/toy_physics \
  --model grok-4.6 --effort high --permission-mode bypassPermissions --output-format plain \
  > …/scratchpad/grok_wl_repair_review.log 2>&1
```
⛔ **These re-review legs are NOT the directive's decision legs** — the repair DIRECTIVE's rule-7 decision gate was
skipped (see the CORRECTION banner); it was run retroactively in `_measurements/S11c_c1_wl_repair_directive_review.md`.

## Deliverable verification (orchestrator, rule 13; vs committed baseline `e139bc61`)
- `git diff --stat e139bc61 13f0bd2c -- <the .wl>` = **+443/−103, engine-scoped (one file)**; the full commit
  `e139bc61..13f0bd2c` is **4 files, +691/−103** (the .wl + `S11c_c1_wl_repair_directive.md` +
  `_legs/S11c_c1_wl_repair_review.md` + this record). Hunks land in `operatorCompositionFromDerivation` (R1),
  `deriveEnergyOperands` + the outgoing-half-space energy construction (R2/R3), the response traction/source-
  equation plumbing (R3/R4a), the real-admissibility test (R4d), the parity-port blocks (R4b), the layer-potential
  operand + `REP_INVARIANCE` emits (R4c). 51 tags unchanged; all tasks ran; blindness scan 0; `mathematica/out/`
  clean; builder isolated-vs-in-repo byte-identical.
- Both legs independently confirmed the protected core byte-identical: `directBoundaryDerivation` sha
  `557ed758…` = baseline (two-leg `DTN_KERNEL`), `DTN_FLAT_SYMBOL`, `DTN_RIGID_SHIFT_*`, `DTN_BY_REGIME_PAIR`,
  `DTN_BY_PARITY`, `DTN_HERMITIAN_PART`/`DTN_REACTIVE_PART`, `DTN_GRAZING_BEHAVIOUR`, `DTN_TERM_ORIGINS`,
  `FACE_RESPONSE_COEFFS`, `NONINVERTIBILITY_CONDITION`, the four non-admissible `DEGENERATE_LOCI_*`, `ZERO_JET_*`,
  `DIMENSIONS`, `HOMOGENEITY_*`, `PERMEABLE_DISSIPATION_VS_OMEGA_TAU`, the T-a..T-i maps, `μ_θ` opacity, reserved
  spellings. Leg 1 proved the differing tags are all R1–R4 + benign `SOURCE_EQUATIONS` propagation (blanking every
  `SOURCE_EQUATIONS` value makes every §5 control byte-identical to baseline).

## The three MUST-FIX — each now BITES (both legs, literal one-sided corruption)
- **R1 `DTN_OPERATOR` carries both legs.** `COMPOSITION_HAS_Q_INPUT`/`_Q_OUTPUT` → `PROVED_TRUE`
  (`qOut(k′)`×2, `qOut(k)`×3 incl. the flat symbol); composition kernel-expands on-shell to `DTN_KERNEL` (residual
  0); **re-freezing the rightmost `N_0` to the output leg drops `qOut(k′)` and the residual goes nonzero** (both
  legs' independent probes; the output-frozen control does NOT match `DTN_KERNEL`). Grok also noted the fix added
  a Fredholm dummy-momentum split so the two-leg operator is not reused as the function-space condition.
- **R2 energy BULK operand = far-field Poynting of `φ` at `|w|→∞`.** `farFieldPhase` fully removed (repaired 0 vs
  baseline 2); bulk operand = `Inactive[Limit][∫ ½Re[δp·v_n*] of the outgoing φ, outwardDistance→∞]` (zeroth-order
  `φ = (V/(I q)) e^{i(k′·x+q w)}`, the half-space outgoing wave — NOT `δp_face·farFieldPhase`). One-sided bulk
  `q_out` branch flip (parameter AND source-level) moves ONLY the bulk operand + `ENERGY_RESIDUAL`, face
  unchanged; NO substitution (`farFieldPhase→1`, releasing the `Limit`) collapses it to a structural zero; face
  vs bulk not `SameQ`. Routes independent.
- **R3 energy FACE operand binds the response `t_s`.** `FACE_HAS_TS_FROM_RESPONSE` → `PROVED_TRUE`,
  `REFERENCED_RESPONSE_FIELDS = {TRACTION_FLAT_VECTOR, TRACTION_FIRST_SHAPE_KERNEL_VECTOR}` (the operand carries
  `LambdaX0`/`LambdaA0`, not a fresh impermeable `−δp n̂`). A one-sided sign flip of the **response** `t_s` at its
  source moves ONLY the face operand + `ENERGY_RESIDUAL`, bulk unchanged. A response `t_s` error is now caught.

## R4 NITs — folded (both legs)
- **R4a `SOURCE_EQUATIONS`** — real re-parseable `Inactive[Equal][FacePressureUnknown[…],…]`; 0 `flatBulk$NNNNN`
  gensyms.
- **R4c `§5a` layer-potential** — decorative `RadiationPreservingLayerPotential[…]` head gone; the second-route
  operand IS the computed `layerMain["KERNEL"]` (data-dependent, `layer−direct=0`), and `REP_INVARIANCE_RESIDUAL`
  references that same computed expression.
- **R4d locus `REAL_ADMISSIBLE`** — genuine `Resolve[Exists[…], Reals]`; the live branch → `UNDECIDED`, ⛔ no
  longer auto-`ADMISSIBLE` (`True→ADMISSIBLE`, `False→EXCLUDED`, non-boolean→`UNDECIDED`).
- **R4b `PERMEABLE_PORT_HERMITIAN`** — the duplicated matrix is gone (both legs: `DELTA_W`≠`ZETA_C`, computed
  nonzero residual + `PARITY_COMBINATION_COMPARISON`).

## Residual NIT (⛔ NOT a MUST-FIX; carried to the c1 step record + T7, ⛔ not a 2nd repair round — rules 10/15)
- **R4b parity LABELLING** (leg 1, per its probe; Grok checked "distinct not duplicated" only, did not examine the
  labelling): under `s→−s`, `ζ_c` is even and `δW` is odd (`zetaForFace = ζ_c + s·δW/2`), so the `DELTA_W`
  diagonal block should be the ODD combination `(plus−minus)` — which **vanishes** at this order — while both
  emitted diagonal blocks carry the scaled EVEN (`ζ_c`) form and the odd combination is routed to the (vanishing)
  off-diagonal coupling. The **content** is correct (even blocks computed, odd = 0 at this order); only the
  `DELTA_W`/`ZETA_C` **labels** are mis-assigned in this emit-only §3a parity diagnostic (`PERMEABLE_PORT_HERMITIAN`
  is NOT exported). Carry to the step record and flag for the T7 comparator (a possible spurious parity
  disagreement vs PY, adjudicated after the run — like the omega/μ_θ representational residuals), ⛔ not a blocker.

## Verdict
Repair-1 CLEARED its own 2 re-review legs (all four re-reviewed controls bite; core byte-identical) — but see the
CORRECTION banner: the repair DIRECTIVE skipped its decision legs, which (run retroactively) + the repair-2
decision gate found two emit-only defects the re-review was structurally blind to (over-protected
`NONINVERTIBILITY_CONDITION`; dead parity axis in `PERMEABLE_DISSIPATION_VS_OMEGA_TAU`). ⭐ **c1's WOLFRAM
per-engine REVIEW is complete; c1 is NOT done** — NEXT = the FULL REMEDIATION (repair-2 → records → both `.out`),
THEN the T7 cross-engine comparator. Plan: `_measurements/S11c_c1_wl_remediation_plan.md`.
