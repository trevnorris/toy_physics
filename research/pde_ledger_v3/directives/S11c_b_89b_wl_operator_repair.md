# S11c-b #89b (WL) — repair round: emitted-operator re-freeze + two broken controls (v2, decision-reviewed)

## 0 · Role and single deliverable

The engine `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` (working tree,
uncommitted) passed its physics on the KERNEL path, but the two build-review legs + orchestrator verification
found three defects (record: `directives/_measurements/S11c_b_89b_wl_operator_build_review.md`). Fix all three.
⛔ Do NOT commit — the orchestrator commits after re-review. ⛔ Change nothing else — the operator un-freeze
mechanism (`constrainedRowsWithLiveEnergyEL`), the tractability activation (the `operatorLive`-based path at
~L2238–2254), the deferral gate (`S11CB_SKIP_HEAVY_CONTROLS`), and every unrelated emit stay exactly as they
are — **EXCEPT** the specific comparison operands named in §1 that MUST be brought to the same form as the
corrected operator to keep their residual like-with-like (they are separate constructors and will otherwise
report a form artifact instead of physics). The full `.out` is NOT regenerated here (needs a ≥64 GB box;
deferred — `DEFERRED_HEAVY_RUNS.md`); do the reduced-scale self-check in §4, then STOP.

Line numbers below are from the current working tree and were confirmed by the decision review; treat them as
signposts, not as a promise the builder may skip re-reading the surrounding code.

## 1 · Defect A (BLOCKER) — the EMITTED operator (and everything compared to it) is re-frozen on the outer divergences

**What is wrong.** The emitted `OPERATOR` (`:1404`) and `SLAB_OPERATOR_TERM_ORIGINS` (`:1405`) apply the final
background reduction to the operator whose OUTER constrained `Inactive[Div]` are STILL HELD (the un-activated
`operatorLive` / `originsLive`). The outer Divs come from `linearVirtualVariation` (`:324–326`), NOT from the
already-live inner energy EL — so they are real and un-activated. The correctly-activated objects are computed
one step earlier (`activatedOperator = activateSpatialDivergences[operatorLive]` `:1376`; `activatedOrigins`
`:1378`) and then DISCARDED. Reducing before activating those outer divergences freezes the background field
inside the held `Div`, so the mixed and higher background-jet terms a live activation generates are permanently
lost. The sibling emits already use the activate-THEN-reduce route (`MU_THETA` `:1377`,`:1403`; faces
`:1379`,`:1406`; divergence-form source `:1380–1381`,`:1407–1408`; θ-variation rule `:1382–1383`,`:1420–1421`).
The L1372–1375 comment claiming the held-`Div` emit feeds the §3c weak split is FALSE — the split consumes
`KERNEL_SOURCE_OPERATOR` (`:1425`, `:2192–2195`, `:1679–1680`), not the emitted `OPERATOR`; correct or delete it.

**The object to produce (name, not recipe).** Every REDUCED emitted operator-family object, AND every operand
differenced or compared against it, must be the **activate-then-reduce** form: outer constrained divergences
activated to derivative-normal form BEFORE the final background reduction, so it retains the full background-jet
tower the live activation generates — identical treatment to `MU_THETA`. ⛔ This does NOT include the objects
intentionally kept un-reduced-live — the two `KERNEL_SOURCE_*` slots AND the `LIVE_DIVERGENCE_FORM_OPERAND` tower
slot (`:2212–2214`, the un-reduced live source captured from `KERNEL_SOURCE_OPERATOR` and emitted at `:2298–2308`);
leave those exactly as they are (see the invariant below). You are responsible for finding ALL the reduced/compared
objects; the decision-reviewed site table is your checklist (⛔ do not treat it as exhaustive — re-grep every
reader of `["OPERATOR"]`, `["ORIGINS"]`, and each per-depth/ablation/limit/cross-model operand):

- **Emit** `OPERATOR` and `SLAB_OPERATOR_TERM_ORIGINS` (`:1404–1405`; the `formOnly` return at `:1409–1412` uses
  the same assignment — one fix covers the parametric/uniform `formOnly` callers at `:2719`,`:2750`,`:2860`).
- **Origins-sum residual** (`:2276–2279`) stays like-with-like ONLY if `OPERATOR` and `ORIGINS` change together.
- **Tower-depth control** (`:2201–2220`): `operatorTruncated`/`operatorExtended` currently reduce un-activated
  `operatorLive` at a retained depth, and the depth itself is measured on un-activated `operatorLive` (`:2201`).
  Activating the outer Divs chain-rules one extra background order, so measure the generated depth on the
  activated operator and reduce THAT object at depth ∓1, so "truncate-one-below moves / extend-one-above does not"
  is still testing the emitted object's tower.
- **Hessian-zero** (`:2578–2588`) is a self-difference of the emit — stays like-with-like once the emit is fixed.
- **Frozen Hessian-witness** (`:2521–2562`, deferred control): `frozenEvaluatedModel` (`:1228–1293`) is a
  SEPARATE freeze-first constructor that never activates. The witness differences the production side (already
  activated) against the frozen side (un-activated) for MULTIPLE slots — `OPERATOR` (`:2523`), `ORIGINS`
  (`:2539`), `MU_THETA` (`:2526`), `DIVERGENCE_FORM_SOURCE` (`:2552`), and the face-substrate fields (`:2557`) —
  and `surfaceDifferenceJetAtoms` (`:1858–1872`) sees jet atoms INSIDE a held `Div`, so any form-only difference
  pollutes the witness (this ALREADY mis-compares MU_THETA/faces on the frozen side, and your fix adds
  OPERATOR/ORIGINS). For EVERY slot the witness differences, name a comparison-only operand: the derivative-normal
  form of the ALREADY-FROZEN object (activate-after-freeze, for the atom difference only). ⛔ Do NOT re-live
  `frozenEvaluatedModel`, and ⛔ do NOT change the frozen replay or the frozen kernel construction — the frozen
  side stays freeze-first (it is the re-freeze reference; its emit `CONTROL_OPERATOR_REFREEZE_REGRESSION_OPERAND`
  at `:2508–2511` must not change).
- **Uniform-limit S11b** (`:2860–2874`): after the fix `uniformLimit[model["OPERATOR"]]` (`:2862`) is
  derivative-normal but `uniformS11b["OPERATOR"]` (`:2807–2839`,`:2871`, built by `constrainedRows` with no
  activation) stays Div-form — so BOTH the `SLAB_OPERATOR` residual (`:2899–2902`) AND the `TRANSVERSE_DISPERSION`
  residual (derived from the SAME S11b operator through `planeWaveCoefficient` `:2845–2850`, which substitutes but
  does not activate) report a Div-vs-expanded artifact. Feed the SAME comparison-only activated S11b operator into
  BOTH the `SLAB_OPERATOR` and the `TRANSVERSE_DISPERSION` operands, ⛔ WITHOUT changing S11b's independent
  reconstruction physics, and ⛔ leave the original S11b object feeding `uniformS11bKernel` unchanged.

⛔ **INVARIANT — do not break the kernel.** `KERNEL_SOURCE_OPERATOR` MUST remain the un-reduced live operator
`operatorLive` (`:1425`), AND `KERNEL_SOURCE_ORIGINS` MUST remain the un-reduced live `originsLive` (`:1426`) —
the kernel origin split (`kernelOriginsFromOrigins[processed["KERNEL_SOURCE_ORIGINS"]]` `:2198–2199` →
`kernelFromOrigin` `:1805–1814`) relies on the held `Inactive[Div]` split the kernel path was verified on. The
`LIVE_DIVERGENCE_FORM_OPERAND` tower slot (`:2212–2214`,`:2298–2308`) is likewise the un-reduced live source and
stays un-changed. Only the REDUCED emitted/compared objects change. Both legs verified the kernel path is correct;
leave it.

**Acceptance — the builder PRINTS, it does not assert the target (rule 5).** ⚠ The independent reference MUST be
rebuilt from the un-reduced live slot: `finalBackgroundReduction[activateSpatialDivergences[KERNEL_SOURCE_OPERATOR
/ operatorLive]]` (activate-then-reduce). ⛔ Do NOT use `activateSpatialDivergences[emitted OPERATOR]` as the
reference — reduce-then-activate is ITSELF the frozen object, so that residual is identically 0 for BOTH the buggy
and the fixed engine (and a leftover-`Inactive[Div]` count of 0 is likewise true of the buggy order, `:1476–1478`),
which is why a bare residual-0 check would silently pass the non-fix. The builder's job ENDS at compute-and-PRINT:
emit, literally, (a) the residual between the fixed emitted operator and that from-`operatorLive` reference, and
(b) the SAME residual with one route swapped to the reduce-then-activate order — and STOP, for the orchestrator to
judge. ⛔ Do NOT state or assert what those residuals should equal, ⛔ do not name the expected jet orders or a
leftover-`Div`/coefficient/rank count, and ⛔ do not iterate toward any target value; print the objects and stop.

## 2 · Defect B (SHOULD-FIX) — the §5.E "independent primitive-atom" dimension walker never attaches

**What is wrong.** `primitiveExpressionDimension` (`:2942–2969`) is defined with composite argument patterns of
the form `Times[factors__]` / `Plus[terms__]`. `Times`/`Plus` are `Flat`+`Orderless`+`OneIdentity`, so those
definitions do not register as usable rules (the product rule is dropped; the sum rule is demoted below the
generic `expression_` catch-all). Every composite invariant (the Kronecker `Times`/`Plus` builds, not `Dot`) falls
through to the `UnassignedPrimitiveDimension` sink, so the emitted `RESIDUAL_DERIVED_MINUS_PRIMITIVE_EXPRESSION`
is non-numeric and does NOT depend on the primitive atom dimensions — the very defect §6 was to remove.

**The object to produce.** A dimension route that actually walks each composite invariant to its primitive atoms
(`uOne`, `D[uOne,x]`, `Derivative[widthBase]/WZero`, …) and returns their combined dimension. **Decisive
acceptance:** mutating one primitive atom's dimension MOVES the emitted residual (it currently does not).
⛔ Do not special-case the invariants or hard-code any dimension.

## 3 · Defect C (SHOULD-FIX) — the independence control is a `base − base = 0` tautology AND an unlike-forms weak split

**What is wrong.** `corrupted = If[branch === First[branches], evaluatedModel[…, corrupted->True], base]`
(`:2466–2468`); `First[branches]` is `"LAB_HELD"` (`:245`). For every later branch (the two MATERIAL_ADVECTED
cases) `corrupted = base = mainModels[key]`, so:
- `SLAB_OPERATOR` and `ADMISSIBILITY_OPERATOR` residuals are `base − base = 0` — a control that cannot fail; and
- the `COUPLING_KERNEL` slot `extractCoupling[corrupted]` (`:2475`) routes through
  `Lookup[model, "KERNEL_SOURCE_OPERATOR", model["OPERATOR"]]` (`:1679–1680`), but `base` has already had
  `KERNEL_SOURCE_OPERATOR`/`_ORIGINS` `KeyDrop`'d (`:2257–2258`), so it falls back to the reduced emitted
  `OPERATOR` — which after Defect A is derivative-normal with no top-level `Inactive[Div]`, so
  `prepareOperatorForWeakSplit` (`:1629–1643`,`:1681`) treats every term as direct. Not a corruption test.

**The object to produce.** Make every branch's independence operand a GENUINE one-sided corruption of that
branch's own input, OR, where the corrupted build is skipped for tractability, emit an EXPLICIT
`"VALIDATED_ON_REPRESENTATIVE_BRANCH" -> <key>` marker that replaces the WHOLE independence package for that
branch (SLAB **and** ADMISSIBILITY **and** COUPLING_KERNEL slots together) — never a silent `base − base` on some
slots with a stale kernel on another. The cleanest genuine fix is a fresh `evaluatedModel[…, corrupted->True]`
per branch (it still carries `KERNEL_SOURCE_OPERATOR`, `:1425`), if it fits the reduced-scale budget. The
`LAB_HELD` branch's live corruption test must remain live. ⛔ A skipped corruption must be self-documenting and
must never read as a passing check.

## 4 · Reduced-scale self-check (⛔ do NOT run the full engine / regenerate the `.out`)

Full-depth is ~16 GB/case and OOMs the 30 GB box. Verify at REDUCED scale only: set
`basisRepresentativeIndices = {16}` (or a small jet-bearing subset), extract the definitions you need into a
`/tmp` harness, wrap every kernel run in `timeout 600`, one kernel at a time, kill any leftover `WolframKernel`
by EXACT pid afterward. ⛔ PRINT the objects below literally and STOP — ⛔ do not assert what any of them should
equal, do not label pass/fail, and do not iterate toward a value; the orchestrator + re-review legs judge.
- **Defect A**: print (a) the fixed emitted operator's leftover `Inactive[Div]` count; (b) the residual between
  the fixed emitted operator and the reference `finalBackgroundReduction[activateSpatialDivergences[operatorLive]]`
  on the U row; (c) the SAME residual with one route swapped to the reduce-then-activate order; (d) the
  `KERNEL_SOURCE_OPERATOR`/`KERNEL_SOURCE_ORIGINS` activation-postcondition and higher-jet-depth objects (to show
  they are unchanged); (e) the frozen-witness and uniform-S11b difference operands (to show they are now the same
  form on both sides).
- **Defect B**: print the §5.E residual, then print it again after mutating one primitive atom's dimension (so
  the orchestrator can see whether it depends on that atom).
- **Defect C**: print a MATERIAL_ADVECTED branch's independence operands (base and corrupted for all three slots),
  or the explicit skip marker if the corruption is deferred there.
Report the harness path, the literal stdout, and a `git diff` of the engine. ⛔ Do not commit.
