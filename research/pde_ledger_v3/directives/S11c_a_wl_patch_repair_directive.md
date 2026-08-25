# S11c-a WL engine PATCH — REPAIR (one blocking defect from review)

You are repairing the blind Wolfram (WL) engine you just patched. An independent review found ONE blocking
defect. Fix ONLY it. ⛔ PRESERVE the three schema/coverage changes (all verified correct) and everything else.

## ⛔⛔ BLIND ENGINE + RULE 5
This engine imports nothing and re-derives from the spec `directives/S11c_a_SHARED_PHYSICS.md`. Cite ONLY the
spec. ⛔ Do not import or reconstruct any sibling engine. ⛔ State no expected value; a residual is PRINTED,
never asserted, and its value is never a target. The three script clauses still bind: PRINT computed objects
(never prose), PRINT the residual (never `assert ==0`), interpretation belongs to the step record.

## THE BLOCKING DEFECT (verbatim from the independent review; orchestrator-confirmed)
> §5a rep-invariance MATERIAL route is UNCOMPUTED for T-c, T-d(traction), T-i.
> `WL_S11CA_REP_INVARIANCE_MATERIAL_OPERAND` and `WL_S11CA_REP_INVARIANCE_RESIDUAL` are the literal
> unevaluated symbol `materialShape["LAB_HELD", 1, "RHO4_CONSTANT", "FLUX"]` (and …"TRACTION", …"CLOSURE")
> for all 16/16 entries each of RELATIVE_FLUX (T-c), TRACTION (T-d traction), CLOSURE_SHAPE_DERIV (T-i) —
> 48 operand + 48 residual entries. VIRTUAL_CONSTRAINT (T-g), EVOLUTION_MASS_BALANCE (T-h), VIRTUAL_WORK
> (T-d work) compute their material route correctly (0 inert).
> Root cause: `materialShape` is defined only at line 595-597 (memoizing `f[x]:=f[x]=…`). It is Cleared in the
> reset blocks (`Clear[…, materialShape]` at ~1456 and ~1973) and is NEVER redefined — unlike `directShape`,
> which has `resetDirectShapeCache`/`resetDirectShapeMemoCache` (582-592) to redefine after clearing. When the
> §5a object loop consumes it (line 1332), the symbol has no definition and returns inert.
> Effect: residual = Eulerian − materialShape[inert]; not the Eulerian-vs-material comparison §5a requires; it
> can never be zero and encodes no route-agreement test — a control that CANNOT FAIL. Corroborated: flipping
> the Eulerian `graphNormalSource` moved `REP_INVARIANCE_EULERIAN_OPERAND` and `_RESIDUAL` but left
> `_MATERIAL_OPERAND` byte-identical.

Orchestrator provenance check (rule 13): PATCH-INTRODUCED. `git diff` shows the patch ADDED the `Clear[…,
materialShape]` reset blocks with only `resetDirectShapeCache[]` after them; the committed baseline `.out`
carries computed material operands and ZERO inert `materialShape[` tokens.

## Spec (why this must compute)
§5a (lines 452-465): "compute T-g, T-h, T-c, T-d, and T-i by **two independent routes**: (1) direct Eulerian …
(2) material-coordinate differentiation using `x=x(X,t)` and the exact face-flattening `w′` … Map route 2 back
to the Eulerian variables before comparison. … For every compared T-object emit the two operators and their
difference under `…EULERIAN_OPERAND`, `…MATERIAL_OPERAND`, `…RESIDUAL`." The material operand must be the
genuinely-computed route-2 object, ⛔ not an inert symbol. §5a (lines 466-474): the one-sided independence
control corrupts one route at a time — which is only meaningful if BOTH routes compute.

## The fix
Ensure `materialShape` is DEFINED (not inert) wherever the §5a controls consume it (line 1332). Model it on
`directShape`'s treatment: `directShape` survives its `Clear` because `resetDirectShapeCache[]` /
`resetDirectShapeMemoCache[]` (582-592) redefine it. Give `materialShape` the same guarantee — either add and
call an analogous `resetMaterialShapeCache[]` (redefining `materialShape` from `materialFaceObjectsSource`
after every `Clear` that names it), or drop `materialShape` from the `Clear` lists if clearing it serves no
purpose (its memo key already carries `(branch,sign,density,name)`). ⛔ Do not change the material route's
construction (`materialFaceObjectsSource` / the w′ face-flattening) — it is the verified §5a route-2; only its
availability at consumption time is broken. Fix nothing else.

## Acceptance — structural, ⛔ no numeric target (rule 5)
1. `REP_INVARIANCE_MATERIAL_OPERAND` and `REP_INVARIANCE_RESIDUAL` for `RELATIVE_FLUX`, `TRACTION`,
   `CLOSURE_SHAPE_DERIV` are COMPUTED expressions containing NO inert `materialShape[` token. Grep the whole
   emitted stream: zero inert `materialShape[` tokens anywhere.
2. The `REP_INVARIANCE_RESIDUAL` for these objects is a genuine difference of two independently-computed
   routes — the control must be ABLE to fail: a one-sided corruption of the Eulerian route alone must move the
   residual, and a one-sided corruption of the material route alone must move it (§5a:466-474). ⛔ Do not state
   what the residual equals; demonstrate the control can move and PRINT it.
3. The three patch changes (projection density axis; kinematic/flux axis drop; §5b/§5c coverage of the 5
   objects) and every other tag are behaviourally unchanged — report the diff is confined to the `materialShape`
   availability fix.
4. Full run evaluates clean (report literal stderr + tag inventory); no new `assert`-before-emit.

## Ops
Codex xhigh, `--sandbox danger-full-access`, background. Wrap kernel runs in a generous `timeout`; ONE kernel
at a time (two seats); watch RSS/orphan kernels. Patch only the one `.wl`. Transcript OUTSIDE the repo.
