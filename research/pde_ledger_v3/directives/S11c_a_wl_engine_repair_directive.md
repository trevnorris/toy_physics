# S11c-a Wolfram engine — full repair directive

## Authority and boundary

Edit `research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` **in place**
(committed baseline `277f3fe7`). Its product is the flushed stdout tag stream; that is its only write. It
remains the **blind** engine: it imports nothing and re-derives every object from the spec's equations. ⛔ Do
not add any import, and ⛔ do not introduce any object, value, or construction from another engine — this
engine's independence from the sibling is the design's one blindness control.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (`2926c71c`) is the sole
physics authority; its §3/§4/§5 define what follows. Preserve the primary Eulerian T-0…T-i derivations
(except T-h below), the `WL_S11CA_<QUANTITY>` tag grammar, and the emit machinery. Add no expected value or
acceptance criterion (`CLAUDE.md` rule 5): a residual is a measurement you compute and emit, ⛔ never a
target you shape a route to hit.

Two independent review legs derived from the spec and located these defects; their statements are quoted
below. Fix the root cause of each; ⛔ never silence a symptom (no `Off`/`Quiet`/`Check` to hide a message,
no deletion of the object that reveals a defect).

## Defect 1 (BLOCKING, physics wrong) — T-h `Σ` drops its thickness and `(1+θ)` via a multi-line parse

> T-h uses `Σ = ρ_4D` instead of the areal density `Σ_E = ρ_4D·W·(1+θ)`. `sigmaEulerianSource` (src
> **481–484**) reads `:= branchedRho4Source[b,d]` newline `(1 + waveOrder θ)` newline
> `graphThicknessSource[b];` — line 482 is already a complete expression, so the newline terminates the `:=`
> RHS and lines 483–484 evaluate as discarded no-op statements. The running `Definition[]` is
> `sigmaEulerianSource := branchedRho4Source[b,d]`. Consequence: `Σ⁰` loses its `W_bg` weighting, `∇Σ⁰ = 0`
> (should be `(rhoBr/W0)∇W_bg`), and the storage term `∂_tΣ ≡ 0` — a mass-conservation law with zero
> storage. T-0 `BACKGROUND_DENSITY_MAP` computes `Σ_E⁰` correctly at line 645 (`Together[ρ_4D·widthAnsatz]`)
> — the same file emits the correct `Σ_E⁰` in T-0 and the wrong `Σ` in T-h.

**Fix.** Make the product one parsed expression (explicit `*` / a trailing operator on each continued line /
wrap the RHS in parentheses / `Times[...]`), so `sigmaEulerianSource` is `ρ_4D·(1+θ)·W` as intended.
⛔⛔ **This is a class of bug, not one site.** Audit **every** multi-line definition in the file whose `:=`
(or `=`) RHS begins with a complete subexpression followed by further implicit-product factors on later
lines, and correct each so it parses as the intended product. Report each site you changed, and confirm by
`Definition[...]` (or an equivalent print of the stored value) that each now carries all its factors.

## Defect 2 (export degenerate) — T-c′ emits a bare `Equal[0,0]`, not two operands

> `WL_S11CA_KINEMATIC_BALANCE` builds `kinematic = n̂·v_bulk − V_s − flux/ρ_m` with `flux` and `V_s` from the
> same objects, so it is algebraically ≡ 0; the emitted `SHAPE_DERIVATIVE` is `0` and `TEST_OBJECT` is
> `Inactive[Equal][0,0]` on all 16 keys. It is the spec's own identity but emitted as a single zero-form
> rather than two independently-computed operands.

**Fix.** Emit T-c′ as **operand A** = the shape derivative of `n̂_s·v_bulk,s`, **operand B** = the shape
derivative of `V_s + J_s/ρ_m` (with `J_s` the single §3b relative-flux object, computed from its own
definition), and their residual A−B — the operand-A/operand-B/residual shape §6 requires of every
comparison. ⛔ Do not emit a bare `Equal[0,0]` and ⛔ do not substitute a literal for either operand: compute
each side from its own definition and print it. Whatever A−B evaluates to is a measurement for the step
record, ⛔ not something the script asserts or is shaped to produce.

## Defect 3 (§5a control) — route 2 is not DOF-keyed

> Route 1's `heightSource` (179–188) carries the DOF bookkeepers through `zetaSource` (129–131:
> `dofCenter·zetaCenter + sign·dofWidth·deltaWidth/2`). Route 2's `flatteningCoordinateSource` (204–213) uses
> the **raw** fields `zetaCenter` and `deltaWidth`, not `dofCenter`/`dofWidth`. So `applyPhysicalDof["DELTA_W"]`
> zeros route 1's `ζ_c` but leaves route 2's `ζ_c`; the DELTA_W representation residual of the normal is not
> a same-DOF comparison (`EULER_N_HAS_ZETACENTER: False`, `MAT_N_HAS_ZETACENTER: True`).

**Fix.** Route 2 must carry the same DOF bookkeepers as route 1, so `applyPhysicalDof` selects the same DOF
in both routes and the representation residual is a same-DOF comparison. Build route 2's center/width from
`dofCenter`/`dofWidth` (as route 1 does), not the raw fields.

## Defect 4 (§5a control) — T-g route 2 is a duplicate of route 1, not a genuine second construction

> `virtualConstraintMaterialSource` (526–549) builds `rhoFour·(1+θ)·denominator·jacobian` and
> `virtualConstraintDirectSource` (490–524) builds `rhoFour·(1+θ)·thickness·jacobian`; for LAB_HELD
> uncorrupted `thickness` and `denominator` are the same expression, so the two are byte-identical
> (`PROBE_PRE_SAMEQ_TG: True`). The material route does not perform the §5a flattening — it calls neither the
> flattening coordinate nor a `Solve` of `w' = s/2` — it is the same formula under a different name.

**Fix.** Per spec §5a, route 2 for T-g must reach the virtual material-constraint object by the genuine
**material-coordinate `w'` face-flattening** construction (the exact `w' = [w−ζ_c]/[W_bg+δW]` flattening,
its Jacobian, mapped back to Eulerian), a construction that does not reuse route 1's `thickness`/mass
product. Whether the two routes then coincide is the measurement `REP_INVARIANCE` records — if the material
mass ratio is representation-invariant they may agree, and that agreement is then a real identity rather than
a duplicate; emit operand A, operand B, and A−B and let it be what it is. ⛔ Do not make route 2 a renamed
copy of route 1.

## Defect 5 (§5a control) — the one-sided independence control mutates a result, not a source

> `WL_S11CA_CONTROL_INDEPENDENCE_*` (1106–1109) sets `corrupted = applyProfileJets[applyPhysicalDof[
> directShape[...], dof], 0, sign==1]`, and `directShape` (419) is `shapeDerivative[...]` — the
> already-computed shape derivative. `applyProfileJets` with the reversal flag (164–165) does
> `result /. w1Jet[1,0,0] -> -w1Jet[1,0,0]`, flipping a jet in the emitted operand. Spec §5a says "mutate
> only one route **at its source**" and "does not authorize editing an already-emitted normal, traction, …".

**Fix.** Per spec §5a, the one-sided independence control must reverse the upper-face `x¹` first jet of
`W_bg` **only in the direct-route source** (the `F_+`/`R_+` face source that route 1 differentiates),
leaving the material-coordinate route and every other source unchanged, and then recompute the object
downstream from that corrupted source. ⛔ Do not apply the reversal to the already-computed
`directShape`/normal/traction; corrupt the source and let the shape derivative recompute.

## Preserve

Everything else is unchanged: the primary Eulerian T-0…T-i derivations and payloads (both legs confirmed
these correct and form-ablation clean, T-h excepted); the orientation law; the §3c shifted trace; T-g's use
of the virtual `δ_vu` and both density representatives; blindness (no import); the `WL_S11CA_<QUANTITY>`
grammar and one-tag-per-object; the emit machinery (`Quit[90/91]` guard tag-name uniqueness only). Re-run to
regenerate the tag stream.

## Confirm

Re-run the finished engine (`wolframscript -file …`, one kernel, `timeout` generously — a full run is
~15 min, ~5 GB RSS; kill a demo at 600 s no-output or 6 GB and report it; runs under scratch, ⛔ never
`mathematica/out/`; `--sandbox danger-full-access`). Verify: (1) stderr carries only the benign
`WolframScript.conf` line — no `Part::pkspec1`, no `AssociateTo::argrx`, no other message; (2) all 40
`WL_S11CA_<QUANTITY>` tags still emit; (3) `sigmaEulerianSource` (and any sibling multi-line site) now carries
all its factors; (4) T-c′ emits two operands + residual; (5) route 2 is DOF-keyed and T-g's route 2 is the
flattening construction; (6) `CONTROL_INDEPENDENCE` corrupts a source. Report the before/after stderr and the
list of changed sites. ⛔ Report any ambiguity or non-computable object in the §8 builder report, not as an
invented tag or a suppressed message.
