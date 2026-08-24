# S11c-a Wolfram engine — repair round 2 (T-f projection + T-0 gradient grading)

## Authority and boundary

Edit `research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` **in place**
(committed baseline `a15bc69c`). Its product is the flushed stdout tag stream; that is its only write. It
remains the **blind** engine — imports nothing, re-derives from the spec's equations. ⛔ Add no import, and
⛔ do not introduce any object, value, or construction from another engine (its independence is the one
blindness control).

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (`2926c71c`) is the sole
physics authority; §2a (`∂_{y_i}W_bg = σ_W ∂_{ξ_i}w_1`, the first-jet map and its bookkeeping), §2b/§4 T-0
(background density map and its in-plane gradient), and §3c/§4 T-f (the dynamic window `Ω^α = 𝒪(G_+,G_-)`
and the projection identity + its shape derivative) define what follows. Add no expected value (`CLAUDE.md`
rule 5): a residual is a computed measurement, ⛔ never a target.

⛔ **This is a surgical repair of exactly two defects.** The prior round's §5a and T-h fixes are verified
good by two legs — ⛔ do not touch route-2, the one-sided source corruption, T-c′, `sigmaEulerianSource`, or
any primary T-a…T-e/T-g/T-h/T-i derivation. Fix the root cause of each defect below; ⛔ never silence a
symptom (no `Off`/`Quiet`/`Check`, no deletion of the object that reveals the defect). Re-run to regenerate
the tag stream.

## Defect 1 (physics wrong) — T-f dynamic projection produces a divergent `∫1 dw`

A review leg derived from the spec and found (quoted):

> `WL_S11CA_PROJECTION_SHAPE_DERIV` / `_DYNAMIC_OPERAND`. `projectionTermsSource` (656–685);
> `shapeDerivative`+`applyProfileJets` (933–972). `D` does penetrate `Inactive[Integrate]`. After
> `Expand`/`Series`, mixed-`η` pieces of the **dynamic** operand are rewritten as w-dependent fields sitting
> **outside** `Inactive[Integrate][1, {normalCoordinate, -∞, ∞}]`. `WL_SLICE_TF_DYNAMIC_SHAPE` contains
> `etaBg*W0*currentWPerturbation[…]*w1Jet[0,0,0]*Inactive[Integrate][1,{normalCoordinate,-Infinity,Infinity}]*Derivative[0,2][windowFunction][…]`.
> The mixed-`η` terms are an ill-defined `∫_{-∞}^{∞} 1 dw` with the `w`-dependence of the window pulled out
> of the integrand. That is a way T-f physics is wrong.

**Root cause.** The base `projectionTermsSource` correctly keeps the window and each field **inside**
`Inactive[Integrate][window·field, {w,-∞,∞}]` (the window localises → the integral converges). The **shape
derivative / `Expand` / `Series`** step then distributes across the held integral and factors a
`w`-dependent window derivative **out**, leaving `Inactive[Integrate][1, {w,-∞,∞}]` (a divergent `∫1 dw`)
multiplied by `windowFunction''` evaluated outside.

**Fix.** The shape derivative of `∫ window·field dw` must remain `∫ ∂(window·field) dw` — every window
factor and its derivatives stay **inside** the `Inactive[Integrate]`, as functions of `w`, so the integrand
still localises and no `Inactive[Integrate][1, …]` (constant integrand over the whole line) ever appears.
⚠ **Two transformation phases operate on the projection object and BOTH must be made integrand-preserving,
not only the first:** (a) `shapeDerivative` takes the `D`/`Expand` step, and (b) `applyProfileJets` then
runs the `etaBg`/`sigmaW` `Series` truncation. Mapping the `D` into the integrand
(`D[Inactive[Integrate][g[w,η], {w,…}], η] → Inactive[Integrate][D[g,η], {w,…}]`) is necessary but **not
sufficient** — the subsequent `Series`/`Expand` on the truncation step can still factor a `w`-dependent
window derivative out of a held integral. Ensure **every** phase that touches these projection integrals
keeps all window factors inside the integrand (e.g. differentiate/`Series` the integrand and re-wrap, or
otherwise prevent distribution across `Inactive[Integrate]`). Recompute the dynamic operand and, **as an
independently computed object (spec §3c/§4 require two operands, not one derived from the other)**, the
static-flat operand, both under the corrected integrand-preserving transformation; then recompute their
residual.

## Defect 2 (grading bookkeeping) — T-0 `RHO4_CONSTANT` gradient graded η, not σ_W

Quoted:

> Engine lines 709–710, emit 720: `gradient = (D[sigmaZero, #] & /@ {y1, y2, y3}) /. etaBg W0/LW -> sigmaW`.
> For `RHO4_CONSTANT`, `Together[ρ_4D W_bg] = rhoBr (1 + η w1)`, so `∇Σ` is `(η ρBr / L_W) ∂w1` and **does
> not contain the product `η W_0 / L_W`** — the replacement is a no-op. The emitted CAS object is graded as
> an `η` zero-jet, not a `σ_W` first-jet. Spec §2a forbids collapsing those bookkeepers.
> `RHOBR_CONSTANT` gradient is the zero vector (correct).

**Root cause.** The gradient is graded through a fragile literal pattern substitution `etaBg W0/LW -> sigmaW`
that only fires when the exact product `etaBg·W0/LW` appears. For `RHO4_CONSTANT` the derivative carries
`etaBg·rhoBr/LW`, so the pattern never matches and `∇Σ` keeps its `η` grade instead of the spec's `σ_W`
first-jet grade.

**Fix.** Grade `∇Σ_E^0` by the spec §2a first-jet map `∂_{y_i}W_bg = σ_W ∂_{ξ_i}w_1` itself, not by a
literal-product rewrite that assumes a particular coefficient. ⚠ **"Retain `W_bg`" is not sufficient on its
own:** `W_bg` must remain an **unexpanded, held factor through the `D`** — a held object, ⛔ not the
substituted product `W_0(1+η w_1)`. `Together` is not the only cancellation that defeats the map:
`Times[(ρ_br/W_0), W_0(1+η w_1)]` **already cancels `W_0`** before differentiation, so `D` never produces a
`∂_{y_i}W_bg` unit and `etaBg W0/LW -> sigmaW` still no-ops. Keep `W_bg` inert through the derivative;
**after `D`, replace each generated `∂_{y_i}W_bg` factor by the §2a RHS `σ_W ∂_{ξ_i}w_1`**, so every such
factor is expressed through `σ_W` and carries the `{ε,η,σ_W}` multigrade the §2a map assigns — independent
of the density representative's prefactor — for **both** representatives. ⛔ Do **not** globally rewrite `η`
in terms of `σ_W` (that collapses the two bookkeepers §2a keeps independent). ⛔ Do not hard-code the answer
and ⛔ do not assume a target grade; derive the gradient, let its grade fall out of the §2a map, and
emit/report the resulting computed expression and its multigrade for each representative. (This directive
supplies no expected gradient value or grade; `CLAUDE.md` rule 5.)

## Preserve

Everything else is unchanged and verified good: route-2 (`flatteningCoordinateSource`,
`parametricNormalSource`, the genuine `w'` flattening), the one-sided source corruption
(`CONTROL_INDEPENDENCE`), T-c′ two operands, `sigmaEulerianSource` and all six repaired multiline sites, and
every primary T-a…T-e/T-g/T-h/T-i derivation. Blindness (no import), the `WL_S11CA_<QUANTITY>` grammar,
one-tag-per-object, and the emit machinery (`Quit[90/91]` guard names only) stay intact.

## Confirm

Re-run the finished engine (`wolframscript -file …`, one kernel, generous `timeout` — a full run is ~20 min
at ~9–10 GB RSS; ⛔ watch for an orphaned kernel if the job dies and `kill -9` it by PID; runs under scratch,
never `mathematica/out/`; `--sandbox danger-full-access`). Verify: (1) stderr carries only the benign
`WolframScript.conf` line; (2) all 40 `WL_S11CA_<QUANTITY>` tags still emit; (3) **every** full-line
projection integral — in `WL_S11CA_PROJECTION_*` **and** in every T-f consumer that carries a projection
integral (the projection entries inside `CONTROL_FORM_*` and `UNIFORM_LIMIT_*`) — retains a localizing
`windowFunction` (or a derivative of it) **inside** its integrand, and satisfies **both**: (i) ⛔ **no**
factor that depends on `normalCoordinate` sits **outside** a held `∫ dw` — this is stronger than "no
`Inactive[Integrate][1, {normalCoordinate, …}]`": a term like `field[w]·Inactive[Integrate][window''[w],
{normalCoordinate, …}]` with the field left **outside** must also be rejected; and (ii) the **dynamic**
operand still depends on `etaBg`, with those `etaBg` factors sitting **inside** the integrand (a fix that
drops the mixed-`η` terms removes the `∫1` but destroys the window's `η` dependence, which §2a forbids); (4) the `RHO4_CONSTANT` and `RHOBR_CONSTANT` `GRADIENT_SIGMA_E_ZERO` payloads are computed via the §2a
map — **report the computed expression and its multigrade for each** (the orchestrator verifies the grade;
this directive supplies no target grade); (5) the §5a/T-h objects are **byte-unchanged** from baseline
`a15bc69c` — byte-compare all fenced payloads (REP_INVARIANCE_*, CONTROL_INDEPENDENCE_*, KINEMATIC_BALANCE,
EVOLUTION_MASS_BALANCE, route-2/T-g/T-c′ objects), not only a spot sample. Report before/after and the
changed sites.
⛔ Report any ambiguity or non-computable object in the §8 builder report, not as an invented tag.
