# S11c-a SymPy engine — §5a route-2 repair directive

## Authority and boundary

Edit `research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py` **in place** (committed
baseline `488c2a65`). Re-run it to regenerate `research/pde_ledger_v3/scripts/S11c_a_exports.py`. Those two
products are the only writes; every other file, including `S11b_exports.py`, is read-only.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (HEAD `2926c71c`) is the sole
physics authority — §5a defines route-2 (the material-coordinate face-flattening derivation) and its
one-sided independence control; §4 T-c′ defines the kinematic-balance object; §6 states the three script
clauses. The build directive `research/pde_ledger_v3/directives/S11c_a_sympy_build_directive.md` (`304fa46f`)
governs all wiring (imports, binds, F9, digests, publish gate, D3, `_RELATIONALS`) unchanged.

**This is a surgical repair of two measured defects, nothing more.** ⛔ Do not touch the primary derivations
T-0…T-i, the export schema, the F9 collision logic, the digests, the D3 round-trip, or the T-h
(`evolution_route`) construction — two review legs confirmed the primary payloads are correctly computed and
form-ablation clean, and that T-h's two routes are already genuinely independent. ⛔ Add no expected value or
acceptance criterion (`CLAUDE.md` rule 5): you compute and emit operand A, operand B and A−B; the residual
is a measurement, **never a target**, and you must not shape any route so that a residual comes out a
particular way.

## The two measured defects — verbatim from the review legs

Two independent legs (a fresh Claude agent and Grok) each derived the objects from the spec with their own
runnable scripts and converged. Their findings, quoted:

### Defect 1 (BLOCKING) — §5a route-2 is a byte-identical alias of route-1 for T-c, T-d, T-i, T-g

> The engine's `route="MATERIAL"` branch does **not** implement the face-flattening `w′` construction. In
> `build_face_source` it computes `tangent_vertical = dx(h_exact,i,scales)` — the *same* `grad_h` as the
> Eulerian branch — and forms `cofactor=(-face*grad_h, face)`, whose normalization is algebraically
> **identical** to the Eulerian normal. The code comment claims "an independent normal construction …
> mapped back to Eulerian jets," but it is the same expression. `virtual_constraint_route` is worse: the
> MATERIAL return `density*thickness*jacobian/background_mass` is character-for-character the same formula
> as the EULERIAN return, with identical `W_anchor`.  … `w′` from §5a is never introduced. … the
> spec-mandated face-flattening `w′` route 2 is not implemented. Reimplement route 2 as the spec's
> material-coordinate face-flattening construction before this control is trusted or committed.

> `reverse_upper_x1` is applied to `scales` *before* the route split. … Shared reverse moved **both** routes
> equally … A residual that stays zero while both operands move is not a two-route test. … corrupting
> **only** the Eulerian source moves EUL but not MAT; corrupting the shared build moves both operands
> identically → residual stays 0. I would not freeze this as having two independent routes until MATERIAL
> actually differentiates on `w′` (and the one-sided control mutates only that direct-route source).

Both legs agree the **only** genuine second route is T-h (`evolution_route`, EULERIAN dilatation+advection
vs MATERIAL material-mass time-derivative through the Jacobian + evaluation-shift). Leave T-h as it is.

### Defect 2 (secondary) — T-c′ `KINEMATIC_BALANCE` is a bare zero that audits its own input

> `kinematic_raw` … Operand B is `V + J/ρ_m` with `J` from `relative_flux_raw` on the same `FaceSource`. …
> The emitted equation is `Eq(<unsimplified 0>, 0)`. That is operand B derived from operand A. It cannot
> fail. … T-c′ `KINEMATIC_BALANCE` is a definitional identity (`n̂·v_bulk = V_s + J_s/ρ_m`) built from the
> same normal/flux/velocity operands, so its residual is structurally 0 — correct physics but a check that
> cannot independently fail.

## What to build

### Fix 1 — implement §5a route-2 as a genuinely independent computation (T-c, T-d, T-i, T-g)

Re-derive the `route="MATERIAL"` operand for the relative flux (T-c), the traction and virtual work (T-d),
the closure (T-i) and the virtual material constraint (T-g) as the **material-coordinate face-flattening
derivation of §5a**: use `x=x(X,t)` and the exact flattening coordinate

```text
w′ = [w−ζ_c(X,t)]/[W_bg(x(X,t))+δW(X,t)]   (LAB_HELD),
w′ = [w−ζ_c(X,t)]/[W_bg(X)+δW(X,t)]        (MATERIAL_ADVECTED),
```

differentiate in the material coordinates, and **map route 2 back to the Eulerian variables before
comparison**, holding the anchoring branch fixed across the two routes (§5a). Route 2 must reach its object
by this construction — ⛔ it may **not** reuse route 1's `grad_h`/normal/measure objects, its
`density*thickness*jacobian` product, or any already-computed route-1 quantity; the whole point of §5a is
that the two routes are **independent derivations of the same object**, so that a source error appears in
only one of them. Whether the mapped-back operands then coincide is the measurement `S11CA_REP_INVARIANCE_*`
records — emit operand A (Eulerian), operand B (material, mapped back) and A−B, and continue (§6).

Then repair the **one-sided independence control** so it does what §5a says: "reverse only the `x¹` first jet
of `W_bg` in the upper-face **direct-route** `F_+/R_+` source, leaving the **material-coordinate route** and
every other source unchanged." The corruption must enter route 1's source alone; route 2 must be rebuilt
from its own uncorrupted source. Emit the uncorrupted-route operand, the corrupted-route operand, and their
residual (`S11CA_CONTROL_INDEPENDENCE_*`), per §5a — ⛔ do not apply the reversal to a shared jet-scale that
both routes read.

Keep the T-g/T-h source mutations for the independence control that §5a already prescribes
(`Σ_E(x(X,t)) → Σ_E(X,t)` for T-g; `∇·(Σv) → Σ∇·v` for T-h) — those enter the **direct route** only, which
is already correct for T-h and must become correct for T-g under the reimplemented route split.

### Fix 2 — T-c′ emits its two operands, not a bare zero

Emit the kinematic-balance object as operand A = the shape derivative of `n̂_s·v_bulk,s`, operand B = the
shape derivative of `V_s + J_s/ρ_m` (with `J_s` the single T-c relative-flux object, computed from its own
§3b definition as in T-c), and their residual A−B — the same operand-A / operand-B / residual pattern §6
requires of every comparison. ⛔ Do not emit a bare `Eq(residual, 0)`, and ⛔ do not substitute a literal for
either operand: compute each operand from its own definition and print it. Whatever A−B evaluates to is a
measurement for the step record, ⛔ not something the script asserts, characterises, or is shaped to produce
(§6 clauses 1–3).

## Preserve

Everything else is unchanged: T-0…T-i primary derivations and their emitted payloads; the export candidate
population and schema; F1/F9 (three-valued, TOTAL over imported row shapes) and the `s11c_a_` F9c prefix;
`BUILD_INPUT_DIGESTS` over the same three inputs; the F6 publish gate (publish only if every primary
T-0…T-i completed); the D3 in-run round-trip; the `_RELATIONALS` reviver; `rho_m → LEDGER['rho_m']` and every
other inherited bind. Re-run the engine so `S11c_a_exports.py` is regenerated from the edited source (its
self-referential digest will update).

## Script clauses (exact, from §6)

1. The script prints computed CAS objects and never states conclusions.
2. For every comparison it emits operand A, operand B and A−B **before** any guard. A physics disagreement
   emits and continues with exit status 0; nonzero exit is reserved for operational failure only.
3. Interpretation belongs to the step record, not the engine.

Emission is never conditional on a payload's value. ⛔ No `assert` precedes the `emit` of the value it
guards. Report under §8 anything you cannot build under one-tag-per-object; ⛔ do not fill a gap with new
physics, an expected result, or a self-review mechanism.
