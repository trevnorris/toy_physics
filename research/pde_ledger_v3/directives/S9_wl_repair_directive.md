# S9 REPAIR — make the dimension block live, and close four smaller holes

You are repairing an existing Mathematica engine. ⛔ **Do not rewrite it.** Make the changes below and
leave everything else alone.

**File (absolute, edit in place):**
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`

Verify by running it:
`cd /var/projects/toy_physics && timeout 600 math -script research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`

⚠ Mathematica has a **two-seat licence** here and seats may be busy. Run **one** kernel at a time; on a
licence error, wait and retry. Iterate until exit 0.

**The three standing clauses still govern every line you touch:** emit **computed objects**, never prose
and never a hand-typed algebraic answer; emit **operands and residual before** any guard; **interpretation
belongs to the step record**, so ⛔ no verdict tag, no `PASS`/`FAIL`, no conclusion in a tag *name*.

⭐ **The structural rule is unchanged and it is what R1 below is about:** the physical symbols may be
combined by hand **only** in constructing the actions and the plane-wave ansatz. Everything else must be
**reached by computation**, and every control must re-enter the chain **at the action**.

---

## R1 — ⭐⭐⭐ THE MAIN REPAIR: the dimension block must be DERIVED FROM THE ACTION

**The defect, measured by ablation.** When the stiffness form in the action was changed from
`μ_R (∇×u)·(∇×u)` to `μ_R (∇·u)²`, the **entire** dimension block emitted **byte-identical** output.
The cause is that the derivative structure is **hand-encoded**: lines of the form

```
timeDerivativeDimension = displacementDimension - timeDimension
curlDimension           = displacementDimension - lengthDimension
```

assert *"the inertia term carries one time derivative of `u`, the stiffness term carries one space
derivative of `u`"*. Nothing reads that off the Lagrangian.

⚠ **For this action the resulting answer is correct** — curl and divergence both carry exactly one space
derivative. ⛔ **But it is invariant because nothing examined it, not because anything checked it.** A
stiffness with a *different* derivative count would give different dimensions and this block would not
notice.

### What to build instead

**Read the derivative structure out of the constructed Lagrangian.** For each additive term of the
Lagrangian, and for each field factor in that term, extract the derivative multi-order by pattern
matching on `Derivative[dt, dx, dy, dz][u_i][t, x, y, z]`. Then the dimension of a term is

```
[coefficient]  +  Σ over field factors of ( [u] − dt·[T] − (dx + dy + dz)·[L] )
```

Require every term of the Lagrangian to equal the energy density on a `D`-dimensional sheet, and
**solve** the resulting linear system for the unknown coefficient dimensions exactly as now.

⛔ Do not hand-write the derivative counts. ⛔ Do not special-case curl or divergence by name. ⭐ The
counts must come from the expression tree, so that changing the action changes them.

Emit, as computed objects:
- the extracted per-term derivative multi-orders,
- the per-term dimension expressions built from them,
- the linear system, and its solution, exactly as now.

### R1b — a control that PROVES the block is now live

Add a control **X7** whose stiffness has a **different derivative count** — flexural stiffness:

```
L = (1/2) ρ_br (∂_t u)·(∂_t u)  −  (1/2) μ_F (∇²u)·(∇²u)
```

with `μ_F` a fresh positive symbol. Run it through the **same** chain as every other control, and emit
its full control package **including its own dimension block**.

⭐ **X7 exists to exercise the dimension block, not the spectrum** — it is the one control whose
stiffness carries a derivative count unlike the others, so it is the case that puts R1 under load.
⛔ Emit `[μ_F]` and the whole block as computed, and ⛔ **do not compare it against anything, in the
script or in your report.** Whether the numbers come out related is not yours to decide and not
something to tune toward — the reading happens elsewhere.

---

## R2 — impose the wavevector assumption instead of asserting it in a comment

`k·k ≠ 0` appears only as a comment; it is imposed nowhere, and the previous version of this file carried
positivity and reality assumptions that the rebuild dropped. `MatrixRank` on a symbolic matrix returns
the **generic** rank, so every rank-derived result is currently a generic-direction statement that the
output does not mark as such.

- Introduce explicit assumptions (`ρ_br > 0`, `μ_R > 0`, `μ_F > 0`, `ρ_z > 0`, `λ_ρ > 0`, `λ_μ > 0`,
  `k·k > 0`, wavevector components real) and thread them through the simplifications.
- ⭐ **Emit the assumption set itself as a CAS object**, so a reader can see what the results are
  conditioned on.

## R3 — emit the polarisation results at SPECIAL DIRECTIONS, not only generic ones

`E2` is the transverse mode count and a later step's entire result rests on it, so a generic-direction
number that the output does not distinguish from a domain-wide one is a real hazard.

For the **main** action and for **X6** (anisotropic inertia), re-run the per-root **E1–E4** with the
wavevector specialised, and emit each as its own tag:

```
k = (kx, ky, kz)   generic          k = (kx, ky, 0)
k = (0, 0, kz)                      k = (kx, 0, 0)          k = (kx, kx, kx)
```

⛔ Emit the computed integers and booleans for every case. ⛔ Do not summarise them, do not state which
directions are special, and ⛔ do not add a tag saying whether they agree — emit the numbers.

## R4 — delete the structurally-forced zero

Delete `WL_S9_X4_ROOT_DIFFERENCE` and the two operand tags feeding it
(`WL_S9_X4_TRANSVERSE_ROOT_OPERAND`, `WL_S9_X4_E3_ROOT_OPERAND`).

⛔ **This comparison was forbidden by the build directive and should not have been emitted.** Under X4
both subsets contain the *same single root*, so the difference is `0` **by construction**, not by
physics. It also calls `First[...]` on a subset that is empty in other controls, which would raise.
⭐ The degeneracy is already carried honestly by `X4_ROOT_MULTISET` and the per-root E1/E3/E4 tests.

## R5 — guard the even-polynomial extraction

`evenPolynomial` sums only the even coefficients in `ω` and **silently discards** any odd power. All
current determinants are even, so nothing is being lost — but the discard is unguarded.

⭐ Emit the **discarded odd part** as its own CAS object, for the main action and every control, and then
apply a hard stop if it is nonzero. ⭐ Emit first, **then** stop.

## R6 — stop emitting operand/residual triples for things that are not comparisons

Two triples are cosmetic and should be reduced to the single object each actually carries:

- `WL_S9_EQUATION_REFERENCE` is the literal zero vector, so `EQUATION_DIFFERENCE` is just the equation of
  motion with extra steps.
- `WL_S9_PLANE_WAVE_ANSATZ_RESIDUAL` differences a field symbol against the ansatz that replaces it,
  which is a definition, not an agreement between two routes.

⭐ **Emit an operand/residual triple only where two genuinely INDEPENDENT routes produce objects that
should agree.** Elsewhere emit the object alone. ⚠ A residual whose operands share a route manufactures
the appearance of a check, which is the pattern this rebuild exists to remove.

⭐ **Keep** the `E1_LEFT`/`E1_RIGHT`/`E1_RESIDUAL` triple — there the operands *are* independent (a
computed matrix rank against a structural dimension).

---

## Report back — under 30 lines

1. Exit status and runtime of your own run; total tag count and confirmation all tags are unique.
2. **For R1: the literal `[μ_F]` dimension vector X7 produced, and the literal `[μ_R]` vector**, quoted
   exactly as emitted. ⛔ Do not interpret them, relate them, or comment on them.
3. For R1, the line numbers where derivative orders are now extracted from the Lagrangian.
4. One line per repair R2–R6: what changed, and the tag names added or removed.
5. Anything you could not do, and what blocked it.
6. ⭐ Anything the run produced that surprised you.

⛔ Do not commit to git. ⛔ Do not modify any other file.
