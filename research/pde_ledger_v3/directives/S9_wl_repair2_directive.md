# S9 REPAIR 2 — six fixes, one of them a wrong value

⛔ **Repair in place. Do not rewrite the script.** Make these changes and leave everything else alone.

**File (absolute):**
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`

Verify by running it:
`cd /var/projects/toy_physics && timeout 600 math -script research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`

⚠ Mathematica has a **two-seat licence**; run **one** kernel at a time and retry on a licence error.
Iterate to exit 0.

**The three standing clauses still govern:** emit **computed objects**, never prose and never a hand-typed
algebraic answer; emit **operands and residual before** any guard; ⛔ no verdict tag, no conclusion in a
tag name. ⭐ **The structural rule holds:** the physical symbols may be combined by hand **only** in the
actions and the ansatz; everything else is reached by computation.

---

## F1 — ⛔⛔ A WRONG VALUE: the implied speed dimension drops the wavevector

`WL_S9_SPEED_SQUARED_IMPLIED_DIMENSION` is currently built by reading the **powers of `rhoBr` and `muR`**
out of the computed speed expression and combining their solved dimension vectors. ⛔ **That silently
discards every other factor**, including explicit wavevector dependence.

⚠ **Measured, and it produces a wrong number today.** For the flexural control the candidate speed
squared carries an **explicit factor of `q`**. The current method reads only the exponents of `rhoBr` and
`muR`, so that factor is dropped and the reported dimension is short by the dimension of `q` — and the
reported difference against squared velocity is correspondingly wrong. ⭐ The main action escapes only
because its speed expression happens to carry no `q`. ⛔ **Verified independently in a second CAS.**

⚠ **It is also flattering, which is worse than merely wrong:** the bogus nonzero difference *looks* like
the script correctly detecting that flexural is dispersive. It is not — it is arithmetic that drops a
factor. ⭐ The tag that genuinely detects dispersion is the homogeneity defect, which is already correct.

### The fix
⭐ **Dimension the WHOLE expression by walking its tree**, exactly as the per-term dimensional analysis in
the Lagrangian block already does — ⛔ do not read parameter powers.

- Build one symbol → dimension-vector map containing: the **solved** coefficient dimensions this script
  derived (`rhoBr`, `muR`, `muF`, `rhoZ`, …), the **definitional** primitives (`kx`,`ky`,`kz` an inverse
  length; `q` an inverse length squared; `omega` an inverse time), and the **dimensionless** symbols
  (`lambdaRho`, `lambdaMu`, `D`).
- Walk the expression: a sum requires **every summand to carry the same dimension** (⛔ otherwise emit the
  summands with their differing dimensions and hard-stop **after** emitting); a product adds dimensions;
  a power with a rational exponent scales them.
- ⛔ **An unknown symbol is a hard stop, named** — ⛔ never assumed dimensionless.

Emit the resulting dimension for **every** candidate speed of **every** action, alongside the squared
velocity reference and their difference, as now. ⭐ Keep the tag names.

---

## F2 — ⛔⛔ A LATENT WRONG VALUE: a term with no derivatives loses its field dimension

The per-term dimensional analysis sums over the `Derivative[…][field][…]` atoms in a term. ⛔ **A term
containing a bare field with no derivative contributes nothing**, so its `[u]` factors are dropped.

⚠ **Measured:** adding a bare-field term to an action makes the script emit, for that term's coefficient,
a dimension **short by exactly the field factors the term contains** — because the term's field content
was never counted. ⛔ **And it does so with exit 0, the full unique tag count, and every gate green**, so
nothing catches it.

⚠ **No S9 action has a bare-field term, so nothing emitted today is wrong.** ⭐ But this block is the
obvious thing to carry into later steps, and a later step **will** have a gap or restoring term.

### The fix
⭐ **Count FIELD FACTORS, not derivative atoms.** For each additive term, enumerate every factor that is a
field — whether it appears bare or under a `Derivative` — and let a bare field carry the multi-order
`(0,0,0,0)`. Each field factor then contributes
```
[u]  −  dt·[T]  −  (dx + dy + dz)·[L]
```
so a bare field contributes exactly `[u]`. ⛔ Do not special-case any operator by name.

⭐ **Prove the fix with a control.** Add **X8**, identical to the main action but with an added bare-field
term `− (1/2) muG fieldVector·fieldVector` and `muG` a fresh positive symbol. Run it through the same
chain as every other control, emitting its **full package including its own dimension block**.
⛔ Do not state anywhere what `[muG]` should be, and ⛔ do not compare it to anything — emit it.

---

## F3 — the assumptions are inert, and the ASSUMPTIONS tag reports a different list than the one imposed

Two separate problems, both measured:

- ⚠ Replacing the whole assumption list with `{True}` leaves the output **byte-identical** — ⛔ **no
  emitted value depends on the declared domain.**
- ⛔⛔ `WL_S9_ASSUMPTIONS` is built from a **different object** than the list actually passed to
  `FullSimplify` / `Assuming`. Gutting the imposed set left the tag unchanged. ⇒ **a reader cannot learn
  the domain the values were computed on**, and the two expressions can drift apart silently.

### The fix
1. ⭐ **One source only.** Emit the *same* object that is passed to the CAS. ⛔ Delete the parallel
   hand-built assumption object entirely; ⛔ do not keep both in sync.
2. ⭐ **Make the domain restriction VISIBLE AS DATA rather than as an inert declaration.** The exclusion
   `k·k ≠ 0` is real physics, so **compute what happens on the excluded locus**: for the main action,
   additionally emit **E1, E2, E3, E4 at `k = (0,0,0)`**, as their own tags.
   ⛔ Do not comment on the result, ⛔ do not compare it to the `k ≠ 0` case, and ⛔ do not add a check
   that it differs. ⭐ Emit the numbers; the difference is then a fact in the output rather than an
   assumption in a comment.

## F4 — a dimension-block failure escapes both the exit code and the tag count

When the dimension system is inconsistent, `Solve` returns `{}` and the script emits unevaluated residue
(`First[{}]`, `… /. First[{}]`, `Select[0, False&]`) **while exiting 0 with the full unique tag count.**
It is caught only because Mathematica writes `Message` lines to stdout.

⭐ **Fix:** after solving the dimension system, emit the solution set **as computed** — an empty list
emitted as an empty list — and then, **if it is empty, hard-stop with a nonzero exit.** ⛔ Never emit an
expression containing `First[{}]` or an unapplied replacement rule. ⭐ Emit first, then stop.

## F5 — eight duplicate payloads under distinct names

`WL_S9_*_FULL_ROOT_MULTISET` emits the identical object as `WL_S9_*_ROOT_MULTISET` in all eight packages,
so the tag count overstates distinct content by 8. ⭐ Delete the `FULL_ROOT_MULTISET` tags.

## F6 — the anisotropic control's mode count is generic in the PARAMETER domain, unmarked

`X6`'s transverse multiplicity is generic in `rhoZ`; on the sub-locus `rhoZ = rhoBr` — which lies inside
the declared positive domain — it takes a different value. The five specialised **directions** sample the
wavevector domain; ⛔ **nothing samples the parameter domain.**

⭐ **Fix:** for `X6`, additionally emit the per-root **E1–E4 with `rhoZ → rhoBr` substituted**, as their
own tags. ⛔ Emit the computed values; ⛔ do not remark on them.

---

## Report back — under 30 lines

1. Exit status, runtime, total tag count, confirmation all tags unique and no untagged output.
2. **F1:** the literal implied-dimension and difference tags for the **main** action and for the
   **flexural** control, quoted exactly. ⛔ Do not interpret them.
3. **F2:** the literal `[muG]` dimension vector X8 produced. ⛔ Do not interpret or compare it.
4. **F3:** the literal E1–E4 values at `k = (0,0,0)`. ⛔ Do not interpret them.
5. **F6:** the literal X6 E1–E4 values at `rhoZ = rhoBr`. ⛔ Do not interpret them.
6. One line each for F4 and F5: what changed, tags added or removed.
7. Anything you could not do, and what blocked it.
8. ⭐ Anything that surprised you. This is wanted.

⛔ Do not commit to git. ⛔ Do not modify any file other than the deliverable.
