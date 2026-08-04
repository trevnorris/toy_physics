# S9 SymPy engine — conformance repair

⛔ **Repair in place. Do not rewrite the script.** Apply the requirements below and change nothing else.

**File (absolute):**
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py`

Verify by running it:
`cd /var/projects/toy_physics/research/pde_ledger_v3 && timeout 600 python3 scripts/S9_light_requires_shear_sympy_audit.py`

Exit 0, well under 600 s. Iterate until it does.

⛔⛔ **DO NOT READ, RUN OR OPEN ANY `.wl` FILE, any other engine under `scripts/`, or the reduction
registry.** ⭐ This engine's entire value is that it was built independently of the Mathematica engine and
can therefore **disagree** with it. A transcription agrees vacuously and is worse than nothing.

**The three standing clauses still govern:** emit **computed objects**, never prose and never a hand-typed
algebraic answer; emit **operands and residual before** any guard; ⛔ no verdict tag, no aggregate
summary, no conclusion in a tag *name*. ⭐ **The structural rule holds:** the physical symbols may be
combined by hand **only** in constructing the actions and the ansatz; everything else must be reached by
computation, and every control re-enters at its **action**.

---

## C1 — dimension EVERY emitted expression by walking its tree

⚠ A dimension obtained by reading the **exponents of the coefficient symbols** out of an expression is
wrong whenever the expression carries any other dimensionful factor — a wavevector, for instance. It
drops that factor silently and reports a dimension that is short by it.

⭐ **Requirement:** wherever this script reports the dimension implied by a **computed** expression, obtain
it by walking the expression tree against a symbol → dimension map, ⛔ **never** by reading coefficient
exponents.

- The map holds: the **solved** coefficient dimensions this script derived; the **definitional**
  primitives (a wavevector component is an inverse length, `q` an inverse length squared, `ω` an inverse
  time, a coordinate and a displacement a length); and the dimensionless symbols.
- A **sum** requires every summand to carry the same dimension. ⛔ If they differ, emit the summands with
  their dimensions and hard-stop **after** emitting. A **product** adds dimensions; a **power** with a
  rational exponent scales them.
- ⛔ **An unknown symbol is a hard stop, named** — ⛔ never assumed dimensionless.

## C2 — count FIELD FACTORS, not derivative atoms

⚠ Per-term dimensional analysis that sums over `Derivative` nodes gives **no contribution at all** for a
term containing a bare, underived field, so that term's field factors vanish from its dimension — and it
does so silently, with a clean exit and a full tag count.

⭐ **Requirement:** for each additive term, enumerate every **field factor**, whether bare or under a
derivative, treating a bare field as multi-order `(0,0,0,0)`. Each contributes
```
[u]  −  dt·[T]  −  (dx + dy + dz)·[L]
```
so a bare field contributes exactly `[u]`. ⛔ Do not special-case any operator by name.

⭐ **Load-test it with a new control X8**: the main action plus a bare-field term
`− (1/2) mu_G u·u`, with `mu_G` a fresh positive symbol. Run X8 through the same chain as every other
control and emit its **full package, including its own dimension block**.
⛔ Do not state anywhere what `[mu_G]` should be, and ⛔ do not compare it to anything — emit it.

## C3 — one assumption object, and put the excluded locus in the output

- ⭐ **Emit exactly the assumption object that is used in the computation.** ⛔ Do not maintain a second,
  parallel object for display — the advertised domain and the imposed domain must be the same expression,
  because otherwise they drift and a reader cannot learn what the values were computed on.
- ⭐ **Make the domain restriction visible as DATA.** The polarisation tests exclude `k·k = 0`. For the
  main action, **additionally emit E1, E2, E3, E4 evaluated at `k = (0,0,0)`**, each as its own tag.
  ⛔ Do not comment on the result, ⛔ do not compare it to the `k ≠ 0` case, and ⛔ do not add a check that
  it differs. ⭐ Emit the numbers.

## C4 — a failed dimension solve must not exit 0

If the dimension system has no solution, emit the solution set **as computed** — an empty result emitted
as an empty result — and then **hard-stop with a nonzero exit**. ⛔ Never emit an expression containing an
unapplied substitution or an indexing of an empty result. ⭐ Emit first, then stop.

## C5 — sample the PARAMETER domain, not only the wavevector domain

The anisotropic control's transverse multiplicity is generic in its independent inertia symbol; on the
sub-locus where that symbol equals `ρ_br` — inside the declared positive domain — it can differ. The
specialised directions sample the **wavevector** domain and nothing samples the **parameter** domain.

⭐ For the anisotropic control, additionally emit the per-root **E1–E4 with `ρ_z → ρ_br` substituted**, as
their own tags. ⛔ Emit the computed values; ⛔ do not remark on them.

## C6 — no duplicate payloads

⛔ Do not emit the same computed object under two different tag names. If two tags currently carry an
identical payload, delete one. ⚠ Duplicate payloads inflate the tag count with content that is not
distinct, and a tag count is quoted as coverage.

---

## Report back — under 30 lines

1. Exit status, runtime, total tag count, confirmation all tags unique and no untagged output.
2. Confirmation you did not open any `.wl`, the registry, or another engine. ⭐ If you did, say so.
3. **C1:** the literal implied-dimension and difference tags for the **main** action and for the
   **flexural** control, quoted exactly. ⛔ Do not interpret them.
4. **C2:** the literal `[mu_G]` dimension vector X8 produced. ⛔ Do not interpret or compare it.
5. **C3:** the literal E1–E4 values at `k = (0,0,0)`. ⛔ Do not interpret them.
6. **C5:** the literal anisotropic E1–E4 values at `ρ_z → ρ_br`. ⛔ Do not interpret them.
7. One line each for C4 and C6: what changed, tags added or removed.
8. Anything you could not do, and what blocked it.
9. ⭐ Anything that surprised you. This is wanted.

⛔ Do not commit to git. ⛔ Do not modify any file other than the deliverable and its `.premises` sidecar.
