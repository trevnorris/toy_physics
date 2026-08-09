# S9/S10 export chain — fix round 4. Decision list.

**You author the directive prose and apply the change.** Read `CLAUDE.md` first.

Scope: both SymPy engines · `scripts/extract_knob_inventory.py` · both generated modules · both
regenerated `scripts/out/*.out` · a new directive under `directives/`.
⛔ Out of scope: every `.wl`, every step record, every `.tex`. ⚠ S10 is slow: wrap runs in `timeout 1800`.

⭐ Four independent legs measured rounds 2 and 3. **Every repair rejects its weakest change, and no
computed physics value moved.** ⛔ Nothing below is a physics defect — each is a way the **ledger's record**
could be wrong without any guard noticing.

⚠ **Read this before starting.** The last two rounds were satisfied *literally*: I named an object, you
made that object true, and the level above it stayed broken — values became immutable while the mapping
stayed mutable; words became `Str` while `Str` stayed unguarded. ⭐ **Each item below names the property
AND the level the previous fix stopped at.** Satisfy the property.

---

## F1 · ⛔⛔ A generated module on disk is a product of the run that last attempted it.

⚠ **The previous fix stopped at run-time failure.** Both engines unlink before writing — but the unlink is
inside `main()` and the upstream import is at **module scope** (`S10:18`, `S9:11`). Measured, with a good
module from an earlier run on disk:

| | |
|---|---|
| a guard fails **during** the run | module removed ⭐ |
| the run fails **at import** | module **SURVIVES, byte-identical** ⛔ |

⇒ **This is the state the chain itself produces:** S9 fails → S9 unlinks its own module → S10 cannot import
it → S10 exits 1 → **the previous run's S10 ledger is still there.** A consumer then binds an S10 ledger
built from an S9 ledger that no longer exists.

**Decision: no failure mode leaves a module that a consumer cannot tell is stale.** ⭐ Two properties, and
you need both: a failed attempt does not leave the previous module behind, **and** a module says what it
was built from. ⚠ Measured: no hash, timestamp, version or upstream marker exists anywhere in either
generated module.

## F2 · ⛔⛔ An authored conclusion cannot become a ledger record — through ANY carrier.

⚠ **The previous fix stopped at `str`.** `contains_authored_text` (`S9:699-706`) tests
`isinstance(value, str)`, so the sympy `Str` **this project introduced last round** returns `False` — and
`S10` has **no such guard at all** (0 occurrences). Measured: emitting
`Str("light_requires_shear_confirmed")` from S10 exits **0** and produces four records classed `DERIVED`,
`step: 'S10'`, with **no field distinguishing them from a computed quantity.**

⇒ ⚠ **This is the defect class the whole rebuild exists to remove**, with a carrier the rebuild created.

**Decision: a record whose whole value is an authored word is distinguishable from a computed quantity,
and the population is emitted so it can be read.**

⚠ **Not every `Str` is a defect and you must not delete the honest ones** — field names inside rows, and
branch outcomes whose operands are emitted alongside (`S10:1200-1208`, `:1555-1583`), are legitimate.
⛔ Nothing currently counts or bounds them. ⛔ If you cannot make the distinction structural, **name that
and stop** — ⛔ do not introduce a maintained list of allowed words.

## F3 · ⛔⛔ The assumptions a consumer inherits are the assumptions the engine reasoned under.

There are **two independent channels**: the sympy assumptions on an exported `Symbol` (what a consumer
inherits) and the `Q.*` predicates the engine refines under (what decides its own results). ⛔ **No residual
compares them.**

⚠ Measured — widen `rho_br` from `positive` to `real` at its single declaration (`S9:28`), which is a real
physics change, a density no longer required positive:

> run exits **0** · every guard passes · module written · **102 of 165** exported values change ·
> all 102 **string-identical** · `display` moves for **0** records · the committed transcript moves by
> **ONE line**, in a control package.

⭐ The channels currently **agree** — 0 disagreements measured. ⇒ this is a **latent gap, ⛔ not an active
error** — and it is load-bearing, because inheriting the upstream assumptions is the entire design intent.
⚠ It is also already a *measured* defect in miniature: `lambdaScale` was declared `real` while the same
engine refined under `Q.positive(lambdaScale)`.

**Decision: a disagreement between the two channels is emitted as operands and a residual.**

## F4 · Two small ones.

- ⭐ **Restore a measure of unit coverage.** The round deleted `COMPUTED_/ABSENT_DIMENSION_COUNT` as dead,
  correctly — but with them went the transcript's **only** statement of how many exported objects carry no
  units. ⛔ Do not restore a counter that cannot vary; emit the coverage of the mechanism that replaced it.
- ⭐ **The mapping the ledger hands a consumer is not editable.** ⚠ Values are immutable; `LEDGER[k] = …`
  and `LEDGER[k]["class"] = …` both still succeed.

---

## Record, ⛔ do not try to fix

⭐ **`symbol_binding_residual` cannot see a quantity under two names.** Measured by two legs
independently: with both spellings bound, residual reads **0**. It counts variants **per name**, and a
quantity→name split is a different question. ⇒ ⭐ **state in your directive that the ledger is clean of
splits by inspection, ⛔ not by instrument**, and that S11 adding names is unpoliced.
⛔ Do not build a check for it — the object it would need is a name→quantity map, and that is the
hand-written pair table whose deletion this design exists to make permanent.

## Constraints

- ⛔ The derivation, action, ansatz and every computed physics value stay untouched. ⭐ Report the complete
  `.out` diff for both engines.
- ⭐ Every S9 record S10 does not overwrite stays identical between the two generated modules.
- ⛔ Do not state anywhere what a count, tally or partition size comes out as. Emit it and let it be read.

## Acceptance

⛔⛔ **A test that passes on a weaker fix is not a test.** For each item construct the weakest change that
should be rejected and show whether it is. ⛔⛔ **A FORM ablation is mandatory.**
⚠ **Verify your mutation is still present in the object you mutated** — a leg lost an injection to
`simplify` and got a false pass.

## Deliverables

The changed files, both literal stdouts, all diffs, and every ablation script with its stdout at named
absolute paths **outside the repository**.
