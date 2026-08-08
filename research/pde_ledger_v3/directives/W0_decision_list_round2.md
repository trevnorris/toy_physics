# W0 round 2 — decision list for repairing the builder and integrator statements

⛔ **The W0 builder statement is NOT safe to issue.** Two independent legs reviewed it; both found a
compliant build that emits a value recombined from the registry rather than derived from the spectrum.
⭐ The orchestrator verified the decisive items. ⇒ this is a **repair round on Codex's own artifact**, ⛔ not
a rewrite by the orchestrator.

⭐ **What round 1 got right, ⛔ do not regress:** the builder/integrator split; the scope table; the
dimension pin that excludes emitting a frequency; the pinned-failure escape; the satisfiable regression
projection. ⚠ Cardinality-one was verified to hold in **all ten in-scope cells**, by both legs, per cell.

---

## ⛔⛔ E0 · SEQUENCING — W0 cannot launch at all right now

⭐ **S10-python does not run.** The uncommitted registry change makes a declared dimension a
`BoundDimensionLaw`, which the engine iterates:

```
TypeError: 'BoundDimensionLaw' object is not iterable      EXIT=1
```

⇒ ⛔ The builder statement's first regression instruction — *"capture the assigned engine's committed
output… compare it with the post-change output"* — is **unexecutable** for that build.
⇒ ⭐ **W0's four builds are blocked on `W3_fix_round2.md` landing.** ⭐ Say so in the integrator statement.

---

## The defects to repair

### ⛔⛔ E1 · The defining residual does not pin `v_T²` to the selected root

⭐ The statement says `v_T² = r_T/K` in **prose**; the only required emitted residual is `L_T = c_T² − v_T²`.
⇒ ⛔ a build may write a value recombined from material coefficients into the `v_T²` slot, and `L_T` is
still `0`, dimensions still `[1,-1,0]`, certificate still satisfied.
⚠ On **every in-scope cell** `r_T/K` is algebraically identical to that recombination, so ⛔ nothing
downstream — including four engines agreeing — can tell the difference.

⇒ ⭐ **Decision: require an emitted residual tying `v_T²` to the selected root and `K`**, with the same
status as `L_T`: both operands and the residual, then guard.
⚠ This is the original S9 defect — *define the object, then assert a residual that is zero by
construction* — relocated one square root later.

### ⛔⛔ E2 · The certificate's predicate is computed and then never used to select

⭐ Five per-root fields are required, including a computed transversality predicate. ⭐ The selected set,
its cardinality and the selected root are required **separately**. ⛔ **No sentence binds the selection to
the predicate.**

⇒ a build selecting by *"the root that is nonzero"* satisfies every sentence, and the certificate is
decoration. ⚠ In scope the two criteria coincide, so ⛔ nothing in scope can expose it — but they do **not**
coincide in general: a control package exists with a **nonzero root whose transverse nullity is `0`**.

⇒ ⭐ **Decision: the selected set is the truth set of the emitted per-root predicate**, and each root's
membership decision is itself an emitted computed object paired with its predicate value.

### ⛔ E3 · "The required cardinality is one" is an expected value handed to the builder

⛔ A builder whose criterion is subtly wrong sees a count of 2, knows the target is 1, and tunes the
criterion until the count matches. ⚠ The discriminating step is exactly the one the number lets them
reverse-engineer.

⇒ ⭐ **Decision: the builder emits the selected set, its computed cardinality, and `v_T²`/`c_T` for each
selected root. ⛔ No cardinality requirement appears in the builder statement.** ⭐ The cardinality verdict
moves to the integrator statement, where verdicts belong.

### ⛔ E4 · The mutation does not separate an honest object from a recombined one

⛔ Injecting a sentinel **into `v_T²`** is downstream of any recombination that writes that symbol — both
paths move. ⇒ the test distinguishes nothing.

⇒ ⭐ **Decision: the sentinel enters upstream of the quotient** — on the selected root or on `K` — with the
production value recomputed from the mutated input and all coefficient emissions held fixed.
⭐ Mutation 2 must hold the **solved spectrum** fixed and change only the selection, ⛔ so it cannot be
satisfied by perturbing an upstream coefficient.

### ⛔⛔ E5 · The mutation's control instructs the builder to construct the target

⭐ A non-negative scalar of dimension `[1,-1,0]` built from the two named material coefficients is
**uniquely determined** — the exponents are forced by the dimension algebra for every `D`. ⇒ ⛔ mandating
that control puts the relation's other operand into the builder's own stdout, where it becomes an
acceptance target to iterate against.

⇒ ⭐ **Decision: name no material coefficients. Hold every coefficient emission fixed.** ⭐ State the
required observable on the production path alone: `c_T` moves with the sentinel and no pre-existing
emitted record does. ⛔ No engine-side reconstruction of the target.

### ⛔ E6 · The boundaries clause forbids what the regression clause requires

⛔ Boundaries forbid changing any *existing emitted name/value/shape*; the regression section explicitly
anticipates that a runtime record moves every run and an emitted tag count must change on any addition.
⇒ read literally, one build is impossible.

⇒ ⭐ **Decision: the boundaries clause excepts exactly the records the regression section names.**

### ⚠ E7 · The manifest is a disclosure, ⛔ not an acceptance criterion

⭐ The bar checks pre-existing records unchanged, additions ≡ manifest, counts consistent. ⛔ All three are
self-consistency against a **builder-authored** manifest; nothing constrains what the additions are.

⇒ ⭐ **Decision: the manifest gives each new record's name AND value**, and the statement says plainly that
it is a review input, ⛔ not something a build can pass by satisfying.

### ⚠ E8 · Placement is unpinned across four builders

⛔ Whether `c_T` sits once per cell or once per root, and how S9's direction specialisations are handled,
is undecided. ⭐ Tag-name divergence the harness can absorb; **placement divergence decides whether a later
cross-engine row can bind the pair at all.** ⚠ This exact class already cost a false disagreement here.

⇒ ⭐ **Decision: pin the placement grain.** ⛔ Leave tag spelling free.

### ⚠ E9 · `L_T = 0` is not evidence, and the integrator currently says it is

⛔ `L_T` is identically zero for any build that defines `c_T` as a square root of `v_T²`. ⭐ The integrator
statement instructs a reviewer to check bindings against it — which will find zeros meaning nothing.

⇒ ⭐ **Decision: the integrator says plainly that `L_T = 0` carries no information**, and names which
residual can actually be nonzero — ⛔ that residual stays out of the builder statement.

---

## Rules for the repair

- ⛔ Do not restate physics that is not measured here, and ⛔ do not add a factual claim you have not
  produced with a command.
- ⭐ Every time you are about to write *how*, write *what* instead.
- ⛔ Do not modify any engine, any `reduction/` file, or any committed `.out`. ⛔ Do not commit.
- ⚠ ⭐ **Verify E0–E9 yourself before repairing.** ⛔ If one is wrong, say so and stop — ⭐ that is worth more
  than the deliverable. ⚠ An earlier decision list handed to you contained a false measurement and you
  correctly refused to build from it.
