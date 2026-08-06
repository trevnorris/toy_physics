# Fix `OpaqueDerivative`'s process-global identity, before any S11 spelling is declared

**Do not commit.** One defect, one file, plus a measurement script. ⭐ **Verify the measurement below
yourself before you rely on it** — it is the orchestrator's and it has not been through a review leg.

---

## The defect

`reduction/engine_output_checks.py:157-171`. `OpaqueDerivative` subclasses `sp.Symbol` and calls
`sp.Symbol.__new__(cls, rendered)`, which is **cached on the symbol name**. The constructor then stamps
`orders`, `function_name` and `variables` onto whatever object the cache returned. ⇒ constructing an
`OpaqueDerivative` with a rendered name that already exists **returns the existing object and overwrites
its identity attributes in place**.

`_map_tree` (`:894-902`) rebuilds the atom under a rename while keeping the old rendered name:

```python
return OpaqueDerivative(str(item), item.orders, function, variables)
```

⇒ the *original* atom is mutated.

### What the orchestrator measured — ⭐ re-run it, ⛔ do not transcribe it

Five probes, run against the committed module:

| probe | rename applied | `mapped is original` | which identity field moved on the original |
|---|---|---|---|
| 1 | key is the differentiated **function name** | `True` | `function_name` |
| 2 | key is a **differentiation variable** | `True` | `variables` |
| 3 | key is neither | `True` | none |
| 4 | none — a **fresh** atom built with the original identity | `True` | it stamped the identity **back** |
| 5 | none — two atoms, **same rendered name, different `(orders, function_name, variables)`** | `True` | the first was overwritten by the second; `==` and `hash` agree |

⇒ two distinct failure channels:

- **Channel A (probe 1, 2, 4).** The attributes are **process-global state keyed on the rendered name**, and
  every construction overwrites them. An atom's identity at any instant is whatever the **most recent**
  construction with that name wrote — regardless of which row, expression or engine is being compared.
  ⚠ S10's record calls this *"LRU-cache dependent"*; the measurement above says it is **last-writer-wins**,
  which is stronger and does not depend on allocation volume.
- **Channel B (probe 5).** Two genuinely different derivatives that render to the same string are the
  **same object**, compare equal, and hash equal. ⚠ **This channel needs no declared rename at all**, so it
  is not gated on the S10 record's *"before any field or coordinate spelling is declared"*.

## What must become true

⭐ State these as invariants and choose your own implementation; ⛔ do not copy a recipe from this file.

1. Constructing an `OpaqueDerivative` **never mutates a previously constructed one.** After building any
   number of atoms, every earlier atom's `(orders, function_name, variables)` is what it was built with.
2. Two `OpaqueDerivative`s are **equal and hash-equal iff their full identity matches** — the rendered
   name **and** `orders` **and** `function_name` **and** `variables`. ⇒ probe 5's two atoms must be
   unequal.
3. ⛔ **`str(atom)` still returns the rendered string, unchanged.** Emission, display and every
   name-keyed code path must be untouched.
4. `CanonicalDerivative` (`:174-191`) already encodes its full identity in its symbol name. ⭐ Check
   whether it has the same defect and say what you found; ⛔ do not change it if it does not.

## What must NOT change

⛔ The comparator's output on the **committed** S9 and S10 configurations. Both are the regression baseline
and neither may move by a single counter.

---

## Deliverables

**A.** The fix in `reduction/engine_output_checks.py`.

**B.** `reduction/measurements/opaque_derivative_identity.py` — a script that prints, for each of the five
probes above, the atom's identity **before**, the identity **after**, whether the rebuilt object **is** the
original, and the equality and hash relations. ⛔ It prints operands and residuals; ⛔ it asserts nothing
and states no conclusion. ⚠ Scripts in that directory resolve the repository root with `parents[2]`.

**C.** A unit test in `reduction/test_engine_output_checks.py` covering invariants 1 and 2. ⛔ It must fail
against the unfixed module — **show that it does**, by reverting the fix, running the test, and pasting
the failure.

---

## Acceptance — run these and paste literal stdout

1. `reduction/measurements/opaque_derivative_identity.py`, whole output.
2. The new test failing on the unfixed module, then passing on the fixed one.
3. The full unit suite from `reduction/`.
4. `harness_ablation.py` — every `ACCEPTANCE` line.
5. `engine_output_checks.py --config reduction/checks_S10.yaml` and `--config reduction/checks_S9.yaml`,
   the counter lines from each, **beside the same lines from before your change**. ⭐ Produce the "before"
   by stashing your edit, not by copying a number out of a document.
6. `git status --short` and `git diff --stat`.

## Constraints

- ⭐ A script may run long **if it is visibly progressing** — print incrementally and flush. ⛔ What is
  forbidden is a long silent stretch, which is indistinguishable from a runaway solve. ⛔ Do not weaken a
  computation to make it faster.
- ⛔ Do not launch Mathematica or `wolframscript`.
- ⛔ No `git add`, no `git commit`, no `git checkout`, no other git write except a stash you restore.
- ⛔ Do not touch `directives/S11_SHARED_PHYSICS.md`, either engine, or any committed `.out`.

## Report back — under 20 lines

1. Deliverables A/B/C: done or not.
2. The acceptance output.
3. ⭐ Whether your own measurement of the five probes agreed with the table above — **including if it did
   not**. Say which probe and what you got.
4. Anything you found that this directive gets wrong.
