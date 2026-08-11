# Independent physics review — the export-key naming rule, and the S11 PY decision list

## Artifacts — ⭐ both, in one round

1. `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` — ⭐ read
   **`F8`** and the amended **`F4`**; `F1`–`F7` are already two-legged and settled.
2. `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_py_decisions_v3.md` — the list a builder
   will be handed to rewrite `scripts/S11_stray_longitudinal_sympy_audit.py`.

Both are **orchestrator-written**. ⚠⚠ **A decision list is the one artifact in this pipeline that is checked
once and then trusted** — everything downstream of it is reviewed twice. Two earlier versions of artifact 2
were blocked by two legs each; ⛔ do not assume this one is closer to right.

## Read in this order — ⛔ the order is load-bearing

1. `/var/projects/toy_physics/CLAUDE.md`.
2. `git show cf4a21a4:research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — the physics spec.
3. `/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md` — ⭐ §7's package
   actions in particular, and `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_exports.py`,
   the committed artifact to be renamed.
4. `/var/projects/toy_physics/research/pde_ledger_v3/S9_REWRITE_PLAN.md` — `D1`–`D12`.
5. ⭐ **Only now** open the two artifacts.

⛔ Do not read them first. Form your own view of what the chain requires, then read what was decided.

## Do not read

- The other leg's output.

## ⚠ For the Codex leg specifically

A design consult was run before `F8` was written, and its recommendation was **partially rejected**. It
proposed a canonicalized `(T,W)` descriptor **hashed** into every key, with a normative **cross-CAS**
serialization. `F8c` rejects the hash and the cross-CAS requirement, on the ground that **only one engine
writes the ledger**, so the rule needs to be stable across **steps**, not across engines.

⭐⭐ **Attack that rejection.** If the one-engine argument is wrong — if anything requires the Wolfram engine
or a second CAS to produce, read or agree on a storage key — then `F8c` is unsound and the hash is needed.
⛔ Verify it from the spec and the engines, ⛔ not from the claim. ⛔ Do not defend the consult's proposal
because it was the consult's; ⭐ show where the artifact is wrong.

## What to check

### ⭐⭐ 1 · `F8a` — is the partition really a measurement?

`F8a` says: mutate the action; if the object does not move it is **universal** and keeps a bare name, and
its cross-step collision is the check working; if it moves, its key must say which action.

- ⭐ Run it. Take S10's committed output, compare `MAIN` against packages that genuinely change the
  stiffness **form**, and produce your own partition. ⛔ Select the controls by reading their actions, ⛔
  never by their tag prefix — the prefix has already been measured to be a false guide.
- ⛔ **Find a false universal**: an object invariant under the mutations available in S10 that would
  nonetheless move under a different action. What happens to it downstream, and does anything catch it?
- ⛔ **Find a false action-dependent**: an object that moves for a reason that is not the action.
- ⭐ Is "does it move" decidable from committed data alone, or does it need runs that do not exist?

### ⭐⭐ 2 · `F8b`/`F8c`/`F8d` — does the naming rule hold up?

- ⭐ `F8c` makes the key a readable **slug** and puts the **canonical descriptor** in its own row, with two
  cross-checks. **Construct the case it misses**: two materials, or two steps, where slug and descriptor
  both agree and the objects still differ — or both differ and the objects are the same.
- ⛔ Is the descriptor computable from the live action **without** the package label and **without** the
  step name? If it cannot be, say so — that is a finding, not a detail.
- ⭐ `F8d` scopes every action-dependent object to the full `(T, W)` pair and records the lost `W`-only
  meeting as a known cost. Is that cost larger than stated? Name the objects that lose a real comparison.
- ⛔ Does `F8` anywhere prevent a collision that **should** happen? That is rule 6, and it is the failure
  that blocked both previous attempts.

### ⭐⭐ 3 · THE RENAME CONTROLS — ⛔ construct the failure, do not reason about it

`F8`'s six controls exist because re-pointing a name is the failure where **every check in the repository
passes and the physics silently moves**.

- ⭐⭐ **Build the derangement.** Swap two destination keys on a copy of the committed artifact — ⭐ choose
  rows whose payloads are equal, zero, or symmetry-invariant, because those are the ones payload checks
  cannot see — and show whether the six controls fail. Report the literal output.
- ⛔ **Then defeat them.** Construct a rename that satisfies all six and still moves physics. If you find
  one, that is the most valuable finding in this round.
- ⭐ Is the map's **totality and injectivity** checkable against the committed artifact today?

### ⭐⭐ 4 · THE PY LIST'S GUARDS — ⛔ each must be shown able to fail

- **`P4`** must be **total** over every value type the export admits and must not throw. ⭐ Take the
  committed export, apply the predicate to every row against its own reconstruction, and report how many
  rows it decides, how many it returns `UNDECIDED` for, and how many it fails on. ⛔ Then find two objects
  it calls the same that differ physically, and two it calls different that are the same.
- **`P5`** compares whole rows against a snapshot and reports **changed minus declared**. ⭐ Construct: an
  undeclared value write; an undeclared metadata write; a declared re-derivation that reproduces its row
  exactly; a declared re-derivation with a **wrong** value. Which are visible, and in which operand?
- **`P6`** requires publication atomic with respect to completeness. ⭐ Construct the incomplete run and show
  what a downstream import sees afterwards.
- **`P7`** compares eligible-emitted against written. ⭐ Omit one computed object and confirm it is the only
  guard that fires.
- **`P2`** reads the row schema from the imported artifact at run time. ⛔ What breaks on the first import
  whose schema differs — and is that better or worse than a typed schema?
- **`P1`** must decide the boundary for the **symbol population**, which is manufactured by traversal, not
  from tags. ⭐ Is it decidable as written?

### ⭐ 5 · THE RULES NEITHER ARTIFACT MAY BREAK

- **Rule 5 — ⛔ never state what anything equals, is expected, or was measured about the physics.** ⚠ The
  previous version was caught handing a builder measured values, counts, membership, and an expected
  vacuity. ⭐ Check both artifacts line by line for a value, count, membership, sign, or outcome a builder
  could converge on. ⚠ A **prohibition** leaks as surely as an assertion.
- **Rule 3 — name the object, ⛔ do not specify the recipe.** Which items specify *how* rather than *what
  must become true*?
- **Rule 6 — ⛔ do not design away the disagreement.**
- **Rule 2 — both operands and the residual, then guard.**
- **Physics filter (standing user instruction):** report an item only if it catches a way the **physics**
  could be wrong. ⛔ If an item is **process rather than physics**, say so — that is a finding, because it
  is cost with no correctness behind it. ⭐ The user's stated goal for this round is **fewer unnecessary
  collisions with less machinery**; an item that adds machinery without catching a physics defect fails it.

### ⭐ 6 · WHAT IS MISSING

⛔ The most expensive defects here have been **absent computations**. What could go wrong in the rename or
the export that **no** item in either artifact would detect? ⚠ Name it concretely, with the mechanism.

## ⛔ Constraints on how you run

- ⛔ Read-only on the working tree. Copy to `/tmp` and modify the copy.
- ⛔ `timeout 600` on every run. ⛔ Never raise it. ⛔ No Mathematica this round.
- ⭐ Save every script and its literal stdout to named absolute paths and report them. ⛔ A prose derivation
  will be discarded — a claim with no computation behind it is the defect class this rebuild removes,
  relocated into the review.

## Report

For each finding: **what is wrong**, **the mechanism by which a wrong result would survive**, **the literal
output that shows it**, and **what must become true instead** — ⛔ not a rewrite.
⭐ Separate blocking from non-blocking. ⛔ "Nothing found" with no ablation behind it is not a result: say
which computation could have found something and did not.
