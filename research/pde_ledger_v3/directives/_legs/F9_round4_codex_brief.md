# Fix brief — `F9`'s guard paragraph. ⛔ ONE defect, ⛔ nothing else.

## Authority

Edit **only** the `F9` section of
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md`
(from `## ⭐⭐ F9 — WHEN A STEP'S EXPORT WRITES A KEY THE IMPORTED LEDGER ALREADY HAS` to the start of
`## ⛔ What this list does not decide`). ⛔ Change nothing else in the file, and ⛔ no other file.

`CLAUDE.md` binds. This is a **decision list**, ⛔ not a build directive: it says **what must become true**,
⛔ never how, and ⛔ never what anything equals, is expected, or was measured
(rule 5 — the previous author breached this twice and both breaches were caught by review).

## What is wrong

`F9` ends its rule with a retraction and **two** requirements on the publish guard:

> the guard sees every write **under its final key**, and it admits an `F9b` whose predecessor is **any**
> imported row.

⭐ Both are **necessary**. ⛔ **Neither, nor both together, is sufficient.** Two independent review legs
built four wrong exports that satisfy them and still publish. Their constructions, abstracted:

```
WRONG_ROUTE_EQUALITY_RESIDUAL=0        WRONG_ROUTE_FINAL_KEY=s11_slot
WRONG_ROUTE_PUBLISHED=True                  # proved-equal pair written PREFIXED

WRONG_PAYLOAD_MOVEMENT_FROM_LIVE_OBJECT=-B*y + B*z
WRONG_PAYLOAD_PUBLISHED=True                # payload is not this step's derived object
```

The four, as classes:

1. ⛔⛔ **A payload that is not this step's derived object publishes.** The writer can put the **imported**
   value back as its own row: the equality residual is `0`, the guard passes, and ⛔ a **form ablation of
   this step's action never moves the published row.** ⚠ `F9`'s existing ablation requirement is aimed at
   the **comparison**, ⛔ not the payload — the comparison genuinely runs, and the payload is borrowed.
2. ⛔ **A proved-equal pair can be written prefixed.** The final key is then absent from the imported
   ledger, so ⛔ nothing looks at it. ⇒ the bare name keeps the predecessor and *"override when it should"*
   is lost silently.
3. ⛔ **Undecided can be collapsed to equal** by handing the comparison two tokens instead of two live
   objects.
4. ⛔ **The check is not closed over what is published** — inspecting only the writer's own records leaves
   an illicit bare write in the merged result invisible.

## What must become TRUE — ⛔ state properties, ⛔ not a procedure

Replace the two-requirement paragraph with what actually has to hold. It must at least close all four
above. ⭐ Say it as properties of the published artifact; ⛔ do not write an algorithm, ⛔ do not name a
routine, and ⛔ do not add a control that lives inside the thing it polices.

⚠ Known traps, each measured:
- ⛔ **Naming these four as exceptions breeds a fifth.** ⭐ State the property; a list of cases is a
  regression waiting for the next one.
- ⛔ **An obligation the shared spec already carries must be POINTED AT, ⛔ never re-worded** — every
  restatement in this workstream came out weaker than the original. ⭐ Search
  `S11_SHARED_PHYSICS.md` first; `§5`'s corollaries are the likely home for at least one of the four.
- ⛔ Anything requiring the two engines to agree does **not** belong here: `F9` is PY-only and the sibling
  engine writes no ledger.

## Boundaries

- ⛔ Do **not** re-open the three-way outcome, the prefix, or the blocked record. Both legs cleared them.
- ⛔ Do **not** compute, quote or imply any S11 spectrum, determinant, root or count. The engine has not
  been built and those answers must not exist yet in any readable artifact.
- ⭐ Keep the retraction visible — the false claim it replaces is part of the record.
- ⭐ Shorter is better. The section is already the length at which this workstream's defect rate rises.

## Deliverable

The edited `F9` section, and a short note of what property you chose for each of the four classes and
which existing obligation you pointed at rather than restated.
