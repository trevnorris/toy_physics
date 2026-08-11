# Build directive — a point is its locus's solve variables, keyed by shared names

⭐ **Target:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`.
⛔ **Edit ONLY that file.** ⛔ Do not commit. ⛔ Do not run either engine.

⚠⚠ **Both of your last report's blocking findings were right, and I have verified both:**
- ⭐ **`D12` does NOT reach a payload's keys.** ⚠ **My directive truncated it at the decisive clause** — the
  full text ends `emit["rho_br_dimension", rhoBrDimension]`, whose payload argument is the raw internal WL
  variable. ⇒ ⛔ **`D12` is not the authority for this and is not cited again.**
- ⭐ **`c_s0` is not a coefficient** by this file's own construction (`:517-520` builds
  `COEFFICIENT_ORDERING` from `L_pkg`'s factors and the kinetic density; `§Q11:997-998` writes
  *"`COEFFICIENT_ORDERING` **together with** `c_s0`"*). ⇒ ⛔ my `P2` wording was wrong.

⭐⭐ **Both items below are rewritten because of your report.** ⛔ Do not invent a mechanism. ⭐ **Reporting
is success** — ⛔ if an item cannot be written, or a fact here proves false, STOP AND REPORT.

---

## ⭐⭐ Q1 · A point assigns exactly its locus's NAMED SOLVE VARIABLES — ⛔ nothing more, nothing less

⛔ Not *"the coefficients and the wavevector"* — ⚠ that is what broke last round.

⭐ **What must become true:** wherever this file requires an exact point — `§5`'s `_REAL_WITNESS`,
`§Q8b`'s `STRATUM<s>_POINT`, and anything that consumes them — the point is an assignment to **exactly the
solve variables that locus already names**, ⛔ and to nothing else.

⭐ **This needs no new category.** Every locus in this file **already names its solve variables**: `§Q8a`'s
three sets, `§Q3`'s coincidence families, `§Q11`'s `k_w² = 0` locus, and `§5`'s `_SOLUTION`, which is
required to be emitted *"with the solve variables named"*.
⇒ ⭐ `c_s0` is carried wherever a locus names it, ⛔ and an object that is never a solve variable is never
in a point.
⚠ ⛔ **Do not write this as a denylist** of what a point excludes. ⭐ State what it assigns.

## ⭐⭐ Q2 · The point payload's KEYS are shared names — ⛔ and this file is the authority, ⛔ not `D12`

⚠ **The problem this closes:** a point keyed by one engine's internal spellings cannot be read by the
other, and translating them needs the hand-written cross-engine pair table this project exists to abolish.

⭐ **What must become true:** wherever this file requires an exact point, its payload's **keys** are the
**standard names this file already uses for those objects**, in both engines. ⛔ An engine's internal
variable spelling never appears as a key. ⭐ Neither engine is obliged to rename an internal variable —
⛔ only what it emits.

⭐ **This file already pins payload FORMS** — the boolean-test record, the failure payload, the orderings.
⇒ ⭐ **pinning a payload's keys is the same kind of rule, and this file may simply state it.**

⛔⛔ **VERIFY AND REPORT, ⛔ do not assume:**
- ⭐ **Is there already a declared list of standard names in this file covering every symbol a point can
  carry** under `Q1`? ⚠ `§Q6r` carries one such list — ⛔ **but you reported it sits inside an
  engine-local section.** ⭐ Say whether it can serve as a shared vocabulary where it stands, whether it
  needs a structural home outside that section, or whether it does not cover everything `Q1` admits.
- ⛔ **If no adequate list exists, STOP AND REPORT** rather than authoring one. ⭐ Tell me exactly which
  names would be missing.

## ⭐ Q3 · Retire the re-spelling instruction — ⛔ only if `Q2` lands

⚠ Two places instruct the orchestrator to re-spell the point with `§Q6r`'s map; ⭐ both review legs
measured that the map cannot do it.
⭐ **Under `Q2` no re-spelling happens at all** ⇒ the instruction goes and is ⛔ **not replaced**.
⛔ **If `Q2` does not land, leave those lines alone and say so** — ⚠ you were right last round that
deleting them without a replacement makes each engine's spelling a silent, unstated binding.

---

## ⛔ Constraints

- ⛔⛔ **THE SPEC SAYS WHAT TO COMPUTE. ⛔ NEVER what anything equals or was measured.** ⛔ No count, symbol
  list, spelling pair, or the fact that some object was found unnecessary, may appear in the file in any
  form a builder could **derive** it from. ⭐ Test: **derivability**.
- ⛔⛔ **`§4`'s quoted block is shared VERBATIM with the other steps' specs — ⛔ do not edit inside it.**
  ⭐ Verify byte-identity against `directives/S10_SHARED_PHYSICS.md:111-113` and report it.
- ⛔ Do not touch `_CANONICAL_LOCUS`. ⛔ Do not delete or weaken any `§5` locus object.
- ⛔ Do not re-open `N1` or `N3`.
- ⛔ No admissibility algorithm, no recursive stratification, ⛔ no rule pinning how an engine **describes**
  a component.
- ⚠⚠ **AFTER EDITING, grep the whole file for sentences quoting a CONSEQUENCE of anything you changed.**
  ⭐⭐ **Five of the six previous rounds left exactly one stale sentence, every one a downstream sentence
  asserting a consequence the edit falsified.**

## Deliverables

1. The edited file.
2. ⭐ Per `Q1`–`Q3`: lines changed and the sentence that now makes it true — ⭐ or why it did not land.
3. ⭐ Your answer to `Q2`'s verification: which declared list you used, where it lives, and whether it
   covers everything `Q1` admits.
4. ⭐ The `§4` byte-identity result.
5. ⛔ Anything you could not write without inventing, or any fact above that proved false ⇒ **STOP AND
   REPORT.**
