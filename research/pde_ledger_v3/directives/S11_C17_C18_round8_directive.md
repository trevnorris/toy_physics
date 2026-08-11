# Build directive — the decision you asked for: declare the vocabulary, and scope the transport to it

⭐ **Target:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`.
⛔ **Edit ONLY that file.** ⛔ Do not commit. ⛔ Do not run either engine.

⚠ **Your `Q2` report is accepted in full.** ⭐ `Q1` stands as written; ⛔ do not re-open it.

⭐⭐ **The decision, and it is mine:** this file **declares a shared payload-key vocabulary**, and the
point transport is **scoped to it**. ⚠ Your point (c) is why the scope clause exists — ⛔ a fixed list
cannot pre-name a parameter a CAS invents, ⭐ so the transport does not pretend to cover that case.

---

## ⭐⭐ R1 · Declare the shared standard-name vocabulary, in a section BOTH engines execute

⭐ **What must become true:** this file declares, in one place, the **standard names** for the objects a
point can be keyed by:
- ⭐ the coefficients, and `c_s0` — ⚠ **the names already listed at `§Q6r`**, ⛔ unchanged in content;
- ⭐ the **wavevector components**, which `§Q6r` omits.

⛔ **It must NOT live inside `§Q6r`** — ⚠ you established that section is engine-local, so a shared
vocabulary cannot be declared there. ⭐ **You choose the structural home** — `§8` governs tag and payload
conventions and is the obvious candidate — and ⭐ `§Q6r` then **references** the declaration instead of
carrying its own copy. ⛔ Do not change what `§Q6r` maps or what it is for.

⚠ ⭐ **The wavevector half is a DECLARATION OF SOMETHING ALREADY TRUE, ⛔ not a new binding** — both
engines already spell those components identically. ⛔ Do not state that fact in the file; ⭐ just declare
the names.

## ⭐ R2 · A point's payload keys are those names — ⛔ never an engine's internal spelling

⭐ **What must become true:** wherever this file requires an exact point, its payload's **keys** are the
`R1` names, in **both** engines. ⛔ An internal variable spelling never appears as a key.
⭐ Neither engine is obliged to rename an internal variable — ⛔ only what it emits.
⚠ This file already pins payload **forms** (the boolean record, the failure payload, the orderings);
⭐ **this is the same kind of rule, on this file's own authority.** ⛔ Do not cite `D12` or `D11`.

## ⭐⭐ R3 · Scope the transport — ⛔ and say plainly what is NOT transported

⚠ **Your point (c), accepted:** a change sub-locus's solve variables come from
`STRATUM<s>_FREE_PARAMETERS`, which this file **deliberately refuses to pin**, and a CAS may introduce a
parameter no declaration can name.

⭐ **What must become true:** `§Q8c` transports a point **only when every solve variable its locus names
is in the `R1` vocabulary.** ⛔ A point whose locus names anything outside it is **not** transported — and
⭐ **that is emitted as an explicit coverage fact**, ⛔ never silently skipped.
⚠ ⭐ This is a **stated limit of the exchange**, ⛔ not a defect to hide: `§9` is where this file records
what it does not establish. ⭐ Put it there too if that is where it belongs.

## ⭐ R4 · Now retire the re-spelling instruction

⭐ Under `R2` no re-spelling happens. ⭐ The two places instructing the orchestrator to re-spell the point
with `§Q6r`'s map (you located them at lines `183` and `822`) go, ⛔ and are **not replaced**.
⚠ ⭐ **This was blocked for three rounds precisely because deleting it without `R2` would leave each
engine's spelling a silent binding.** ⭐ `R2` is that replacement.

---

## ⛔ Constraints

- ⛔⛔ **THE SPEC SAYS WHAT TO COMPUTE. ⛔ NEVER what anything equals or was measured.** ⛔ No count,
  spelling pair, or measured fact may be **derivable** from the file. ⚠ In particular ⛔ do not write that
  the two engines already agree on the wavevector spellings, ⛔ or that anything was found unnecessary.
- ⛔⛔ **`§4`'s quoted block is shared VERBATIM — ⛔ do not edit inside it.** ⭐ Verify byte-identity against
  `directives/S10_SHARED_PHYSICS.md:111-113` and report it.
- ⛔ Do not touch `_CANONICAL_LOCUS`. ⛔ Do not delete or weaken any `§5` locus object.
- ⛔ Do not re-open `N1`, `N3`, or `Q1`.
- ⛔ No admissibility algorithm, no recursive stratification, ⛔ no rule pinning how an engine **describes**
  a component. ⚠ ⛔ `R1` declares **names**; it does ⛔ **not** pin `STRATUM<s>_FREE_PARAMETERS`.
- ⚠⚠ **AFTER EDITING, grep for sentences quoting a CONSEQUENCE of anything you changed.** ⭐ Your `Q1`
  sweep was the first clean one in seven rounds — ⭐ do that again, including the `§Q6r` reference you
  create and anything that reads the transported point.

## Deliverables

1. The edited file.
2. ⭐ Per `R1`–`R4`: lines changed and the sentence that now makes it true.
3. ⭐ Where you put the `R1` declaration and why; and what `§Q6r` now says.
4. ⭐ The `§4` byte-identity result, and your stale-sentence sweep.
5. ⛔ Anything you could not write without inventing, or any fact above that proved false ⇒ **STOP AND
   REPORT.** ⭐ Reporting is success.
