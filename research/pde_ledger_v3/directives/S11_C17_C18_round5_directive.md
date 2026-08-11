# Build directive — three measured items on `S11_SHARED_PHYSICS.md`

⭐ **Target:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (1228
lines). ⛔ **Edit ONLY that file.** ⛔ Do not commit. ⛔ Do not run either physics engine.

⚠⚠ **Five revisions; the first four bred thirteen defects, every one inside a mechanism its author
INVENTED to make two engines agree.** ⭐ The fifth was restricted to measured items and told to **report**
anything it could not write without inventing — ⭐ **it reported two, and both were right.** ⭐ That is the
behaviour this round wants.

⇒ ⛔⛔ **DO NOT INVENT A MECHANISM.** ⭐ Each item below is settled by a measurement given with it.
⛔ **If an item cannot be written without inventing a construction, or a fact it rests on turns out false,
STOP AND REPORT IT.** ⭐ Reporting is a success here, ⛔ not a failure.

---

## ⭐ N1 · *"reduced modulo `STRATUM<s>_DEFINING_EQUATIONS`"* — ⛔ it never says WHOSE

**Where:** `§Q8b`, the component-scoped comparison paragraph.

⚠ Each engine emits **its own** defining equations, in **its own** chart. ⛔ The text names the tag and
binds it to no engine.

⭐ **MEASURED — on the committed `MAIN, D = 2` component the ambiguity does not bite, and the reason is
structural:** after name alignment both engines emit the **same principal ideal**, and generator
orientation (`mu_R − B_comp` vs `B_comp − mu_R`) does not change an ideal ⇒ all readings must agree.
⭐ **MEASURED — and here is a case where it does bite**, from a review leg:

| payloads | ideal A | ideal B | modulo A | modulo B |
|---|---|---|---|---|
| `p_A = y`, `p_B = y + z` | `{x}` | `{x, z}` | **`−z`** | **`0`** |

⇒ ⛔ **whenever the two paired strata do not cut the same variety, the reading decides between AGREEMENT
and a nonzero residual.**

⭐ **What must become true:** ⛔ do not pick one engine's equations. ⭐ **Both reductions are emitted** —
the difference reduced modulo each engine's own emitted defining equations — ⭐ **as two operands, with the
comparison stated on the pair.** ⚠ Where the two reductions disagree, ⛔ that is a **finding**, ⛔ never a
silent choice of the more convenient one.
⚠ ⭐ This is `CLAUDE.md` rule 2 — *emit both operands and the residual, then guard* — applied to the
comparison itself.

---

## ⭐⭐ N2 · The point cannot be transported BY NAME — ⛔ transport it POSITIONALLY

**Where:** `§4`'s witness-input block and `§Q8c` step 2.

**Now reads:** the orchestrator re-spells the point into the receiver's own coordinate names *"with the
closed map at `§Q6r`"*.

⛔⛔ **That map cannot do it.** ⭐ **MEASURED, by both review legs:**
- `§Q6r`'s map is **this file's notation → `LEDGER` key names**, ⚠ engine-local to the **importing** engine
  — ⛔ it is not an engine-to-engine map.
- Its right-hand side is **one ASCII spelling** (`rho_br`, `mu_R`, `B_comp`, …). ⛔ It contains **no**
  Wolfram spelling, and ⚠ **Wolfram symbols cannot contain `_`** — the committed Wolfram point uses
  `rhoBr`, `muR`, `bComp`, `muBr`.
- It covers the **coefficients only**. ⚠ The committed Wolfram `STRATUM<s>_POINT` payload carries
  `a1…a5`, `bComp`, `beta`, `dimsslot1–3`, `k1…k5`, `muBr`, `muR`, `rhoBr`, `s` ⇒ ⛔ **17 of 19 symbols are
  not on that map at all.**

⭐⭐ **What must become true — ⛔ and the point is to stop transporting NAMES:**
1. ⭐ **A point is an assignment to the coefficients in `COEFFICIENT_ORDERING` and to the wavevector
   components, and to nothing else.** ⚠ The committed points carry amplitudes and an engine-local
   `dimsslot*`; ⛔ those have no place in a point on a locus in coefficient-and-wavevector space.
2. ⭐ **It is delivered POSITIONALLY, against orderings this file ALREADY PINS** — ⛔ no symbol name
   crosses the boundary, so ⛔ no name map is needed and ⛔ none may be authored.
3. ⭐ The receiver reads each value into **its own** symbol at that position.

⛔⛔ **VERIFY BEFORE WRITING, and report if any of it is false:**
- ⭐ that **both** engines emit `COEFFICIENT_ORDERING`, and that this file already requires it per
  package-and-dimension;
- ⭐ that the wavevector components have a pinned index order both engines share;
- ⚠ that the orchestrator's existing ordering-alignment work is enough to make positions correspond —
  ⛔ **if the two `COEFFICIENT_ORDERING`s can differ and nothing aligns them, SAY SO** rather than
  assuming.

---

## ⭐ N3 · The fifth stale sentence — ⛔ entailed by the last round's coverage widening

**Where:** `§Q8c`, immediately below the witness block.

**Now reads:** *"If `WITNESS<w>_POINT_COVERAGE` is `INCOMPLETE_POINT`, the evaluated objects remain partly
symbolic and their counts are generic."*

⛔ **No longer true.** ⭐ Coverage is now decided over the **union** of every evaluated locus's symbols
**and** `M`'s ⇒ ⚠ the token can read `INCOMPLETE_POINT` because a **locus** equation lacks a symbol while
`M` is fully assigned and its spectrum is genuinely **point-local**.

⭐ **MEASURED**, on the committed Wolfram `MAIN, D = 2` point:

```
M after point:  Matrix([[omegaSquared/2 - 1, 0], [0, omegaSquared/2 - 1]])
remaining M symbols: [omegaSquared]
roots with multiplicity: {2: 2}
locus residual: -2*(cs0 - 1)*(cs0 + 1)/cs0**2
```

⇒ ⛔ aggregate coverage is `INCOMPLETE_POINT`, ⭐ yet the `M`-derived counts are **not** generic.
⚠ As written, valid spectrum evidence would be discarded or mis-scoped.

⭐ **What must become true:** the sentence says what is actually true of each part — ⛔ it must not
attribute genericity to objects the point fully determined.

---

## ⛔ Constraints

- ⛔⛔ **THE SPEC SAYS WHAT TO COMPUTE. ⛔ NEVER what anything equals, is expected, or was measured.**
  ⚠⚠ **Every measurement above is FOR YOU. ⛔ None of it may appear in the file, in any form a builder
  could DERIVE it from** — ⛔ no count, value, sign, coefficient, symbol list or named locus.
  ⭐ The test is **derivability**, ⛔ not literal presence.
- ⛔⛔ **`§4`'s quoted block is shared VERBATIM with the other steps' specs — ⛔ do not edit inside it.**
  ⭐ Verify `sed -n '156,158p'` of the target stays byte-identical to `sed -n '111,113p'` of
  `directives/S10_SHARED_PHYSICS.md`, and report the result.
- ⛔ Do not touch `_CANONICAL_LOCUS` — ⚠ it is deliberately inert on non-polynomial systems.
- ⛔ Do not delete or weaken any `§5` locus object.
- ⛔ Do not add an admissibility algorithm, recursive stratification, or any rule pinning **how an engine
  describes a component**. ⚠ A previous round pinned an elimination and it **deleted a branch**.
- ⛔ Do not re-open `M1`, `M2`'s rule, or `M4` from the last round — ⭐ they are discharged.
- ⚠⚠ **AFTER EDITING, RE-READ EACH EDITED SECTION AND ITS INTRODUCTION**, and grep the whole file for
  sentences that quote a **consequence** of anything you changed. ⭐⭐ **All five previous rounds left
  exactly one stale sentence, and every one was a downstream sentence asserting a consequence that the
  edit falsified.** ⛔ `N3` is the current one; ⛔ do not create the sixth.

## Deliverables

1. The edited file.
2. ⭐ Per `N1`–`N3`: the lines changed and the sentence that now makes it true.
3. ⭐ The `§4` byte-identity result, and the outcome of `N2`'s three verifications.
4. ⛔ Anything you could not do without inventing, any fact above that proved false, or any conflict ⇒
   **STOP AND REPORT.** ⛔ Do not invent a compromise.
