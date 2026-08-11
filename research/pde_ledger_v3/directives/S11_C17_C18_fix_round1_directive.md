# Fix round 1 — `S11_SHARED_PHYSICS.md`, after two review legs

⭐ **Target:** `directives/S11_SHARED_PHYSICS.md` (1118 lines). ⛔ **Edit ONLY that file.**
⚠ The decisions remain `S11_C17_C18_spec_repair_decisions_v2.md` (`T1`–`T7`). ⛔ Do not reopen them.

⛔⛔ **`T3` IS NOT DISCHARGED, and it is the one decision that closes the measured case.** ⭐ Everything else
below is smaller. ⚠ Two legs and the orchestrator verified each item at source.

## ⛔⛔⛔ R1 · Q8c CONTRADICTS §4 OF THIS SAME FILE — ⭐ the highest-severity item

⚠ `§4:156-158` is labelled **verbatim, non-negotiable**:
> *"The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE ACTION and the
> ANSATZ. Every other expression involving them must be REACHED BY COMPUTATION."*

⚠ and corollary 1's test (`:160-162`): *"if you deleted the computation, would this tag change?"*
⚠ and `§5`'s live-read exemptions (`:205-212`) are *"**closed** and are **exactly these four entries**."*

⛔ `WITNESS<w>_RECEIVED_POINT` is an assignment over the physical symbols that the receiving engine did
**not** reach by computation, and it is **not** in that closed list. ⇒ ⛔⛔ **A builder obeying `§4` will not
build `Q8c` at all.**

⭐ **Make the received point a NAMED, PERMITTED INPUT** — extend the closed exemption list, or `§4`'s
"action and ansatz" clause, whichever fits the file's structure. ⚠ Say plainly that it is an **input**,
⛔ not a computed object, and that ⛔ **no other object may be received.**

## ⛔⛔ R2 · Q8c CANNOT EXECUTE, AND THE FILE SANCTIONS THE EMPTY RESULT

⛔ `Q8c` requires knowing that *"one engine emitted a point and the other emitted none"* ⇒ that comparison
cannot exist on a first pass. ⚠ And `:747-749` explicitly blesses the empty outcome:
> *"An absent counterpart point produces no indexed witness object."*

⇒ ⛔⛔ **Both engines may emit `WITNESS_ORDERING: {}` and be fully compliant**, and the case `T3` exists to
close reverts to **silence** — ⛔ which `T5` names unacceptable.

⭐ **State the SHAPE plainly:** both engines run and emit natively; the **orchestrator** aligns candidate
components across engines and supplies each engine the counterpart point; each engine is then **run again**
and emits its witness objects. ⛔ Name the orchestrator as the agent. ⛔ No engine reads the other's output.

⭐ Two gaps inside the same block, both concrete:
- ⛔ `WITNESS<w>_SOURCE_LOCUS_TAGS` carries `ROOT<r>_…` identities whose `<r>` indexes an **engine-local**
  `ROOT_ORDERING` (`§8:1039-1041`) ⇒ ⚠ something must translate root indices across engines. Say what.
- ⛔ Nothing requires the received point to cover the **receiver's** symbols ⇒ ⚠ a partially symbolic
  `WITNESS<w>_OWN_M_EVALUATED` with silently generic counts and **no token for it.**

## ⛔ R3 · The count tag's own payload under `VARIES` is unspecified — ⚠ this is the `C17` case

⭐ `:702-704` pins the payload for *cannot-compute* (`NOT_COMPUTED_COMPONENT_WIDE`); ⛔ `:683-698` pins
nothing for `VARIES`, while forbidding a bare integer. ⇒ ⚠ three builders, three **types** for
`STRATUM<s>_ROOT_COUNT_DISTINCT`. ⭐ Give `<COUNT>` **one type in all three statuses.**

## ⛔ R4 · The component's PARAMETERISATION is unpinned — ⚠ measured, ⛔ not hypothetical

⛔ `:677-681` says restrict and leave the free parameters free, but ⛔ never says **which variable is
eliminated**. ⚠ **Measured in the committed outputs, same component:**
`WL_S11_MAIN_D3_ROOT_COINCIDENCE_COEFF_SOLUTION: … {{muR -> bComp}}` vs
`PY_S11_MAIN_D3_ROOT_COINCIDENCE_COEFF_SOLUTION: … ({B_comp: mu_R},)`.
⇒ ⛔ every component-scoped **symbolic** payload differs for a reason that is not physics. ⭐ Counts are
unaffected. ⭐ **Pin the elimination by a rule computable from the objects**, or state that the
component-scoped symbolic payloads are **inspection-only**, as `:420-422` already does for `N6_BASIS`.

## ⛔ R5 · `_CANONICAL_LOCUS` is inert exactly where it is needed

⚠ **Measured in BOTH CASes** on `MAIN`'s coincidence residual
`((k1²+k2²+k3²)(bComp−muR))/rhoBr`: not polynomial in all symbols ⇒ `NOT_APPLICABLE` — ⛔ **only because of
a premise-positive denominator.** ⭐ `ρ_br > 0` is a `§3` premise, so clearing it is equivalence-preserving
on the admissible region, and the **numerator is polynomial in both CASes.**
⭐ Apply the predicate to the numerator over premise-positive denominators, ⛔ or state that the object is
inert wherever such a denominator appears.

## ⛔ R6 · *"whose payload is an integer count"* is a VALUE-INSPECTED predicate

⛔ `:683-684` defines the companion families by inspecting the payload, with a **non-exhaustive** list.
⚠ Genuinely ambiguous: `ROOT_DEGREE_RESIDUAL`, `N4_NULLITY_DIFFERENCE`, `N7_RESIDUAL` — each an integer,
each a **difference**. ⇒ ⛔ the two engines emit **different tag sets**, against `§8:1025`'s *"both engines
must produce PARALLEL tag sets."* ⭐ **List the qualifying quantity names**, as `§8` does everywhere else.

## ⛔ R7 · Q10 still speaks the pre-repair dialect — ⚠ found by both legs AND the orchestrator

⛔ `:838` `RECOMPUTED_ROOTS: <the stratum's live recomputed ROOT_DISTINCT object>` — ⚠ there are now **two**
such objects (component-scoped and `POINT_EVIDENCE_`), and this pinned field does not say which.
⛔ `:822` `STRATUM<s>_ROOT_COEFFICIENT_JACOBIAN_RESTRICTED` is *"evaluated at `STRATUM<s>_POINT`"* while
carrying the **bare component scope** ⇒ ⚠ the `C17` shape, in a second place `T2` did not reach.

## ⛔ R8 · The new "never a native boolean" sentence discriminates nothing

⚠ **Measured:** `str(sympy.S.false)` = `'False'`, `srepr(sympy.S.false)` = `'false'`, and `§8:1056-1057`
permits **either**. ⇒ ⛔ `False`/`False` fires the comparator defect; ⛔ `False`/`false` reports a mismatch
that is not one. ⭐ The real protection is `STATUS_TOKEN`. ⭐ **Pin the `TEST_OBJECT` serialisation for
boolean-valued objects, ⛔ or state that `STATUS_TOKEN` is the comparison field and `TEST_OBJECT` is
inspection-only.**

---

## ⛔ Constraints — unchanged from the build directive

- ⛔⛔ **The spec says what to compute. ⛔ NEVER what anything equals, is expected to be, or was measured.**
  ⚠ `R4` and `R5` above quote **measured payloads**; ⛔ **none of that may enter the file.** ⭐ The test is
  **derivability**, ⛔ not literal presence.
- ⛔ No admissibility algorithm. ⛔ No recursive stratification. ⛔ No engine imports, agrees with, or
  adjusts toward the other — ⚠ `R1`/`R2` make a point an **input**; ⛔ they do not relax that.
- ⛔ Do not delete or weaken any `§5` locus object. ⛔ Do not renumber sections.
- ⚠ **After editing, re-read each edited section AND its introduction** for a sentence the edit made false.

## Deliverables

1. The edited file.
2. ⭐ Per `R1`–`R8`: the lines changed and the sentence that now makes it true. ⛔ If one needs no change,
   say so and why.
3. ⛔ Any conflict between two items, or an item that cannot be met without stating a value ⇒ **STOP and
   report.** ⛔ Do not invent a compromise.

⛔ Do not commit. ⛔ Do not edit any other file. ⛔ Do not run either engine.
