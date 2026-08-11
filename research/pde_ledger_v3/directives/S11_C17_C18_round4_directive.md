# Build directive — four measured items on `S11_SHARED_PHYSICS.md`

⭐ **Target:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (1199
lines). ⛔ **Edit ONLY that file.** ⛔ Do not commit. ⛔ Do not run either engine.

⚠⚠ **This file has had FOUR revisions and every one bred defects — thirteen in total.** ⭐ **Every single
one lived in a mechanism its author INVENTED to make two engines describe a locus the same way.**
⇒ ⛔⛔ **DO NOT INVENT A MECHANISM. Each of the four items below is already settled by a measurement, and
the measurement is given with it.** ⭐ Your job is to write what the measurement says, ⛔ not to design.

⚠ If you find yourself reasoning about whether some construction is well-defined, **stop and report it** —
that is the exact failure mode this round exists to avoid.

---

## ⭐ M1 · `POINT_COVERAGE` is under-scoped — ⛔ it can say COMPLETE while residuals are symbolic

**Where:** `§Q8c`'s witness block, the `WITNESS<w>_POINT_COVERAGE` line.

**Now reads:** coverage is tested over *"every symbol this engine's own `M` depends on"*.

⛔ **But `WITNESS<w>_OWN_LOCUS_RESIDUALS` evaluates EVERY locus this engine emitted**, and a locus can
depend on symbols `M` does not.

⭐ **MEASURED** — a review leg supplied the counterpart point to a reconstruction of the SymPy engine's own
`XFORM_EXTRA, D = 2` loci and got, literally:

```
POINT_COVERAGE_BY_SPEC COMPLETE_POINT
ZERO_LOCUS_COUNT 8
ROOT1_KW_ZERO_LOCUS = -2*(c_s0 - 1)*(c_s0 + 1)/c_s0**2
ROOT2_KW_ZERO_LOCUS = -2*(c_s0 - 1)*(c_s0 + 1)/c_s0**2
SYMBOLIC_LOCUS_COUNT 2
```

⇒ `COMPLETE_POINT` with **two** residuals still symbolic in `c_s0`.

⭐ **What must become true:** coverage is decided over the **union of every symbol appearing in everything
that witness slot actually evaluates** — every emitted locus's equations **and** `M`. ⛔ Not `M` alone.
⚠ An alternative that also satisfies this is a coverage token **per locus**; ⭐ either is acceptable,
⛔ pick one and be consistent.

---

## ⭐⭐ M2 · Restore the symbolic component comparison — ⛔ as a REDUCTION, ⛔ never a normalisation

**Where:** `§Q8b`, the paragraph beginning *"THE TWO ENGINES MAY DESCRIBE ONE COMPONENT IN DIFFERENT
VARIABLES"*.

**Now reads:** the component-scoped **symbolic** payloads are *"INSPECTION-ONLY … not cross-engine
comparison rows"*, and only counts and statuses are compared.

⛔⛔ **That gives up the only check that catches a wrong dispersion coefficient on a stratum.**

⭐ **MEASURED.** On the component `B_comp = mu_R`, reducing the difference modulo the defining equations:

| operands | reduced difference |
|---|---|
| two engines in **different charts** (`B_comp·k²/ρ` vs `mu_R·k²/ρ`) | **`0`** |
| a dispersion relation **off by a factor of two** | **`−k²·mu_R`** |

⭐ And the **counts are identical in both cases** — `distinct roots 1 · rank 0 · nullity 3`.
⇒ ⛔ counts cannot see the factor-of-two error; ⭐ the reduction can.

⭐ **What must become true:** a component-scoped symbolic payload **is** a cross-engine comparison row, and
what is compared is **the difference of the two payloads REDUCED MODULO `STRATUM<s>_DEFINING_EQUATIONS`**.

⛔⛔ **CRITICAL — this must not become a normalisation:**
- ⛔ **Neither engine is told which variable to eliminate.** ⚠ A previous round pinned that and it **deletes
  a branch**: on `x·y = 0` the pinned rule returns `y = 0` and loses `x = 0` with `y` free.
- ⭐ Each engine keeps its own chart and emits its own payload, unchanged.
- ⭐ **The reduction is the comparison's operation, ⛔ not an engine's obligation.** Say so plainly.
- ⭐ Where the reduction does not settle, the outcome is an explicit **undecided**, ⛔ never a silent pass.

---

## ⭐ M3 · The point cannot be consumed — the two engines spell its coordinates differently

**Where:** `§4`'s witness-input block and/or `§Q8c` step 2.

⭐ **MEASURED, from the committed outputs:** the same point is emitted as
`rhoBr`, `muR`, `bComp` (Wolfram, a rules list) and `rho_br`, `mu_R`, `B_comp` (SymPy, a `Dict`). ⇒ ⛔ a
receiver handed the raw counterpart point does not recognise the symbols, and every witness would come back
`INCOMPLETE_POINT` for a reason that is **not physics**.

⭐ **What must become true:** the spec says the **orchestrator** delivers the point in the **receiver's own
coordinate names**, and names the map it uses. ⚠ ⭐ **A closed map already exists in this file at `§Q6r`** —
⭐ use it; ⛔ do not author a second one. ⛔ If it does not cover every symbol a point can carry, **report
that** rather than extending it.

---

## ⭐ M4 · One stale sentence — ⛔ found by BOTH review legs

**Where:** `§4`, in the paragraph after the witness-input block.

**Now reads:** *"it names **these two fields**, in this pass, and nothing else."*

⛔ The input became **one field** (the point) in the last revision; ⚠ this sentence still counts two.
⭐ **What must become true:** it names **one** field. ⚠ A builder reading "two" can reintroduce the removed
selector, which undoes the whole design.

---

## ⛔ Constraints

- ⛔⛔ **THE SPEC SAYS WHAT TO COMPUTE. ⛔ NEVER what anything equals, is expected to be, or was measured.**
  ⚠⚠ **The measurements above are for YOU, ⛔ not for the file.** ⛔ No count, value, sign, coefficient,
  named locus or example from them may appear in the spec, in any form a builder could **derive** them
  from. ⭐ The test is **derivability**, ⛔ not literal presence.
- ⛔ Do not renumber sections. ⛔ Do not reflow untouched prose. ⭐ Match the file's existing voice.
- ⛔⛔ **`§4`'s quoted block is shared VERBATIM with the other steps' specs — ⛔ do not edit inside it.**
  ⭐ Verify: `sed -n '156,158p'` of the target must stay byte-identical to `sed -n '111,113p'` of
  `directives/S10_SHARED_PHYSICS.md`.
- ⛔ Do not touch `_CANONICAL_LOCUS`. ⚠ It is deliberately inert on non-polynomial systems; ⭐ a previous
  round made it fire and it **divided the critical component out of its own canonical object**.
- ⛔ Do not delete or weaken any `§5` locus object.
- ⛔ Do not add an admissibility algorithm, recursive stratification, or any rule pinning how an engine
  presents a component.
- ⚠⚠ **AFTER EDITING, RE-READ EACH EDITED SECTION AND ITS INTRODUCTION**, and grep the whole file for
  sentences referring to anything you changed. ⭐ **All four previous rounds left exactly one stale
  sentence; `M4` is the current one.** ⛔ Do not make the fifth.

## Deliverables

1. The edited file.
2. ⭐ Per `M1`–`M4`: the lines changed and the sentence that now makes it true.
3. ⭐ The result of the `§4` byte-identity check above.
4. ⛔ Anything you could not do without inventing a construction, or any conflict between two items ⇒
   **STOP and report it.** ⛔ Do not invent a compromise.
