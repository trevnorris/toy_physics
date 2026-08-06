# Loosen `S11_SHARED_PHYSICS.md` — ⛔ do NOT try to complete it

**Do not commit.** One file: `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`.

---

## ⭐⭐⭐ THE GOVERNING INSTRUCTION — read this before any edit

> **Prefer REMOVING or LOOSENING a requirement to adding one.** Where a review finding can be answered
> either by tightening the file or by deleting the thing that made it tightenable, ⭐ **delete.**

**Why, and it is measured, ⛔ not a preference.** This file has been through two review rounds. Round one
returned thirteen findings, all genuine design problems, and folding them closed all four of the ones both
legs agreed on. Round two then returned seven **new** findings — and almost every one was created by the
fold itself, in material that had just been rewritten to be correct. The pattern is that each attempt to
make the file airtight introduces the next round's defects.

⇒ ⭐ **The engines are the better instrument.** §10 already requires the builder to report any named object
it could not emit. A build surfaces a gap concretely in one run; a review can only predict one.

⛔ **Therefore: do not attempt to make this specification complete, closed, or self-consistent in every
clause.** ⛔ Do not add new tags, new rules, new guards, or new checks beyond the five changes below.
⛔ Do not "improve" a section you were not asked to touch.

⚠ **A signal, ⛔ not a gate:** if your diff makes the file substantially *longer*, you have probably added
where you should have cut. Report the net line change.

---

## R1 · ⭐⭐ THE ONE REAL FIX — pin the monomial ordering

**This is the only change here that affects physics.** §Q9 requires each basis emitted in reduced row
echelon form, and §7 builds one package's action from that basis. ⛔ **RREF is canonical only relative to a
coordinate ordering, and the file pins none.**

Measured by a review leg: at `D = 2`, two orderings of the quadratic monomials give `P_D` and `−P_D`.
⚠ That is a **relative sign between two stiffness terms**, ⛔ not an overall scale of `L` — so the two
engines would build **genuinely different actions** for that package, and `β`'s free sign does not absorb
it (`β → −β` maps the family, not the expression).

⭐ **Fix it by pinning one explicit ordering in the file**, stated so both engines construct the identical
coordinate list. ⛔ Do not leave it to an engine's "own stated rule". ⭐ Choose the ordering yourself; any
deterministic, fully specified one will do. ⛔ State no property of the resulting basis.

## R2 · LOOSEN — §8's closed-vocabulary mandate

§8 asserts that `<QUANTITY>` is **exactly** the set of code-font names in §6 and §7. ⛔ **That claim is
false in its own file** — a review enumerated 114 name patterns and found required emissions with no name:
the Euler–Lagrange system, `M_A − M_B`, the entry ratio, the coding-consistency tag, every premise tag,
§Q6r's registry objects, and the `_LOCAL_` inventory tag.

⛔ **Do NOT fix this by naming all of them.** ⭐ **Delete the closure claim.** What the architecture
actually needs is a shared **decomposition**, ⛔ not shared spelling — the cross-engine pairing is declared
downstream in the harness configuration, which can map two spellings onto one row. What it cannot repair is
one engine bundling two objects into a single payload where the other emits them separately.

⇒ Replace the mandate with the three requirements that carry the weight:

- emit **one tag per named object**, at the scope §8 gives it;
- ⛔ **do not bundle two named objects into one payload, and do not split one across several** — ⚠ that,
  and not naming, is what left the two engines being replaced sharing a single tag suffix;
- where you emit something this file names no name for, choose one and **list it in your §10 report**.

## R3 · REMOVE — `PD_BASIS_DEPENDENT`

Delete the tag and the sentence introducing it. It fires only when `V6_DIM` exceeds one, which is emission
conditional on a computed value and violates the file's own corollary 4. ⭐ Once `R1` pins the ordering it
has nothing left to check.

## R4 · LOOSEN — §Q8's "irreducible component"

A review leg measured that primary decomposition is not stable across CAS implementations, so two engines
need not return the same components or the same representative points.

⭐ Delete the word **irreducible** and the requirement that the enumeration be one. ⭐ Keep: emit the
components your CAS returns, each with its defining equations and the chosen point, plus the ordering. ⭐ Add
one sentence saying stratum indices are aligned by the orchestrator from the emitted defining equations,
⛔ never assumed to correspond.

## R5 · FIX — two small scope and wording errors

- `KW_ZERO_LOCUS_*` sits inside §Q11's per-root block but carries no `ROOT<r>` scope, so per-root emission
  collides on one tag name. Give it `ROOT<r>` scope.
- §Q3 calls the sign test **three-way** and then lists **four** tokens. Make the wording match the token
  list; ⛔ do not change the token list.

## R6 · REMOVE — the bare integers in §Q9's cost note

The note quotes `10` and `325`. They are `dim Quad(D)`, but they sit beside the census and a leg flagged
that they read as expected output. ⭐ Delete the numbers; keep the statement that the solve grows with `D`
and that narrowing the sweep and substituting a floating-point solve are both forbidden.

---

## ⛔ Do not touch

The three clauses and five corollaries · the no-verdict rule · §5's five-object locus protocol · the pinned
sign and marker token sets · §Q7's five tags · §Q11's bulk field, ansatz, and the ban on the
back-substitution residual · the package list and the `D` sweep · §Q6d's vacuity explanation, which is
inherited deliberately from the previous step · §7's runtime rule.

⚠ §Q6d's explanation has been flagged as a leak by a reviewer and that finding was **considered and
rejected** — it states a general fact about dimensional analysis, and removing it leaves a check that looks
meaningful and is not. ⛔ Do not remove it.

## Acceptance

1. `git diff --stat` on the file, and the **net line change**.
2. `grep -n` showing the pinned ordering from `R1`, quoted in full.
3. `grep -n "PD_BASIS_DEPENDENT\|irreducible"` — expect no hits.
4. `grep -n "KW_ZERO_LOCUS"` with its surrounding scope line.
5. The paragraph you wrote for `R2`, quoted in full.

## Constraints

- ⛔ No `git add`, no `git commit`, no other git write.
- ⛔ Touch no other file. ⛔ Do not create a new file.
- ⛔ Do not run either engine; neither exists yet.

## Report back — under 15 lines

1. `R1`–`R6`: done or not, and the net line change.
2. The acceptance output.
3. ⭐ Anything you removed that you think the file needed, and anything in `R1`–`R6` you judge wrong.
   ⛔ Do not list further defects you noticed and were not asked to fix — ⭐ note only that they exist and
   how many.
