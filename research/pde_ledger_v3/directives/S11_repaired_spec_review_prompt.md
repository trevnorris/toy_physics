# Independent review — the repaired S11 specification, before either engine is built

## Artifact under review

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — **999 lines**, as it
now stands in the working tree. It was **914 lines** at `HEAD`; the repair is uncommitted, so
`git diff` shows exactly what changed.

It is **Codex-written**, applied against a decision list that had already passed two review legs. You are
one of two independent legs; the other is not visible to you.

⚠⚠ **This file is the shared physics for TWO engines that are rebuilt from it independently — one SymPy,
one Wolfram — and then compared.** ⇒ **An error here makes both engines agree on the same wrong thing**,
and a comparison cannot detect it. That is why the spec is reviewed before either engine exists.

## Required reading

1. `research/pde_ledger_v3/directives/S11_spec_repair_decisions_v2.md` — the decision list the repair was
   applied against. **This is the contract.**
2. `git diff research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — what actually changed.
3. The repaired file **in full**. ⛔ Reviewing only the diff will miss the most important failure mode
   below.

## What to check

### 1. ⭐⭐ IS EVERY ITEM APPLIED, AND APPLIED WHERE IT HAD TO BE?

Walk items 1–5 of the decision list. For each: is it in the file, and is it in **every** place it had to
reach? ⚠ Several items change a definition **and** the rules that consume it — a package now has two
members, a coefficient inventory clause was singular, a classification rule read only one member of the
package. ⭐ **An item applied in one place and missed in another is the characteristic failure**, and it is
invisible if you read only the diff hunks.

### 2. ⛔⛔ DID THE REPAIR BREAK SOMETHING IT DID NOT TOUCH? — the most important question

⚠ **The repair inserted 124 lines into a file whose sections cross-reference each other heavily.**
⭐ **Read the whole file and find every passage that is now wrong, stale, unreachable, or contradicted** —
including passages the diff does not show, because the sentence that broke may be a hundred lines from the
edit that broke it. Specifically worth tracing: anything that names a package, the action, the coefficient
inventory, the kinetic or stiffness terms, tag scopes, or the locus protocol.

### 3. ⛔ DOES THE FILE NOW STATE, OR IMPLY, AN ANSWER?

⛔ This specification may **never** state what anything comes out to be — no value, count, sign, rank,
dimension, or expected effect of any package.
⚠⚠ **The test is DERIVABILITY, not literal presence.** A sentence that prints no number but from which a
builder can derive one is a leak. ⭐ This exact defect was caught in the decision list: a *justification*
for a modelling choice let the reader eliminate two equations and fix a dimension vector. **Look for
justifications, motivations, and "because" clauses**, not just numerals.
⚠ `DEFECT_REGISTER.md`'s `C16`, `C17`, `C18` entries contain **measured values** — root lists, mode counts,
stratum equations. ⭐ Check that **none** of them, and nothing from which they follow, reached this file.

### 4. ⛔ DID THE REPAIR DO MORE THAN THE LIST AUTHORISED?

⭐ The decision list has a *"What this round does not do"* section and a *"Registered separately"* section.
Those are **deliberate exclusions**. ⚠ Did the repair implement any of them anyway, or make any change the
list did not ask for? ⛔ Name it, and name what it damages. A repair that grows past the defect it names
breeds defects in the material it changes.

### 5. ⭐⭐ IS THE RESULT BUILDABLE BY TWO ENGINES THAT NEVER SPEAK?

⭐ For each repaired passage: if a SymPy builder and a Wolfram builder both implement it faithfully, can
they produce **different objects**? Name the specific freedom.
⚠ Pay particular attention to anything that pins a **payload shape**, a **tag name**, a **term
decomposition**, a **count**, or a **name mapping** — those exist precisely because an earlier wording let
two engines diverge. ⭐ Check that what the list pinned was transcribed exactly, not paraphrased into
something looser.
⚠ The file's own `§8` requires **parallel tag sets**. Does every newly named object have a name both
engines will emit identically?

### 6. ⭐ IS ANYTHING NOW UNSATISFIABLE, OR SATISFIABLE ONLY BY INVENTING?

⚠ Is any requirement impossible to meet from what `§1`–`§3` supply? ⚠ Does any requirement force a builder
to manufacture an object merely to make a declaration appear wired? ⭐ One item deliberately narrows an
exemption to a **single field** of a structured object whose other fields must stay live-read — check that
the narrowing is exact in both directions.

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, or a way this file would cause a
**wrong engine, or two divergent engines**, to be built. ⛔ Not style, not formatting, not "a builder might
misread this" absent a concrete reading that produces a wrong build.

## Method

- ⭐⭐ **Quote both sides for every finding**: the repaired file's text, and the decision-list or spec text
  it fails against. A finding without both quotations is not usable.
- ⭐ For any claim about a file, a line, or a count: **state what you opened or ran, and what you found.**
  ⛔ A path that resolves is not a source — open it.
- ⭐ Where a question is settled by computation, **write a script, run it, and paste its literal stdout.**
  ⛔ A prose derivation is not accepted for a contested claim.
- ⛔ Do **not** edit any file. Read-only. ⛔ Do not write engine code.
- ⭐ State explicitly **what you checked that could have failed and did not.**
