# ⭐⭐ AUTHOR a MINIMAL amendment to the shared spec — ⛔ three items, ⛔ no physics

⚠ **You are the author.** `CLAUDE.md` **rule 15** is standing for this material. ⛔ You will not review it —
a fresh Claude agent and Grok will.

## ⛔⛔ What you are touching, and why that is serious

`research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (**1149 lines, `cf4a21a4`**) is the artifact
**both engines read**. ⚠⚠ **An error in it makes both engines agree on the same wrong thing** — rule 7.
⭐ It took **nine rounds** to close. The user has authorised reopening it for **these three items only**.

⛔⛔ **MINIMAL DIFF.** ⛔ Do not restructure, ⛔ do not improve prose, ⛔ do not revisit a closed question,
⛔ do not touch `§4`'s quoted structural rule (`:156-158`, byte-identical to S10's and verified every round).
⭐ Every line you change must be traceable to one of the three items.

## ⭐ Why these three, and ⛔ why they cannot be fixed anywhere else

⭐ Two independent legs blocked the S11 SymPy build directive. ⭐ Full reports:
`directives/_legs/round5_grok_build_directive.txt` and `directives/_legs/round5_claude_build_directive.md`.
⚠ ⛔ **`B2` in the Claude leg is RETRACTED by the leg itself — ⛔ do not act on it.**

⭐ Each item below is a **joint property of two engines**. ⛔ A single-engine directive cannot establish one,
and the Wolfram builder will never read the SymPy directive. ⇒ ⭐ they belong here or nowhere.

**Item 1 · The shared objects this file ORDERS but never NAMES.**
⭐ `§Q1` says *"Emit `L` expanded"* and *"Emit the full system"*; `§Q2` orders the two-route difference and
the entry ratio. ⛔ None has a `<QUANTITY>` spelling anywhere in this file, while `M_A`, `M_B`,
`KINETIC_TERMS` and `STIFFNESS_TERMS` do. ⚠ A leg walked **every** "Emit" instruction and measured that
these are **exactly** the shared objects left unnamed.
⇒ ⛔ **Measured consequence of leaving them unnamed:** with the relative-weight error `§3:103-106` warns
about injected, matching names give **1 comparison row, residual nonzero**; differing names give
**0 comparison rows and two orphan tags** ⇒ the physics error is **never compared** and both engines exit 0.
⭐ Name them. ⛔ Do not name anything else.

**Item 2 · The decomposition rule, stated where both engines read it.**
⚠ **The previously measured catastrophe:** the two engines being replaced were built from directives that
**decomposed the work differently** — one bundled a root's value, nullity and orientation into a single
payload where the other emitted three — and **shared one tag suffix between them.** `§8:1074-1084` states
one-tag-per-named-object; ⛔ a leg found the *positive* rule that makes it checkable — an ordered family
named by one tag stays **one ordered payload**, and only scopes this file **explicitly indexes** become
separate tags — stated only in the SymPy directive. ⇒ ⭐ state it here.

**Item 3 · A boolean payload must be rendered so a CAS boolean is distinguishable from a host boolean.**
⛔⛔ **Measured: the compliant and the forbidden implementation are BYTE-IDENTICAL** under the `str()`
rendering `§8:1086-1087` currently permits, and differ under `srepr`. ⇒ ⭐ this is `§4:162`'s own test
(*"if you deleted the computation, would this tag change?"*) failing on the objects `§5` calls load-bearing
— `_IDENTICALLY_SATISFIED`, `_INCONSISTENT`, `_REAL_ADMISSIBLE`'s test object, the homogeneity booleans.

⚠⚠ **THE WOLFRAM HALF IS AN OPEN QUESTION AND MAY NOT BE ASSUMED.** `InputForm[True]` is `True` whether it
came from a symbolic decision or a host literal. ⇒ ⛔ **Do not invent a Wolfram mechanism.** ⭐ Either state
a requirement **both** CASes can meet, ⭐ or state the requirement and record explicitly that the Wolfram
side is unresolved. ⚠ **Reporting is success here** — every mechanism invented to close a gap in this
workstream has been blocked by two legs.

## ⛔ Constraints

1. **`CLAUDE.md`** in full.
2. ⛔⛔ **RULE 5 — ⛔ THE SPEC SAYS WHAT TO COMPUTE, ⛔ NEVER WHAT ANYTHING EQUALS, IS EXPECTED, OR WAS
   MEASURED.** ⛔ No value, count, membership, sign, spectrum or outcome. ⚠ Two artifacts leaked here this
   session, the second one round after the first was caught.
3. ⭐ **Rule 3** — name the object; ⛔ do not specify a derivation path. ⚠ Item 3 is a **rendering** rule, ⛔
   not a recipe for deciding a boolean.
4. ⭐ **Rule 6** — ⛔ do not make divergence impossible. ⚠ Thirteen of sixteen bred defects in this file's
   history lived in machinery invented to stop two engines describing something differently. ⭐ These three
   items give the engines a **shared vocabulary**; ⛔ they must not narrow what either engine may compute.
5. ⭐ **PHYSICS, ⛔ NOT PROCESS** (standing user instruction).

## Report back — ⛔ under 12 lines

1. The diff: lines added, removed, changed, and the new total. ⭐ Confirm `§4:156-158` is untouched.
2. ⭐ For each of the three items, the section you changed.
3. ⭐ What you left **open** — ⭐ item 3's Wolfram half in particular — and why.
4. ⛔ Do not summarise the physics; ⛔ do not claim the amendment is correct.

⛔ Edit only `S11_SHARED_PHYSICS.md`. ⛔ Read-only otherwise.
