# Independent physics review — the S11 shared-spec amendment

## Artifact

The **working-tree** `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (1159 lines),
against its committed baseline `git show cf4a21a4:research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`
(1149 lines). ⭐ **The amendment is the diff: 13 lines added, 3 changed.** Codex-authored.

⚠⚠ **This file is read by BOTH engines. An error in it makes both engines agree on the same wrong thing**
(rule 7). ⭐ It took **nine rounds** to close, and it was reopened for **three items only**:

1. name the shared objects the file **orders but never named** — `§Q1`'s expanded `L`, `§Q1`'s
   Euler–Lagrange system, `§Q2`'s two-route difference and entry ratio;
2. state, where **both** engines read it, the rule that an ordered family named by one tag stays **one
   ordered payload**;
3. require a boolean payload to be rendered so a **live CAS boolean** is distinguishable from a
   **host-language native boolean**.

## Read in this order — ⛔ load-bearing

1. `/var/projects/toy_physics/CLAUDE.md`.
2. `git show cf4a21a4:research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — **the baseline**. ⭐ Form
   your own view of what these three items would have to say to work.
3. `directives/_legs/round5_grok_build_directive.txt` and `directives/_legs/round5_claude_build_directive.md`
   — the two legs whose findings drove the amendment. ⚠ ⛔ **`B2` in the Claude leg is RETRACTED by the leg
   itself.**
4. ⭐ **Only now** read the working-tree file and `git diff` it.

## Do not read

- The other leg's output.

## What to check

### ⭐⭐ 1 · IS THE DIFF MINIMAL, AND DOES IT ADD PHYSICS?

- ⭐ Verify `§4`'s quoted structural rule is **byte-identical** to the baseline and to
  `S10_SHARED_PHYSICS.md`'s corresponding block. ⛔ Report the hashes.
- ⛔ Is any changed line **not** traceable to one of the three items? ⭐ Name it.
- ⛔⛔ **Does any added line state what a computation will produce** — a value, count, membership, sign,
  spectrum, dimension or outcome? ⚠ Rule 5, and two artifacts leaked here this session.
- ⛔ Does any added line **narrow what either engine may compute**, as opposed to naming what it emits?
  ⚠ Rule 6 — thirteen of sixteen bred defects in this file's history lived in machinery invented to stop
  two engines describing something differently.

### ⭐⭐ 2 · ITEM 1 — THE NAMES

- ⭐ Walk **every** "Emit" instruction in the amended file. Is the set of shared objects that are **ordered
  but unnamed** now empty? ⛔ Name any that remain.
- ⛔ Does any new name **collide** with an existing `<QUANTITY>` in this file, or with a name used for a
  different object anywhere the engines read? ⭐ Check, do not assume.
- ⛔ Does naming them break any existing cross-reference in the file?

### ⭐⭐ 3 · ITEM 2 — DOES THE DECOMPOSITION RULE ACTUALLY BITE?

⚠ **The measured catastrophe it exists to prevent:** the two engines being replaced were built from
directives that decomposed the work differently — one bundled a root's value, nullity and orientation into
a single payload where the other emitted three — and **shared one tag suffix between them**.

- ⭐⭐ **Construct that case against the amended rule.** Would a builder following it produce the bundled
  decomposition or the split one? Show the two tag sets and whether the rule separates them.
- ⛔ Is the rule **decidable** for the objects the file actually orders? ⭐ Walk the hardest ones — `§5`'s
  locus protocol families, `§Q8b`'s stratum-scoped counts, `§Q9`'s census, `§Q11`'s `C1`–`C4` — and say for
  each whether the rule determines one tag or several. ⛔ Any object where two builders could still differ
  is a finding.

### ⭐⭐ 4 · ITEM 3 — MEASURE WHAT EACH ENGINE ACTUALLY EMITS

⛔ **Do not reason about this. Run it.**

- ⭐ Take each boolean-valued object the rule names. Write **both** a compliant implementation and one that
  substitutes a host-language boolean for the live CAS decision, render each under the rule, and report the
  **literal** output of both. ⛔ Is the forbidden one distinguishable in every case, or only in some?
- ⛔⛔ **Try to defeat the rule.** Is there an implementation that emits a host boolean and still renders
  identically to the compliant one? ⭐ If there is, the rule is weaker than it reads.
- ⭐ **Then look across the two engines.** The rule pins a rendering for one CAS and records the other as
  unresolved. ⛔ **What does the cross-engine comparator see on a boolean row when both engines are
  correct?** ⚠ Is that a physics disagreement, a vocabulary difference, or indistinguishable from one?
- ⭐ Is recording the Wolfram half as **unresolved** the right disposition, or does an unresolved half make
  the pinned half harmful? ⛔ Argue it with the measurement.

### ⭐ 5 · WHAT THE AMENDMENT DID NOT FIX

⭐ The two legs raised findings that were routed to the **build directive** rather than here. ⛔ Check that
routing: is any of them in fact a **shared** obligation that only this file can carry? ⚠ Name it and say
why the other level cannot establish it.

### ⭐ 6 · WHAT IS MISSING

⛔ With these three items in place, what could still make the two engines produce **uncomparable** tag sets
for the same physics? ⚠ Name it concretely, with the mechanism.

## ⛔ Constraints on how you run

- ⛔ Read-only on the working tree; copy to `/tmp` to modify. ⛔ `timeout 600` on every run, never raised.
- ⛔ **No Mathematica** — the licence has two seats. ⭐ Reason about Wolfram rendering from the
  documentation and from the committed `.wl` outputs under `mathematica/out/`, ⛔ not by running a kernel.
- ⭐ Save every script and its literal stdout to named absolute paths and report them. ⛔ A prose derivation
  will be discarded.

## Report

For each finding: **what is wrong**, **the mechanism by which a wrong result would survive**, **the literal
output that shows it**, and **what must become true instead** — ⛔ not a rewrite.
⭐ Separate blocking from non-blocking. ⛔ "Nothing found" with no computation behind it is not a result.
