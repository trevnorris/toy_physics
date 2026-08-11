# Independent physics review — the S11 SymPy export decision list

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_sympy_export_decisions_codex.md` — 101
lines, **Codex-authored**.

⚠⚠ **This is the list a builder will be handed, and it is the one artifact in this pipeline checked once
and then trusted.** Everything downstream of it is reviewed twice. ⭐ Three earlier attempts at this same
material were each blocked by two independent legs; ⛔ do not assume this one is closer to right because it
is shorter.

## Read in this order — ⛔ the order is load-bearing

1. `/var/projects/toy_physics/CLAUDE.md`.
2. `git show cf4a21a4:research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — the physics spec.
3. `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_exports.py` — the committed artifact this
   engine imports — and `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py`,
   which wrote it.
4. `/var/projects/toy_physics/research/pde_ledger_v3/S9_REWRITE_PLAN.md` — `D1`–`D12` of the plan.
   ⚠ ⛔ Note the name clash: the plan's `D`-numbers are **not** the artifact's `D`-numbers.
5. `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` —
   `F1`–`F7`. ⚠ `F8` is a **blocked record**.
6. ⭐ **Only now** open the artifact.

⛔ Do not open the artifact first. Form your own view of what the chain requires, then read what was
decided.

## Do not read

- The other leg's output.

## ⭐⭐ What this artifact changed, so you can attack the right thing

It **supersedes `F1`** — previously settled and two-legged — which said storage keys are flat and that a
later step re-deriving an object writes **the same key**, because *the collision is the measurement*.

In its place: a key is an **opaque locator** with no meaning; identity moves into a typed per-row
`object_claim`; the merge joins on **claims** rather than keys. Nothing is renamed.

⭐ **The three prior attempts all tried to make a key string carry identity, and all three were blocked.**
Their leg reports are in `directives/_legs/` — ⭐ read them, so you do not spend this round re-finding what
is already known.

## What to check

### ⭐⭐ 1 · WHO WRITES A CLAIM, AND WHAT CATCHES A WRONG ONE?

`D2` makes `object_claim` the load-bearing object in the whole design: it decides what meets what. It is
composed per object from the specification's named role, its scope, and references to the claims of its
defining inputs and premises.

- ⛔ **A claim is authored, not computed.** How many claims does a full `MAIN` sweep require? Count them
  from the spec and the declared sweep, and report the number.
- ⛔⛔ **What in this list detects a WRONG claim** — one that references the wrong defining inputs, omits a
  premise the object depends on, or gives two genuinely different objects the same claim tree? `D8` checks
  the **population**; does anything check **correctness**? ⭐ If nothing does, state the mechanism by which
  wrong physics then survives.
- ⚠ Compare against what it replaced. A hand-authored identity is what `D11` of the plan calls out as
  having cost 3,121 lines of hand-written name pairs. ⭐ **Is `object_claim` that defect returning in a new
  shape, or is it genuinely different?** Argue it either way, with the mechanism.

### ⭐⭐ 2 · DOES THE CLAIM JOIN ACTUALLY SEPARATE AND UNITE THE RIGHT OBJECTS?

⭐ Construct both cases against the committed artifact and report literal output:

- Two objects of the **same kind** from **different actions** — they must land in different rows, and
  neither may overwrite the other.
- Two genuinely **identical** objects derived by different steps — they must meet and be compared. ⭐ Find
  such a pair in the committed data yourself; ⛔ do not take one from a leg report.
- ⛔ **Find a pair the join gets wrong** in either direction. A false meeting or a lost meeting is the
  highest-value finding available in this round.

### ⭐⭐ 3 · WHAT THE SUPERSESSION OF `F1` COSTS

`F1`'s stated reason was that two steps deriving one object must be able to **meet**, and that *the
collision is the measurement*.

- ⭐ Under `D1`–`D3`, is that measurement preserved, moved, or lost? Show it.
- ⛔ `D6` makes imported physics **append-only**: a verified re-derivation appends a comparison record and
  does **not** replace the imported value. ⚠ The plan's `D5` says each step *"overwrites what it derives"*
  and that the `step` field *"records who last set the entry"*. ⭐ Is that superseded silently? What does a
  consumer of the merged export now read for an object two steps derived?
- ⛔ Does anything in the artifact prevent a collision that **should** happen (rule 6)?

### ⭐⭐ 4 · CAN EACH GUARD FAIL? ⛔ CONSTRUCT IT; DO NOT REASON ABOUT IT

- **`D5`** requires a structure-preserving residual per `value_kind`, with `AGREE`/`DISAGREE`/`UNDECIDED`,
  and says agreement requires the residual to establish no difference. ⭐ Apply it to **every row** of the
  committed export against its own reconstruction and report how many land in each status, and how many
  the rule cannot handle. ⚠ Most committed values do not support subtraction, and a zero matrix is not
  `== 0`.
- **`D6`** compares the imported reconstruction against the merged reconstruction using an independently
  formed claim join. ⭐ Construct an undeclared value write, an undeclared metadata write, and a redirected
  reference; show which operand each appears in.
- **`D8`** derives the required population from the spec and the declared sweep, "independently of the
  emitter and writer". ⭐ **Is that independence real, or does the required population come from the same
  place the emitter does?** Construct an object omitted *before* emission and show whether it is visible.
- **`D9`** permits publication only when all same-claim comparisons are `AGREE`. ⛔ So a genuine cross-step
  **disagreement blocks publication**. ⭐ Is that right? A disagreement is a finding, and the record is
  supposed to carry findings. Argue it and say what must become true.
- **`D4`**: a reference carries both locator and target claim. ⭐ Construct a `dimension_key` redirect and
  confirm which check catches it.

### ⭐⭐ 5 · `D10` — THE UPSTREAM REACH

`D10` says S10 must be regenerated before S11 may publish, and that the regeneration *"reaches upstream as
far as the imported population requires"*.

- ⭐ **How far is that, measured?** Does it reach S9? Report what the imported population actually requires.
- ⛔ Is the reach **bounded**, or does it recurse to the first step in the chain? If unbounded, say so.
- ⭐ `D10` allows S11 to compute and report but not publish until that regeneration passes its own legs. ⛔ Is
  a computed-but-unpublished S11 run useful, or does something downstream need the export to exist?

### ⭐ 6 · THE RULES

- **Rule 5 — ⛔ no measured physics in a builder-facing artifact.** ⭐ Scan every line for a value, count,
  membership, sign, spectrum, movement, or expected outcome a builder could converge on. ⚠ Two earlier
  attempts leaked here.
- **Rule 3 — name the object, ⛔ do not specify the recipe.**
- **Rule 2 — both operands and the residual, then guard**; ⛔ no check that audits its own input.
- **Rule 11 — cost is never a reason** to drop a control.
- ⭐⭐ **PHYSICS, ⛔ NOT PROCESS (standing user instruction).** Report an item only if it catches a way the
  physics could be wrong. ⛔ **If an item is ceremony, say so — that is a finding.** ⚠ The user's stated
  goal: *"simple here is usually better than more process"*. ⭐ Price `object_claim` against everything it
  removes and say whether the design is net simpler or net heavier.

### ⭐ 7 · WHAT IS MISSING

⛔ The most expensive defects in this project have been **absent computations**. What could go wrong in this
export chain that **no** item in the list would detect? ⚠ Name it concretely, with the mechanism.

## ⛔ Constraints on how you run

- ⛔ Read-only on the working tree. Copy to `/tmp` and modify the copy.
- ⛔ `timeout 600` on every run; ⛔ never raise it. ⛔ No Mathematica this round.
- ⭐ Save every script and its literal stdout to named absolute paths and report them. ⛔ A prose derivation
  will be discarded.

## Report

For each finding: **what is wrong**, **the mechanism by which a wrong result would survive**, **the literal
output that shows it**, and **what must become true instead** — ⛔ not a rewrite of the list.
⭐ Separate blocking from non-blocking. ⛔ "Nothing found" with no ablation behind it is not a result: say
which computation could have found something and did not.
