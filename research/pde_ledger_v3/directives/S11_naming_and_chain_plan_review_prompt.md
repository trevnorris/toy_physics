# Independent review — the naming-and-chain plan, before any of it is executed

## Artifact under review

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_naming_and_chain_plan.md`

It is **orchestrator-written**. It proposes a route through three interacting defects and states four
alternatives. **Nothing has been executed.** You are one of two independent legs; the other is not visible
to you.

⚠ The plan decides **what gets rebuilt and in what order**. A wrong choice here costs engine rebuilds and
can silently corrupt an export chain that later steps consume.

## ⭐⭐ REQUIRED READING ORDER — ⛔ do not read the plan first

Reading the plan first anchors you to its framing, which is the thing under test.

1. `research/pde_ledger_v3/DEFECT_REGISTER.md`, entries `C19` and `C20`.
2. `research/pde_ledger_v3/S9_REWRITE_PLAN.md` — the export-chain design, and `D1`–`D12`.
3. `research/pde_ledger_v3/scripts/S10_exports.py` — ⭐ **import and inspect it programmatically**, do not
   read it; it is large.
4. `directives/S10_SHARED_PHYSICS.md` and `directives/S11_SHARED_PHYSICS.md` — enough of `§6`/`§8` in each
   to see what each names its objects.
5. **Write down, before step 6:** how *you* would sequence these fixes, and what you would step-scope.
   ⭐ Keep that list.
6. **Only now** read the plan.

## What to check

### 1. ⭐⭐ IS THE CENTRAL ARGUMENT TRUE? — everything else follows from it

The plan claims: **unifying the emitted vocabulary across steps would turn one collision into many**,
because S10's `MAIN` and S11's `MAIN` are different actions running a similar question list.

⭐ **Test it, do not reason about it.** Take S11's spec quantity names and S10's emitted `MAIN` quantities,
apply the plan's own key transform, and report how many keys would collide **if both steps emitted the
spec's names**. Paste the literal stdout.
⛔ If the number is small, the plan's core justification fails and Alternative A becomes correct. ⭐ Say so.

### 2. ⭐⭐ IS THE TWO-NAMESPACE SPLIT RIGHT?

The plan separates the **emitted tag name** (py↔wl, within a step) from the **`LEDGER` export key** (across
steps), and says conflating them caused the earlier critical defect.
⭐ Is that split real, or is it a distinction without a difference given the comparator joins on tag names
and the chain joins on keys? ⚠ Does anything in `S9_REWRITE_PLAN.md` or either spec **require** them to be
the same string? ⛔ Quote it if so.

### 3. ⭐⭐ IS DECISION 2's BOUNDARY DECIDABLE AT ALL?

The plan says step-specific results get step-scoped keys while "objects genuinely carried between steps"
keep flat keys — and then explicitly **declines to draw that line**, deferring it.
⚠ **Is that line drawable by a mechanical rule?** ⭐ Try to state one, and test it against
`S10_exports.LEDGER`: which of its 617 entries would be flat, which step-scoped, and does the rule need a
human judgement per object? ⛔ If it does, the plan defers its hardest problem into a decision list that has
already failed once on exactly this.
⚠ Check what it costs: `S9 → S10` reportedly corroborated through 3 overwritten-and-agreeing rows.
⭐ **Verify that number**, and report what the plan's rule would do to those rows.

### 4. ⭐ IS THE ORDERING JUSTIFICATION SOUND?

The plan defers S10's retrofit on the grounds that **S10's computed result is not wrong** — the joined
objects are the right objects, and only the name is unpoliced.
⭐ Is that true? ⚠ Can you construct a case where a name that both engines apply to the same object, but
which the spec did not authorise, produces a **wrong claim** in a step record or a **false agreement** in
the comparator? ⛔ If you can, the deferral is wrong and the retrofit is urgent.

### 5. ⭐ ARE THE ALTERNATIVES FAIRLY STATED, AND IS ONE MISSING?

⭐ Each alternative should be recognisable to someone who favours it. ⚠ Is any straw-manned, or is its cost
overstated to make the proposal look better? ⛔ Name it.
⭐ **Is there a fifth option the plan does not list?** In particular, the plan raises amending
`S10_SHARED_PHYSICS.md` to match its engines and declines to propose it — ⭐ say whether that is the right
call and why.

### 6. ⛔ DOES THE PLAN LEAK, OR COMMISSION DAMAGE?

⛔ Does anything in it state or imply what a computation will produce? ⚠ **The test is derivability, not
literal presence.**
⚠ Does any proposed step touch a correct artifact more broadly than the defect requires? ⭐ Name what would
be damaged.

### 7. ⭐ WHAT IS MISSING?

⭐ Compare your step-5 list. What must be settled that the plan leaves open? ⚠ The plan's *"does not
decide"* section is a deliberate exclusion list — ⛔ is anything in it actually load-bearing for step 4 and
therefore wrongly deferred?

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, a way a **claim would be
unsupported**, or a way this plan would cause **wrong or corrupted artifacts** to be built. ⛔ Not style.

## Method

- ⭐⭐ **Quote both sides for every finding**: the plan's text, and the source it fails against.
- ⭐ For any claim about a file, a line, or a count: **state what you opened or ran, and what you found.**
  ⛔ A path that resolves is not a source.
- ⭐ Checks 1, 3 and 4 are settled by **running code**. Write the script, run it, paste literal stdout.
  ⛔ Prose is not accepted for those.
- ⛔ Do **not** edit anything and ⛔ do not execute any part of the plan. Read-only.
- ⭐ End with: which of your step-5 items the plan handles, which it misses, **which alternative you would
  choose and why**, and what you checked that could have failed and did not.
