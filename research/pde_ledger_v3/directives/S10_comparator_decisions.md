# The cross-engine comparator (F-2). Decision list.

**You author the directive prose and build it.** Read `CLAUDE.md` first.

⭐ Build this **independently**. ⛔ Do not ask for, and do not be shown, any prior comparison output — an
orchestrator prototype exists and is deliberately withheld. **If your tally differs from it, that
disagreement is a finding and it is why you are building this separately.**

---

## What it is

D11's claim is that both engines emit **the standard name of the object**, so cross-engine comparison is a
**join on the name**. ⛔ **Nothing in the repository performs that join.** Every engine binds a standard
name to an object through its own table, and a review leg re-pointed one entry so the light cone became an
`ω²` instead of a speed squared — **engine exit 0, roundtrip residual 0, mapped diff exit 0, export
regenerated with the same entry count. Every check in the repository passed.** ⚠ The Wolfram side has the
identical hole.

**Decision: the join exists as one committed script.**

It reads **exactly two `.out` files** given as arguments. ⛔ No config file, ⛔ no YAML, ⛔ no per-step file,
⛔ no manifest, ⛔ no runner, ⛔ no hand-written name→name pair table — the pair table is the 3,121-line
artifact whose deletion this design exists to make permanent.

## What it must produce

⭐ **For every name both engines emit: both operands and the residual, then a guard** — rule 2. ⛔ A residual
asserted zero always prints `0` and carries no information.

⭐ It must also emit **what it could not compare, and why**, as its own category. ⚠ A row silently dropped
because the comparator could not parse it reads exactly like a row that agreed, and that is the failure
mode this instrument exists to remove. ⛔ **Never let an unparsed row count as agreement.**

⭐ Names only one engine emits are **inventory, not failure** — report the counts on each side.

## The correspondences it may apply — and the ones it may not

⭐ **Allowed, because they are surface syntax and not claims about the objects:** the two engines' symbol
spellings are related by a **mechanical** rule (snake_case ↔ lowerCamel), and their sequence and power
syntaxes differ. ⭐ Apply those mechanically.

⛔⛔ **Not allowed: absorbing a symbol-name difference that is not mechanical.** Where the two engines use
genuinely different names for the same object, that is a **D12 naming defect and the comparator's job is to
surface it** — ⛔ not to paper over it with a special case. **Report those pairs as output**; they are the
work list for the naming pass that follows.

⚠ **A non-canonical object must not be compared as though it were canonical.** A nullspace basis is one:
S10 previously produced 11 rows reporting disagreement on representation alone. Decide what the comparator
does with such an object and **state the decision**; ⛔ do not let representation masquerade as physics.

## ⛔⛔ It must be demonstrated able to fail

⭐ **Re-point one standard name in one engine to a neighbouring object and show the residual move.** Save
the ablation and its literal stdout. ⛔ Without that demonstration the comparator is decoration — and a
comparator that never compares has already passed this project's acceptance criterion once.

⛔ Do not state anywhere what any count, tally or agreement total comes out as, in the script or the
directive. **Emit it and let it be read.** ⚠ An acceptance criterion naming an expected agreement count is
a target a builder converges on; that is why none is given here.

## Scope

New script under `research/pde_ledger_v3/scripts/`, plus its committed stdout under `scripts/out/`, plus
your directive under `directives/`. ⛔ Do not modify either engine — this round measures them.
⛔ Do not start `wolframscript`; both transcripts are committed.
