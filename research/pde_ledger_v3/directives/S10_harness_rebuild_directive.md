# Rebuild the cross-engine harness — ⛔ a REBUILD, ⛔ not a re-land

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/engine_output_checks.py`
**Config:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/checks_S10.yaml`, `checks_S9.yaml`
**Tests:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/test_engine_output_checks.py`

⛔ **Do not commit.** ⛔ Do not edit anything under `steps/`, `paper/`, or `REBUILD_HANDOFF.md`.
⚠ A previous attempt violated exactly that and landed a step record asserting conclusions.

⚠ **Prior art, ⛔ not a baseline:** an earlier attempt is reachable at
`git show wip-2026-08-05-unreviewed` (commit `92461853`). ⭐ **Read it for what it tried.** ⛔ **Do not
assume any of it is correct** — two independent reviews found it not re-landable, and its own headline
numbers cannot be reproduced from the tree.

---

## ⭐⭐⭐ WHAT THIS INSTRUMENT IS FOR — read this before any line of code

Two engines compute the same physics independently and emit tagged lines. ⭐ **This harness decides whether
they agree, and whether each engine's output actually depends on a computation.**

⛔⛔ **THE GOVERNING FAILURE MODE: a layer that reports success without having compared anything.** ⚠ Every
counter in this file has read as **health** while measuring **nothing** at some point in its history.
⇒ ⭐ **For every layer you touch, answer in a comment: what does this print when its input is empty,
mis-keyed, or filtered to nothing, and is that distinguishable from a clean pass?**

⛔ **A checker that mis-parses two expressions into agreement is worse than no checker**, because it
manufactures confidence.

---

## ⭐⭐⭐ THE ONE STRUCTURAL CHANGE — ⭐ it closes three defects at once

⛔⛔ **Coverage is currently measured against WHAT THE ENGINE HAPPENED TO EMIT.** ⇒ three consequences,
**all three measured**:

- ⛔ **Deleting output improves the score.** Removing all but one main tag drove the reported uncovered
  fraction to zero and turned the guard **green**.
- ⛔ **A declared control absent from an entire dimension is silently dropped** — deleting every tag of one
  control package at one dimension left **every counter byte-identical**, with nothing reported missing.
- ⛔ **A "packages matched" count reads as full coverage** while a package is inert over most of the sweep.

⭐⭐ **THE FIX: the config DECLARES the expected control matrix, and coverage is measured against the
DECLARATION.**

```yaml
# the (package, dimension) pairs this step EXPECTS to exist
control_matrix:
  MAIN:      [2, 3, 4, 5]
  <control>: [3, 4]
  ...
```

⇒ ⭐ **An expected cell that emitted nothing is a REPORTED GAP, ⛔ never an absence.** ⭐ Deleting output
then makes the score **worse**, which is the only orientation a coverage measure may have.

⚠ **Where the declaration and reality disagree, that is the finding** — ⛔ do not "fix" it by trimming the
declaration to match what is emitted. ⭐ Report it.

---

## ⛔ THE DEFECTS — ⭐ each one measured, ⛔ none hypothetical

### ⛔⛔ D1 — the whole guard set is switched off for most configs

`coverage_required=main_package is not None`. ⚠ Only one config declares `main_package`, so ⛔ **every new
guard silently does nothing for the others** — and destroying the control naming in such a config produced
a coverage fraction that read **better than the healthy run**.
⭐ **Fix:** ⛔ no guard may be conditional on the presence of an optional config field. ⭐ If a config does
not declare what it expects, ⛔ **that is itself a failure**, ⛔ not a licence to skip the check.

### ⚠ D2 — the boolean/number false agreement: ⛔ THE FIX IS DEFERRED, ⭐ but add a TRIPWIRE

⛔⛔ **DO NOT FIX `symbolic_equal` HERE** (user decision, 2026-08-05: this is post-S10 work).
⇒ `DEFECT_REGISTER.md#f7`.

⭐ **Why deferring is safe, and it was MEASURED rather than assumed:** the comparison kernel routes to an
equivalence test when **either** side is boolean, so a boolean compares equal to any nonzero number.
⭐ **Both configured steps have ZERO cross-engine exposure** — not one configured pair on either step has a
boolean on either side. ⇒ ⛔ **no agreement number is contaminated.** ⚠ The residual effect is confined to
the ablation layer, where it can only turn RESPONSIVE into INVARIANT ⇒ it **hides** evidence of
computation and ⛔ can never manufacture it.

⛔⛔ **BUT THAT SAFETY IS A PROPERTY OF TODAY'S CONFIG, ⛔ NOT OF THE CODE.** ⚠ Adding one boolean-valued
cross-engine pair later would silently start reporting false agreement, and ⛔ nothing would say so.
⭐ **So add the TRIPWIRE, ⛔ not the fix:** ⭐ **count the cross-engine comparisons in which either side is
boolean, and make a nonzero count an OPERATIONAL FAILURE naming the pairs.**
⚠ Today that count is **zero**, so the tripwire is silent — ⭐ **that is the point.** ⇒ the deferral stays
safe **by construction** instead of by one measurement that ages.

### ⛔⛔ D3 — one stray CAS message deletes the entire report

The parser treats any word followed by `:` as a tag, so a repeated solver message becomes a duplicate tag;
and the newer variant **raises** on a non-grammar line, which `main` swallows into a single stderr line.
⇒ ⛔ **the run prints nothing at all** — ⛔ not a partial report.
⭐ **Fix:** ⭐ **a non-grammar line is REPORTED AND SKIPPED, ⛔ never fatal, and never a tag.** ⭐ Emit the
count of skipped lines and the first few verbatim. ⚠ A run that silently ate output is the defect; a run
that says *"I ignored N lines, here they are"* is not.

### ⛔⛔ D4 — the dimension layer cannot fail on empty work

An all-integer ablation gave `compared=0` with **no** operational failure. ⚠ And the buckets **do not
partition** the tag set — a few hundred tags land in neither `compared` nor `vacuous` and are invisible on
the summary line.
⭐ **Fix:** ⭐ zero comparisons is a **failure**, and ⭐ **the buckets must sum to the tag total.** ⛔ Print
the arithmetic so the reader can check it.

### ⛔ D5 — `compared` is inflated by a check true by construction

The dimension walker counts a comparison at a power site even when the exponent is a **literal integer**,
where "the exponent is dimensionless" holds by construction. ⚠ Measured: a substantial fraction of
`compared` tags had **no** rule, relational, or additive comparison at all.
⭐ **Fix:** ⛔ a site that cannot fail does not increment `compared`. ⭐ Report the increment counts **per
site kind**, so the composition of the number is visible rather than asserted.

### ⛔⛔ D6 — the dimension table: ⭐ DECLARE the shape, ⛔ do not guess it

⚠ **This and `D8` were previously deferred to a "part B". ⛔ The split was wrong** — ⭐ the instrument
cannot be exercised **at all** without them, so they are in scope here. ⭐ The material below survived two
review legs; ⭐ treat it as specification.

**D6a — the payload is a LIST of per-coefficient dimension vectors** while the reader requires a single
triple ⇒ ⛔ it **raises**, which is why the run produces nothing. ⭐ Collapse a family to the vector its
members share, and ⛔ **FAIL if they do not share one** — the config maps the family to **one** symbol, so
a disagreement means the mapping is ill-posed. ⛔ **Do not take the first entry.**

**D6b — ⛔⛔ THE COLLAPSE MUST NOT BE A SHAPE HEURISTIC.** ⚠ **Measured by a review leg:** a square matrix
whose rows happen to agree collapses *successfully* and installs a **wrong dimension with no complaint**;
an all-zero one collapses to a zero vector. ⚠ **And a real payload in this step has exactly that shape**,
so this fires on live data. ⇒ ⭐ **the config declares which it is:**

```yaml
derived_dimensions:
  <symbol>: {tag: <TAG>, shape: family}
  <symbol>: {tag: <TAG>, shape: vector}
```

⭐ Keep the plain `name: TAG` form working as `shape: vector`, so the other config is untouched.
⛔ **Collapse only when the config said `family`.** ⭐ Reject a genuinely 2-D matrix **before** consulting
entry shapes. ⚠ A payload whose shape contradicts its declaration is an **error**, ⛔ not a coercion.

**D6c — point the config at the SYMBOLIC tags**, whose components are expressions in the brane dimension,
⛔ **not** the specialised ones. ⚠ The specialised payloads are numeric at a single dimension, so the
current config **silently applies one dimension's table to every package in the sweep.**
⛔⛔ **Do not restate this as "correct at every `D`".** ⭐ What is true: a review leg verified the two tables
give **byte-identical verdicts across every tag**, while also constructing a case where the numeric table
passes an accidental coincidence the symbolic one catches. ⇒ ⭐ symbolic is **strictly stronger**,
⛔ not infallible.

**D6d — the table is built from ONE package and applied to all.** ⭐ **Mandatory:** emit a report line
whenever another package's corresponding tag disagrees with the one used, naming both packages and both
vectors. ⛔ Not fatal — a form control may legitimately change a coefficient's dimension — ⛔ but ⛔ not
silent either.

### ⛔⛔ D8 — the reader cannot parse the ACTION. ⭐ This is the item that matters most.

⚠ **Measured: a few hundred payloads are UNPARSED**, and they include the Lagrangian and the
Euler–Lagrange system. ⇒ ⛔ **until this lands, the dimension and ablation layers never see the action at
all**, and the cross-engine layer cannot compare the objects everything else is derived from.

⭐ Six constructs the reader rejects, **all measured in committed output**:

1. ⭐⭐ **A derivative head of ARBITRARY arity applied to an arbitrary argument list.** The reader hardcodes
   the fixed shape of an earlier step, which swept a single dimension; ⚠ this step sweeps several, so both
   the arity and the argument list vary. ⭐ **This is the bulk of the unparsed payloads and it is what
   blocks the action.** ⇒ ⭐ Fix it **generally**: any arity, any argument list, ⛔ no hardcoded coordinate
   names or ordering.
2. **A set-membership assertion over the integers.** ⚠ The reader accepts only the reals. ⛔ **Do not "fix"
   conjunctions** — membership inside a boolean chain parses correctly today; ⚠ four forms were tested.
3. **Association syntax with a string key** — the tokenizer rejects the key.
4. **A marker head applied to a LIST argument.**
5. **A piecewise payload** — rejected as "not a scalar expression". ⚠ This one is the **other** config's
   single unparsed payload; ⭐ fixing it makes that run clean. ⛔ Change nothing else about that config.
6. ⚠ **A conditional-expression head**, which currently parses as a **generic function** ⇒ ⛔ it silently
   carries a condition into a dimension walk. ⭐ Handle it explicitly: take the expression, and ⭐ **emit
   the condition somewhere visible** — ⛔ never drop it silently.

⭐ **Acceptance for `D8`:** report the unparsed count **before and after**, and ⭐ confirm the Lagrangian and
the Euler–Lagrange system now parse **in both engines and at every dimension**.
⛔ Where a construct genuinely cannot be parsed, ⭐ the payload stays UNPARSED **with its reason** —
⛔ never silently dropped, ⛔ **never coerced into something that parses.**

### ⛔ D7 — the ablation instrument is broken in the direction of over-trust

`reduction/harness_ablation.py` reads counters with an **integer-only** regex, so a fractional coverage
value is truncated to `0` ⇒ ⛔ **it cannot read the number the coverage guard rests on.** It also uses a
budget above the project limit, and seeds corruption with a **per-process randomised** hash ⇒ ⛔ not
reproducible.
⭐ **Fix all three**, and ⭐ make it take the harness, config and outputs as **arguments** rather than
hardcoding tree paths, so it can be pointed at a candidate build.

---

## ⭐ ACCEPTANCE — ⛔ run these and paste literal output; ⛔ do not assert them

⚠ **Judge by the reported counters, ⛔ not by exit code.**

1. ⭐ **Both configs produce a full report.** ⛔ Not a partial one, ⛔ not an exception. Paste both.
2. ⭐⭐ **THE ABLATION BATTERY — ⛔ each must be DISTINGUISHABLE from a clean pass.** For each, paste the
   counters **before and after**:
   - delete every control tag at one dimension for one package
   - delete all but one main tag
   - rename the main package
   - rename one control package
   - corrupt one invariant payload to a different value
   - corrupt a **boolean** payload to a nonzero integer
   - replace every dimension payload with a bare integer
   - feed a file containing a stray solver message
   ⛔ **Any ablation that leaves the counters unchanged, or improves them, is a FAILURE of this rebuild.**
3. ⭐ **Coverage orientation:** show that **removing** output makes the coverage number **worse**.
4. ⭐ **Bucket arithmetic:** show the dimension buckets summing to the tag total, on both configs.
5. ⭐ **The unit suite**, with a one-line cause for **every** failure — real regression, already-wrong test,
   or a test this change should have updated. ⛔ Do not re-baseline a frozen expectation silently; ⭐ if you
   update one, say which and why.

## ⛔ WHAT NOT TO DO

⛔ Do not make a physics disagreement exit non-zero — ⭐ only **operational** failure does. ⚠ Otherwise a
builder iterating to exit 0 can make a disagreement disappear.
⛔ Do not trim a declaration to match emission. ⛔ Do not delete a test to make the suite green.
⛔ Do not touch the engines or their committed output.

## Report back — ⛔ under 40 lines, plus the acceptance output

1. One line per `D1`–`D8`: fixed / partially / not, with line numbers.
2. The ablation battery table.
3. ⭐ Anything in this directive you believe is **wrong**. ⭐ This is wanted.
4. ⛔ Do not report what any engine's values came out to be.
