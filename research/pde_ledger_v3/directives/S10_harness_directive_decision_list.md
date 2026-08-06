# Decision list — revise the harness rebuild directive

**You are the author of the revision.** The orchestrator wrote the current directive and it failed two
review legs twice. The orchestrator's job is this decision list; yours is the document.

**Edit in place:** `research/pde_ledger_v3/directives/S10_harness_rebuild_directive.md`.
Do not commit. Do not touch anything else under `research/`.

Every decision below is a *requirement on the revised text*, not a patch to apply literally. Where a
decision says "name X", the revised directive must state X — not instruct a future builder to name it.

## Context you need

The directive tells a builder how to rebuild `reduction/engine_output_checks.py` and its two configs so
that two independently-constructed engines can be compared. Round 1 found the directive had been drafted
against reverted code (`92461853`, tag `wip-2026-08-05-unreviewed`, reverted by `ab77f25d`); it was
re-anchored to HEAD. Round 2 found the items below. Both rounds' full reports are in
`_scratch/S10_harness_directive_review_fold.md` and the leg logs referenced there.

Governing failure mode throughout: **a layer that reports success without having compared anything.**

---

## D1 — A1 states a substring count as an object count

`78 / 130 LAGRANGIAN tags` are substring hits; most are Q6 dimension metadata *about* the Lagrangian.
Measured truth: the action Lagrangian is **13 per engine**, under different names —
WL `*_Q1_LAGRANGIAN`, PY `*_Q1_LAGRANGIAN_EXPANDED`. A separate `*_AVERAGED_LAGRANGIAN` family also
exists, 13 per engine.

The revision must name the tag pairs (or a declaration schema that names them) and the expected
cardinality, so that configuring a companion metadata tag as "the action" is not conforming. State
explicitly whether the averaged family is in scope; if it is, it gets its own declared pairs.

## D2 — A2 tells the builder to name the object instead of naming it

"Name the object both engines are to be compared on" is an instruction to the builder to make the
decision. Make it: the compared object is the **ordered variational residual vector `δL/δu_i`** — `lhs −
rhs` for an equation, the residual directly where an engine already emits one — with **order and scale
preserved**.

Rationale to carry into the text: a zero-locus/solution-set convention would make `x` and `2x` equal, which
would erase a wrong overall factor in the action and report agreement.

"Do not make the two sides equal by construction" is unfalsifiable advice. Replace it with an acceptance
test: a same-shape **unequal** residual pair must still report `DISAGREE`.

## D3 — A3 is factually wrong

`SHAPE_MISMATCH` is **already** an operational failure (`engine_output_checks.py:242`). The claim that it
"absorbs this silently" is false. What is actually silent is that **no configured action row ever reaches
a comparison verdict**. Require the report to distinguish configured action rows that reached a verdict
from those that did not, and to fail when the former is zero.

## D4 — B3 and B5 contradict each other

B3 requires per-engine control and dimension checks. B5's schema supplies **one tag per symbol**, and both
currently configured tags are WL-only, so a PY run would silently have no derived table. The declaration
schema must be **engine- and package-specific**, and a configured source that is missing or empty must
fail rather than default.

## D5 — B6 orders a warning while continuing with the wrong table

Reporting a cross-package dimension disagreement and then applying the other package's table anyway can
produce wrong dimensional verdicts precisely when the warning fires. Use each package's own vector, or
mark that package unassessable. A warning is not a disposition.

## D6 — B8 is stale, carried over from the reverted build

`self.comparisons` occurs **zero times** on HEAD; there is no per-site comparison counter to fix. Replace
B8 with the real HEAD defect it was groping at: `checked` / `homogeneous` increment once per successfully
walked tag, so payloads such as `2`, `x`, and `Element[x, Reals]` report `checked=1 homogeneous=1` having
had no homogeneity comparison performed at all. The requirement is that a tag on which no comparison was
possible must not be counted as one that was checked.

## D7 — B9 has three inaccuracies

- **Bare string payloads: there are none.** Zero top-level bare strings; 71 unparsed *containers* contain
  quoted strings. Restate accordingly.
- **ConditionalExpression already preserves and compares its condition**, and its dimension walk raises a
  visible error. The current requirement 7 describes a defect that is not there. Either drop it or restate
  it as whatever remains true after measurement.
- **The one real literal** (`WL_S10_RUNTIME_SECONDS`) is mentioned in prose but has no numbered
  requirement, so a builder can leave it unparsed and still claim B9 complete.

## D8 — B10's "boolean" is undefined, and its acceptance proves the wrong path

Define the term. The defect class is truthiness matching between Python booleans and numbers
(`symbolic_equal(True, 2) = True`, `symbolic_equal(False, 2) = False`); `_is_boolean_value` also treats
relationals as boolean, which after B9 and A2 covers a large and growing population. State which
population each of the two counts covers.

The acceptance demonstration must use a **nested-only** synthetic pair — selected-operand count zero, tree
count nonzero — because the top-level path is shape-gated and an implementation with
`anywhere_count == selected_count` would otherwise pass while the real false-`AGREE` path survives.
Measured example of that path: `symbolic_equal({a: True}, {a: 1}) = True`.

## D9 — the form ablation as written cannot be performed, and does not select form controls

The directive freezes the engines and their committed output, so the only available mutation is payload
substitution — a coefficient control under a form-shaped name.

It does not need an engine change. Genuine form controls are already in committed output
(`scripts/S10_brane_mode_spectrum_sympy_audit.py:51-55`): the stiffness form is the third `Package`
argument, and only `XFORM_FULLGRAD` ("fullgrad") and `XFORM_DIVONLY` ("divonly") differ from `MAIN`
("curl"). `XFORM_SIGNFLIP` and `XFORM_ANISO` are **`curl` with a changed coefficient**, despite the
`XFORM_` prefix.

Restate the ablation as: the action comparison must respond across packages whose **stiffness functional**
differs, and must not be moved merely by a coefficient change. Select by the stiffness-form argument,
never by tag prefix — and say so, because the prefix is the trap.

## D10 — the acceptance battery has no oracle

"Any counter changed" and "improves" are not oracles: total-tag or kind-bucket movement can satisfy them
while coverage, parsing or observation stays dead, and "improves" is undefined when counters move in
different directions. **Each ablation names its target counter or status and the required direction, or the
new operational finding it must produce.** The coverage item must show its formula and its declared
denominator.

Acceptance must also force the items it currently leaves untested — D4's per-engine checks, D5's
per-package tables, and D6's replacement — since a build can otherwise ship green with all three surviving.

## D11 — duplicate and empty comparisons inflate the counts (unnamed in both rounds)

`N` identical configured rows produce `N` `AGREE` results from **one** distinct tag pair, and empty lists,
mappings and multisets report `AGREE` without comparing anything. Require distinct-pair counting, and a
per-object nonempty/cardinality declaration so an empty payload cannot agree.

## D12 — `registry_residual` is the same defect as B1, unnamed

It is permanently `configured=0` on both configs and this is never treated as a declaration gap. Bring it
under the same rule as B1.

## D13 — record what may not be reachable, and do not promise agreement

A review leg measured that even after an `Eq → residual` normalisation, the `MAIN D3` Euler–Lagrange
residuals do **not** match: one engine emits `OpaqueDerivative[...][x1,x2,x3,t]`, the other
`Derivative(u(t,x1,...))`. Representation alignment is a further problem, and B9's "no hardcoded coordinate
ordering" constrains how it may be solved.

The revised directive must **not** read as promising that the action will agree. A disagreement is a
finding. State that Part A's requirement is that the action be *compared*, and that an honest persistent
`DISAGREE` is a valid outcome to report, not a failure to iterate away.

---

## Constraints on the revision

- Keep it anchored to HEAD. Every factual claim must be one you have measured against the working tree or
  committed output — round 1 failed precisely because claims were inherited rather than measured. Where you
  cannot measure a claim, remove it.
- **Do not state what any counter will read after the fix.** Present-state diagnostics are fine; post-fix
  expected values are not. Report-what-you-read is the standing instruction to the builder.
- Prefer required semantics over prescribed mechanism, except where a mechanism is itself the decision
  (the residual convention in D2, the shape declaration in D4).
- Undecorated prose. No emphasis glyphs.
- The document is for a builder who has not read this list.

## Report back — under 30 lines

1. One line per `D1`–`D13`: how the revision satisfies it, with the section it lands in.
2. Any decision above you believe is **wrong or unmeasurable** — say so rather than drafting around it.
3. Anything you measured that contradicts this list.
