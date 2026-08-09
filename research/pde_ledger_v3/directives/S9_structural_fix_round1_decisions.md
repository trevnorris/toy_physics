# S9 structural export — fix round 1. Decision list.

**You author the directive prose and apply the change.** Read `CLAUDE.md` first.

Files: `scripts/S9_light_requires_shear_sympy_audit.py` · `scripts/S9_exports.py` (generated) ·
`scripts/out/S9_light_requires_shear_sympy_audit.out` (regenerated) · your directive at
`directives/S9_export_structural_directive.md` (update it). ⛔ Nothing else. ⛔ Not the `.wl`.

⭐⭐ **This round REMOVES code. It adds no mechanism.** If you find yourself adding a check, a guard or a
helper, stop and say why.

---

## F1 · Delete the new in-run check battery.

Two independent legs measured every check added in the previous round and none can fail:

- Reverting the class gate to its previous single class — **undoing the entire change** — leaves every new
  assertion satisfied and the run at exit 0, because an `all(...)` over an empty collection is vacuously
  true and a sum over an empty generator is zero.
- Doubling the shear modulus in the main action moves hundreds of export lines, and the residual named for
  preserving existing entries still reads zero. Flipping a previously exported entry's class: also zero.
- The selection residual compares a list against an identity copy of itself, so it fires only for a filter
  inserted at one line nobody would write it on. Placed where anyone would — inside the classification —
  an entry silently leaves the ledger at exit 0.

**Decision: remove those emissions, their residuals and their asserts.** ⛔ Do not repair them and ⛔ do not
replace them with anything.

⚠ **The reason is structural, not a quality complaint.** The export writer's operands are its own inputs,
so no check inside the run can establish that the export is right. What establishes it is a **diff against
the committed baseline**, from outside the run. That is where the claim lives, and the plan already says so
for chain integrity.

⭐ **Keep** the pre-existing round-trip, count and class-tally emissions that predate the previous round.

## F2 · `identity3` must not reach the ledger, and its component count must be read.

A 3×3 identity is a construction at the engine's component count, ⛔ not a parameter a consumer can use —
a consumer working at another count needs its own, and a fixed-count matrix in a ledger that reads as
general is a false record.

Measured: its construction **types** the component count where the computation already has it, and reading
it instead leaves the output byte-identical.

**Decision: read the count rather than typing it, and give the object a classification that reflects what
it is** — a constructed object, not a declared parameter. ⛔ Do not add a name-based exclusion to the
export path. **Report the classification you chose and why.**

## F3 · Engine machinery is `CONTROL`.

Decided by the user, 2026-08-08. The tag-uniqueness state, the name-substitution table, the standard-name
collection, the class-assignment map, the production-metadata key and the action package table are all
`CONTROL`. ⚠ One declaration in the previous round was tagged `DERIVED` instead and is a Python string, not
a computed physics object — **fix that one to match.**

⚠ **State plainly in your directive that `CONTROL` now carries two meanings in this engine** — an ablation
control coefficient, and engine machinery — so a reader of the knob inventory is not misled.

---

## Constraints

- ⭐ Every entry that survived the previous round keeps its key, value, dim, class and step, **except**
  whatever F2 removes.
- ⭐ The exact-`srepr` round trip still covers every exported record.
- ⛔ The derivation, action, ansatz, assumptions and every computed value stay untouched.
- ⭐ **Establish preservation by diffing the generated export against the committed one** and report that
  diff. ⛔ Not by an emission.
- ⭐ Report the complete `.out` diff against the committed baseline.
- ⭐ **Report how many lines this round removed and how many it added.**

## Acceptance

⛔⛔ **A test that passes on a weaker fix is not a test.** Show what the deletions do under a change that
should be caught, and say plainly what is now caught by the run and what is only caught by the diff.

⛔ Do not state in the directive or a comment what any count or tally comes out as. Emit it and let it be
read.

## Report, do not fix

If removing a check leaves something genuinely unpoliced that you think matters, **name it and stop** —
that is a finding for the step record, ⛔ not a reason to build a replacement.
