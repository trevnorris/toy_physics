# Independent review — W3, the registry's bound dimension law

⚠ **Warnings — read these before you open anything.**
- `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.
- `research/pde_ledger_v3/directives/W3_registry_D_laws.md` (the build directive) and
  `research/pde_ledger_v3/reduction/W3_REGISTRY_D_LAWS_REPORT.md` (the builder's report) **both state the
  dimension laws this review is supposed to adjudicate.** ⭐ Derive yours **first** and save it to a named
  absolute path; ⛔ only then open either file. An artifact can satisfy its directive and still be wrong.

## Artifact — Codex-written, first round, uncommitted

Repository: `/var/projects/toy_physics`. Everything below is under
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/`.

| file | state |
|---|---|
| `dimension_laws.py` | new |
| `dimension_law_check.py` · `dimension_law_able_to_fail.py` · `test_dimension_laws.py` | new |
| `registry_read.py` · `registry_schema.yaml` · `quantities.yaml` | modified |
| `dimensional_homogeneity_gate.py` · `show_reduced.py` · `README.md` | modified |
| `relations.yaml` | modified — ⚠ the builder claims **version bump only** |

Some changes are staged and some are not ⇒ use `git diff HEAD -- research/pde_ledger_v3/reduction/` and
`git status`, ⛔ not bare `git diff`.

## What it is meant to do

`reduction/quantities.yaml` declares each registry quantity's dimension as an **LTM exponent vector**
(length, time, mass). Three quantities live on a brane whose dimension `D` is itself a registry quantity
(`Q.brane.D_brane`), so their true exponents are **functions of `D`** — but the schema could previously
only hold three bare integers, and recorded the `D = 3` evaluation with nothing saying so. W3 introduces a
**bound dimension law**: an exponent component may be an expression in a symbol bound to a structural
quantity, and a consumer may ask for the law or for its evaluation at a stated `D`.

## ⭐⭐ ① The physics — derive it yourself, in a script, before reading anything

⛔ **Do not accept the repository's answer, the directive's answer, or the report's answer.**

Work out from the model what the dimensions of the brane inertia density, the brane first-gradient modulus,
and the brane compression modulus **must be on a `D`-dimensional brane**, and what the two speeds built
from them must be. Then compare against what `quantities.yaml` now declares and what the bindings point at.

⭐ Sources of truth to read first: `reduction/relations.yaml` (5 relations, 14 quantities, 2 sectors),
`reduction/quantities.yaml` at HEAD (`git show HEAD:research/pde_ledger_v3/reduction/quantities.yaml`), and
the S9/S10 step records under `research/pde_ledger_v3/steps/`.

⛔⛔ **A PROSE DERIVATION IS WORTH NOTHING.** Write your own derivation script **before** opening the
artifact, and save **both the script and its literal stdout** to named absolute paths. Without these your
derivation claims will be discarded.

⭐ Two specific things to compute rather than assume:
- Does `D` genuinely **cancel** in each speed, so the speeds are `D`-independent — or does it cancel in one
  and not the other? Check it against **R1, R2.h0, R2.xi_h, R4 and R5**, not just the pair it is claimed on.
- Is every quantity that *should* carry `D`-dependence now declaring it? ⭐ Enumerate the brane sector
  yourself and say which entries you would have bound and which the artifact bound. ⛔ A missing binding
  looks exactly like a correct constant.

If your answer is *"nobody is wrong"*, say precisely what is missing and where it belongs.

## ⭐⭐ ② Is the law load-bearing, or decoration?

⚠ `quantities.yaml` appears to retain a plain exponent triple **alongside** the law. The builder states the
triple is only the checked evaluation at reference values and is **never** a fallback.

⛔ **Test that, don't read it.** Make the law and the retained triple **disagree**, load the registry, and
report which one every consumer actually uses — the gate, `show_reduced.py`, `engine_output_checks.py`,
`acceptance_check.py`. ⭐ If a disagreement loads silently, or if the triple wins anywhere, the law is
decoration and every downstream green is uninformative.

⭐ Same question for evaluation: can a consumer obtain a dimension **without** stating `D`? Try it. If an
unstated binding resolves to a default — a reference value, a fallback, a cached triple — that is the exact
silent `D = 3` specialisation this change exists to remove, reintroduced one layer down.

⭐ And: what happens when a binding names a quantity that does not exist, is not an integer, or is itself
`D`-dependent? Each must be `UNDETERMINED` or a hard error, ⛔ never a presumed value.

## ⭐⭐ ③ The regression bar — reproduce it, ⛔ do not take it

The acceptance criterion is: **every result that exists today at `D = 3` is unchanged.** The orchestrator
reports having observed, after the change:

- `dimensional_homogeneity_gate.py` → `POPULATION_COUNTS: HOMOGENEOUS=5 INHOMOGENEOUS=0 UNDETERMINED=0`,
  `PASS`, exit 0
- `engine_output_checks.py --config checks_S9.yaml` → `CROSS_ENGINE: agree=12`,
  `REGISTRY_RESIDUAL: nonzero=1`
- `relations.yaml`'s diff is a version bump and nothing else

⛔ **These are the orchestrator's observations, not findings.** Reproduce each yourself from a pristine
`HEAD` tree and report agreement or contradiction with literal stdout on both sides.

⭐ To get the baseline: `git archive HEAD | tar -x -C <your-temp-dir>` — the checks resolve the registry
root from the config path, so a bare copy of `reduction/` alone will fail; extract the tree.
⛔⛔ **NEVER `git checkout`, `git stash`, or `git restore` anything.** The working tree holds uncommitted
work that a revert destroys.

⭐⭐ **And then ask the question that matters more than the diff:** *does the baseline command reach the new
code at all?* Byte-identical output is worth **nothing** if the D=3 path reads a cached triple and never
evaluates a law. ⛔ Prove reachability — instrument, break, or delete the law-evaluation path and show what
moves. If nothing moves, say so; that is the finding.

## ⭐ ④ Able-to-fail — run the cases, then defeat them

`dimension_law_able_to_fail.py` ships four cases (wrong law, wrong binding, unbound, unresolvable). Run
each and paste literal stdout with exit codes.

⭐ Then attack them:
- ⛔ **A test that "covers" an invariant demonstrates it on one example.** Build the **weaker** wrong laws —
  a wrong law that still cancels in the relations, a wrong binding that happens to resolve to the same
  integer — and show whether each case still fires. In this repository a half-fix has three times passed
  the new test, the whole suite, the full battery **and** produced byte-identical output.
- ⭐ Is a wrong law caught because it is wrong, or only because it breaks a relation? A quantity that
  appears in **no** relation carries a law nothing constrains — check whether one exists.
- ⭐ Name a defect class none of the four cases would catch.

⚠ **A blindness already exists and W3 does not claim to fix it:** a common-mode shift applied to all three
brane constituents leaves all five relations homogeneous. ⭐ Verify that claim is stated honestly in
`README.md` and that nothing elsewhere invites a green gate to be read as evidence about **absolute** brane
dimensions. ⛔ Do not treat the blindness itself as a finding; ⭐ do report any place the artifact overstates.

## ⛔⛔ ⑤ FORM ablation — mandatory, not optional

⭐ A **coefficient** rescale tests arithmetic; only a **FORM** change tests physics, because scaling never
leaves the family. Change the **structure**: collapse two independently-bound symbols into one; flip the
sign of the `D` coefficient in exactly one law; make two quantities share a binding that should be distinct.
Re-run the gate and the S9/S10 configs and report the **literal** diff.

⚠ **Measured in this repository:** a script emitted a physics conclusion as a typed sentence with no
computation behind it and the ablation produced **byte-identical output**; eight fidelity legs missed it.
⇒ ⭐ **Ask of every claim in this artifact: WHICH LINE COMPUTED THIS?** Give the line number, or report the
claim as uncomputed.
⚠ An `assert` placed **before** the value it guards hides this — a perturbation strong enough to flip the
check kills the process, so you see only PASS-or-crash. ⭐ Report every such `assert`.

## ⭐ ⑥ Do the new tests bind, and does the version bump?

- ⚠ A sibling artifact in this same tree shipped tests that reported `Ran 0 tests` under `unittest` because
  they were free functions. ⭐ Run `test_dimension_laws.py` under **both** `pytest` and
  `python -m unittest` and report the collected count from each.
- ⭐ Mutate `dimension_laws.py` and confirm the suite goes red. A suite that stays green under mutation is
  not a suite.
- ⭐ The schema/convention version was advanced. Is the new version **enforced**? Load a HEAD-form
  `quantities.yaml` under the new loader and a new-form one under the old rules; report what happens. An
  unenforced version bump is a comment.

## ⭐ ⑦ Rule 2

Does every script print **computed objects** — both operands and the residual — and then guard? ⛔ Any
status typed rather than computed? ⛔ Any conclusion emitted as an unconditional literal? ⛔ Any check whose
expected value lives inside the artifact it checks?

## Method and operational constraints

⭐ Order: (1) read the sources of truth and write your own derivation script; (2) read the diff and the
artifact; (3) read the directive and the builder's report **last**.

- ⭐ Work under a temp directory of your own choosing, outside the repository and outside
  `/tmp/claude-1000/`. ⛔ **Do not modify the working tree.** Report the absolute paths you used.
- ⛔ Wrap every long run in `timeout 600`. A 600 s hit is a failed check — report it and move on;
  ⛔ never raise the timeout. ⚠ `engine_output_checks.py --config checks_S10.yaml` is the slow one.
- ⭐ Save **every** script you write and its **literal stdout** to named absolute paths, and report those
  paths. Claims without them are discarded.
- ⚠ Exit codes carry no signal on the S9/S10 configs — both exit 2 in the healthy state. **Judge by
  counters.**

## Physics filter

> Report a finding only if it changes what the project computes, what it may claim, or what method it
> adopts. ⛔ Do not report "the script would be wrong on a different input."

⭐ Worth most: *"the law is never read — the triple still decides"* · *"a dimension can still be obtained
without stating `D`"* · *"this quantity should be `D`-bound and is not"* · *"the D=3 regression is
uninformative because the baseline never reaches the new path"* · *"the able-to-fail case passes on a
weaker wrong law."*

⚠ A leg returning *"nothing survives the filter"* is weak evidence on its own. State what you checked, what
you could not, and what would have had to be true for you to find something. ⛔ Do not manufacture findings.

Rank most-severe first: **claim · evidence (quotation with file:line, or a literal command and its output)
· what must change.**
